
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_MFParallelForC.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuComplex.H>
#include<AMReX_Print.H>
#include<AMReX_TableData.H>
#include <AMReX_TinyProfiler.H>
#include <AMReX_GpuUtility.H>

#include <cmath>
#include <iomanip>
#include <math.h> 
using namespace amrex;
using MatrixDType = amrex::GpuComplex<amrex::Real>;
using Matrix1D = TableData<MatrixDType, 1>;
using Matrix2D = TableData<MatrixDType, 2>;

template<typename T>
class TD;

namespace MathConst
{
    static constexpr amrex::Real pi = static_cast<amrex::Real>(3.14159265358979323846);
}

namespace PhysConst
{
    static constexpr auto q_e   = static_cast<amrex::Real>( 1.602176634e-19 );
    static constexpr auto ep0   = static_cast<amrex::Real>( 8.8541878128e-12 );
//    static constexpr auto hbar  = static_cast<amrex::Real>( 1.054571817e-34 );
}

MatrixDType get_beta(amrex::Real gamma, int M, int J) 
{
   amrex::GpuComplex arg(0., -MathConst::pi*J/M);
   return 2 * gamma * cos(-1*arg.imag()) * exp(arg); 
}

MatrixDType conjugate(MatrixDType a) 
{
   amrex::GpuComplex a_conj(a.real(), -1.*a.imag());
   return a_conj;
}

void PrintTable(amrex::Vector< Matrix2D >& G_blkvec)
{
    int num_blocks = G_blkvec.size();

    for (int b = 0; b < num_blocks; ++b) /*loop over blocks*/
    {
        auto tlo = G_blkvec[b].lo();
        auto thi = G_blkvec[b].hi();
	

        auto const& G = G_blkvec[b].table();

        amrex::Print() << "block: " << b << "\n";
	amrex::Print() << "tlo: " << tlo[0] << " " << tlo[1] 
		       << " thi: " << thi[0] << " " << thi[1] << "\n";

        for (int i = tlo[0]; i < thi[0]; ++i) 
	{ 
            for (int j = tlo[1]; j < thi[1]; ++j) //slow moving index. printing slow
	    { 
                amrex::Print() << std::setw(12) << std::setprecision(2) << G(i,j);
            }
            amrex::Print() << "\n";
        }
    }
}

template<std::size_t N>
void MatInv_BlockTriDiagonal(const amrex::Vector< Matrix1D >& A_blkvec, 
                             const Matrix1D& B_data, 
                             const Matrix1D& C_data,
                             amrex::Vector< Matrix2D >& G_blkvec)
{
    auto const& B = B_data.const_table();
    auto const& C = C_data.const_table();

    int num_blocks = A_blkvec.size();
    amrex::Vector< int > cumulative_blk_size(num_blocks);

    amrex::Vector< Matrix1D > Xtil_blkvec(num_blocks);
    amrex::Vector< Matrix1D > Ytil_blkvec(num_blocks);
    amrex::Vector< Matrix1D > X_blkvec(num_blocks);
    amrex::Vector< Matrix1D > Y_blkvec(num_blocks);

    Matrix1D Xtil_gl_data({0},{N});
    Matrix1D Ytil_gl_data({0},{N});

    auto const& Ytil_gl = Ytil_gl_data.table();
    auto const& Xtil_gl = Xtil_gl_data.table();

    cumulative_blk_size[0] = 0;
    for(int b=1; b < num_blocks; ++b)
    {
	int prev_blk_size = A_blkvec[b-1].hi()[0];    
        cumulative_blk_size[b] = cumulative_blk_size[b-1] + prev_blk_size;
    }
    for(int b=0; b < num_blocks; ++b)
    {
	int blk_size = A_blkvec[b].hi()[0];     
        Xtil_blkvec[b].resize({0},{blk_size}); 
        Ytil_blkvec[b].resize({0},{blk_size}); 
           X_blkvec[b].resize({0},{blk_size}); 
           Y_blkvec[b].resize({0},{blk_size}); 
    }

    MatrixDType A_prev, Y_prev;

    for (int b = 0; b < num_blocks; ++b) /*loop over blocks*/
    {
        auto const& A    =    A_blkvec[b].const_table();
        //auto const& Ytil = Ytil_blkvec[b].table();
        auto const& Y    =    Y_blkvec[b].table();
        int blk_size = A_blkvec[b].hi()[0];


	if(b > 0)
	{
            int n = 0 + cumulative_blk_size[b]; /*id if array was not partitioned in blocks*/
	    int p = (n-1)%2;

            Ytil_gl(n) = C(p) / ( A_prev - Y_prev );
            Y(0)    = B(p) * Ytil_gl(n);

            //Ytil(0) = C(p) / ( A_prev - Y_prev );
            //Y(0)    = B(p) * Ytil(0);
	} 
	else 
	{
	    Ytil_gl(0) = 0; 
	    //Ytil(0) = 0; 

	    Y(0)    = 0; 
	}

        for (int e = 1; e < blk_size; ++e) /*loop over elements of a block*/
	{
            int n = e + cumulative_blk_size[b]; /*id if array was not partitioned in blocks*/		 int p = (n-1)%2;
            
	    Ytil_gl(n) = C(p) / ( A(e-1) - Y(e-1) );
            Y(e)       = B(p) * Ytil_gl(n);

            //Ytil(e) = C(p) / ( A(e-1) - Y(e-1) );
            //Y(e)    = B(p) * Ytil(e);
        } 
        Y_prev = Y(blk_size-1); 
        A_prev = A(blk_size-1); 
    }

    MatrixDType A_next, X_next;

    for (int b = num_blocks-1; b > -1; --b) /*loop over blocks*/
    {
        auto const& A    =    A_blkvec[b].const_table();
        //auto const& Xtil = Xtil_blkvec[b].table();
        auto const& X    =    X_blkvec[b].table();
        int blk_size = A_blkvec[b].hi()[0];


	if(b < num_blocks - 1)
	{
            int n = (blk_size-1) + cumulative_blk_size[b]; /*id if array was not partitioned in blocks*/
	    int p = n%2;

            Xtil_gl(n)    = B(p) / ( A_next - X_next );
            X(blk_size-1) = C(p) * Xtil_gl(n);

            //Xtil(blk_size-1) = B(p) / ( A_next - X_next );
            //X(blk_size-1)    = C(p) * Xtil(blk_size-1);
	} 
	else 
	{
            int n = (blk_size-1) + cumulative_blk_size[b]; 
	    Xtil_gl(n) = 0; 
	    //Xtil(blk_size-1) = 0; 
	    
	    X(blk_size-1)    = 0; 
	}

        for (int e = blk_size-2; e > -1; --e) /*loop over elements of a block*/
	{
            int n = e + cumulative_blk_size[b]; /*id if array was not partitioned in blocks*/
	    int p = n%2;

            Xtil_gl(n) = B(p) / ( A(e+1) - X(e+1) );
            X(e)    = C(p) * Xtil_gl(n);

            //Xtil(e) = B(p) / ( A(e+1) - X(e+1) );
            //X(e)    = C(p) * Xtil(e);
        } 
        X_next = X(0); 
        A_next = A(0); 
    }

    //amrex::Print() << "\nYtil & Y: \n";
    //for (int b = 0; b < num_blocks; ++b) /*loop over blocks*/
    //{
    //    auto const& A    =    A_blkvec[b].const_table();
    //    //auto const& Ytil = Ytil_blkvec[b].table();
    //    auto const& Y    =    Y_blkvec[b].table();
    //    int blk_size = A_blkvec[b].hi()[0];

    //    for (int e = 0; e < blk_size; ++e) /*loop over elements of a block*/
    //    {
    //        int n = e + cumulative_blk_size[b]; 
    //        amrex::Print() << std::setw(25) << Ytil_gl(n) << std::setw(25) << Y(e)<< "\n";
    //    } 
    //}

    //amrex::Print() << "\nXtil & X: \n";
    //for (int b = num_blocks-1; b > -1; --b) /*loop over blocks*/
    //{
    //    auto const& A    =    A_blkvec[b].const_table();
    //    //auto const& Xtil = Xtil_blkvec[b].table();
    //    auto const& X    =    X_blkvec[b].table();
    //    int blk_size = A_blkvec[b].hi()[0];

    //    for (int e = blk_size-1; e > -1; --e) /*loop over elements of a block*/
    //    {
    //        int n = e + cumulative_blk_size[b]; 
    //        amrex::Print() << std::setw(25) << Xtil_gl(n) << std::setw(25) << X(e)<< "\n";
    //    } 
    //}

    for (int b = 0; b < num_blocks; ++b) /*loop over blocks*/
    {
        auto const& A    =    A_blkvec[b].const_table();
        //auto const& Ytil = Ytil_blkvec[b].const_table();
        auto const& Y    =    Y_blkvec[b].const_table();
        //auto const& Xtil = Xtil_blkvec[b].const_table();
        auto const& X    =    X_blkvec[b].const_table();
        int num_columns = A_blkvec[b].hi()[0];

        auto const& G = G_blkvec[b].table();

	int cumulative_columns = cumulative_blk_size[b];

        amrex::ParallelFor(num_columns, [=] AMREX_GPU_DEVICE (int n) noexcept
        {
	    int n_gl = n + cumulative_columns; /*global column number*/

            G(n_gl,n) =  1./(A(n) - X(n) - Y(n)); 

            for (int m = n_gl; m > 0; m--)
            {   
                G(m-1,n) =  -Ytil_gl(m)*G(m,n);
            }
            for (int m = n_gl; m < N-1; ++m)
            {   
                G(m+1,n) = -Xtil_gl(m)*G(m,n);
            }
        });
    }
    amrex::Gpu::streamSynchronize();

}

int main (int argc, char* argv[])
{

    amrex::Initialize(argc,argv);

    //Construct Tridiagonal Dummy Hamiltonian
    std::array<amrex::Real,3> point_charge_loc {0., 0., 1e-9};

    const int N_total = 60000;
    const int N = 20000;
    const int num_blocks = ceil(static_cast<amrex::Real>(N_total)/N);
    amrex::Print() << "number of blocks: " << num_blocks << "\n";

    amrex::Vector<int> diag_block_size(num_blocks);

    /*following alg for determining diagonal block size can be different*/
    /*block sizes can be different for some intricate material structures*/
    for(int p=0; p < num_blocks-1; ++p)
    {
        diag_block_size[p] = N;
    }
    diag_block_size[num_blocks-1] = N_total - (num_blocks-1)*N;


    /*setting A, the diagonal 1D vector, which is partitioned in the form of some blocks*/
    /*in general A can be 2D*/
    amrex::Vector< Matrix1D > A_blkvec(num_blocks);
    for(int p=0; p < num_blocks; ++p)
    {
        A_blkvec[p].resize({0},{diag_block_size[p]}); 
    }

    /*R this is the size of repeated B and C vectors in Hamiltonian in the mode-space approximation for CNT*/
    const int R = 2; 

    /*In general, B and C are 2D blocks of different sizes 2D blocks*/
    Matrix1D B_data({0},{R});
    Matrix1D C_data({0},{R});
    auto const& B = B_data.table();
    auto const& C = C_data.table();


    /*specifying A block vector for CNT with point charge potential*/
    int counter = 0;
    for(auto& A_blk: A_blkvec) 
    {
        auto const& A = A_blk.table();
        auto num_elem = A_blk.hi();
        for (std::size_t i = 0; i < num_elem[0]; ++i)
        {
            amrex::Real layer_loc_y = -10e-9 
		    + (static_cast<amrex::Real>(counter)/(N_total-1))*20e-9;        
            
            amrex::Real r = pow( (  pow((0.          - point_charge_loc[0]),2)
                                  + pow((layer_loc_y - point_charge_loc[1]),2)
                                  + pow((0.          - point_charge_loc[2]),2) ), 0.5);

            A(i) = -(PhysConst::q_e)/(4. * MathConst::pi * PhysConst::ep0 * r);

	    counter += 1;
	    //amrex::Print() <<std::setw(5) << counter << std::setw(15) << std::setprecision(2) << A(i) << "\n";
        }
    }	
 
    /*below are some parameters for (17,0) carbon nanotube*/
    amrex::Real gamma = 2.5; //eV
    int M=17;
    MatrixDType beta = get_beta(gamma,M,6);

    /*specifying B and C vectors for the special case of mode-space Hamiltonian of CNT*/
    for (std::size_t i = 0; i < R; ++i)
    {
       if(i%2 == 0) {
          B(i) = conjugate(beta);
          C(i) = beta;
       }
       else {
          B(i) = gamma;
          C(i) = gamma;
       } 
    }
  
    /*allocating G in the vector of blocks format,size: N-1 x size of 1D diagonal block*/
    amrex::Vector< Matrix2D > G_blkvec(num_blocks);
    for(int p=0; p < num_blocks; ++p)
    {
        G_blkvec[p].resize({0,0},{N_total,diag_block_size[p]}); 
    }

    amrex::Real mat_inv_beg_time = amrex::second();

    MatInv_BlockTriDiagonal<N_total>(A_blkvec, B_data, C_data, G_blkvec);   

    amrex::Real mat_inv_time = amrex::second() - mat_inv_beg_time;
    amrex::Print() << "Matrix inversion time: " << std::setw(15) << mat_inv_time << "\n";

    //amrex::Print() << "G:\n";
    //PrintTable(G_blkvec);
 
    amrex::Finalize();

}
