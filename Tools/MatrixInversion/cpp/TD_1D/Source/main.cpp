
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

amrex::GpuComplex<amrex::Real> get_beta(amrex::Real gamma, int M, int J) 
{
   amrex::GpuComplex arg(0., -MathConst::pi*J/M);
   return 2 * gamma * cos(-1*arg.imag()) * exp(arg); 
}

amrex::GpuComplex<amrex::Real> conjugate(amrex::GpuComplex<amrex::Real> a) 
{
   amrex::GpuComplex a_conj(a.real(), -1.*a.imag());
   return a_conj;
}

void PrintTable(const TableData<MatrixDType, 2>& G)
{
    auto const& h_table = G.table();
    auto tlo = G.lo();
    auto thi = G.hi();
    for (int i = tlo[0]; i <= thi[0]; ++i) { 
        for (int j = tlo[1]; j <= thi[1]; ++j) { //slow moving index. printing slow
            amrex::Print() <<std::setw(12) << std::setprecision(2) << h_table(i,j);
        }
        amrex::Print() << "\n";
    }
}

template<std::size_t N>
void MatInv_BlockTriDiagonal(TableData<MatrixDType, 1>& A_data, 
                             TableData<MatrixDType, 1>& B_data, 
                             TableData<MatrixDType, 1>& C_data,
                             TableData<MatrixDType, 2>& G_data)
{
    TableData<MatrixDType, 1> X_tilde_data({0},{N});
    TableData<MatrixDType, 1> Y_tilde_data({0},{N});
    TableData<MatrixDType, 1> X_data({0},{N});
    TableData<MatrixDType, 1> Y_data({0},{N});

    auto const& X_tilde = X_tilde_data.table();
    auto const& Y_tilde = Y_tilde_data.table();
    auto const& X = X_data.table();
    auto const& Y = Y_data.table();
    auto const& A = A_data.const_table();
    auto const& B = B_data.const_table();
    auto const& C = C_data.const_table();

    Y(0) = 0;
    for (int n = 1; n < N; ++n)
    {   
	int p = (n-1)%2;    
        Y_tilde(n) = C(p) / ( A(n-1) - Y(n-1) );
        Y(n) = B(p) * Y_tilde(n);
    }

    X(N-1) = 0;
    for (int n = N-2; n > -1; n--)
    {   
	int p = n%2;    
        X_tilde(n) = B(p)/(A(n+1) - X(n+1));
        X(n) = C(p)*X_tilde(n);
    }

    amrex::Print() << "\nYtil & Y: \n";
    for (std::size_t n = 0; n < N; ++n)
    {   
        amrex::Print() << std::setw(25)<< Y_tilde(n) << std::setw(25) << Y(n)<< "\n";
    }
    amrex::Print() << "\nXtil & X: \n";
    for (int n = N-1; n > -1; n--)
    {   
        amrex::Print() << std::setw(25)<< X_tilde(n) << std::setw(25) << X(n)<< "\n";
    }

    auto const& G = G_data.table();
    amrex::ParallelFor(N, [=] AMREX_GPU_DEVICE (int n) noexcept
    {
        G(n,n) =  1./(A(n) - X(n) - Y(n)); 

        for (int m = n; m > 0; m--)
        {   
            G(m-1,n) =  -Y_tilde(m)*G(m,n);
        }
        for (int m = n; m < N-1; ++m)
        {   
            G(m+1,n) = -X_tilde(m)*G(m,n);
        }
    });
    amrex::Gpu::streamSynchronize();

}

int main (int argc, char* argv[])
{

    amrex::Initialize(argc,argv);

    //Construct Tridiagonal Dummy Hamiltonian
    std::array<amrex::Real,3> point_charge_loc {0., 0., 1e-9};

    const int N = 5;
    const int R = 2;
    amrex::Real gamma = 2.5; //eV
    int M=17;

    TableData<MatrixDType, 1> A_data({0},{N});
    TableData<MatrixDType, 1> B_data({0},{R});
    TableData<MatrixDType, 1> C_data({0},{R});
    auto const& A = A_data.table();
    auto const& B = B_data.table();
    auto const& C = C_data.table();

    for (std::size_t i = 0; i < N; ++i)
    {
        amrex::Real layer_loc_y = -10e-9 + (static_cast<amrex::Real>(i)/(N-1))*20e-9;        
        
        amrex::Real r = pow( (  pow((0.          - point_charge_loc[0]),2)
                              + pow((layer_loc_y - point_charge_loc[1]),2)
                              + pow((0.          - point_charge_loc[2]),2) ), 0.5);

        A(i) = -(PhysConst::q_e)/(4. * MathConst::pi * PhysConst::ep0 * r);
    }

    amrex::GpuComplex<amrex::Real> beta = get_beta(gamma,M,6);

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
 
    Array<int,2> tlo{0,0};
    Array<int,2> thi{N-1,N-1};

    TableData<MatrixDType, 2> G_data(tlo,thi);

    amrex::Real mat_inv_beg_time = amrex::second();

    MatInv_BlockTriDiagonal<N>(A_data,B_data,C_data,G_data);   

    amrex::Real mat_inv_time = amrex::second() - mat_inv_beg_time;
    amrex::Print() << "Matrix inversion time: " << std::setw(15) << mat_inv_time << "\n";

    amrex::Print() << "G:\n";
    PrintTable(G_data);
 
    amrex::Finalize();

}
