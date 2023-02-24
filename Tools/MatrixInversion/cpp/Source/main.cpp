
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
    static constexpr auto hbar  = static_cast<amrex::Real>( 1.054571817e-34 );
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
void MatInv_BlockTriDiagonal(amrex::GpuArray<amrex::GpuComplex<amrex::Real>, N> A, 
                             amrex::GpuArray<amrex::GpuComplex<amrex::Real>, N-1> B, 
                             amrex::GpuArray<amrex::GpuComplex<amrex::Real>, N-1> C,
                             TableData<MatrixDType, 2>& G)
{
    amrex::GpuArray<amrex::GpuComplex<amrex::Real> , N> X_tilde, Y_tilde;
    amrex::GpuArray<amrex::GpuComplex<amrex::Real>, N> X, Y;


    Y[0] = 0;
    for (int n = 1; n < N; ++n)
    {   
        Y_tilde[n] = C[n-1]/(A[n-1] - Y[n-1]);
        Y[n] = B[n-1]*Y_tilde[n];
    }
    //amrex::Print() << "\nY: \n";
    //for (std::size_t n = 1; n < N; ++n)
    //{   
    //    amrex::Print() << Y_tilde[n] << "\n";
    //}

    X[N-1] = 0;
    for (int n = N-2; n > -1; n--)
    {   
        X_tilde[n] = B[n]/(A[n+1] - X[n+1]);
        X[n] = C[n]*X_tilde[n];
    }
    //amrex::Print() << "\nX: \n";
    //for (int n = N-2; n > -1; n--)
    //{   
    //    amrex::Print() << X_tilde[n] << "\n";
    //}

    //auto const& h_table = G.table();
    //for (int n = 0; n < N; ++n) /*This loop can be parallelized*/
    //{
    //    h_table(n,n) = 1./(A[n] - X[n] - Y[n]);
    //}
    //for (int n = 0; n < N; ++n) /*This loop can be parallelized*/
    //{   
    //    for (int m = n; m > 0; m--)
    //    {   
    //        h_table(m-1,n) = -Y_tilde[m]*h_table(m,n);
    //    }
    //    for (int m = n; m < N-1; ++m)
    //    {   
    //        h_table(m+1,n) = -X_tilde[m]*h_table(m,n);
    //    }
    //}

    auto const& table = G.table();
    amrex::ParallelFor(N, [table=table,A=A, X=X, Y=Y,Y_tilde=Y_tilde,X_tilde=X_tilde] AMREX_GPU_DEVICE (int n) noexcept
    {
        table(n,n) =  1./(A[n] - X[n] - Y[n]); 

        for (int m = n; m > 0; m--)
        {   
            table(m-1,n) =  -Y_tilde[m]*table(m,n);
        }
        for (int m = n; m < N-1; ++m)
        {   
            table(m+1,n) = -X_tilde[m]*table(m,n);
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
    amrex::GpuArray<amrex::GpuComplex<amrex::Real>, N  > A; //diagonal
    amrex::GpuArray<amrex::GpuComplex<amrex::Real>, N-1> B; //+1 diagonal
    amrex::GpuArray<amrex::GpuComplex<amrex::Real>, N-1> C; //-1 diagonal

    //amrex::Print() << "\nA & r: \n";
    for (std::size_t i = 0; i < N; ++i)
    {
        amrex::Real layer_loc_y = -10e-9 + (static_cast<amrex::Real>(i)/(N-1))*20e-9;        
        
        amrex::Real r = pow( (  pow((0.          - point_charge_loc[0]),2)
                              + pow((layer_loc_y - point_charge_loc[1]),2)
                              + pow((0.          - point_charge_loc[2]),2) ), 0.5);

        A[i] = -(PhysConst::q_e)/(4. * MathConst::pi * PhysConst::ep0 * r);
    //    amrex::Print() << i << std::setw(25) << A[i] << std::setw(25) << layer_loc_y << "\n"; 
    }

    amrex::Real gamma = 2.5; //eV
    int M=17;
    amrex::GpuComplex<amrex::Real> beta = get_beta(gamma,M,6);

    //amrex::Print() << "\nB & C: \n";
    for (std::size_t i = 0; i < N-1; ++i)
    {
       if(i%2 == 0) {
          B[i] = conjugate(beta);
          C[i] = beta;
       }
       else {
          B[i] = gamma;
          C[i] = gamma;
       } 
    //   amrex::Print() << i << std::setw(25) << B[i] << std::setw(25) << C[i] << "\n"; 
    }
 
    Array<int,2> tlo{0,0};
    Array<int,2> thi{N-1,N-1};

    //#ifdef AMREX_USE_GPU
    //    TableData<MatrixDType,2> G_inv(tlo, thi, The_Pinned_Arena());
    //#else
    //    TableData<MatrixDType,2> G_inv(tlo, thi);
    //#endif
    //auto const& G_inv_table = G_inv.table();

    //for (int j = tlo[1]; j <= thi[1]; ++j) { //running over columns
    //for (int i = tlo[0]; i <= thi[0]; ++i) { //running over rows
    //    G_inv_table(i,j) = 0.; //i + j*1e3;
    //}}

    //for (int i = tlo[0]; i <= thi[0]; ++i) {
    //    G_inv_table(i,i)   = A[i];

    //    if(i < thi[0]) 
    //    {
    //        G_inv_table(i,i+1) = C[i];
    //        G_inv_table(i+1,i) = B[i];
    //    }  
    //}
    ////TD<decltype(h_table)> h_table_type; //const amrex::Table2D<amrex::Real>&    
    //amrex::Print() << "G_inv:\n";
    //PrintTable(G_inv);


    TableData<MatrixDType, 2> G(tlo,thi);

    amrex::Real mat_inv_beg_time = amrex::second();

    MatInv_BlockTriDiagonal<N>(A,B,C,G);   

    amrex::Real mat_inv_time = amrex::second() - mat_inv_beg_time;
    amrex::Print() << "Matrix inversion time: " << std::setw(15) << mat_inv_time << "\n";

    amrex::Print() << "G:\n";
    PrintTable(G);
 
    amrex::Finalize();

}
