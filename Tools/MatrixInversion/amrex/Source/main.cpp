
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_GpuComplex.H>
#include<AMReX_Print.H>
#include <AMReX_TinyProfiler.H>

#include <cmath>
#include <iomanip>
#include <math.h> 
using namespace amrex;

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

template<std::size_t N>
void MatInv_BlockTriDiagonal(amrex::GpuArray<amrex::GpuComplex<amrex::Real>, N> A, 
                             amrex::GpuArray<amrex::GpuComplex<amrex::Real>, N-1> B, 
                             amrex::GpuArray<amrex::GpuComplex<amrex::Real>, N-1> C)
{
    amrex::Print() << "\nA: \n";
    for (std::size_t i = 0; i < N; ++i)
    {
        amrex::Print() << i << std::setw(25) << A[i] << "\n"; 
    }

    amrex::Print() << "\nB & C: \n";
    for (std::size_t i = 0; i < N-1; ++i)
    {
        amrex::Print() << i << std::setw(25) << B[i] << std::setw(25) << C[i] << "\n"; 
    }

    amrex::GpuArray<amrex::GpuComplex<amrex::Real> , N> X_tilde, Y_tilde;
    amrex::GpuArray<amrex::GpuComplex<amrex::Real>, N> X, Y;

    amrex::Array<amrex::Array<amrex::GpuComplex<amrex::Real>,N>,N> G_inv;
    Y[0] = 0;
    for (int n = 1; n < N; ++n)
    {   
        Y_tilde[n] = C[n-1]/(A[n-1] - Y[n-1]);
        Y[n] = B[n-1]*Y_tilde[n];
    }
    //amrex::Print() << "\nY: \n";
    //for (std::size_t n = 1; n < N; ++n)
    //{   
    //    amrex::Print() << Y[n] << "\n";
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
    //    amrex::Print() << X[n] << "\n";
    //}

    for (int n = 0; n < N; ++n) /*This loop can be parallelized*/
    {   
        G_inv[n][n] = 1./(A[n] - X[n] - Y[n]);
        amrex::Print() << G_inv[n][n] << "\n";  
    }

    for (int n = 0; n < N; ++n) /*This loop can be parallelized*/
    {   
        for (int m = n; m > 0; m--)
        {   
            G_inv[m-1][n] = -Y_tilde[m]*G_inv[m][n];
        }
        for (int m = n; m < N-1; ++m)
        {   
            G_inv[m+1][n] = -X_tilde[m]*G_inv[m][n];
        }
    }

    amrex::Print() << "G_inverse (real): " << "\n";
    for (std::size_t m = 0; m < N; ++m)
    {   
        amrex::Print() << "ROW: "<< m << "\n";
        for (std::size_t n = 0; n < N; ++n) 
        {   
            amrex::Print() <<std::setw(15) << std::setprecision(2) << G_inv[m][n].real();
        }
        amrex::Print() << "\n";
    }

    amrex::Print() << "G_inverse (imag): " << "\n";
    for (std::size_t m = 0; m < N; ++m) 
    {   
        amrex::Print() << "ROW: "<< m << "\n";
        for (std::size_t n = 0; n < N; ++n) 
        {   
            amrex::Print() <<std::setw(15) << std::setprecision(2) << G_inv[m][n].imag();
        }
        amrex::Print() << "\n";
    }
   
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

    amrex::Print() << "\nA & r: \n";
    for (std::size_t i = 0; i < N; ++i)
    {
        amrex::Real layer_loc_y = -10e-9 + (static_cast<amrex::Real>(i)/(N-1))*20e-9;        
        
        amrex::Real r = pow( (  pow((0.          - point_charge_loc[0]),2)
                              + pow((layer_loc_y - point_charge_loc[1]),2)
                              + pow((0.          - point_charge_loc[2]),2) ), 0.5);

        A[i] = -(PhysConst::q_e)/(4. * MathConst::pi * PhysConst::ep0 * r);

        amrex::Print() << i << std::setw(25) << A[i] << std::setw(25) << layer_loc_y << "\n"; 
    }

    amrex::Real gamma = 2.5; //eV
    int M=17;
    amrex::GpuComplex<amrex::Real> beta = get_beta(gamma,M,6);

    amrex::Print() << "\nB & C: \n";
    for (std::size_t i = 0; i < N-1; ++i)
    {
       if(i%2 == 0) {
          B[i] = beta;
          C[i] = beta;
       }
       else {
          B[i] = gamma;
          C[i] = gamma;
       } 
        amrex::Print() << i << std::setw(25) << B[i] << std::setw(25) << C[i] << "\n"; 
    }
 
    MatInv_BlockTriDiagonal<N>(A,B,C);   
          
    amrex::Finalize();

}
