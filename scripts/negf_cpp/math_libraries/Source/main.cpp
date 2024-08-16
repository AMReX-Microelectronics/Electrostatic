#include "MathLib.H"
#include "MatrixDef.H"

#include <AMReX_ParmParse.H>
using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    int num_proc = ParallelDescriptor::NProcs(); 
    int my_rank = ParallelDescriptor::MyProc();
    amrex::Print() << "total number of procs: " << num_proc << "\n";

    int N_total = 4; /*matrix size*/
    amrex::ParmParse pp;
    pp.query("N_total", N_total);

    int print_matrix_flag = false;
    pp.query("print_matrix", print_matrix_flag);
    
    ComplexType num1(1.,2.);
    ComplexType num2(-5.,3.);

    Matrix2D h_A_data({0,0},{N_total,2*N_total},The_Pinned_Arena());
    SetVal_Table2D(h_A_data, num1);

    Matrix2D d_A_data({0,0},{N_total,2*N_total},The_Arena());
    h_A_data.copy(d_A_data);

    Matrix2D h_B_data({0,0},{2*N_total,N_total},The_Pinned_Arena());
    SetVal_Table2D(h_A_data, num2);

    Matrix2D d_B_data({0,0},{2*N_total,N_total},The_Arena());
    h_A_data.copy(d_A_data);

    Matrix2D h_C_data({0,0},{N_total,N_total},The_Pinned_Arena());
    Matrix2D d_C_data({0,0},{N_total,N_total},The_Arena());

    const auto& dim_A = d_A_data.hi();
    const auto& dim_B = d_B_data.hi();
    const auto& d_A = d_A_data.const_table();
    const auto& d_B = d_B_data.const_table();
    const auto& d_C = d_C_data.table();

    MathLib::MatrixMatrixMultiply(d_C.p, d_A.p, d_B.p, dim_A[0], dim_A[1], dim_B[0]);
    d_A_data.copy(h_A_data);
    Gpu::streamSynchronize();

    amrex::Finalize();
}
