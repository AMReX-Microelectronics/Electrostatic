#include <AMReX_ParmParse.H>

#include "MathLib.H"
#include "MatrixDef.H"
#include "GlobalFuncs.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    int num_proc = ParallelDescriptor::NProcs(); 
    int my_rank = ParallelDescriptor::MyProc();
    amrex::Print() << "total number of procs: " << num_proc << "\n";

    amrex::ParmParse pp;

    int A_rows=4, A_cols=4, B_cols = 4; //assume, B_rows = A_cols
    pp.query("A_rows", A_rows);
    pp.query("A_cols", A_cols);
    pp.query("B_cols", B_cols);

    int print_matrix_flag = false;
    pp.query("print_matrix", print_matrix_flag);
    
    //Create A, B, C matrices as 2D tables. We want to perform C = A * B.
    Matrix2D h_A_data({0,0},{A_rows, A_cols},The_Pinned_Arena());
    Matrix2D d_A_data({0,0},{A_rows, A_cols},The_Arena());

    Matrix2D h_B_data({0,0},{A_cols, B_cols},The_Pinned_Arena());
    Matrix2D d_B_data({0,0},{A_cols, B_cols},The_Arena());

    Matrix2D h_C_data({0,0},{A_rows, B_cols},The_Pinned_Arena());
    Matrix2D d_C_data({0,0},{A_rows, B_cols},The_Arena());

    //define matrices A & B
    ComplexType num1(1.,2.);
    ComplexType num2(-5.,3.);
    ComplexType zero(0.,0.);
    Define_Table2D(h_A_data, num1);
    Define_Table2D(h_B_data, num2);
    SetVal_Table2D(h_C_data, zero);

    //copy to A & B to device 
    d_A_data.copy(h_A_data);
    d_B_data.copy(h_B_data);
    d_C_data.copy(h_C_data);

    //get references to tables
    const auto& dim_A = d_A_data.hi();
    const auto& dim_B = d_B_data.hi();
    const auto& d_A = d_A_data.const_table();
    const auto& d_B = d_B_data.const_table();
    const auto& d_C = d_C_data.table();

    amrex::Print() << "dim_A (rows/cols): " << dim_A[0] << " " << dim_A[1] << "\n";
    amrex::Print() << "dim_B (rows/cols): " << dim_B[0] << " " << dim_B[1] << "\n";

    //print A & B
    Print_Table2D(h_A_data, "A");
    const auto& h_A = h_A_data.const_table();
    const auto& h_B = h_B_data.const_table();

    amrex::Print() << "\nPrinting h_A using h_A.p\n";

    //Usage Error: In the forloop we should be using (i < dim_A[0]*dim_A[1])
    //But currently we need to add a buffer of 1 unit size to print properly!
    for(int i=0; i<(dim_A[0]+1)*dim_A[1]; ++i) 
    { 
        amrex::Print() << i << " "<< *(h_A.p+i) << "\n";
    }

    Print_Table2D(h_B_data, "B");

    //Perform C = A * B
    MathLib::MatrixMatrixMultiply(d_C.p, d_A.p, d_B.p, dim_A[0], dim_A[1], dim_B[1]);

    //copy C from device to host and print
    h_C_data.copy(d_C_data);
    Gpu::streamSynchronize();

    Print_Table2D(h_C_data, "C");

    h_A_data.clear();
    h_B_data.clear();
    h_C_data.clear();
    d_A_data.clear();
    d_B_data.clear();
    d_C_data.clear();
    
    d_A_data.clear();
    d_B_data.clear();
    d_C_data.clear();

    amrex::Finalize();
}
