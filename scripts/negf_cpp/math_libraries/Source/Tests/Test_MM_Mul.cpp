#include "Test_MM_Mul.H"
#include "GlobalFuncs.H"
#include "MathLib.H"
#include <AMReX_Gpu.H>

using namespace amrex;

Test_MM_Mul::Test_MM_Mul(int a_rows, int a_cols, int b_cols, const std::vector<int>& flags):
        A_rows(a_rows), A_cols(a_cols), B_cols(b_cols), print_flags(flags) {}


void
Test_MM_Mul:: Define() 
{
    h_A_data.resize({0,0},{A_rows-1, A_cols-1},The_Pinned_Arena());
    d_A_data.resize({0,0},{A_rows-1, A_cols-1},The_Arena());
    h_B_data.resize({0,0},{A_cols-1, B_cols-1},The_Pinned_Arena());
    d_B_data.resize({0,0},{A_cols-1, B_cols-1},The_Arena());
    h_C_data.resize({0,0},{A_rows-1, B_cols-1},The_Pinned_Arena());
    d_C_data.resize({0,0},{A_rows-1, B_cols-1},The_Arena());
}


void
Test_MM_Mul:: Initialize() 
{
    ComplexType zero(0., 0.);

    Define_Table2D(h_A_data, A_init_val);
    Define_Table2D(h_B_data, B_init_val);
    SetVal_Table2D(h_C_data, zero);

    d_A_data.copy(h_A_data);
    d_B_data.copy(h_B_data);
    d_C_data.copy(h_C_data);
}


void
Test_MM_Mul:: Print_Input() 
{
    if(print_flags[0]) {
        Print_Table2D(h_A_data, "A");
        const auto& h_A = h_A_data.const_table();
        const auto& h_B = h_B_data.const_table();

        amrex::Print() << "\nPrinting h_A using h_A.p\n";

        //Usage Error: In the forloop we should be using (i < dim_A[0]*dim_A[1])
        //But currently we need to add a buffer of 1 unit size to print properly!
        for(int i=0; i< A_rows*A_cols; ++i)
        {
            amrex::Print() << i << " "<< *(h_A.p+i) << "\n";
        }

        Print_Table2D(h_B_data, "B");
    }
}


void
Test_MM_Mul:: Perform_Test() 
{
    const auto& d_A = d_A_data.const_table();
    const auto& d_B = d_B_data.const_table();
    const auto& d_C = d_C_data.table();

    MathLib::MatrixMatrixMultiply(d_C.p, d_A.p, d_B.p, A_rows, A_cols, B_cols);
}


void
Test_MM_Mul:: Print_Output() 
{
    if(print_flags[1]) 
    {
        //copy C from device to host and print
        h_C_data.copy(d_C_data);
        Gpu::streamSynchronize();

        Print_Table2D(h_C_data, "C");
    }
}
