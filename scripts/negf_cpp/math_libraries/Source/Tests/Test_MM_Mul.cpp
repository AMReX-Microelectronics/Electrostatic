#include "Test_MM_Mul.H"
#include "GlobalFuncs.H"
#include "MathLib.H"
#include <AMReX_Gpu.H>

using namespace amrex;

Test_MM_Mul::Test_MM_Mul(int a_rows, int a_cols, int b_cols):
        A_rows(a_rows), A_cols(a_cols), B_cols(b_cols) 
{
    amrex::Print() << "\n##### Testing Matrix Matrix Multiplication #####\n";
}


void
Test_MM_Mul:: Define() 
{
    h_A_data.resize({0,0},{A_rows-1, A_cols-1},The_Pinned_Arena());
    d_A_data.resize({0,0},{A_rows-1, A_cols-1},The_Arena());
    h_B_data.resize({0,0},{A_cols-1, B_cols-1},The_Pinned_Arena());
    d_B_data.resize({0,0},{A_cols-1, B_cols-1},The_Arena());
    h_C_data.resize({0,0},{A_rows-1, B_cols-1},The_Pinned_Arena());
    d_C_data.resize({0,0},{A_rows-1, B_cols-1},The_Arena());
    h_ANS_data.resize({0,0},{A_rows-1, B_cols-1},The_Pinned_Arena());
}


void
Test_MM_Mul:: Initialize() 
{
    ComplexType zero(0., 0.);

    Define_Table2D(h_A_data, A_init_val);
    Define_Table2D(h_B_data, B_init_val);
    SetVal_Table2D(h_C_data, zero);
    SetVal_Table2D(h_ANS_data, zero);

    d_A_data.copy(h_A_data);
    d_B_data.copy(h_B_data);
    d_C_data.copy(h_C_data);
}


void
Test_MM_Mul:: Print_Input() 
{
    Print_Table2D(h_A_data, "A");
    const auto& h_A = h_A_data.const_table();
    const auto& h_B = h_B_data.const_table();

    Print_Table2D(h_B_data, "B");
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
    if(!output_copied_to_host) Copy_Soln_To_Host();

    Print_Table2D(h_C_data, "C = A * B");
}


void
Test_MM_Mul:: Generate_Answer() 
{
    const auto& h_A = h_A_data.const_table();
    const auto& h_B = h_B_data.const_table();
    const auto& h_ANS = h_ANS_data.table();

    for (int i = 0; i < A_rows; ++i)
    {
        for (int j = 0; j < B_cols; ++j) //slow access
        {
            ComplexType sum(0.,0.);
            for (int k = 0; k < A_cols; ++k) //slow access
            {
                sum += h_A(i,k) * h_B(k,j);
            }
            h_ANS(i,j) = sum;
        }
    }
}

void
Test_MM_Mul:: Print_Answer()
{
    Print_Table2D(h_ANS_data, "ANS = A * B");
}


void
Test_MM_Mul:: Copy_Soln_To_Host() 
{
    //copy C from device to host and print
    h_C_data.copy(d_C_data);
    Gpu::streamSynchronize();
    output_copied_to_host = true;
}


bool
Test_MM_Mul:: Check_Answer() 
{
    const auto& h_C = h_C_data.const_table();
    const auto& h_ANS = h_ANS_data.const_table();

    bool test_passed = true;
    for (int i = 0; i < A_rows; ++i)
    {
        for (int j = 0; j < B_cols; ++j) //slow access
        {
            if((fabs(h_ANS(i,j).real() - h_C(i,j).real()) > 1e-8) or
               (fabs(h_ANS(i,j).imag() - h_C(i,j).imag()) > 1e-8)) {
                return false;
            }
        }
    }
    return test_passed;
}


void
Test_MM_Mul:: Verify(bool flag_print_answer)
{
    if(!output_copied_to_host) Copy_Soln_To_Host();

    Generate_Answer();

    if(flag_print_answer) Print_Answer();

    bool test_passed = Check_Answer();
    if (test_passed)
        amrex::Print() << "\nMatrix Matrix Multiplication Test Passed!\n"; 
    else
        amrex::Print() << "\nMatrix Matrix Multiplication Test Failed!\n"; 
}

void
Test_MM_Mul::Clear()
{
    h_A_data.clear();
    h_B_data.clear();
    h_C_data.clear();
    h_ANS_data.clear();

    d_A_data.clear();
    d_B_data.clear();
    d_C_data.clear();
}
