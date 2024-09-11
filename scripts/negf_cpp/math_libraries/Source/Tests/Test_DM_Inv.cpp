#include "Test_DM_Inv.H"
#include "GlobalFuncs.H"
#include "MathLib.H"
#include <AMReX_Gpu.H>

using namespace amrex;

Test_DM_Inv::Test_DM_Inv(int a_rows, int a_cols):
        A_rows(a_rows), A_cols(a_cols) {}


void
Test_DM_Inv:: Define() 
{
    h_A_data.resize({0,0},{A_rows-1, A_cols-1},The_Pinned_Arena());
    d_A_data.resize({0,0},{A_rows-1, A_cols-1},The_Arena());
    h_Ainv_data.resize({0,0},{A_rows-1, A_cols-1},The_Pinned_Arena());
    d_Ainv_data.resize({0,0},{A_rows-1, A_cols-1},The_Arena());
    h_ANS_data.resize({0,0} ,{A_rows-1, A_cols-1},The_Pinned_Arena());
}


void
Test_DM_Inv:: Initialize() 
{
    ComplexType zero(0., 0.);

    Define_Table2D(h_A_data, A_init_val);
    SetVal_Table2D(h_Ainv_data, zero);
    SetVal_Table2D(h_ANS_data, zero);

    d_A_data.copy(h_A_data);
    d_Ainv_data.copy(h_Ainv_data);
}


void
Test_DM_Inv:: Print_Input() 
{
    Print_Table2D(h_A_data, "A");
}


void
Test_DM_Inv:: Perform_Test() 
{
    const auto& d_A = d_A_data.table();
    const auto& d_Ainv = d_Ainv_data.table();

    MathLib::DenseMatrixInversion(d_Ainv.p, d_A.p, A_rows, A_cols);
}


void
Test_DM_Inv:: Print_Output() 
{
    if(!output_copied_to_host) Copy_Soln_To_Host();

    amrex::Print() << "Printing output: " << "\n";
    Print_Table2D(h_Ainv_data, "Ainv");
}


void
Test_DM_Inv:: Generate_Answer() 
{
    const auto& h_A = h_A_data.const_table();
    const auto& h_Ainv = h_Ainv_data.const_table();
    const auto& h_ANS = h_ANS_data.table();

    //We store ANS = A*Ainv
    for (int i = 0; i < A_rows; ++i)
    {
        for (int j = 0; j < A_cols; ++j) //slow access
        {
            ComplexType sum(0.,0.);
            for (int k = 0; k < A_cols; ++k) //slow access
            {
                sum += h_A(i,k) * h_Ainv(k,j);
            }
            h_ANS(i,j) = sum;
        }
    }
}


void
Test_DM_Inv:: Copy_Soln_To_Host() 
{
    //copy Ainv from device to host and print
    h_Ainv_data.copy(d_Ainv_data);
    Gpu::streamSynchronize();
    output_copied_to_host = true;
}


bool
Test_DM_Inv:: Check_Answer() 
{
    const auto& h_Ainv = h_Ainv_data.const_table();
    const auto& h_ANS = h_ANS_data.const_table();

    bool test_passed = true;
    //ANS should be an identity matrix.
    for (int i = 0; i < A_rows; ++i)
    {
        for (int j = 0; j < A_cols; ++j) //slow access
        {
            if(i==j) {
                if((fabs(h_ANS(i,j).real() - 1.) > 1e-8) &&
                   (fabs(h_ANS(i,j).imag()) > 1e-8)) {
                    return false;
                }
            } else if (i != j) {
                if((fabs(h_ANS(i,j).real()) > 1e-8) &&
                   (fabs(h_ANS(i,j).imag()) > 1e-8)) {
                    return false;
                }
            }
        }
    }
    return test_passed;
}


void
Test_DM_Inv:: Verify()
{
    if(!output_copied_to_host) Copy_Soln_To_Host();

    Generate_Answer();
    bool test_passed = Check_Answer();
    if (test_passed)
        amrex::Print() << "\nDense Matrix Inversion Test Passed!\n"; 
    else
        amrex::Print() << "\nDense Matrix Inversion Test Failed!\n"; 
}
