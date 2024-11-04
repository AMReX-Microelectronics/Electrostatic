#include "Test_ConjT.H"

#include <AMReX_Gpu.H>

#include "GlobalFuncs.H"
#include "MathLib.H"

using namespace amrex;

Test_ConjT::Test_ConjT(int a_rows, int a_cols) : A_rows(a_rows), A_cols(a_cols)
{
}

void Test_ConjT::Define()
{
    h_A_data.resize({0, 0}, {A_rows - 1, A_cols - 1}, The_Pinned_Arena());
    d_A_data.resize({0, 0}, {A_rows - 1, A_cols - 1}, The_Arena());
    h_AconjT_data.resize({0, 0}, {A_cols - 1, A_rows - 1}, The_Pinned_Arena());
    d_AconjT_data.resize({0, 0}, {A_cols - 1, A_rows - 1}, The_Arena());
    h_ANS_data.resize({0, 0}, {A_cols - 1, A_rows - 1}, The_Pinned_Arena());
}

void Test_ConjT::Initialize()
{
    ComplexType zero(0., 0.);

    Define_Table2D(h_A_data, A_init_val);
    SetVal_Table2D(h_AconjT_data, zero);
    SetVal_Table2D(h_ANS_data, zero);

    d_A_data.copy(h_A_data);
    d_AconjT_data.copy(h_AconjT_data);
}

void Test_ConjT::Print_Input() { Print_Table2D(h_A_data, "A"); }

void Test_ConjT::Perform_Test()
{
    const auto& d_A = d_A_data.const_table();
    const auto& d_AconjT = d_AconjT_data.table();

    MathLib::ConjugateTranspose(d_AconjT.p, d_A.p, A_rows, A_cols);
}

void Test_ConjT::Print_Output()
{
    if (!output_copied_to_host) Copy_Soln_To_Host();

    Print_Table2D(h_AconjT_data, "AconjT");
}

void Test_ConjT::Generate_Answer()
{
    const auto& h_A = h_A_data.const_table();
    const auto& h_AconjT = h_AconjT_data.const_table();
    const auto& h_ANS = h_ANS_data.table();

    // We store ANS = A*AconjT
    for (int i = 0; i < A_rows; ++i)
    {
        for (int j = 0; j < A_cols; ++j)  // slow access
        {
            ComplexType temp(h_A(i, j).real(), -1. * h_A(i, j).imag());
            h_ANS(j, i) = temp;
        }
    }
}

void Test_ConjT::Print_Answer() { Print_Table2D(h_ANS_data, "ANS = AconjT"); }

void Test_ConjT::Copy_Soln_To_Host()
{
    // copy AconjT from device to host and print
    h_AconjT_data.copy(d_AconjT_data);
    Gpu::streamSynchronize();
    output_copied_to_host = true;
}

bool Test_ConjT::Check_Answer()
{
    const auto& h_AconjT = h_AconjT_data.const_table();
    const auto& h_ANS = h_ANS_data.const_table();

    bool test_passed = true;
    // ANS should be equal to the conjugate transpose.
    for (int i = 0; i < A_rows; ++i)
    {
        for (int j = 0; j < A_cols; ++j)  // slow access
        {
            if ((fabs(h_ANS(i, j).real() - h_AconjT(i, j).real()) > 1e-8) or
                (fabs(h_ANS(i, j).imag() - h_AconjT(i, j).imag()) > 1e-8))
            {
                return false;
            }
        }
    }
    return test_passed;
}

void Test_ConjT::Verify(bool flag_print_answer)
{
    if (!output_copied_to_host) Copy_Soln_To_Host();

    Generate_Answer();

    if (flag_print_answer) Print_Answer();

    bool test_passed = Check_Answer();
    if (test_passed)
        amrex::Print() << "\nConjugate Transpose Test Passed!\n";
    else
        amrex::Print() << "\nConjugate Transpose Test Failed!\n";
}

void Test_ConjT::Clear()
{
    h_A_data.clear();
    h_AconjT_data.clear();
    h_ANS_data.clear();

    d_A_data.clear();
    d_AconjT_data.clear();
}
