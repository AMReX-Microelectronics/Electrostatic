#include "Test_DM_Eig.H"

#include <AMReX_Gpu.H>

#include "GlobalFuncs.H"
#include "MathLib.H"

using namespace amrex;

Test_DM_Eig::Test_DM_Eig(int a_rows, int a_cols)
    : A_rows(a_rows), A_cols(a_cols)
{
}

void Test_DM_Eig::Define()
{
    h_A_data.resize({0, 0}, {A_rows - 1, A_cols - 1}, The_Pinned_Arena());
    d_A_data.resize({0, 0}, {A_rows - 1, A_cols - 1}, The_Arena());
    h_U_data.resize({0, 0}, {A_rows - 1, A_cols - 1}, The_Pinned_Arena());
    d_U_data.resize({0, 0}, {A_rows - 1, A_cols - 1}, The_Arena());
    h_Lambda_data.resize({0}, {A_rows - 1}, The_Pinned_Arena());
    d_Lambda_data.resize({0}, {A_rows - 1}, The_Arena());
    h_ANS_data.resize({0, 0}, {A_rows - 1, A_cols - 1}, The_Pinned_Arena());
}

void Test_DM_Eig::Initialize()
{
    ComplexType zero(0., 0.);

    Define_Table2D(h_A_data, A_init_val);
    SetVal_Table2D(h_U_data, zero);
    SetVal_Table1D(h_Lambda_data, zero);
    SetVal_Table2D(h_ANS_data, zero);

    d_A_data.copy(h_A_data);
    d_U_data.copy(h_A_data);  // copy A into U for doing eigendecomposition of
                              // A.
    // U will store eigenvectors.

    d_Lambda_data.copy(h_Lambda_data);
}

void Test_DM_Eig::Print_Input() { Print_Table2D(h_A_data, "A"); }

void Test_DM_Eig::Perform_Test()
{
    const auto& d_A = d_A_data.const_table();
    const auto& d_U = d_U_data.table();
    const auto& d_Lambda = d_Lambda_data.table();

    MathLib::DenseMatrixEigendecomposition(d_U.p, d_Lambda.p, A_rows, A_cols);
}

void Test_DM_Eig::Print_Output()
{
    if (!output_copied_to_host) Copy_Soln_To_Host();

    Print_Table2D(h_U_data, "U");
    Print_Table1D(h_Lambda_data, "Lambda");
}

void Test_DM_Eig::Generate_Answer()
{
    const auto& h_A = h_A_data.const_table();
    const auto& h_U = h_U_data.const_table();
    const auto& h_ANS = h_ANS_data.table();

    // We store ANS = U*Lambda*U^dagger
    // for (int i = 0; i < A_rows; ++i)
    //{
    //     for (int j = 0; j < A_cols; ++j) //slow access
    //     {
    //         ComplexType sum(0.,0.);
    //         for (int k = 0; k < A_cols; ++k) //slow access
    //         {
    //             sum += h_A(i,k) * h_U(k,j);
    //         }
    //         h_ANS(i,j) = sum;
    //     }
    // }
}

Test_MM_Mul::Print_Answer() { Print_Table2D(h_ANS_data, "ANS ="); }

void Test_DM_Eig::Copy_Soln_To_Host()
{
    // copy U, Lambda from device to host and print
    h_U_data.copy(d_U_data);
    h_Lambda_data.copy(d_Lambda_data);
    Gpu::streamSynchronize();
    output_copied_to_host = true;
}

bool Test_DM_Eig::Check_Answer()
{
    const auto& h_A = h_A_data.const_table();
    const auto& h_ANS = h_ANS_data.const_table();

    bool test_passed = true;
    // ANS should be equal to the original matrix.
    for (int i = 0; i < A_rows; ++i)
    {
        for (int j = 0; j < A_cols; ++j)  // slow access
        {
            if ((fabs(h_ANS(i, j).real() - h_A(i, j).real()) > 1e-8) or
                (fabs(h_ANS(i, j).imag() - h_A(i, j).imag()) > 1e-8))
            {
                return false;
            }
        }
    }
    return test_passed;
}

void Test_DM_Eig::Verify()
{
    if (!output_copied_to_host) Copy_Soln_To_Host();

    Generate_Answer();
    if (flag_print_answer) Print_Answer();

    bool test_passed = Check_Answer();
    if (test_passed)
        amrex::Print() << "\nDense Matrix Eigendecomposition Test Passed!\n";
    else
        amrex::Print() << "\nDense Matrix Eigendecomposition Test Failed!\n";
}

void Test_MM_Mul::Clear()
{
    h_A_data.clear();
    h_U_data.clear();
    h_Lambda_data.clear();
    h_ANS_data.clear();

    d_A_data.clear();
    d_U_data.clear();
    d_Lambda_data.clear();
}
