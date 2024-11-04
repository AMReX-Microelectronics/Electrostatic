#include "MatrixDef.H"

class Test_DM_Eig
{
   private:
    int A_rows, A_cols;
    bool output_copied_to_host=false;

    Matrix2D h_A_data, h_U_data, h_ANS_data;
    Matrix2D d_A_data, d_U_data;
    Matrix1D h_Lambda_data, d_Lambda_data;

    ComplexType A_init_val= ComplexType(1., 2.);

    void Generate_Answer();
    bool Check_Answer();
    void Copy_Soln_To_Host();

   public:
    Test_DM_Eig(int a_rows, int a_cols);

    void Define();
    void Initialize();
    void Perform_Test();
    void Print_Input();
    void Print_Output();
    void Verify(bool flag_print_answer=false);
    void Clear();

};
