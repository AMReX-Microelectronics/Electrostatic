#include "MatrixDef.H"

class Test_MM_Mul
{
   private:
    int A_rows, A_cols, B_cols;
    bool output_copied_to_host=false;

    Matrix2D h_A_data, h_B_data, h_C_data, h_ANS_data;
    Matrix2D d_A_data, d_B_data, d_C_data;

    ComplexType A_init_val= ComplexType(1., 2.);
    ComplexType B_init_val= ComplexType(-5., 3.);

    void Generate_Answer();
    void Print_Answer();
    bool Check_Answer();
    void Copy_Soln_To_Host();

   public:
    Test_MM_Mul(int a_rows, int a_cols, int b_cols);

    void Define();
    void Initialize();
    void Perform_Test();
    void Print_Input();
    void Print_Output();
    void Verify(bool flag_print_answer=false);
    void Clear();

};
