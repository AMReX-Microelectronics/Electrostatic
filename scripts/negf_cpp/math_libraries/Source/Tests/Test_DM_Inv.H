#include "MatrixDef.H"

class Test_DM_Inv
{
   private:
    int A_rows, A_cols;
    bool output_copied_to_host=false;

    Matrix2D h_A_data, h_Ainv_data, h_ANS_data;
    Matrix2D d_A_data, d_Ainv_data;

    ComplexType A_init_val= ComplexType(1., 2.);

    void Generate_Answer();
    void Print_Answer();
    bool Check_Answer();
    void Copy_Soln_To_Host();

   public:
    Test_DM_Inv(int a_rows, int a_cols);

    void Define();
    void Initialize();
    void Perform_Test();
    void Print_Input();
    void Print_Output();
    void Verify(bool flag_print_answer=false);
    void Clear();
};
