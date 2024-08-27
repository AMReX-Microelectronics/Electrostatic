#include "MatrixDef.H"

class Test_MM_Mul
{
   private:
    int A_rows, A_cols, B_cols;
    std::vector<int> print_flags;
    
    Matrix2D h_A_data, h_B_data, h_C_data;
    Matrix2D d_A_data, d_B_data, d_C_data;

    ComplexType A_init_val= ComplexType(1., 2.);
    ComplexType B_init_val= ComplexType(-5., 3.);

   public:
    Test_MM_Mul(int a_rows, int a_cols, int b_cols, 
                const std::vector<int>& flags);

    void Define();
    void Initialize();
    void Perform_Test();
    void Print_Input();
    void Print_Output();

};
