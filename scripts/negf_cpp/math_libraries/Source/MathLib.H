#include "MatrixDef.H"

namespace MathLib 
{
    void MatrixMatrixMultiply(ComplexType* d_C,
                              const ComplexType* d_A,  
                              const ComplexType* d_B,
                              unsigned int A_rows,
                              unsigned int A_cols,
                              unsigned int B_cols);

    void DenseMatrixInversion(ComplexType* d_Anv,
                              ComplexType* d_A,  
                              unsigned int A_rows,
                              unsigned int A_cols);

    void ConjugateTranspose(ComplexType* d_AconjT,
                            const ComplexType* d_A,  
                            unsigned int A_rows,
                            unsigned int A_cols);

    void DenseMatrixEigendecomposition(ComplexType* d_U,
                                       ComplexType* d_Lambda,
                                       unsigned int A_rows,
                                       unsigned int A_cols);


}
