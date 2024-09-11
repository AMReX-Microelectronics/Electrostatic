#include "MatrixDef.H"

namespace MathLib 
{
    /* MatrixMatrixMultiply */
    void MatrixMatrixMultiply(ComplexType* d_C,
                              const ComplexType* d_A,  
                              const ComplexType* d_B,
                              unsigned int A_rows,
                              unsigned int A_cols,
                              unsigned int B_cols);

    /* DenseMatrixInversion */
    void DenseMatrixInversion(ComplexType* d_Anv,
                              ComplexType* d_A,  
                              unsigned int A_rows,
                              unsigned int A_cols);

    /* DenseMatrixEigendecomposition */
    void DenseMatrixEigendecomposition(ComplexType* d_U,
                                       ComplexType* d_Lambda,
                                       const ComplexType* d_A,  
                                       unsigned int A_rows,
                                       unsigned int A_cols);


}
