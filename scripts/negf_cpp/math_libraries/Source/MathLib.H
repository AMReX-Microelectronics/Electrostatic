#include "MatrixDef.H"

namespace MathLib 
{
    /* MatrixMatrixMultiply */
    void MatrixMatrixMultiply(ComplexType* d_C,
                              const ComplexType* d_A,  
                              const ComplexType* d_B,
                              unsigned int wA,
                              unsigned int hA,
                              unsigned int wB);

}
