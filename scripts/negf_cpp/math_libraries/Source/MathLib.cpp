#include "MathLib.H"
#include "cudaErrorCheck.H"

#ifdef AMREX_USE_CUDA
    // CUDA runtime
    #include <cuda_runtime.h>
    
    // CUSOLVER
    #include "cusolverSp.h"
    #include "cusolverSp_LOWLEVEL_PREVIEW.h"
    //#include "helper_cuda.h"
    //#include "helper_cusolver.h"
    
    // CUBLAS
    #include <cublas_v2.h>

    //#include <helper_functions.h>
    //#include <helper_cuda.h>
#endif

void MathLib::MatrixMatrixMultiply(ComplexType* d_C,
                                   const ComplexType* d_A, 
                                   const ComplexType* d_B,
                                   unsigned int wA, 
                                   unsigned int hA, 
                                   unsigned int wB) 
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_CUDA
    cublasHandle_t handle;
    const cuDoubleComplex alpha = make_cuDoubleComplex(1.0, 0.0);
    const cuDoubleComplex beta = make_cuDoubleComplex(0.0, 0.0);
    checkCudaErrors(cublasCreate(&handle));

    // Perform matrix multiplication: C = alpha * A * B + beta * C
    checkCudaErrors(cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, wB, hA, wA, &alpha, 
                    reinterpret_cast<const cuDoubleComplex*>(d_B), wB, 
                    reinterpret_cast<const cuDoubleComplex*>(d_A), wA, &beta, 
                    reinterpret_cast<cuDoubleComplex*>(d_C), wB));

    checkCudaErrors(cublasDestroy(handle));
#elif AMREX_USE_HIP
    /*HIP version*/
#elif AMREX_USE_SYCL
    /*SYCL version*/
#endif   
#else
    /*LAPACK version*/
#endif
}
