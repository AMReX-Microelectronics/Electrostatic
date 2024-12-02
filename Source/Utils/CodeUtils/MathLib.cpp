#include <AMReX_Print.H>

#ifdef AMREX_USE_CUDA
#include <cublas_v2.h>     //CUBLAS
#include <cuda_runtime.h>  //CUDA runtime
#include <cusolverDn.h>    //CUSOLVER
#endif

#include "MathLib.H"
#include "cudaErrorCheck.H"

void MathLib::MatrixMatrixMultiply(ComplexType* d_C, const ComplexType* d_A,
                                   const ComplexType* d_B, unsigned int A_rows,
                                   unsigned int A_cols, unsigned int B_cols)
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_CUDA
    cublasHandle_t handle;
    const cuDoubleComplex alpha = make_cuDoubleComplex(1.0, 0.0);
    const cuDoubleComplex beta = make_cuDoubleComplex(0.0, 0.0);
    checkCudaErrors(cublasCreate(&handle));

    // Perform matrix multiplication: C = alpha * A * B + beta * C
    // see: https://docs.nvidia.com/cuda/archive/10.0/cublas/index.html
    cublasStatus_t status =
        cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, A_rows, B_cols, A_cols,
                    &alpha, reinterpret_cast<const cuDoubleComplex*>(d_A),
                    A_rows, reinterpret_cast<const cuDoubleComplex*>(d_B),
                    A_cols, &beta, reinterpret_cast<cuDoubleComplex*>(d_C),
                    A_rows);
    checkCudaErrors(status);

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

void MathLib::DenseMatrixInversion(ComplexType* d_Ainv, ComplexType* d_A,
                                   unsigned int A_rows, unsigned int A_cols)
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_CUDA
    cublasHandle_t cublasH = nullptr;
    cusolverDnHandle_t cusolverH = nullptr;

    checkCudaErrors(cusolverDnCreate(&cusolverH));
    checkCudaErrors(cublasCreate(&cublasH));

    int* d_info = nullptr;
    checkCudaErrors(cudaMalloc(&d_info, sizeof(int)));

    int* d_pivot = nullptr;
    checkCudaErrors(cudaMalloc(&d_pivot, A_rows * sizeof(int)));

    int lwork = 0;
    cuDoubleComplex* d_work = nullptr;

    // Compute A_inverse from A using LU decomposition.
    checkCudaErrors(
        cusolverDnZgetrf_bufferSize(cusolverH, A_rows, A_cols,
                                    reinterpret_cast<cuDoubleComplex*>(d_A),
                                    A_rows, &lwork));

    checkCudaErrors(cudaMalloc(&d_work, lwork * sizeof(cuDoubleComplex)));

    // Perform LU decomposition using cuSolver.
    // P A = L U (P is permutation matrix)
    //  see: https://docs.nvidia.com/cuda/cusolver/index.html
    checkCudaErrors(cusolverDnZgetrf(cusolverH, A_rows, A_cols,
                                     reinterpret_cast<cuDoubleComplex*>(d_A),
                                     A_rows, d_work, d_pivot, d_info));

    // Check if LU decomposition was successful
    int info = 0;
    checkCudaErrors(
        cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (info != 0)
    {
        std::cerr << "LU decomposition failed with info = " << info << "\n";
        cudaDeviceReset();
        exit(1);
    }

    // Solve for Ainv using eq: A * A_inv = I, where A is factorized.
    // op(A) * X = B
    checkCudaErrors(
        cusolverDnZgetrs(cusolverH, CUBLAS_OP_N, A_rows, A_rows,
                         reinterpret_cast<const cuDoubleComplex*>(d_A), A_rows,
                         d_pivot, reinterpret_cast<cuDoubleComplex*>(d_Ainv),
                         A_rows, d_info));

    checkCudaErrors(cudaFree(d_info));
    checkCudaErrors(cudaFree(d_pivot));
    checkCudaErrors(cudaFree(d_work));

    checkCudaErrors(cusolverDnDestroy(cusolverH));
    checkCudaErrors(cublasDestroy(cublasH));

#elif AMREX_USE_HIP
    /*HIP version*/
#elif AMREX_USE_SYCL
    /*SYCL version*/
#endif
#else
    /*LAPACK version*/
#endif
}

void MathLib::ConjugateTranspose(ComplexType* d_AconjT, const ComplexType* d_A,
                                 unsigned int A_rows, unsigned int A_cols)
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_CUDA
    cublasHandle_t handle;
    const cuDoubleComplex alpha = make_cuDoubleComplex(1.0, 0.0);
    const cuDoubleComplex beta = make_cuDoubleComplex(0.0, 0.0);
    checkCudaErrors(cublasCreate(&handle));

    // Perform matrix conjugate transpose: d_AconjT = conj(d_A)
    // https://docs.nvidia.com/cuda/cublas/
    cublasStatus_t status =
        cublasZgeam(handle, CUBLAS_OP_C, CUBLAS_OP_N, A_cols, A_rows, &alpha,
                    reinterpret_cast<const cuDoubleComplex*>(d_A), A_rows,
                    &beta, reinterpret_cast<const cuDoubleComplex*>(d_A),
                    A_rows, reinterpret_cast<cuDoubleComplex*>(d_AconjT),
                    A_cols);
    checkCudaErrors(status);

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

void MathLib::DenseMatrixEigendecomposition(ComplexType* d_U,
                                            ComplexType* d_Lambda,
                                            unsigned int A_rows,
                                            unsigned int A_cols)
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_CUDA
    cusolverDnHandle_t cusolverH;
    checkCudaErrors(cusolverDnCreate(&cusolverH));

    // Copy the matrix to be decomposed, say A, into d_U before calling this.
    // so d_U will hold the eigenvectors upon completion

    // Workspace and info setup for cuSolver
    int work_size = 0;
    int* devInfo;
    checkCudaErrors(cudaMalloc((void**)&devInfo, sizeof(int)));

    //// Query workspace size for eigendecomposition of a complex Hermitian
    /// matrix
    // checkCudaErrors(cusolverDnZheevd_bufferSize(
    //     cusolverH, CUBLAS_FILL_MODE_UPPER, A_rows,
    //     reinterpret_cast<cuDoubleComplex*>(d_U), A_rows,
    //     reinterpret_cast<cuDoubleComplex*>(d_Lambda), &work_size));

    //// Allocate workspace
    // ComplexType* d_work;
    // checkCudaErrors(cudaMalloc((void**)&d_work, work_size *
    // sizeof(ComplexType)));

    //// Perform the eigendecomposition
    // checkCudaErrors(cusolverDnZheevd(
    //     cusolverH, CUBLAS_FILL_MODE_UPPER, A_rows,
    //     reinterpret_cast<cuDoubleComplex*>(d_U), A_rows,
    //     reinterpret_cast<cuDoubleComplex*>(d_Lambda),
    //     reinterpret_cast<cuDoubleComplex*>(d_work), work_size, devInfo));

    // Clean up
    // checkCudaErrors(cudaFree(d_work));
    // checkCudaErrors(cudaFree(devInfo));
    checkCudaErrors(cusolverDnDestroy(cusolverH));
#elif AMREX_USE_HIP
    /*HIP version*/
#elif AMREX_USE_SYCL
    /*SYCL version*/
#endif
#else
    /*LAPACK version*/
#endif
}
