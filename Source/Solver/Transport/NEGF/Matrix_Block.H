#ifndef MATRIX_BLOCK_H_
#define MATRIX_BLOCK_H_

#include <AMReX_BoxArray.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <math.h>
#include <stdlib.h>

#include <cmath>

#include "../../../Code_Definitions.H"

using namespace amrex;
using ComplexType = amrex::GpuComplex<amrex::Real>;

template <typename T>
struct MatrixBlock
{
   public:
    T block;

    // Non-const .data() for 1D array
    template <typename U = T>
    auto data() ->
        typename std::enable_if<std::rank<U>::value == 1,
                                typename std::remove_extent<U>::type *>::type
    {
        return &block[0];
    }

    // Const .data() for 1D array
    template <typename U = T>
    auto const_data() const -> typename std::enable_if<
        std::rank<U>::value == 1,
        const typename std::remove_extent<U>::type *>::type
    {
        return &block[0];
    }

    // Non-const .data() for 2D array
    template <typename U = T>
    auto data() -> typename std::enable_if<
        std::rank<U>::value == 2,
        typename std::remove_all_extents<U>::type *>::type
    {
        return reinterpret_cast<typename std::remove_all_extents<U>::type *>(
            &block[0][0]);
    }

    // Const .data() for 2D array
    template <typename U = T>
    auto const_data() const -> typename std::enable_if<
        std::rank<U>::value == 2,
        const typename std::remove_all_extents<U>::type *>::type
    {
        return reinterpret_cast<
            const typename std::remove_all_extents<U>::type *>(&block[0][0]);
    }

    AMREX_GPU_HOST_DEVICE MatrixBlock<T> operator=(const amrex::Real c);
    AMREX_GPU_HOST_DEVICE MatrixBlock<T> operator=(const ComplexType c);

    AMREX_GPU_HOST_DEVICE MatrixBlock<T> operator*(const ComplexType c);
    AMREX_GPU_HOST_DEVICE MatrixBlock<T> operator*(const MatrixBlock<T> &rhs);
    template <typename U>
    friend AMREX_GPU_HOST_DEVICE MatrixBlock<U> operator*(
        const ComplexType c, const MatrixBlock<U> &rhs);
    template <typename U>
    friend AMREX_GPU_HOST_DEVICE MatrixBlock<U> operator*(
        const amrex::Real c, const MatrixBlock<U> &rhs);

    AMREX_GPU_HOST_DEVICE MatrixBlock<T> operator+(const ComplexType c);
    AMREX_GPU_HOST_DEVICE MatrixBlock<T> operator+(const MatrixBlock<T> &rhs);
    template <typename U>
    friend AMREX_GPU_HOST_DEVICE MatrixBlock<U> operator+(
        const ComplexType c, const MatrixBlock<U> &rhs);

    AMREX_GPU_HOST_DEVICE MatrixBlock<T> operator-(const amrex::Real c);
    AMREX_GPU_HOST_DEVICE MatrixBlock<T> operator-(const ComplexType c);
    AMREX_GPU_HOST_DEVICE MatrixBlock<T> operator-(const MatrixBlock<T> &rhs);
    template <typename U>
    friend AMREX_GPU_HOST_DEVICE MatrixBlock<U> operator-(
        const ComplexType c, const MatrixBlock<U> &rhs);
    template <typename U>
    friend AMREX_GPU_HOST_DEVICE MatrixBlock<U> operator-(
        const MatrixBlock<U> &B, const MatrixBlock<U> &C);

    AMREX_GPU_HOST_DEVICE MatrixBlock<T> operator/(const ComplexType c);
    AMREX_GPU_HOST_DEVICE MatrixBlock<T> operator/(const MatrixBlock<T> &rhs);

    template <typename U>
    friend AMREX_GPU_HOST_DEVICE MatrixBlock<U> operator/(
        const ComplexType c, const MatrixBlock<U> &C);

    template <typename U>
    friend std::ostream &operator<<(std::ostream &stream,
                                    const MatrixBlock<U> &rhs);

    AMREX_GPU_HOST_DEVICE MatrixBlock<T> Conj() const;   /*Conjugate*/
    AMREX_GPU_HOST_DEVICE MatrixBlock<T> Tran() const;   /*Transpose*/
    AMREX_GPU_HOST_DEVICE MatrixBlock<T> Dagger() const; /*Conjugate-Transpose*/
    AMREX_GPU_HOST_DEVICE MatrixBlock<T> Inverse() const; /*Inverse*/
    AMREX_GPU_HOST_DEVICE ComplexType DiagSum() const;    /*Trace*/
    AMREX_GPU_HOST_DEVICE ComplexType FrobeniusNorm() const;
    AMREX_GPU_HOST_DEVICE void SetDiag(
        const ComplexType c); /*set diagonal element*/

    template <typename U>
    AMREX_GPU_HOST_DEVICE ComplexType
    DiagDotSum(U &vec) const; /*SUM of { BLOCK(n,n)*vec[n] }*/

    template <typename U>
    AMREX_GPU_HOST_DEVICE MatrixBlock<T> DiagMult(
        U &vec) const; /*{ BLOCK(n,n)*vec[n] }*/
    // MatrixBlock<T> AtomicAdd(MatrixBlock<T>& B) const; /*Atomic add for the
    // block*/ template<typename U, std::size_t N> MatrixBlock<amrex::Real>
    // Real() const;
    //
};

#endif
