#ifndef MATRIX_BLOCK_H_
#define MATRIX_BLOCK_H_

#include "../../../Code_Definitions.H"

#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Print.H>
#include <AMReX_GpuComplex.H>

#include <cmath>
#include <math.h>
#include<stdlib.h>

using namespace amrex;
using ComplexType = amrex::GpuComplex<amrex::Real>;


template<typename T>
struct MatrixBlock
{
   public:
   T block;

   //AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
   //template<typename U>
   //friend MatrixBlock<U> real(const MatrixBlock<T>& rhs);

   MatrixBlock<T> operator=(const amrex::Real c);
   MatrixBlock<T> operator=(const ComplexType c);

   MatrixBlock<T> operator*(const amrex::Real c);
   MatrixBlock<T> operator*(const MatrixBlock<T>& rhs);
   template<typename U>
   friend MatrixBlock<U> operator*(const amrex::Real c, const MatrixBlock<U>& rhs);

   MatrixBlock<T> operator+(const ComplexType c);
   MatrixBlock<T> operator+(const MatrixBlock<T>& rhs);
   template<typename U>
   friend MatrixBlock<U> operator+(const ComplexType c, const MatrixBlock<U>& rhs);

   MatrixBlock<T> operator-(const amrex::Real c);
   MatrixBlock<T> operator-(const ComplexType c);
   MatrixBlock<T> operator-(const MatrixBlock<T>& rhs);
   template<typename U>
   friend MatrixBlock<U> operator-(const ComplexType c, const MatrixBlock<U>& rhs);

   MatrixBlock<T> operator/(const ComplexType c);
   MatrixBlock<T> operator/(const MatrixBlock<T>& rhs);

   template<typename U>
   friend std::ostream& operator<<(std::ostream& stream, const MatrixBlock<U>& rhs);
};

#endif
