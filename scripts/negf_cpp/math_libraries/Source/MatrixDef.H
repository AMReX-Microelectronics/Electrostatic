#ifndef MATRIX_DEF_H_
#define MATRIX_DEF_H_

#include <AMReX_GpuComplex.H>
#include<AMReX_TableData.H>
#include <AMReX_REAL.H>

using ComplexType = amrex::GpuComplex<amrex::Real>;
using Matrix2D = amrex::TableData<ComplexType, 2>;

#endif
