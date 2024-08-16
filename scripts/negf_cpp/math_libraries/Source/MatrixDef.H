#ifndef MATRIX_DEF_H_
#define MATRIX_DEF_H_
#include <AMReX_GpuComplex.H>
#include <AMReX_REAL.H>
#include<AMReX_TableData.H>
#include <AMReX_GpuUtility.H>

using ComplexType = amrex::GpuComplex<amrex::Real>;
using Matrix2D = amrex::TableData<ComplexType, 2>;

template<typename U, typename V>
void
SetVal_Table2D (U& Tab2D_data, V val)
{
    auto tlo = Tab2D_data.lo();
    auto thi = Tab2D_data.hi();

    auto const& Tab2D = Tab2D_data.table();

    for (int i = tlo[0]; i < thi[0]; ++i)
    {
        for (int j = tlo[1]; j < thi[1]; ++j) //slow moving index. printing slow
        {
            Tab2D(i,j) = val;
        }
    }
}
#endif
