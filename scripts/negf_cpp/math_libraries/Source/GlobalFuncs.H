#ifndef GLOBAL_FUNCS_H_
#define GLOBAL_FUNCS_H_

#include <AMReX_GpuUtility.H>

#include <iomanip>
#include <string>

#include "MatrixDef.H"

template<typename U, typename V>
void
Define_Table2D (U& Tab2D_data, V val)
{
    auto tlo = Tab2D_data.lo();
    auto thi = Tab2D_data.hi();

    auto const& Tab2D = Tab2D_data.table();

    for (int i = tlo[0]; i <= thi[0]; ++i)
    {
        for (int j = tlo[1]; j <= thi[1]; ++j) //slow access
        {
            ComplexType new_val(val.real()*(j+1), val.imag()*(j+1));

            Tab2D(i,j) = new_val; 
        }
    }
}


template<typename U, typename V>
void
SetVal_Table2D (U& Tab2D_data, V val)
{
    auto tlo = Tab2D_data.lo();
    auto thi = Tab2D_data.hi();

    auto const& Tab2D = Tab2D_data.table();

    for (int i = tlo[0]; i <= thi[0]; ++i)
    {
        for (int j = tlo[1]; j <= thi[1]; ++j) //slow access
        {
            Tab2D(i,j) = val; 
        }
    }
}

template<typename U>
void Print_Table2D(const U& Tab2D_data, const std::string tablename="")
{
    auto tlo = Tab2D_data.lo();
    auto thi = Tab2D_data.hi();

    auto const& Tab2D = Tab2D_data.table();

    std::cout << "\nPrinting Table: " << tablename << "\n";
    for (int i = tlo[0]; i <= thi[0]; ++i)
    {
        for (int j = tlo[1]; j <= thi[1]; ++j) //slow access
        {
            std::cout << std::setw(12) << std::setprecision(6) << std::fixed
                      << Tab2D(i,j).real() << " + "
                      << Tab2D(i,j).imag() << "i";

            if (j < thi[1] - 1) std::cout << ", "; 
        }
        std::cout << "\n";
    }
}
#endif
