#ifndef BROYDEN_NAMESPACE_H_
#define BROYDEN_NAMESPACE_H_

#include <AMReX_Geometry.H>

namespace Broyden
{
    extern AMREX_GPU_MANAGED amrex::Real Broyden_Norm;
    extern AMREX_GPU_MANAGED amrex::Real Broyden_NormSum_Curr;
    extern AMREX_GPU_MANAGED amrex::Real Broyden_Denom;
    extern AMREX_GPU_MANAGED int Broyden_Threshold_MaxStep;
}

#endif
