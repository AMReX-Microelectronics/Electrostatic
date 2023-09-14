#include "Broyden_Namespace.H"

#include <AMReX.H>
#include <AMReX_GpuContainers.H>


AMREX_GPU_MANAGED amrex::Real Broyden::Broyden_Norm;
AMREX_GPU_MANAGED amrex::Real Broyden::Broyden_NormSum_Curr;
AMREX_GPU_MANAGED amrex::Real Broyden::Broyden_Denom;
AMREX_GPU_MANAGED         int Broyden::Broyden_Threshold_MaxStep;
