#ifndef MATH_FUNCTIONS_H_
#define MATH_FUNCTIONS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuContainers.H>

namespace MathFunctions
{
    AMREX_GPU_HOST_DEVICE
    amrex::Real Sigmoid(amrex::Real x) {
        return 1.0 / (1.0 + exp(-x));
    }
}
#endif
