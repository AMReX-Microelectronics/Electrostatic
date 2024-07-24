#ifndef CLOUD_IN_CELL_H_
#define CLOUD_IN_CELL_H_

#include <AMReX_Real.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuContainers.H>

namespace CloudInCell
{
    AMREX_GPU_DEVICE amrex::Real 
        Gather(const amrex::Real& pos, const amrex::Real& plo, 
               const amrex::Real& dx,  const amrex::MultiFab& phi);
}
#endif
