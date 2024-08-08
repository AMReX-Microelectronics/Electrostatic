#ifndef CLOUD_IN_CELL_H_
#define CLOUD_IN_CELL_H_

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuContainers.H>


namespace CloudInCell
{
    AMREX_GPU_HOST_DEVICE amrex::Real 
    Gather_Trilinear(const amrex::RealVect& pos, 
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo, 
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,  
                     const amrex::Array4<amrex::Real>& phi);

    AMREX_GPU_HOST_DEVICE void
    Deposit_Trilinear(const amrex::RealVect& pos, 
                      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo, 
                      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,  
                      const amrex::Array4<amrex::Real>& rho,
                      const amrex::Real qp=0.);
}
#endif
