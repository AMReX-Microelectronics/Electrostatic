#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "Global.H"


using namespace amrex;

void InitializeBoxArrayAndGeom(BoxArray& ba, 
                               Geometry& geom, 
                               GpuArray<Real, AMREX_SPACEDIM> const& prob_lo,
                               GpuArray<Real, AMREX_SPACEDIM> const& prob_hi,
                               GpuArray<int, AMREX_SPACEDIM> const& n_cell,
                               int max_grid_size
                               );

void InitializeCharge(MultiFab& rho,
                   GpuArray<Real, AMREX_SPACEDIM> const& prob_lo,
                   GpuArray<Real, AMREX_SPACEDIM> const& prob_hi,
                   const Geometry& geom);

void InitializePermittivity(MultiFab& eps_cc,
                   GpuArray<Real, AMREX_SPACEDIM> const& prob_lo,
                   GpuArray<Real, AMREX_SPACEDIM> const& prob_hi,
                   const Geometry& geom,
                   Real eps_0,
                   Real eps_r);
#endif
