#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "Global.H"


using namespace amrex;

void InitializeCharge(MultiFab& rho,
                   GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                   GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                   const       Geometry& geom);

void InitializePermittivity(MultiFab& eps_cc,
                   GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                   GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                   const       Geometry& geom,
                   Real eps_0,
                   Real eps_r);

void AveragePermittivityToCellFaces(MultiFab& cc_arr,
                   std::array< MultiFab, AMREX_SPACEDIM >& face_arr);
#endif
