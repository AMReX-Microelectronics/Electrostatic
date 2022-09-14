#ifndef INPUT_H_
#define INPUT_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
using namespace amrex;

void ReadBasicInput(GpuArray<int, AMREX_SPACEDIM>& n_cell, 
                    int& max_grid_size, 
                    GpuArray<Real, AMREX_SPACEDIM>& prob_lo, 
                    GpuArray<Real, AMREX_SPACEDIM>& prob_hi, 
                    int& nstep, 
                    int& plot_int, 
                    Real& dt);

void ReadMaterialInput(Real& eps_0,
                       Real& eps_r);

void ReadMLMGInput(Real& alpha, Real& beta, int& set_verbose_param);

#endif
