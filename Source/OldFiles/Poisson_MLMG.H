#ifndef POISSON_MLMG_H_
#define POISSON_MLMG_H_
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include "DataAveragingRoutines.H"

using namespace amrex;

void setupMLABecForPoissonProblem(MLABecLaplacian& mlabec,
                                  Real mlmg_ascalar,
                                  Real mlmg_bscalar,
                                  int mlmg_max_order,
                                  const BoxArray& ba,
                                  const DistributionMapping& dm,
                                  MultiFab& PoissonPhi,
                                  MultiFab& eps_cc,
                                  MultiFab& rho);

void GetEFieldUsingMLMG(MLMG& mlmg,
                        Array< MultiFab, AMREX_SPACEDIM >& E,
                        const BoxArray& ba,
                        const DistributionMapping& dm);


void GetFluxesUsingMLMG(MLMG& mlmg,
                        amrex::Vector<Array< MultiFab, AMREX_SPACEDIM >>& minus_epsGradPhi,
                        const BoxArray& ba,
                        const DistributionMapping& dm);
#endif
