#ifndef DATA_AVERAGING_ROUTINES_H_
#define DATA_AVERAGING_ROUTINES_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void AverageCellCenteredMultiFabToCellFaces(const MultiFab& cc_arr,
                   std::array< MultiFab, AMREX_SPACEDIM >& face_arr);

//void ReverseSignOfValuesOfMultiFab(MultiFab& mFab);

void ComputeSumOfFluxes(amrex::Vector<Array< MultiFab, AMREX_SPACEDIM >>& minus_epsGradPhi);

#endif
