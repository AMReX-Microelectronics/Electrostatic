#ifndef GEOMETRY_PROPERTIES_H_
#define GEOMETRY_PROPERTIES_H_

#include "GeometryProperties_fwd.H"

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>

//#include <AMReX_BaseFwd.H>



#include <string>

class 
c_GeometryProperties
{
public: 
    c_GeometryProperties();

    amrex::GpuArray<int, AMREX_SPACEDIM>  n_cell;     // number of cells in each dimension
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>  prob_lo; // physical lo coordinate
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>  prob_hi; // physical hi coordinate
    amrex::IntVect max_grid_size;
    amrex::IntVect blocking_factor;
    amrex::BoxArray ba; // a list of boxes that cover the domain
    amrex::Geometry geom;
    amrex::DistributionMapping dm;
//    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx;

    void ReadData();
    void InitData();

private:
    void ParseBasicDomainInput();
    void InitializeBoxArrayAndDistributionMap ();

};

#endif
