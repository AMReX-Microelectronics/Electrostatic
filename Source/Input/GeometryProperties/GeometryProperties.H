#ifndef GEOMETRY_PROPERTIES_H_
#define GEOMETRY_PROPERTIES_H_

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include "GeometryProperties_fwd.H"

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBSupport.H>

#include "EmbeddedBoundaries.H"
#endif

#include <string>

class c_GeometryProperties
{
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi;
    amrex::GpuArray<int, AMREX_SPACEDIM>
        n_cell;  // number of cells in each dimension
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;
    amrex::IntVect max_grid_size;
    amrex::IntVect blocking_factor;
    std::array<int, AMREX_SPACEDIM> is_periodic;
    amrex::CoordSys::CoordType coord_sys;
    int embedded_boundary_flag;

    void ParseBasicDomainInput();
    void InitializeBoxArrayAndDistributionMap();

   public:
    amrex::BoxArray ba;  // a list of boxes that cover the domain
    amrex::DistributionMapping dm;
    amrex::Geometry geom;

#ifdef AMREX_USE_EB
    std::unique_ptr<c_EmbeddedBoundaries> pEB;
#endif

    c_GeometryProperties();
    void ReadData();
    void InitData();

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> &get_ProbLo() const
    {
        return prob_lo;
    }
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> &get_ProbHi() const
    {
        return prob_hi;
    }
    const amrex::GpuArray<int, AMREX_SPACEDIM> &get_NumCells() const
    {
        return n_cell;
    }
    const amrex::Array<amrex::Real, AMREX_SPACEDIM> &get_CellSize() const
    {
        return dx;
    }
    const amrex::IntVect &get_MaxGridSize() const { return max_grid_size; }
    const amrex::IntVect &get_BlockingFactor() const { return blocking_factor; }
    const std::array<int, AMREX_SPACEDIM> &get_PeriodicityArray() const
    {
        return is_periodic;
    }
    const amrex::CoordSys::CoordType &get_CoordinateSystem() const
    {
        return coord_sys;
    }
    bool is_eb_enabled() const { return embedded_boundary_flag; }

    bool Is_Point_Inside_Physical_Domain(const amrex::Real *pos);

    bool Is_Position_Inside_BoxArray(const amrex::Real *pos);

    int Get_BoxID(const amrex::Real *pos);
};
#endif
