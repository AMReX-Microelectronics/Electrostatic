#ifndef GEOMETRY_PROPERTIES_H_
#define GEOMETRY_PROPERTIES_H_

#include "GeometryProperties_fwd.H"

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>

#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBSupport.H>


#include <string>
class
c_EmbeddedBoundary
{
public:
    c_EmbeddedBoundary();

    amrex::EBSupport support;
    std::unique_ptr<amrex::EBFArrayBoxFactory> pFactory;
    int required_coarsening_level;
    int max_coarsening_level;

    void ReadGeometry();
private:

    std::map<std::string,amrex::EBSupport>
    map_eb_support = {
                     {"full", amrex::EBSupport::full},
             		 {"basic", amrex::EBSupport::basic},
        	    	 {"volume",amrex::EBSupport::volume},
                     {"Full", amrex::EBSupport::full},
        		     {"Basic", amrex::EBSupport::basic},
            		 {"Volume",amrex::EBSupport::volume},
                     };
};


class 
c_GeometryProperties:public c_EmbeddedBoundary
{
public: 
    c_GeometryProperties();
    ~c_GeometryProperties();

    amrex::GpuArray<int, AMREX_SPACEDIM>  n_cell;     // number of cells in each dimension
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>  prob_lo; // physical lo coordinate
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>  prob_hi; // physical hi coordinate
    amrex::IntVect max_grid_size;
    amrex::IntVect blocking_factor;
    std::array<int, AMREX_SPACEDIM> is_periodic;
    amrex::CoordSys::CoordType coord_sys;

    amrex::BoxArray ba; // a list of boxes that cover the domain
    amrex::Geometry geom;
    amrex::DistributionMapping dm;

    int embedded_boundary_flag;
    c_EmbeddedBoundary eb;

    void ReadData();
    void InitData();

    void CreateEmbeddedBoundary();

private:
    void ParseBasicDomainInput();
    void InitializeBoxArrayAndDistributionMap ();

};

#endif
