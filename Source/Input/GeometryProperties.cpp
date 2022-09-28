#include "GeometryProperties.H"

#include "../Utils/WarpXUtil.H"
//#include "Code.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_RealBox.H>

using namespace amrex;

c_GeometryProperties::c_GeometryProperties ()
{
    ReadData();
} 

void 
c_GeometryProperties::ReadData()
{    
     ParseBasicDomainInput();
}


void 
c_GeometryProperties::InitData()
{
    InitializeBoxArrayAndDistributionMap();
}


void
c_GeometryProperties::ParseBasicDomainInput()
{

    amrex::Vector<int> num_cell;
    amrex::Vector<amrex::Real> prob_min(AMREX_SPACEDIM);
    amrex::Vector<amrex::Real> prob_max(AMREX_SPACEDIM);
    amrex::Vector<amrex::Real> mg(AMREX_SPACEDIM);
    amrex::Vector<amrex::Real> bf(AMREX_SPACEDIM);

    amrex::ParmParse pp_domain("domain");

    getArrWithParser(pp_domain, "prob_lo", prob_min, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_lo.size() == AMREX_SPACEDIM);

    getArrWithParser(pp_domain, "prob_hi", prob_max, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_hi.size() == AMREX_SPACEDIM);

    pp_domain.getarr("n_cell", num_cell, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(n_cell.size() == AMREX_SPACEDIM);

    pp_domain.queryarr("max_grid_size", mg);

    pp_domain.queryarr("blocking_factor", bf);
    bf.resize(std::max(static_cast<int>(bf.size()),1),8);

    pp_domain.addarr("n_cell", num_cell);
    pp_domain.addarr("prob_lo", prob_min);
    pp_domain.addarr("prob_hi", prob_max);
    pp_domain.addarr("max_grid_size", mg);
    pp_domain.addarr("blocking_factor", bf);

    for (int i=0; i<AMREX_SPACEDIM; ++i) 
    {
        n_cell[i] = num_cell[i];
        prob_lo[i] = prob_min[i]; //Converting vector to GpuArray
        prob_hi[i] = prob_max[i]; 
        max_grid_size[i] = mg[i];  //Converting Vector to IntVect
        blocking_factor[i] = bf[i]; 
    }

}

void 
c_GeometryProperties::InitializeBoxArrayAndDistributionMap()
{

    amrex::IntVect dom_lo(AMREX_D_DECL(       0,        0,        0)); // domain low
    amrex::IntVect dom_hi(AMREX_D_DECL(n_cell[0]-1, n_cell[1]-1, n_cell[2]-1)); // domain high

    amrex::Box domain(dom_lo, dom_hi); // Make a single box that is the entire domain

    ba.define(domain); // initialize the boxarray 'ba' from the single box 'domain'

    ba.maxSize(max_grid_size); // break up ba into chunks no larger than 'max_grid_size' along a direction

    amrex::RealBox real_box({AMREX_D_DECL( prob_lo[0], prob_lo[1], prob_lo[2])},
                    {AMREX_D_DECL( prob_hi[0], prob_hi[1], prob_hi[2])});  //physical domain

    amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)}; // 0: not periodic, 1: periodic

    geom.define(domain, real_box, CoordSys::cartesian, is_periodic); //define the geom object

    dm.define(ba);

}
