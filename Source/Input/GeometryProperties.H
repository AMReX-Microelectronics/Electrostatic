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

//struct s_Domain
//{
//
//    amrex::GpuArray<int, AMREX_SPACEDIM> n_cell;     // number of cells in each dimension
//    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo; // physical lo coordinate
//    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi; // physical hi coordinate
//    amrex::Vector<int> max_grid_size;
//    amrex::Vector<int> blocking_factor;
//
//};


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
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx;

    void ReadData();
    void InitData();
//    amrex::MultiFab& get_epsilon_mf() {return (*m_p_epsilon_mf);}

private:
    void ParseBasicDomainInput();
    void InitializeBoxArrayAndDistributionMap ();
//    void ReadPermittivity();
//    void InitializePermittivity();
//    void InitializeMacroMultiFabUsingParser (amrex::MultiFab *macro_mf,
//                                             amrex::ParserExecutor<3> const& macro_parser,
//                                             const int lev);

//    std::string m_epsilon_s = "constant";
//    amrex::Real m_epsilon = PhysConst::ep0;
 //   amrex::Real m_epsilon = 8.8541878128e-12;
 //   std::unique_ptr<amrex::MultiFab> m_p_epsilon_mf; //Pointer to multifab of permittivity, eps
 //   std::unique_ptr<amrex::Parser> m_p_epsilon_parser; //Pointer to parser object to store function parser
 //   std::string m_str_epsilon_function; //string for storing epsilon function
};

#endif
