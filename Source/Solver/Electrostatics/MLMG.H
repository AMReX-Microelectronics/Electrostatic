#ifndef MLMG_H_
#define MLMG_H_

#include "MLMG_fwd.H"
#include "Input/BoundaryConditions/BoundaryConditions_fwd.H"

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>

//#include <AMReX_BaseFwd.H>



#include <string>


/**
 * Multi-Level Multi-Grid solver is used which solves for phi in (A*alpha_cc - B * div (beta_fc grad)) phi = rho
 * For Poisson solver, A and alpha_cc are set to zero.
 * see AMReX_MLLinOp.H for the defaults, accessors, and mutators
 */

template<typename T>
class TD;

class 
c_MLMGSolver
{
public: 
    amrex::Real ascalar, bscalar;
    int set_verbose;
    int max_order;
    amrex::Real  relative_tolerance, absolute_tolerance;

    amrex::LPInfo info;
    amrex::MLABecLaplacian mlabec;
    std::unique_ptr<amrex::MLMG> pMLMG;
    std::array< amrex::MultiFab, AMREX_SPACEDIM > beta_fc; 
    c_MLMGSolver();
    ~c_MLMGSolver();


    void InitData();

    void Solve_PoissonEqn();

    void Compute_vecField(std::array<amrex::MultiFab, AMREX_SPACEDIM>& E);

    void Compute_vecFlux(std::array<amrex::MultiFab, AMREX_SPACEDIM>& flux);

    //amrex::MultiFab& get_mf  (std::string macro_str) {auto macro_num = param_map[macro_str]; return (*m_p_mf[macro_num]);}
    //amrex::MultiFab* get_p_mf (std::string macro_str) {auto macro_num = param_map[macro_str]; return m_p_mf[macro_num].get();}

    std::array<std::array<amrex::LinOpBCType,AMREX_SPACEDIM>,2> LinOpBCType_2d;


private:
    amrex::MultiFab* alpha=nullptr;
    amrex::MultiFab* beta=nullptr;
    amrex::MultiFab* soln=nullptr;
    amrex::MultiFab*  rhs=nullptr;

    amrex::MultiFab* robin_a=nullptr;
    amrex::MultiFab* robin_b=nullptr;
    amrex::MultiFab* robin_f=nullptr;

    std::string alpha_str;
    std::string beta_str;
    std::string soln_str;
    std::string rhs_str;
    std::string robin_a_str;
    std::string robin_b_str;
    std::string robin_f_str;


    void ReadData();

    void Set_MLMG_CellCentered_Multifabs();

    void Setup_MLABecLaplacian_ForPoissonEqn();

    void AssignLinOpBCTypeToBoundaries ();

    bool all_homogeneous_boundaries;
    bool some_functionbased_inhomogeneous_boundaries;
    bool some_constant_inhomogeneous_boundaries;

    void Fill_Constant_Inhomogeneous_Boundaries();
    void Fill_FunctionBased_Inhomogeneous_Boundaries();

};

#endif
