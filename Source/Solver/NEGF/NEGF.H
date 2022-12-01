#ifndef NEGF_H_
#define NEGF_H_

#include "NEGF_fwd.H"
#include "../../Code_definitions.H"
//#include "Input/BoundaryConditions/BoundaryConditions_fwd.H"
//
//#include <AMReX_REAL.H>
//#include <AMReX_GpuContainers.H>
//#include <AMReX_Vector.H>
//#include <AMReX_IntVect.H>
//#include <AMReX_BoxArray.H>
//#include <AMReX_Geometry.H>
//#include <AMReX_MLLinOp.H>
//#include <AMReX_MLABecLaplacian.H>
//#include <AMReX_MLEBABecLap.H>
//#include <AMReX_MLMG.H>
//
////#include <AMReX_BaseFwd.H>

#include <string>

class 
c_NEGFSolver:public c_CNT
{
public: 
    c_AtomicStructure cnt;

//    amrex::Real ascalar, bscalar;
//    int set_verbose;
//    int max_order;
//    amrex::Real  relative_tolerance, absolute_tolerance;
//
//    amrex::LPInfo info;
//    std::unique_ptr<amrex::MLMG> pMLMG;
//    std::array< amrex::MultiFab, AMREX_SPACEDIM > beta_fc; 
//
//#ifdef AMREX_USE_EB
//    std::unique_ptr<amrex::MLEBABecLap> p_mlebabec;
//#endif
//    std::unique_ptr<amrex::MLABecLaplacian> p_mlabec;
//
//    c_MLMGSolver();
//    ~c_MLMGSolver();
//
//
//    void InitData();
//
//    amrex::Real Solve_PoissonEqn();
//
//    void Compute_vecField(std::array<amrex::MultiFab, AMREX_SPACEDIM>& E);
//
//    void Compute_vecFlux(std::array<amrex::MultiFab, AMREX_SPACEDIM>& flux);
//
//    //amrex::MultiFab& get_mf  (std::string macro_str) {auto macro_num = param_map[macro_str]; return (*m_p_mf[macro_num]);}
//    //amrex::MultiFab* get_p_mf (std::string macro_str) {auto macro_num = param_map[macro_str]; return m_p_mf[macro_num].get();}
//
//    std::array<std::array<amrex::LinOpBCType,AMREX_SPACEDIM>,2> LinOpBCType_2d;
//
//    void UpdateBoundaryConditions();
//
//    /*The following is public for GPUs*/
//    void Fill_Constant_Inhomogeneous_Boundaries(); 
//    void Fill_FunctionBased_Inhomogeneous_Boundaries(); 
//    std::string soln_str;
//    amrex::MultiFab* soln=nullptr;
//
private:
//    amrex::MultiFab* alpha=nullptr;
//    amrex::MultiFab* beta=nullptr;
//    amrex::MultiFab*  rhs=nullptr;
//
//    amrex::MultiFab* robin_a=nullptr;
//    amrex::MultiFab* robin_b=nullptr;
//    amrex::MultiFab* robin_f=nullptr;
//
//    std::string alpha_str;
//    std::string beta_str;
//    std::string rhs_str;
//    std::string robin_a_str;
//    std::string robin_b_str;
//    std::string robin_f_str;
//
//
//    void ReadData();
//
//    void Set_MLMG_CellCentered_Multifabs();
//
//    void Setup_MLABecLaplacian_ForPoissonEqn();
//
//#ifdef AMREX_USE_EB 
//    void Setup_MLEBABecLaplacian_ForPoissonEqn();
//#endif
//
//    void AssignLinOpBCTypeToBoundaries ();
//
//    bool all_homogeneous_boundaries;
//    bool some_functionbased_inhomogeneous_boundaries;
//    bool some_constant_inhomogeneous_boundaries;
//
//
};

#endif
