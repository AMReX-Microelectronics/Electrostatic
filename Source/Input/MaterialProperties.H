#ifndef MATERIAL_PROPERTIES_H_
#define MATERIAL_PROPERTIES_H_

#include "MaterialProperties_fwd.H"

#include "../Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>

#include <string>
//#include <AMReX_Extension.H>
//#include <AMReX_GpuQualifiers.H>

class 
c_MaterialProperties
{
public: 
    c_MaterialProperties();
    void ReadData();
    void InitData();
    amrex::MultiFab& get_epsilon_mf() {return (*m_p_epsilon_mf);}

private:
    void ReadPermittivity();
    void InitializePermittivity();
//    void InitializeMacroMultiFabUsingParser (amrex::MultiFab *macro_mf,
//                                             amrex::ParserExecutor<3> const& macro_parser,
//                                             const int lev);

    std::string m_epsilon_s = "constant";
//    amrex::Real m_epsilon = PhysConst::ep0;
    amrex::Real m_epsilon = 8.8541878128e-12;
    std::unique_ptr<amrex::MultiFab> m_p_epsilon_mf; //Pointer to multifab of permittivity, eps
    std::unique_ptr<amrex::Parser> m_p_epsilon_parser; //Pointer to parser object to store function parser
    std::string m_str_epsilon_function; //string for storing epsilon function
};

#endif
