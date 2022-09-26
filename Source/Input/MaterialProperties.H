#ifndef MACROSCOPIC_PROPERTIES_H_
#define MACROSCOPIC_PROPERTIES_H_

#include "MacroscopicProperties_fwd.H"

#include "../Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>

#include <string>
//#include <AMReX_Extension.H>
//#include <AMReX_GpuQualifiers.H>

class 
c_MacroscopicProperties
{
public: 

    c_MacroscopicProperties();
    ~c_MacroscopicProperties();

    void ReadData();
    void InitData();

    amrex::MultiFab& get_epsilon_mf  () {return (*m_p_epsilon_mf);}
    amrex::MultiFab* get_p_epsilon_mf () {return m_p_epsilon_mf.get();}


private:
    const int Ncomp1=1;
    const int Nghost1=1;
    const int Nghost0=0;

    void ReadPermittivity();
    void InitializePermittivity();
    void InitializeMacroMultiFabUsingParser (amrex::MultiFab *mf,
                                             amrex::ParserExecutor<3> const& parser);

    int DefineParameterMap();
    map<string,int> param_map;
    int num_params;

    void DefineMacroVariableVectorSizes();

    std::Vector<std::unique_ptr<amrex::Parser>> m_p_macro_parser; //pointer to parser object
    std::Vector<std::string> m_macro_str_function;
    std::Vector<std::string> m_macro_type;
    std::Vector<amrex::Real> m_macro_value;

    std::unique_ptr<amrex::Parser> m_p_epsilon_parser; //Pointer to parser object to store function parser
    std::string m_str_epsilon_function; //string for storing epsilon function
    std::string m_epsilon_s = "constant";
//    amrex::Real m_epsilon = PhysConst::ep0;
    amrex::Real m_epsilon = 8.8541878128e-12;

    std::unique_ptr<amrex::MultiFab> m_p_epsilon_mf; //Pointer to multifab of permittivity, eps
};

#endif
