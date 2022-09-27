#ifndef POST_PROCESSOR_H_
#define POST_PROCESSOR_H_

#include "PostProcessor_fwd.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

#include <string>


struct s_PostProcessMacroName {
    enum macro_name { vecE, vecFlux };
};

class 
c_PostProcessor
{
public: 

    c_PostProcessor();
    ~c_PostProcessor();

    int ReadData();
    void InitData();
    void Compute();

    amrex::MultiFab& get_mf  (std::string macro_name) {auto m = map_param_mf[macro_name]; return (*m_p_mf[m]);}
    amrex::MultiFab* get_p_mf (std::string macro_name) {auto m = map_param_mf[macro_name]; return m_p_mf[m].get();}

    std::array<amrex::MultiFab, AMREX_SPACEDIM>&
    get_array_mf (std::string macro_name) 
    {
        auto m = map_param_arraymf[macro_name];
        return *m_p_array_mf[m];
    }

    std::array<amrex::MultiFab, AMREX_SPACEDIM>*
    get_p_array_mf (std::string macro_name) 
    {
        auto m = map_param_arraymf[macro_name];
        return m_p_array_mf[m].get();
    }

    amrex::MultiFab&
    get_array_mf_component (std::string macro_name, int component) 
    {
        auto m = map_param_arraymf[macro_name];
        return (*m_p_array_mf[m])[component];
    }

    amrex::MultiFab*
    get_p_array_mf_component (std::string macro_name, int component) 
    {
        auto m = map_param_arraymf[macro_name];
        return &(*m_p_array_mf[m])[component];
    }

    std::map<std::string,s_PostProcessMacroName::macro_name> map_macro_name;

private:
    void DefineParameterNameMap()
    {
        map_macro_name["vecE"] = s_PostProcessMacroName::vecE;
        map_macro_name["vecFlux"] = s_PostProcessMacroName::vecFlux;
    };
    int SortParameterType();

    std::map<std::string,int> map_param_all;
    std::map<std::string,int> map_param_mf;
    std::map<std::string,int> map_param_arraymf;
    int num_all_params;
    int num_mf_params;
    int num_arraymf_params;


    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_p_mf;
    amrex::Vector<  std::unique_ptr < std::array<amrex::MultiFab, AMREX_SPACEDIM> >  > m_p_array_mf;

};

#endif
