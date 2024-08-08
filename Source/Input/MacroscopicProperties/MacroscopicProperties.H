#ifndef MACROSCOPIC_PROPERTIES_H_
#define MACROSCOPIC_PROPERTIES_H_

#include "MacroscopicProperties_fwd.H"
#include "ChargeDensitySource.H"

#include "../Utils/SelectWarpXUtils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>

#include <string>
#include <any>


struct s_MacroscopicPropertiesMacroName {
    enum macro_name { epsilon, charge_density, phi };
};


class 
c_MacroscopicProperties
{

private:

    std::map<std::string,amrex::Real> map_default_value;
    std::map<std::string,int> map_num_ghostcell;

    std::unique_ptr<PointChargeSource> p_PointChargeSource = nullptr;

    void Define_ExternalChargeDensitySources();

#ifdef TIME_DEPENDENT
    amrex::Vector<std::string> m_parser_varname_vector = {"x", "y", "z", "t"};
#else
    amrex::Vector<std::string> m_parser_varname_vector = {"x", "y", "z"};
#endif

    void DefineParameterNameMap()
    {
        map_macro_name["epsilon"] = s_MacroscopicPropertiesMacroName::epsilon;
        map_macro_name["charge_density"] = s_MacroscopicPropertiesMacroName::charge_density;
        map_macro_name["phi"] = s_MacroscopicPropertiesMacroName::phi;
    };

    void DefineDefaultValueMap()
    {
        map_default_value["epsilon"] = PhysConst::ep0;
        map_default_value["charge_density"] = 0.0;
        map_default_value["phi"] = 0.0;
    };

    void ReadPermittivity();
    void InitializePermittivity();

    int ReadParameterMapAndNumberOfGhostCells();

    void DefineMacroVariableVectorSizes();

    amrex::Vector<std::unique_ptr<amrex::Parser>> m_p_macro_parser; //pointer to parser object
    amrex::Vector<std::string> m_macro_str_function;
    amrex::Vector<std::string> m_macro_type;
    amrex::Vector<amrex::Real> m_macro_value;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_p_mf; //Pointer to multifab of permittivity, eps

    template<class T> 
        void ReadMacroparam(std::string, T default_value);

    void DefineAndInitializeMacroparam(std::string macro_str,
                                       int macro_num,
                                       amrex::BoxArray& ba,
                                       amrex::DistributionMapping& dm,
                                       amrex::Geometry& geom,
                                       int Ncomp, int Nghost);

public: 

    c_MacroscopicProperties();

    void ReadData();
    void InitData();

    amrex::MultiFab& get_mf  (std::string macro_str) 
    {   
        auto macro_num = map_param_all[macro_str]; 
        return (*m_p_mf[macro_num]);
    }
    amrex::MultiFab* get_p_mf (std::string macro_str) 
    {    
        if(map_param_all.find(macro_str) != map_param_all.end()) {
            auto macro_num = map_param_all[macro_str]; 
            return m_p_mf[macro_num].get();
        }
        return nullptr;
    }

    std::map<std::string,s_MacroscopicPropertiesMacroName::macro_name> map_macro_name;
    std::map<std::string,int> map_param_all;

    int num_params;
    void ReInitializeMacroparam(std::string macro_str);

    void Deposit_AllExternalChargeDensitySources();

};
#endif
