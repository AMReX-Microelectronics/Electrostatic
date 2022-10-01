#ifndef BOUNDARY_CONDITIONS_H_
#define BOUNDARY_CONDITIONS_H_

#include "BoundaryConditions_fwd.H"

#include "../Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>

#include <string>
#include <vector>
#include <any>
//#include <AMReX_Extension.H>
//#include <AMReX_GpuQualifiers.H>


struct s_BoundaryConditions {
    enum boundary_type { dir, neu, per, rob, ref };
};


class 
c_BoundaryConditions
{
public: 

    c_BoundaryConditions();
    ~c_BoundaryConditions();

    void ReadData();
    void InitData();

    amrex::MultiFab& get_mf  (std::string macro_str) {auto macro_num = map_function_parser_name[macro_str]; return (*m_p_mf[macro_num]);}
    amrex::MultiFab* get_p_mf (std::string macro_str) {auto macro_num = map_function_parser_name[macro_str]; return m_p_mf[macro_num].get();}

    std::map<std::string, s_BoundaryConditions::boundary_type> map_boundary_type;

    /*The following 2d arrays are for lo and hi sides each with AMREX_SPACEDIM */
    std::array<std::array<std::any,AMREX_SPACEDIM>,2> bcAny_2d; 
    std::array<std::array<std::string,AMREX_SPACEDIM>,2> bcType_2d;
    std::map<int,std::string> map_bcAny_2d[2];
    std::map<std::string, int> map_function_parser_name;

private:
    void DefineBoundaryTypeMap()
    {
        map_boundary_type["dir"] = s_BoundaryConditions::dir;
        map_boundary_type["neu"] = s_BoundaryConditions::neu;
        map_boundary_type["per"] = s_BoundaryConditions::per;
        map_boundary_type["rob"] = s_BoundaryConditions::rob;
        map_boundary_type["ref"] = s_BoundaryConditions::ref;
    };

    void ReadBoundaryConditionsType();

    void SortBoundaryTypeArrayString(const amrex::Vector<std::string>& bc_str, 
                                     std::array< std::string, AMREX_SPACEDIM >& bcType, 
                                     std::array< std::any, AMREX_SPACEDIM >& bcAny, 
                                     std::map<int,std::string>& map_bcAny);

    std::array<std::string,3> robin_substrings { "_a", "_b", "_f"};

    int num_function_parsers;
    void DefineMacroVariableVectorSizes();

    amrex::Vector<std::unique_ptr<amrex::Parser>> m_p_macro_parser; //pointer to parser object
    amrex::Vector<std::string> m_macro_str_function;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_p_mf; //Pointer to multifab of permittivity, eps

    void ReadBoundaryConditionsParser(std::string macro_str, int macro_num);
#ifdef TIME_DEPENDENT
    amrex::Vector<std::string> m_parser_varname_vector = {"x", "y", "z", "t"};
#else
    amrex::Vector<std::string> m_parser_varname_vector = {"x", "y", "z"};
#endif


};

#endif
