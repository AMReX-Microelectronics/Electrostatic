#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "Output_fwd.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

#include <string>

struct s_Subscript {
    enum subscript_name { x, y, z };
};


class 
c_Output
{
public: 

    c_Output();
    ~c_Output();
//
    int ReadData();
    void InitData();
    void WriteSingleLevelPlotFile(int step, amrex::Real time);
    void AssimilateDataPointers();

    std::map<std::string,s_Subscript::subscript_name> map_subscript_name;

private:
  void DefineSubscriptNameMap()
    {
        map_subscript_name["_x"] = s_Subscript::x;
        map_subscript_name["_y"] = s_Subscript::y;
        map_subscript_name["_z"] = s_Subscript::z;
    };

    std::unique_ptr<amrex::MultiFab> m_p_mf_all;
    std::string plot_file_name;
    const int plt_name_digits=5;

    std::map<std::string,int> map_param_all;
    int num_all_params;

    amrex::Vector< amrex::MultiFab* > m_p_mf;
    amrex::Vector< std::string > m_p_name_str;

    int Evaluate_TypeOf_MacroStr(std::string macro_str);

};

#endif
