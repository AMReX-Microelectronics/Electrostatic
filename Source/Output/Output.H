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


//class
//c_Write {
//
//}


class 
c_Output
{
public: 

    c_Output();
    ~c_Output();

    int ReadData();
    void InitData();
    void WriteOutput(int step, amrex::Real time);
    void WriteOutput_AfterInit();

    int m_write_after_init;
    const int m_plt_name_digits=4;

    int m_num_params_plot_single_level;

    std::string foldername_str;

private:

    std::map<std::string,s_Subscript::subscript_name> m_map_subscript_name;

    void DefineSubscriptNameMap()
    {
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_Output::DefineSubscriptNameMap()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

        m_map_subscript_name["_x"] = s_Subscript::x;
        m_map_subscript_name["_y"] = s_Subscript::y;
        m_map_subscript_name["_z"] = s_Subscript::z;

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_Output::DefineSubscriptNameMap()************************\n";
#endif
    };

    amrex::Vector< int > m_output_option; //1 or 2 specified in the input file or 0 if not specified.

    std::string m_filename_prefix_str;
    std::string m_plot_file_name;
    amrex::Vector<std::string> m_extra_dirs;
    int m_write_interval;
    int m_rawfield_write_interval;

    int m_raw_fields_to_plot;

    std::map<std::string,int> m_map_field_to_plot_after_init;

    amrex::Vector< amrex::MultiFab* > m_p_mf;
    std::unique_ptr<amrex::MultiFab> m_p_mf_all_init;

    void AssimilateDataPointers();

    const std::string m_default_level_prefix {"Level_"};
    std::map<std::string,int> m_map_param_all;
    std::unique_ptr<amrex::MultiFab> m_p_mf_all;

    void WriteSingleLevelPlotFile(int step, amrex::Real time, std::unique_ptr<amrex::MultiFab>& p_all_mf, std::map<std::string,int>& map_all_mf);
    void WriteRawFields(std::map<std::string,int>& map_all_mf);
};

#endif
