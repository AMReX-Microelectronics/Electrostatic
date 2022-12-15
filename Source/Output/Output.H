#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "Output_fwd.H"
#include "Field_Plot_Essentials.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

#include <string>


class 
c_Output
       : public c_Field_Plot_Essentials
{
public: 

    c_Output();
    ~c_Output();

    int ReadData();
    void InitData();
    void WriteOutput(int step, amrex::Real time);
    void WriteOutput_AfterInit();

    int m_write_after_init;

    int m_num_params_plot_single_level;

private:

    amrex::Vector< amrex::MultiFab* > m_p_mf;

    std::map<std::string,int> m_map_param_all;
    std::map<std::string,int> m_map_field_to_plot_after_init;

    std::unique_ptr<amrex::MultiFab> m_p_mf_all;
    std::unique_ptr<amrex::MultiFab> m_p_mf_all_init;

    void AssimilateDataPointers();

};

#endif
