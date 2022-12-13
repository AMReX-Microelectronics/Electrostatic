#ifndef FIELD_PLOT_ESSENTIALS_H_
#define FIELD_PLOT_ESSENTIALS_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

#include <string>

struct s_Subscript {
    enum subscript_name { x, y, z };
};

class
c_Field_Plot_Essentials 
{
public:
    c_Field_Plot_Essentials();
   ~c_Field_Plot_Essentials(); 
    std::string get_folder_name() { return _foldername_str; };

protected:

    std::string _plt_str = "/plt";
    const int _plt_name_digits=4;
    amrex::Vector<std::string> _extra_dirs;
    int _write_interval;
    int _rawfield_write_interval;
    std::string _filename_prefix_str;
    std::string _plot_file_name;
    std::string _foldername_str;
    amrex::Vector< int > _output_option; 
    //1 or 2 specified in the input file or 0 if not specified.
    int _raw_fields_to_plot_flag;

    const std::string _default_level_prefix {"Level_"};

    void Init_Plot_Field_Essentials(const amrex::Geometry& geom, int flag);

    void WriteSingleLevelPlotFile(int step, 
                                  amrex::Real time, 
                                  const amrex::Vector<amrex::MultiFab*>& m_p_mf,
                                  std::unique_ptr<amrex::MultiFab>& p_all_mf, 
                                  std::map<std::string,int>& map_all_mf);

    void WriteRawFields(const amrex::Vector<amrex::MultiFab*>& m_p_mf, 
                        std::map<std::string,int>& map_all_mf);

    int SpecifyOutputOption
        (amrex::Vector<std::string> fields_to_plot_withGhost_str, 
         std::map<std::string,int>& m_map_param_all);

    std::map<std::string,s_Subscript::subscript_name>
    _map_subscript_name =
    {
        {"_x", s_Subscript::x },
        {"_y", s_Subscript::y },
        {"_z", s_Subscript::z }
    };

private:

    const amrex::Geometry * _geom;
    int _embedded_boundary_flag;

};
#endif
