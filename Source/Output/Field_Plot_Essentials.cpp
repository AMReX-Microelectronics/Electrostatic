#include "Field_Plot_Essentials.H"

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

using namespace amrex;

c_Field_Plot_Essentials::c_Field_Plot_Essentials()
{
    _raw_fields_to_plot_flag = 0;
}

c_Field_Plot_Essentials::~c_Field_Plot_Essentials()
{
    _output_option.clear();
    _raw_fields_to_plot_flag = 0;
    _extra_dirs.clear();
}

void c_Field_Plot_Essentials::Init_Plot_Field_Essentials(
    const amrex::Geometry &geom, int flag)
{
    _geom = &geom;
    _embedded_boundary_flag = flag;
}

void c_Field_Plot_Essentials::WriteSingleLevelPlotFile(
    int step, amrex::Real time, const amrex::Vector<amrex::MultiFab *> &m_p_mf,
    std::unique_ptr<amrex::MultiFab> &p_all_mf,
    std::map<std::string, int> &map_all_mf)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Output::"
                      "WriteSingleLevelPlotFile()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__
                   << "\n";
    std::string prt = "\t\t\t";
#endif

    int Ncomp1 = 1;
    int Nghost0 = 0;

    amrex::Vector<std::string> m_p_name_str;

    int c = 0;
    for (auto it : map_all_mf)
    {
        auto field_number = it.second;
        auto field_name = it.first;

        if (_output_option[field_number] == 0 or
            _output_option[field_number] == 1)
        {
#ifdef PRINT_HIGH
            amrex::Print() << prt << " copying parameter: " << field_name
                           << " counter : " << c
                           << " field_number: " << field_number
                           << " _output_option: "
                           << _output_option[field_number] << "\n";
#endif

            m_p_name_str.push_back(field_name);

            amrex::MultiFab::Copy(*p_all_mf, *m_p_mf[field_number], 0, c,
                                  Ncomp1, Nghost0);
            ++c;
        }
    }

#ifdef AMREX_USE_EB
    if (_embedded_boundary_flag)
    {
        amrex::EB_WriteSingleLevelPlotfile(_plot_file_name, *p_all_mf,
                                           m_p_name_str, *_geom, time, step,
                                           "HyperCLaw-V1.1",
                                           _default_level_prefix, "Cell",
                                           _extra_dirs);
    }
#endif
    if (!_embedded_boundary_flag)
    {
        amrex::WriteSingleLevelPlotfile(_plot_file_name, *p_all_mf,
                                        m_p_name_str, *_geom, time, step,
                                        "HyperCLaw-V1.1", _default_level_prefix,
                                        "Cell", _extra_dirs);
    }

    m_p_name_str.clear();

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_Output::"
                      "WriteSingleLevelPlotFile()************************\n";
#endif
}

void c_Field_Plot_Essentials::WriteRawFields(
    const amrex::Vector<amrex::MultiFab *> &m_p_mf,
    std::map<std::string, int> &map_all_mf)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Output::"
                      "WriteRawFields()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__
                   << "\n";
    std::string prt = "\t\t\t";
#endif

    for (auto it : map_all_mf)
    {
        auto field_number = it.second;
        auto field_name = it.first;

        if (_output_option[field_number] == 1 or
            _output_option[field_number] == 2)
        {
            int lev = 0;
            const std::string raw_pltname = _plot_file_name + "/raw_fields";

            std::string prefix =
                amrex::MultiFabFileFullPrefix(lev, raw_pltname,
                                              _default_level_prefix,
                                              field_name);
            VisMF::Write(*m_p_mf[field_number], prefix);
        }
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_Output::WriteRawFields("
                      ")************************\n";
#endif
}

int c_Field_Plot_Essentials::SpecifyOutputOption(
    amrex::Vector<std::string> fields_to_plot_withGhost_str,
    std::map<std::string, int> &map_param_all)
{
    int num_params_without_option_2 = 0;
    const int varsize = fields_to_plot_withGhost_str.size();

    amrex::Vector<std::string> fields_to_plot_str;
    amrex::Vector<std::string> extra_str;

    fields_to_plot_str.resize(varsize);
    extra_str.resize(varsize);

    int c = 0;
    for (auto str : fields_to_plot_withGhost_str)
    {
        std::string second_last_str = str.substr(str.length() - 2, 1);
        if (second_last_str == ".")
        {
            std::string last_str = str.substr(str.length() - 1, 1);

            if (last_str == "1" or last_str == "2")
            {
                std::string string_without_ghost =
                    str.substr(0, str.length() - 2);
                extra_str[c] = str.substr(str.length() - 1, 1);
                fields_to_plot_str[c] = string_without_ghost;
                _raw_fields_to_plot_flag = true;
            }
            else
            {
                // assert unknown option.
            }
        }
        else if (second_last_str != ".")
        {
            extra_str[c] = "0";
            fields_to_plot_str[c] = str;
        }
        ++c;
    }

    std::map<std::string, int>::iterator it_map_param_all;

    c = 0;
    for (std::size_t i = 0; i < fields_to_plot_str.size(); ++i)
    {
        std::string str = fields_to_plot_str[i];
        std::string first_three_char_str = str.substr(0, 3);

        if (first_three_char_str == "vec")
        {
            std::string rest_of_string = str.substr(3);
            for (auto subscript : _map_subscript_name)
            {
                std::string component = rest_of_string + subscript.first;

                it_map_param_all = map_param_all.find(component);

                if (it_map_param_all == map_param_all.end())
                {
                    if (extra_str[i] == "0")
                    {
                        _output_option.push_back(0);
                        ++num_params_without_option_2;
                    }
                    else if (extra_str[i] == "1")
                    {
                        _output_option.push_back(1);
                        ++num_params_without_option_2;
                    }
                    else if (extra_str[i] == "2")
                    {
                        _output_option.push_back(2);
                    }
                    map_param_all[component] = c;
                    ++c;
                }
            }
        }
        else
        {
            it_map_param_all = map_param_all.find(str);

            if (it_map_param_all == map_param_all.end())
            {
                if (extra_str[i] == "0")
                {
                    _output_option.push_back(0);
                    ++num_params_without_option_2;
                }
                else if (extra_str[i] == "1")
                {
                    _output_option.push_back(1);
                    ++num_params_without_option_2;
                }
                else if (extra_str[i] == "2")
                {
                    _output_option.push_back(2);
                }
                map_param_all[str] = c;
                ++c;
            }
        }
    }
    fields_to_plot_str.clear();
    extra_str.clear();

    return num_params_without_option_2;
}
