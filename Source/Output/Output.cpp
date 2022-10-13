#include "Output.H"

#include "Code.H"
#include "Input/GeometryProperties/GeometryProperties.H"
#include "Input/MacroscopicProperties/MacroscopicProperties.H"
#include "PostProcessor/PostProcessor.H"

#include <AMReX_ParmParse.H>
#include "AMReX_PlotFileUtil.H"
#include <AMReX_VisMF.H>


using namespace amrex;


c_Output::c_Output ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Output Constructor()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    DefineSubscriptNameMap();
    num_all_params = ReadData();

    m_p_mf.resize(num_all_params);
    m_p_name_str.resize(num_all_params);

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_Output Constructor()************************\n";
#endif
} 


c_Output::~c_Output ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Output Destructor()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    map_param_all.clear();
    m_p_mf.clear();
    num_ghost.clear();
    m_p_name_str.clear();

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_Output Destructor()************************\n";
#endif
} 


//void
//c_Output::Reset ()
//{
//#ifdef PRINT_NAME
//    amrex::Print() << "\n\n\t\t\t{************************c_Output::Reset()************************\n";
//    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
//#endif
//   
//    num_all_params = 0;
//    map_param_all.clear();
//    m_p_mf.clear();
//    num_ghost.clear();
//    m_p_name_str.clear();
//
//#ifdef PRINT_NAME
//    amrex::Print() << "\t\t\t}************************c_Output::Reset()************************\n";
//#endif
//} 


int 
c_Output::ReadData()
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_Output::ReadData()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t\t\t";
#endif

    amrex::Vector< std::string > fields_to_plot_withGhost_str;
    amrex::Vector< std::string > fields_to_plot_str;
    amrex::Vector< std::string > ghost_str;
    amrex::ParmParse pp_plot("plot");
    
    filename_str = "plt";

    pp_plot.query("filename", filename_str);

    amrex::ParmParse pp_plot_file(filename_str);
    bool varnames_specified = pp_plot_file.queryarr("fields_to_plot", fields_to_plot_withGhost_str);

    const int varsize = fields_to_plot_withGhost_str.size();

    fields_to_plot_str.resize(varsize);
    ghost_str.resize(varsize);

    int c=0;
    for (auto str: fields_to_plot_withGhost_str) 
    {
        std::string second_last_str = str.substr(str.length()-2,1);

        if(second_last_str == ".") 
        { 
             std::string string_without_ghost = str.substr(0,str.length()-2);
             ghost_str[c] = str.substr(str.length()-1,1);
             fields_to_plot_str[c] = string_without_ghost;
        }
        else {
             fields_to_plot_str[c] = str;
        }
        ++c;
    }

    std::map<std::string,int>::iterator it_map_param_all;

    c=0;
    for (std::size_t i = 0; i < fields_to_plot_str.size(); ++i)
    {
        std::string str = fields_to_plot_str[i];
        std::string first_three_char_str = str.substr(0, 3);

        if(first_three_char_str == "vec") 
        {
            std::string rest_of_string = str.substr(3);
            for(auto subscript : map_subscript_name) 
            {
                std::string component = rest_of_string + subscript.first;

                it_map_param_all = map_param_all.find(component);

                if (it_map_param_all == map_param_all.end()) {
                    map_param_all[component] = c;
                    if(ghost_str[i] == "1") {
                       num_ghost.push_back(1); 
                    }
                    else {
                       num_ghost.push_back(0); 
                    }
                    ++c;
                }  
            }   
        }
        else {

            it_map_param_all = map_param_all.find(str);

            if (it_map_param_all == map_param_all.end()) {
                map_param_all[str] = c;
                if(ghost_str[i] == "1") {
                   num_ghost.push_back(1); 
                }
                else {
                   num_ghost.push_back(0); 
                }
                ++c;
            }
        }
    }
#ifdef PRINT_HIGH
    for (auto it: map_param_all) 
    {
        amrex::Print() << prt << "field: " << it.first << " counter " << it.second << " ghost cells: " << num_ghost[it.second] << "\n";
    }
#endif
    fields_to_plot_withGhost_str.clear();
    fields_to_plot_str.clear();
    ghost_str.clear();

#ifdef PRINT_LOW
    amrex::Print() << prt <<  "fields to plot: \n";
    for (auto it: map_param_all) 
    {
        amrex::Print() << prt <<  it.first << "   " << it.second << "\n";
    }
    amrex::Print() << prt << "total parameters to plot (final): " << map_param_all.size() << "\n\n";
#endif 

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_Output::ReadData()************************\n";
#endif

    return map_param_all.size();

}


void 
c_Output::InitData()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t{************************c_Output::InitData()************************\n";
    amrex::Print() << "\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif
    for (auto it : map_param_all) {
        m_p_name_str[it.second] = it.first;
    }

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;
    int Nghost0=0;

    m_p_mf_all = std::make_unique<amrex::MultiFab>(ba, dm, num_all_params, Nghost0); //cell-centered multifab

#ifdef PRINT_NAME
    amrex::Print() << "\t\t}************************c_Output::InitData()************************\n";
#endif
}


void
c_Output::WriteSingleLevelPlotFile(int step, amrex::Real time)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t{************************c_Output::WriteSingleLevelPlotFile()************************\n";
    amrex::Print() << "\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t";
#endif

    AssimilateDataPointers();

    plot_file_name = amrex::Concatenate(filename_str, step, plt_name_digits);

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& geom = rGprop.geom;
    int Ncomp1=1;
    int Nghost0=0;

    for (auto it : m_p_name_str) {
        int field_number = map_param_all[it];

        #ifdef PRINT_HIGH
        amrex::Print() << prt << "copying parameter: " << it << " field_number: " << field_number << " num_ghost: " << num_ghost[field_number] << "\n";
        #endif
        amrex::MultiFab::Copy(*m_p_mf_all, *m_p_mf[field_number], 0, field_number, Ncomp1, Nghost0);

    }

    int iteration=0;
    Vector<std::string> extra_dirs;
    if(num_ghost.size() > 0)  extra_dirs.emplace_back("raw_fields");

    amrex::WriteSingleLevelPlotfile( plot_file_name, 
                                    *m_p_mf_all, m_p_name_str, 
                                     geom, 
                                     time, iteration, 
                                     "HyperCLaw-V1.1", default_level_prefix, "Cell",
                                     extra_dirs);


    /*Writing raw fields with ghost cells*/
    for (auto it : m_p_name_str) 
    {
        int field_number = map_param_all[it];
        if(num_ghost[field_number] == 1) {
            int lev = 0;
            const std::string raw_pltname = plot_file_name + "/raw_fields";
            std::string prefix = amrex::MultiFabFileFullPrefix(lev,
                                raw_pltname, default_level_prefix, m_p_name_str[field_number]);
            VisMF::Write(*m_p_mf[field_number], prefix);
        }
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t}************************c_Output::WriteSingleLevelPlotFile()************************\n";
#endif
}

    
void
c_Output::AssimilateDataPointers()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Output::AssimilateDataPointers()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t\t";
#endif

    auto& rCode = c_Code::GetInstance();
    auto& rPost = rCode.get_PostProcessor(); 
    auto& rMprop = rCode.get_MacroscopicProperties();

    for (auto it : map_param_all) {
        std::string str = it.first;
        int field_counter = it.second;

        std::string second_last_str = str.substr(str.length()-2, 1);
        std::string last_str = str.substr(str.length()-1);
        #ifdef PRINT_HIGH
        amrex::Print() << "\n" << prt << "macro_str " << str << " second_last_str " << second_last_str << "\n";
        #endif

        if( second_last_str != "_" or (second_last_str == "_" and last_str != "x" and last_str != "y" and last_str != "z"))
        {
            #ifdef PRINT_HIGH
            amrex::Print() << prt << str << " doesn't have a subscript \n";
            #endif
            switch ( Evaluate_TypeOf_MacroStr(str) ) 
            {
                case 0:
                {
                    #ifdef PRINT_HIGH
                    amrex::Print() << prt << "field counter: " << field_counter << " fetching pointer for : " << str << "\n";
                    #endif
                    m_p_mf[field_counter] = rMprop.get_p_mf(str);
                    break;
                }
                case 1:
                {
                    #ifdef PRINT_HIGH
                    amrex::Print() << prt << "field counter: " << field_counter << " fetching pointer for : " << str  << "\n";
                    #endif
                    m_p_mf[field_counter] = rPost.get_p_mf(str);
                    break;
                }
                default :
                {
                    #ifdef PRINT_HIGH
                    amrex::Print() << prt << "field counter: " << field_counter << " fetching pointer for: " << str<< "\n";
                    #endif
                    amrex::Print() << "macro " << str << " is not defined. \n";
                    break;
                }
            }
        }
        else 
        {
            std::string subscript = str.substr(str.length()-2);
    //        amrex::Print() << "macro_str " << str << " subscript " << subscript << " map_subscript_name: " << map_subscript_name[subscript] << "\n";
            switch (map_subscript_name[subscript])
            {
                case s_Subscript::x :
                {
                    std::string str_without_subscript = str.substr(0,str.length()-2);
                    std::string vector_str = "vec" + str_without_subscript;
                    #ifdef PRINT_HIGH
                    amrex::Print() << prt << "field counter: " << field_counter << "\n";
                    amrex::Print() << prt << "Inside c_Subscript::x \n";
                    amrex::Print() << prt << "Getting x component of vector: "<< vector_str << "\n";
                    #endif

                    m_p_mf[field_counter] = rPost.get_p_array_mf_component(vector_str, 0);
                    break;
                }
                case s_Subscript::y :
                {
                    std::string str_without_subscript = str.substr(0,str.length()-2);
                    std::string vector_str = "vec" + str_without_subscript;
                    #ifdef PRINT_HIGH
                    amrex::Print() << prt << "field counter: " << field_counter << "\n";
                    amrex::Print() << prt << "Inside c_Subscript::y \n";
                    amrex::Print() << prt << "Getting y component of vector: "<< vector_str << "\n";
                    #endif

                    m_p_mf[field_counter] = rPost.get_p_array_mf_component(vector_str, 1);
                    break;

                }
                case s_Subscript::z :
                {
                    std::string str_without_subscript = str.substr(0,str.length()-2);
                    std::string vector_str = "vec" + str_without_subscript;
                    #ifdef PRINT_HIGH
                    amrex::Print() << prt << "field counter: " << field_counter << "\n";
                    amrex::Print() << prt << "Inside c_Subscript::z \n";
                    amrex::Print() << prt << "Getting z component of vector: "<< vector_str << "\n";
                    #endif

                    m_p_mf[field_counter] = rPost.get_p_array_mf_component(vector_str, 2);
                    break;
                }
            }
        }
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_Output::AssimilateDataPointers()************************\n";
#endif
}

int c_Output::Evaluate_TypeOf_MacroStr(std::string macro_str) 
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_Output::Evaluate_TypeOf_MacroStr(*)************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t\t\t";
#endif
    /** if macro_str is defined in c_MacroscopicProperties then return is 0.
     *  if macro_str is defined in c_PostProcessor then return is 1.
     **/
    int return_type = -1;
    auto& rCode = c_Code::GetInstance();
    auto& rPost = rCode.get_PostProcessor();
    auto& rMprop = rCode.get_MacroscopicProperties();

    std::map<std::string,int>::iterator it_Mprop;
    std::map<std::string,s_PostProcessMacroName::macro_name>::iterator it_Post;

    it_Mprop = rMprop.map_param_all.find(macro_str);
    it_Post = rPost.map_macro_name.find(macro_str);

    if(it_Mprop != rMprop.map_param_all.end())
    { 
        #ifdef PRINT_HIGH
        amrex::Print() << prt << "macro_string " << macro_str << " is a part of c_MacroscopicProperties. \n";
        #endif 
        return_type = 0;
    }
    else if(it_Post != rPost.map_macro_name.end())
    { 
        #ifdef PRINT_HIGH
        amrex::Print() << prt << "macro_string " << macro_str << " is a part of c_PostProcessor. \n";
        #endif 
        return_type = 1;
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_Output::Evaluate_TypeOf_MacroStr(*)************************\n";
#endif
    return return_type;
}
