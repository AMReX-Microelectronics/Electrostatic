#include "Output.H"

#include "Code.H"
#include "Input/GeometryProperties/GeometryProperties.H"
#include "Input/MacroscopicProperties/MacroscopicProperties.H"
#include "PostProcessor/PostProcessor.H"

#include <AMReX_ParmParse.H>
#include "AMReX_PlotFileUtil.H"


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

//    m_p_mf_all.release();

} 


int 
c_Output::ReadData()
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_Output::ReadData()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    amrex::Vector< std::string > m_fields_to_plot;
    amrex::ParmParse pp_plot("plot");
    bool varnames_specified = pp_plot.queryarr("fields_to_plot", m_fields_to_plot);

    std::map<std::string,int>::iterator it_map_param_all;


    int c=0;
    for (auto it: m_fields_to_plot) 
    {
    //    amrex::Print() << "reading plot field parameter: " << it  << "\n";
        std::string first_three_char_str = it.substr(0, 3);

        if(first_three_char_str == "vec") 
        {
            std::string rest_of_string = it.substr(3);
            for(auto subscript : map_subscript_name) 
            {
                std::string component = rest_of_string + subscript.first;

                it_map_param_all = map_param_all.find(component);

                if (it_map_param_all == map_param_all.end()) {
                    map_param_all[component] = c;
                    ++c;
                }  
            }   
        }
        else {

            it_map_param_all = map_param_all.find(it);

            if (it_map_param_all == map_param_all.end()) {
                map_param_all[it] = c;
                ++c;
            }
        }
                                           
    }
    m_fields_to_plot.clear();

    //amrex::Print() <<  " map_param_all: \n";
    //for (auto it: map_param_all) 
    //{
    //    amrex::Print() <<  it.first << "   " << it.second << "\n";
    //}
    //amrex::Print() << "total parameters to plot (final): " << map_param_all.size() << "\n\n";

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
#endif

    AssimilateDataPointers();

    plot_file_name = amrex::Concatenate("plt", step, plt_name_digits);

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& geom = rGprop.geom;
    int Ncomp1=1;
    int Nghost0=0;

    int c=0;
    for (auto& it : map_param_all) {
        amrex::MultiFab::Copy(*m_p_mf_all, *m_p_mf[it.second], 0, c, Ncomp1, Nghost0);
        ++c;
    }
    amrex::WriteSingleLevelPlotfile(plot_file_name, *m_p_mf_all, m_p_name_str, geom, time, 0);

#ifdef PRINT_NAME
    amrex::Print() << "\t\t}************************c_Output::WriteSingleLevelPlotFile()************************\n";
#endif
}


template<typename T>
class TD;

void
c_Output::AssimilateDataPointers()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Output::AssimilateDataPointers()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    auto& rCode = c_Code::GetInstance();
    auto& rPost = rCode.get_PostProcessor(); 
    auto& rMprop = rCode.get_MacroscopicProperties();

    int c=0;
    for (auto it : map_param_all) {
        std::string str = it.first;


        std::string second_last_str = str.substr(str.length()-2, 1);
    //     amrex::Print() << "\nmacro_str " << str << " second_last_str " << second_last_str << "\n";

        if( second_last_str != "_" )
        {
    //        amrex::Print() << str << " doesn't have a subscript \n";
            switch ( Evaluate_TypeOf_MacroStr(str) ) 
            {
                case 0:
                {
    //               amrex::Print() << "counter, c: " << c << "\n";
                   m_p_mf[c] = rMprop.get_p_mf(str);
                   break;
                }
                case 1:
                {
    //                amrex::Print() << "counter, c: " << c << "\n";
                    m_p_mf[c] = rPost.get_p_mf(str);
                    break;
                }
                default :
                {
    //                amrex::Print() << "counter, c: " << c << "\n";
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
    //                amrex::Print() << "counter, c: " << c << "\n";
    //                amrex::Print() << "Inside c_Subscript::x \n";
                    std::string str_without_subscript = str.substr(0,str.length()-2);
                    std::string vector_str = "vec" + str_without_subscript;
    //                amrex::Print() << "Getting x component of vector: "<< vector_str << "\n";

                    m_p_mf[c] = rPost.get_p_array_mf_component(vector_str, 0);
                    break;
                }
                case s_Subscript::y :
                {
    //                amrex::Print() << "counter, c: " << c << "\n";
    //                amrex::Print() << "Inside c_Subscript::y \n";
                    std::string str_without_subscript = str.substr(0,str.length()-2);
                    std::string vector_str = "vec" + str_without_subscript;
    //                amrex::Print() << "Getting y component of vector: "<< vector_str << "\n";

                    m_p_mf[c] = rPost.get_p_array_mf_component(vector_str, 1);
                    break;

                }
                case s_Subscript::z :
                {
    //                amrex::Print() << "counter, c: " << c << "\n";
    //                amrex::Print() << "Inside c_Subscript::z \n";
                    std::string str_without_subscript = str.substr(0,str.length()-2);
                    std::string vector_str = "vec" + str_without_subscript;
    //                amrex::Print() << "Getting z component of vector: "<< vector_str << "\n";

                    m_p_mf[c] = rPost.get_p_array_mf_component(vector_str, 2);
                    break;
                }
            }
        }
        ++c;
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
    //   amrex::Print() << "macro_string " << macro_str << " is a part of c_MacroscopicProperties. \n";
        return_type = 0;
    }
    else if(it_Post != rPost.map_macro_name.end())
    { 
    //    amrex::Print() << "macro_string " << macro_str << " is a part of c_PostProcessor. \n";
        return_type = 1;
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_Output::Evaluate_TypeOf_MacroStr(*)************************\n";
#endif
    return return_type;
}
