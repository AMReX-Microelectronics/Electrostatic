#include "Output.H"

#include "Code.H"
#include "Input/GeometryProperties.H"
#include "Input/MacroscopicProperties.H"
#include "PostProcessor/PostProcessor.H"


#include "AMReX_PlotFileUtil.H"


using namespace amrex;


c_Output::c_Output ()
{

    DefineSubscriptNameMap();
    num_all_params = ReadData();

    m_p_mf.resize(num_all_params);
    m_p_name_str.resize(num_all_params);

} 


c_Output::~c_Output ()
{

//    m_p_mf_all.release();

} 


int 
c_Output::ReadData()
{ 

    map_param_all["E_x"] = 0;
    map_param_all["E_y"] = 1;
    map_param_all["E_z"] = 2;
    map_param_all["Flux_x"] = 3;
    map_param_all["Flux_y"] = 4;
    map_param_all["Flux_z"] = 5;
    map_param_all["epsilon"] = 6;
    map_param_all["charge_density"] = 7;
    map_param_all["phi"] = 8;
//        std::string first_three_char_str = str.substr(0, 3);

    return map_param_all.size();

}


void 
c_Output::InitData()
{
    for (auto it : map_param_all) {
        m_p_name_str[it.second] = it.first;
    }

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;
    int Nghost0=0;

    m_p_mf_all = std::make_unique<amrex::MultiFab>(ba, dm, num_all_params, Nghost0); //cell-centered multifab

}


void
c_Output::WriteSingleLevelPlotFile(int step, amrex::Real time)
{

    AssimilateDataPointers();

    plot_file_name = amrex::Concatenate("plt", step, plt_name_digits);

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& geom = rGprop.geom;
    int Ncomp1=1;
    int Nghost0=0;

    for (auto& it : map_param_all) {
        amrex::MultiFab::Copy(*m_p_mf_all, *m_p_mf[it.second], 0, 0, Ncomp1, Nghost0);
    }
    amrex::WriteSingleLevelPlotfile(plot_file_name, *m_p_mf_all, m_p_name_str, geom, time, 0);

}


template<typename T>
class TD;

void
c_Output::AssimilateDataPointers()
{

    auto& rCode = c_Code::GetInstance();
    auto& rPost = rCode.get_PostProcessor(); 
    auto& rMprop = rCode.get_MacroscopicProperties();

    int c=0;
    for (auto it : map_param_all) {
        std::string str = it.first;


        std::string second_last_str = str.substr(str.length()-2, 1);
         amrex::Print() << "\nmacro_str " << str << " second_last_str " << second_last_str << "\n";

        if( second_last_str != "_" )
        {
            amrex::Print() << str << " doesn't have a subscript \n";
            switch ( Evaluate_TypeOf_MacroStr(str) ) 
            {
                case 0:
                {
                   m_p_mf[c] = rMprop.get_p_mf(str);
                   break;
                }
                case 1:
                {
                    m_p_mf[c] = rPost.get_p_mf(str);
                    break;
                }
                default :
                {
                    amrex::Print() << "macro " << str << " is not defined. \n";
                    break;
                }
            }
        }
        else 
        {
            std::string subscript = str.substr(str.length()-2);
            amrex::Print() << "macro_str " << str << " subscript " << subscript << " map_subscript_name: " << map_subscript_name[subscript] << "\n";
            switch (map_subscript_name[subscript])
            {
                case s_Subscript::x :
                {
                    amrex::Print() << "Inside c_Subscript::x \n";
                    std::string str_without_subscript = str.substr(0,str.length()-2);
                    std::string vector_str = "vec" + str_without_subscript;

                    m_p_mf[c] = rPost.get_p_array_mf_component(vector_str, 0);
                    break;
                }
                case s_Subscript::y :
                {
                    amrex::Print() << "Inside c_Subscript::y \n";
                    std::string str_without_subscript = str.substr(0,str.length()-2);
                    std::string vector_str = "vec" + str_without_subscript;

                    m_p_mf[c] = rPost.get_p_array_mf_component(vector_str, 1);
                    break;

                }
                case s_Subscript::z :
                {
                    amrex::Print() << "Inside c_Subscript::z \n";
                    std::string str_without_subscript = str.substr(0,str.length()-2);
                    std::string vector_str = "vec" + str_without_subscript;

                    m_p_mf[c] = rPost.get_p_array_mf_component(vector_str, 2);
                    break;
                }
            }
        }
        ++c;
    }

}

int c_Output::Evaluate_TypeOf_MacroStr(std::string macro_str) 
{
    /** if macro_str is defined in c_MacroscopicProperties then return is 0.
     *  if macro_str is defined in c_PostProcessor then return is 1.
     **/
    int return_type = -1;
    auto& rCode = c_Code::GetInstance();
    auto& rPost = rCode.get_PostProcessor();
    auto& rMprop = rCode.get_MacroscopicProperties();

    std::map<std::string,s_MacroscopicPropertiesMacroName::macro_name>::iterator it_Mprop;
    std::map<std::string,s_PostProcessMacroName::macro_name>::iterator it_Post;

    it_Mprop = rMprop.map_macro_name.find(macro_str);
    it_Post = rPost.map_macro_name.find(macro_str);

    if(it_Mprop != rMprop.map_macro_name.end())
    { 
        amrex::Print() << "macro_string " << macro_str << " is a part of c_MacroscopicProperties. \n";
        return_type = 0;
    }
    else if(it_Post != rPost.map_macro_name.end())
    { 
        amrex::Print() << "macro_string " << macro_str << " is a part of c_PostProcessor. \n";
        return_type = 1;
    }

    return return_type;
}
