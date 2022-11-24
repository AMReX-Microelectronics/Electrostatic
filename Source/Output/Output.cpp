#include "Output.H"

#include "Code.H"
#include "../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../Utils/CodeUtils/CodeUtil.H"
#include "Input/GeometryProperties/GeometryProperties.H"
#include "Input/MacroscopicProperties/MacroscopicProperties.H"
#include "PostProcessor/PostProcessor.H"

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>


using namespace amrex;


c_Output::c_Output ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Output Constructor()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t\t";
#endif

    DefineSubscriptNameMap();
    m_num_params_plot_single_level = 0;
    m_raw_fields_to_plot = 0;
    m_write_after_init = 0;

    m_num_params_plot_single_level = ReadData();

    if(m_raw_fields_to_plot)  m_extra_dirs.emplace_back("raw_fields");

    m_p_mf.resize(m_map_param_all.size());

#ifdef PRINT_LOW
    amrex::Print() << prt << "number of all parameters: " << m_map_param_all.size() << "\n"; 
    amrex::Print() << prt << "m_raw_fields_to_plot: " << m_raw_fields_to_plot << "\n\n"; 
    amrex::Print() <<"\n" << prt << "number of parameters to plot in a single level plot file: " << m_num_params_plot_single_level << "\n";
    amrex::Print() << prt << "m_write_after_init: " << m_write_after_init << "\n\n"; 
#endif

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

    m_output_option.clear();
    m_map_param_all.clear();
    m_map_field_to_plot_after_init.clear();     

    m_num_params_plot_single_level = 0;
    m_raw_fields_to_plot = 0;
    m_write_after_init = 0;
    m_p_mf.clear();
    m_extra_dirs.clear();
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
//    m_output_option.clear();
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
    amrex::Vector< std::string > extra_str;
    
    int num_params_without_option_2 = 0;

    amrex::ParmParse pp_plot("plot");
    
    foldername_str = "output";
    pp_plot.query("folder_name", foldername_str);

    m_filename_prefix_str = foldername_str + "/plt";

    pp_plot.query("write_after_init", m_write_after_init);

    m_write_interval = 1; //default
    m_rawfield_write_interval = 1; //default
    queryWithParser(pp_plot, "write_interval", m_write_interval);
    queryWithParser(pp_plot, "rawfield_write_interval", m_rawfield_write_interval);

//    amrex::ParmParse pp_plot_file(m_filename_prefix_str);
    bool varnames_specified = pp_plot.queryarr("fields_to_plot", fields_to_plot_withGhost_str);

    const int varsize = fields_to_plot_withGhost_str.size();

    fields_to_plot_str.resize(varsize);
    extra_str.resize(varsize);

    int c=0;
    for (auto str: fields_to_plot_withGhost_str) 
    {
        std::string second_last_str = str.substr(str.length()-2,1);
        if(second_last_str == ".")
        {
            std::string last_str = str.substr(str.length()-1,1);

            if (last_str == "1" or last_str == "2") 
            { 
                std::string string_without_ghost = str.substr(0,str.length()-2);
                extra_str[c] = str.substr(str.length()-1,1);
                fields_to_plot_str[c] = string_without_ghost;
                m_raw_fields_to_plot = true;
            }
            else 
            {
               //assert unknown option.
            }
        }
        else if (second_last_str != ".")
        {
             extra_str[c] = "0";
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
            for(auto subscript : m_map_subscript_name) 
            {
                std::string component = rest_of_string + subscript.first;

                it_map_param_all = m_map_param_all.find(component);

                if (it_map_param_all == m_map_param_all.end()) 
                {
                    if(extra_str[i] == "0")
                    {
                       m_output_option.push_back(0); 
                       ++num_params_without_option_2;
                    }
                    else if(extra_str[i] == "1")
                    {
                       m_output_option.push_back(1); 
                       ++num_params_without_option_2;
                    }
                    else if(extra_str[i] == "2")
                    {
                       m_output_option.push_back(2); 
                    }
                    m_map_param_all[component] = c;
                    ++c;
                }  
            }   
        }
        else 
        {
            it_map_param_all = m_map_param_all.find(str);

            if (it_map_param_all == m_map_param_all.end()) 
            {
                if(extra_str[i] == "0")
                {
                   m_output_option.push_back(0); 
                   ++num_params_without_option_2;
                }
                else if(extra_str[i] == "1")
                {
                   m_output_option.push_back(1); 
                   ++num_params_without_option_2;
                }
                else if(extra_str[i] == "2")
                {
                   m_output_option.push_back(2); 
                }
                m_map_param_all[str] = c;
                ++c;
            }
        }
    }
    fields_to_plot_withGhost_str.clear();
    fields_to_plot_str.clear();
    extra_str.clear();

    amrex::Print() <<  "\n##### OUTPUT #####\n\n";
    amrex::Print() <<  "##### file name: " << m_filename_prefix_str << "\n";
    amrex::Print() <<  "##### fields_to_plot:\n";
    amrex::Print() <<  "##### " << std::setw(20)  << "name" << std::setw(10) << "number" << std::setw(14) << "output_option\n";
    for (auto it: m_map_param_all) 
    {
        amrex::Print() << "##### " << std::setw(20) << it.first << std::setw(10) << it.second << std::setw(10) << m_output_option[it.second] << "\n";
    }
    amrex::Print() <<  "\n##### write_after_init?: " << m_write_after_init << "\n";

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_Output::ReadData()************************\n";
#endif

    return num_params_without_option_2;

}


void 
c_Output::InitData()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t{************************c_Output::InitData()************************\n";
    amrex::Print() << "\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t";
#endif

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;
    int Nghost0=0;

#ifdef AMREX_USE_EB
    if(rGprop.embedded_boundary_flag) {
        m_p_mf_all = std::make_unique<amrex::MultiFab>(ba, dm, m_num_params_plot_single_level, Nghost0, MFInfo(), *rGprop.pEB->p_factory_union); 
    }
#endif
    if(!rGprop.embedded_boundary_flag) {
         m_p_mf_all = std::make_unique<amrex::MultiFab>(ba, dm, m_num_params_plot_single_level, Nghost0);  
    }

    if(m_write_after_init) 
    {
        int counter=0;
        for (auto it : m_map_param_all) 
        {
            std::string field_name = it.first;
            int field_number = it.second;

            if ( Evaluate_TypeOf_MacroStr(field_name) == 0 ) /*it is a part of MacroscopicProperties so the multifab pointer is available after initialization*/
            {
                m_map_field_to_plot_after_init[field_name] = field_number;    

                if(m_output_option[field_number] == 0 or m_output_option[field_number] == 1) 
                {
                    ++counter;
                }
            }
        }
#ifdef AMREX_USE_EB
        if(rGprop.embedded_boundary_flag) {
            m_p_mf_all_init = std::make_unique<amrex::MultiFab>(ba, dm, counter, Nghost0, MFInfo(), *rGprop.pEB->p_factory_union);
        }
#endif
        if(!rGprop.embedded_boundary_flag) {
            m_p_mf_all_init = std::make_unique<amrex::MultiFab>(ba, dm, counter, Nghost0);
        }

#ifdef PRINT_HIGH
        amrex::Print() << "\n" << prt <<  "m_write_after_init is true! \n";
        for (auto it : m_map_field_to_plot_after_init) 
        {
            std::string field_name = it.first;
            int field_number = it.second;
            amrex::Print() << prt  << "field to plot after init, name " << field_name << ", field_number " << field_number << ", output option " << m_output_option[field_number] << "\n";
        }
        amrex::Print() << prt <<  "size of m_p_mf_all_init: " << counter << "\n";
#endif
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t}************************c_Output::InitData()************************\n";
#endif
}


void
c_Output::WriteOutput(int step, amrex::Real time)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t{************************c_Output::WriteOutput()************************\n";
    amrex::Print() << "\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t";
#endif

    if((step+1)%m_write_interval == 0) {
        AssimilateDataPointers();
        m_plot_file_name = amrex::Concatenate(m_filename_prefix_str, step, m_plt_name_digits);

        WriteSingleLevelPlotFile(step, time, m_p_mf_all, m_map_param_all);
    }

    if((m_raw_fields_to_plot == 1) and ((step+1)%m_rawfield_write_interval == 0)) WriteRawFields(m_map_param_all); 

#ifdef PRINT_NAME
    amrex::Print() << "\t\t}************************c_Output::WriteOutput()************************\n";
#endif
}

void
c_Output::WriteOutput_AfterInit()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t{************************c_Output::WriteOutput_AfterInit()************************\n";
    amrex::Print() << "\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t";
#endif

    AssimilateDataPointers();
    int step = 0;
    amrex::Real time = 0;
    m_plot_file_name = amrex::Concatenate(m_filename_prefix_str+"_init", step, 0);
    WriteSingleLevelPlotFile(step, time, m_p_mf_all_init, m_map_field_to_plot_after_init);

    if(m_raw_fields_to_plot) WriteRawFields(m_map_field_to_plot_after_init); 

#ifdef PRINT_NAME
    amrex::Print() << "\t\t}************************c_Output::WriteOutput_AfterInit()************************\n";
#endif
}


void
c_Output::WriteSingleLevelPlotFile(int step, amrex::Real time, std::unique_ptr<amrex::MultiFab>& p_all_mf, std::map<std::string,int>& map_all_mf)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Output::WriteSingleLevelPlotFile()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t\t";
#endif

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& geom = rGprop.geom;
    int Ncomp1=1;
    int Nghost0=0;

    amrex::Vector< std::string > m_p_name_str;

    int c=0;
    for (auto it : map_all_mf) 
    {
        auto field_number = it.second;
        auto field_name = it.first;

        if(m_output_option[field_number] == 0 or m_output_option[field_number] == 1) 
        {
            #ifdef PRINT_HIGH
            amrex::Print() << prt << " copying parameter: " << field_name 
                                  << " counter : " << c 
                                  << " field_number: " << field_number 
                                  << " m_output_option: " << m_output_option[field_number] << "\n";
            #endif

            m_p_name_str.push_back(field_name);

            amrex::MultiFab::Copy(*p_all_mf, *m_p_mf[field_number], 0, c, Ncomp1, Nghost0);
            ++c;
        }
    }

#ifdef AMREX_USE_EB
    if(rGprop.embedded_boundary_flag) {
        amrex::EB_WriteSingleLevelPlotfile( m_plot_file_name, 
                                        *p_all_mf, m_p_name_str, 
                                         geom, 
                                         time, step, 
                                         "HyperCLaw-V1.1", m_default_level_prefix, "Cell",
                                         m_extra_dirs);
   }
#endif
    if(!rGprop.embedded_boundary_flag) {
        amrex::WriteSingleLevelPlotfile( m_plot_file_name, 
                                        *p_all_mf, m_p_name_str, 
                                         geom, 
                                         time, step, 
                                         "HyperCLaw-V1.1", m_default_level_prefix, "Cell",
                                         m_extra_dirs);
    }

    m_p_name_str.clear();

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_Output::WriteSingleLevelPlotFile()************************\n";
#endif
}

 
//void
//c_Output::WriteSingleLevelPlotFile_AfterInit(int step, amrex::Real time)
//{
//#ifdef PRINT_NAME
//    amrex::Print() << "\n\n\t\t\t{************************c_Output::WriteSingleLevelPlotFile_AfterInit(*)************************\n";
//    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
//    std::string prt = "\t\t\t";
//#endif
//
//    auto& rCode = c_Code::GetInstance();
//    auto& rGprop = rCode.get_GeometryProperties();
//    auto& geom = rGprop.geom;
//    int Ncomp1=1;
//    int Nghost0=0;
//
//    amrex::Vector< std::string > m_p_name_str;
//
//    int c=0;
//    for (auto it : m_map_field_to_plot_after_init) 
//    {
//        auto field_number = it.second;
//        auto field_name = it.first;
//
//        if(m_output_option[field_number] == 0 or m_output_option[field_number] == 1) 
//        {
//            m_p_name_str.push_back(field_name);
//
//            #ifdef PRINT_HIGH
//            amrex::Print() << prt << " copying parameter: " << field_name 
//                                  << " counter : " << c 
//                                  << " field_number: " << field_number 
//                                  << " m_output_option: " << m_output_option[field_number] << "\n";
//            #endif
//
//            amrex::MultiFab::Copy(*m_p_mf_all_init, *m_p_mf[field_number], 0, c, Ncomp1, Nghost0);
//            ++c;
//        }
//    }
//
//    amrex::WriteSingleLevelPlotfile( m_plot_file_name, 
//                                    *m_p_mf_all_init, m_p_name_str, 
//                                     geom, 
//                                     time, step, 
//                                     "HyperCLaw-V1.1", m_default_level_prefix, "Cell",
//                                     m_extra_dirs);
//
//    m_p_name_str.clear();
//
//#ifdef PRINT_NAME
//    amrex::Print() << "\t\t\t}************************c_Output::WriteSingleLevelPlotFile_AfterInit(*)************************\n";
//#endif
//}

void
c_Output::WriteRawFields(std::map<std::string,int>& map_all_mf)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Output::WriteRawFields()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t\t";
#endif

    for (auto it : map_all_mf) 
    {
        auto field_number = it.second;
        auto field_name = it.first;

        if(m_output_option[field_number] == 1 or m_output_option[field_number] == 2) 
        {
            int lev = 0;
            const std::string raw_pltname = m_plot_file_name + "/raw_fields";

            std::string prefix = amrex::MultiFabFileFullPrefix( lev, raw_pltname, 
                                                                m_default_level_prefix, field_name );
            VisMF::Write(*m_p_mf[field_number], prefix);
        }
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_Output::WriteRawFields()************************\n";
#endif
}


//void
//c_Output::WriteRawFields_AfterInit()
//{
//#ifdef PRINT_NAME
//    amrex::Print() << "\n\n\t\t\t{************************c_Output::WriteRawFields_AfterInit()************************\n";
//    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
//    std::string prt = "\t\t\t";
//#endif
//
//    for (auto it : m_map_field_to_plot_after_init) 
//    {
//        auto field_number = it.second;
//        auto field_name = it.first;
//
//        if(m_output_option[field_number] == 1 or m_output_option[field_number] == 2) 
//        {
//
//            int lev = 0;
//            const std::string raw_pltname = m_plot_file_name + "/raw_fields";
//
//            std::string prefix = amrex::MultiFabFileFullPrefix( lev, raw_pltname, 
//                                                                m_default_level_prefix, field_name );
//            VisMF::Write(*m_p_mf[field_number], prefix);
//
//            amrex::Print() << "\n" << prt << " writing raw field after init: " << field_name << ", number: " << field_number << " output_option: " << m_output_option[field_number] << ", prefix: " << prefix << "\n";
//        }
//    }
//
//#ifdef PRINT_NAME
//    amrex::Print() << "\t\t\t}************************c_Output::WriteRawFields_AfterInit()************************\n";
//#endif
//}
//

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

    for (auto it : m_map_param_all) 
    {
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
    //        amrex::Print() << "macro_str " << str << " subscript " << subscript << " m_map_subscript_name: " << m_map_subscript_name[subscript] << "\n";
            switch (m_map_subscript_name[subscript])
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

//int c_Output::Evaluate_TypeOf_MacroStr(std::string macro_str) 
//{
//#ifdef PRINT_NAME
//    amrex::Print() << "\n\n\t\t\t\t{************************c_Output::Evaluate_TypeOf_MacroStr(*)************************\n";
//    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
//    std::string prt = "\t\t\t\t";
//#endif
//    /** if macro_str is defined in c_MacroscopicProperties then return is 0.
//     *  if macro_str is defined in c_PostProcessor then return is 1.
//     **/
//    int return_type = -1;
//    auto& rCode = c_Code::GetInstance();
//    auto& rPost = rCode.get_PostProcessor();
//    auto& rMprop = rCode.get_MacroscopicProperties();
//
//    std::map<std::string,int>::iterator it_Mprop;
//    std::map<std::string,s_PostProcessMacroName::macro_name>::iterator it_Post;
//
//    it_Mprop = rMprop.map_param_all.find(macro_str);
//    it_Post = rPost.map_macro_name.find(macro_str);
//
//    if(it_Mprop != rMprop.map_param_all.end())
//    { 
//        #ifdef PRINT_HIGH
//        amrex::Print() << prt << "macro_string " << macro_str << " is a part of c_MacroscopicProperties. \n";
//        #endif 
//        return_type = 0;
//    }
//    else if(it_Post != rPost.map_macro_name.end())
//    { 
//        #ifdef PRINT_HIGH
//        amrex::Print() << prt << "macro_string " << macro_str << " is a part of c_PostProcessor. \n";
//        #endif 
//        return_type = 1;
//    }
//
//#ifdef PRINT_NAME
//    amrex::Print() << "\t\t\t\t}************************c_Output::Evaluate_TypeOf_MacroStr(*)************************\n";
//#endif
//    return return_type;
//}
