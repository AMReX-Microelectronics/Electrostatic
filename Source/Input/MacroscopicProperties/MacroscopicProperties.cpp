#include "MacroscopicProperties.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/CodeUtils/CodeUtil.H"

#include "Code.H"
#include "GeometryProperties.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <ctype.h>

using namespace amrex;

c_MacroscopicProperties::c_MacroscopicProperties ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_MacroscopicProperties Consructor()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    DefineParameterNameMap();
    DefineDefaultValueMap();
    ReadData();

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_MacroscopicProperties Consructor()************************\n";
#endif
} 


c_MacroscopicProperties::~c_MacroscopicProperties ()
{
//    for (auto& elem : m_p_mf) {
//        elem.release();
//    }
} 


void 
c_MacroscopicProperties::ReadData()
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_MacroscopicProperties::ReadData()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    num_params = ReadParameterMapAndNumberOfGhostCells();
    
    DefineMacroVariableVectorSizes();
    std::map<std::string,amrex::Real>::iterator it_default;

    for (auto it: map_param_all)
    {
        amrex::Real default_val;

        it_default = map_default_value.find(it.first);

        if (it_default == map_default_value.end()) {
            default_val = 0.0;
        }
        else {
            default_val = map_default_value[it.first];
        }
        ReadMacroparam(it.first, default_val);
    }
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_MacroscopicProperties::ReadData()************************\n";
#endif
}


int 
c_MacroscopicProperties::ReadParameterMapAndNumberOfGhostCells()
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_MacroscopicProperties::ReadParameterMapAndNumberOfGhostCells()************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t\t\t\t";
#endif


    amrex::Vector< std::string > fields_to_define;
    amrex::Vector< std::string > ghostcells_for_fields;

    amrex::ParmParse pp_macro("macroscopic");

    bool varnames_specified = pp_macro.queryarr("fields_to_define", fields_to_define);
    bool ghostcells_specified = pp_macro.queryarr("ghostcells_for_fields", ghostcells_for_fields);

    std::map<std::string,int>::iterator it_map_param_all;

    int c=0;
    for (auto it: fields_to_define)
    {

        it_map_param_all = map_param_all.find(it);

        if (it_map_param_all == map_param_all.end()) {
                map_param_all[it] = c;
                ++c;
        }
    }
    fields_to_define.clear();

#ifdef PRINT_LOW 
    amrex::Print() <<  "\n" << prt << "map_param_all:\n";
    for (auto it: map_param_all)
    {
        amrex::Print() << prt <<  it.first << "   " << it.second << "\n";
    }
    amrex::Print() << prt << "Total parameters to define (final): " << map_param_all.size() << "\n\n";
#endif

    for (auto it: map_param_all)
    {
        std::string str = it.first;
        bool comparison_true=false;
        std::string str_with_ghost;
        for (auto str_ghost: ghostcells_for_fields) {
             int compare = strncmp(str.c_str(), str_ghost.c_str(), str.length());
             if(compare == 0){
               comparison_true = true; 
               str_with_ghost = str_ghost;
               break;
             }   
        }
        if (comparison_true) {
            std::string num_ghost_str = str_with_ghost.substr(str.length()+1);

            if(std::isdigit(*num_ghost_str.c_str())) {
                amrex::Real num_ghost = std::stod(num_ghost_str);
                map_num_ghostcell[str] = num_ghost;
            } 
            else {
                map_num_ghostcell[str] = 0;
                std::stringstream warnMsg;
                warnMsg << "Macroscopic specification: '" << str_with_ghost  << "' does not contain number of ghost cells.\n"
                << "Default ghost cell value of 0 is chosen\n";
                c_Code::GetInstance().RecordWarning("Macroscopic properties", warnMsg.str());
            }
        }
        else {
            map_num_ghostcell[str] = 0;
            std::stringstream warnMsg;
            warnMsg << "For macroscopic property: '" << str  << "' ghost cell value is not specified.\n"
            << "Default ghost cell value of 0 is chosen\n";
            c_Code::GetInstance().RecordWarning("Macroscopic properties", warnMsg.str());
        }
        
    }
    ghostcells_for_fields.clear();

#ifdef PRINT_LOW 
    amrex::Print() << prt <<  "map_num_ghost_cell: \n";
    for (auto it: map_num_ghostcell)
    {
        amrex::Print() << prt <<  it.first << "   " << it.second << "\n";
    }
#endif
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_MacroscopicProperties::ReadParameterMapAndNumberOfGhostCells()************************\n";
#endif

    return map_param_all.size();

}


void 
c_MacroscopicProperties::DefineMacroVariableVectorSizes()
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_MacroscopicProperties::DefineMacroVariableVectorSizes()************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif
    
    m_macro_type.resize(num_params, "constant"); //initialized to constant
    m_macro_value.resize(num_params);
    m_macro_str_function.resize(num_params);
    m_p_macro_parser.resize(num_params);
    m_p_mf.resize(num_params);

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_MacroscopicProperties::DefineMacroVariableVectorSizes()************************\n";
#endif
}


void 
c_MacroscopicProperties::InitData()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t{************************c_MacroscopicProperties::InitData()************************\n";
    amrex::Print() << "\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    const int Ncomp1=1;

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;
    auto& geom = rGprop.geom;

    for (auto it: map_param_all)
    {
        auto macro_str = it.first;
        auto macro_num = it.second;
        DefineAndInitializeMacroparam(macro_str, macro_num, ba, dm, geom, Ncomp1, map_num_ghostcell[macro_str]);
    }

    //auto& eps = get_mf("epsilon");
    //const auto& eps_arr = eps[0].array();
    //amrex::Print() << "eps_0,0,0 :  " << eps_arr(0,0,0) << "\n";
    //amrex::Print() << "eps_15,49,49:  " << eps_arr(15,49,49) << "\n";
    //amrex::Print() << "eps_24,49,49:  " << eps_arr(24,49,49) << "\n";
    //amrex::Print() << "eps_25,49,49:  " << eps_arr(25,49,49) << "\n";
    //amrex::Print() << "eps_49,49,49:  " << eps_arr(49,49,49) << "\n";
    //amrex::Print() << "eps_74,49,49:  " << eps_arr(74,49,49) << "\n";
    //amrex::Print() << "eps_75,49,49:  " << eps_arr(75,49,49) << "\n";
    //amrex::Print() << "eps_85,49,49:  " << eps_arr(85,49,49) << "\n";
    //auto& rho = get_mf("charge_density");
    //const auto& rho_arr = rho[0].array();
    //amrex::Print() << "rho_0,0,0 :  " << rho_arr(0,0,0) << "\n";
    //amrex::Print() << "rho_15,49,49:  " << rho_arr(15,49,49) << "\n";
    //amrex::Print() << "rho_24,49,49:  " << rho_arr(24,49,49) << "\n";
    //amrex::Print() << "rho_25,49,49:  " << rho_arr(25,49,49) << "\n";
    //amrex::Print() << "rho_49,49,49:  " << rho_arr(49,49,49) << "\n";
    //amrex::Print() << "rho_74,49,49:  " << rho_arr(74,49,49) << "\n";
    //amrex::Print() << "rho_75,49,49:  " << rho_arr(75,49,49) << "\n";
    //amrex::Print() << "rho_85,49,49:  " << rho_arr(85,49,49) << "\n";
#ifdef PRINT_NAME
    amrex::Print() << "\t\t}************************c_MacroscopicProperties::InitData()************************\n";
#endif
}


template < class T >
void 
c_MacroscopicProperties::ReadMacroparam(std::string macro_str, 
                                        T default_value)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_MacroscopicProperties::ReadMacroparam()************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    auto macro_num = map_param_all[macro_str];
    m_macro_value[macro_num] = default_value;

    ParmParse pp_macroscopic("macroscopic");

    bool specified = false; 
    std::string macro_functionXYZ = macro_str+"_function(x,y,z)";
    if (queryWithParser(pp_macroscopic, macro_str.c_str() , m_macro_value[macro_num]) ) {
        m_macro_type[macro_num] = "constant";
        specified = true;
    }
    if (pp_macroscopic.query( macro_functionXYZ.c_str(), m_macro_str_function[macro_num]) ) {
        m_macro_type[macro_num] = "parse_" + macro_str + "_function";
        specified = true;
    }
    if (!specified) {
        std::stringstream warnMsg;
        warnMsg << "Macroscopic parameter '" << macro_str << "' is not specified in the input file. The default value of " 
                <<  m_macro_value[macro_num]
                << " is used.";
        c_Code::GetInstance().RecordWarning("Macroscopic properties", warnMsg.str());
    }

    if (m_macro_type[macro_num] == "parse_" + macro_str + "_function") 
    { 
        Store_parserString(pp_macroscopic, macro_functionXYZ.c_str(),  m_macro_str_function[macro_num]);

        m_p_macro_parser[macro_num] = std::make_unique<amrex::Parser>(
                                           makeParser( m_macro_str_function[macro_num], {"x","y","z"}));
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_MacroscopicProperties::ReadMacroparam()************************\n";
#endif
}


void 
c_MacroscopicProperties::DefineAndInitializeMacroparam(std::string macro_str, 
                                                       int num,
                                                       amrex::BoxArray& ba, 
                                                       amrex::DistributionMapping& dm, 
                                                       amrex::Geometry& geom, 
                                                       int Ncomp, 
                                                       int Nghost)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_MacroscopicProperties::DefineAndInitializeMacroparam()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    auto macro_num = map_param_all[macro_str];
    //amrex::Print()  << " Initializing macro_str: " << macro_str << " macro_num: " << macro_num << " macro_type: " << m_macro_type[macro_num] << "\n";

    m_p_mf[macro_num] = std::make_unique<amrex::MultiFab>(ba, dm, Ncomp, Nghost); //cell-centered multifab

    if (m_macro_type[macro_num] == "constant") {
    //    amrex::Print()  << macro_num << " parse function is constant and set to : " << m_macro_value[macro_num] << "\n";
        m_p_mf[macro_num] -> setVal(m_macro_value[macro_num]);

    } else if (m_macro_type[macro_num] == "parse_" + macro_str + "_function") {
    //    amrex::Print() << macro_num << " parse function is used with name: " << m_macro_type[macro_num] << "\n";

        Multifab_Manipulation::InitializeMacroMultiFabUsingParser_3vars(m_p_mf[macro_num].get(), m_p_macro_parser[macro_num]->compile<3>(), geom);

    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_MacroscopicProperties::DefineAndInitializeMacroparam()************************\n";
#endif
}
