#include "BoundaryConditions.H"

#include "../Utils/WarpXUtil.H"
#include "Utils/WarpXConst.H"

#include "Code.H"
#include "GeometryProperties.H"
#include "../../Utils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <ctype.h>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

using namespace amrex;

c_BoundaryConditions::c_BoundaryConditions ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_BoundaryConditions Constructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    DefineBoundaryTypeMap();
    ReadData();

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_BoundaryConditions Constructor************************\n";
#endif
} 


c_BoundaryConditions::~c_BoundaryConditions ()
{
//    for (auto& elem : m_p_mf) {
//        elem.release();
//    }
}


void 
c_BoundaryConditions::ReadData()
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_BoundaryConditions::ReadData()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    ReadBoundaryConditionsType();

    DefineMacroVariableVectorSizes();

    for (auto it: map_function_parser_name)
    {
        ReadBoundaryConditionsParser(it.first, it.second);
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_BoundaryConditions::ReadData()************************\n";
#endif
}


void
c_BoundaryConditions::SortBoundaryTypeArrayString(const amrex::Vector<std::string>& bc_str, std::array< std::string, AMREX_SPACEDIM >& bcType, std::array< std::any, AMREX_SPACEDIM >& bcAny,  std::map<int,std::string>& map_bcAny)
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t\t{************************c_BoundaryConditions::SortBoundaryTypeArrayString(*)************************\n";
    amrex::Print() << "\t\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    int c=0;
    for (auto str: bc_str)
    { 

        std::string first_three_letters = str.substr(0,3);
        bcType[c] = first_three_letters;
#ifdef PRINT_HIGH
        amrex::Print() << "\nstring: " << str  << " first_three_letters: " << first_three_letters<< "\n";
#endif

        if(bcType[c] == "dir" or bcType[c] == "neu" or bcType[c] == "rob")
        { 
            std::string fourth_char = str.substr(3,1);
            std::string last_char = str.substr(str.length()-1);
#ifdef PRINT_HIGH
            amrex::Print() << "fourth_char: " << fourth_char << "\n";
            amrex::Print() << "last_char: " << last_char << "\n";
#endif

            if(fourth_char == "(" or last_char == ")" )
            {
                if(fourth_char == "(" and last_char == ")" )
                {
                }
                if(fourth_char == "(" and last_char != ")") //throw warning
                {
                    std::stringstream warnMsg;
                    warnMsg << "In the input file, specification of '" << str  << "' is missing a closed bracket\")\".\n"
                    << "Interpreting it as \n" << str.insert(str.length(),")") << "\n";
                    c_Code::GetInstance().RecordWarning("Boundary Conditions", warnMsg.str());

                }
                else if(fourth_char != "(" and last_char == ")") //throw warning
                {
                    std::stringstream warnMsg;
                    warnMsg << "In the input file, specification of '" << str  << "' is missing an open bracket\"(\".\n"
                    << "Interpreting it as \n" << str.insert(3,"(")  << "\n";
                    c_Code::GetInstance().RecordWarning("Boundary Conditions", warnMsg.str());

                }
                std::string bracketed_str = str.substr(4,str.length()-5);
               
#ifdef PRINT_HIGH
                amrex::Print() << "bracketed_str: " << bracketed_str << "\n";
#endif
                std::string stripped_bracketed_str;
                if(bracketed_str.substr(0,1) == "-" or bracketed_str.substr(0,1) == "+") {
                   stripped_bracketed_str = bracketed_str.substr(1);
                }
                else {
                   stripped_bracketed_str = bracketed_str;
                }
                if(std::isdigit( *stripped_bracketed_str.c_str()) ) 
                {   
                    map_bcAny[c] = "inhomogeneous_constant"; 
                    bcAny[c] = std::stod(bracketed_str);
#ifdef PRINT_HIGH
                    amrex::Print() << "inhomo constant: " << bracketed_str << "\n";
#endif
                }
                else 
                {
                    map_bcAny[c] = "inhomogeneous_function"; 
                    bcAny[c] = bracketed_str;
#ifdef PRINT_HIGH
                    amrex::Print() << "inhomo function with parameter name: " << bracketed_str << "\n";
#endif
                }
            }
            else if(fourth_char != "(" and last_char != ")") 
            {
                map_bcAny[c] = "homogeneous"; 
                bcAny[c] = 0.0;
#ifdef PRINT_HIGH
                amrex::Print() << "homo constant 0.0 " << "\n";
#endif
            }
        }
        else if(bcType[c] == "per") {
            map_bcAny[c] = "periodic"; 
        }
        else if(bcType[c] == "ref") {
            map_bcAny[c] = "reflect";
        }
        ++c;
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t\t}************************c_BoundaryConditions::SortBoundaryTypeArrayString(*)************************\n";
#endif
}


void 
c_BoundaryConditions::ReadBoundaryConditionsType()
{ 

#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_BoundaryConditions::ReadBoundaryConditionsType()************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    amrex::Vector<amrex::Vector<std::string>> bc_str_2d(2);

    amrex::ParmParse pp_boundary("boundary");
    bool bc_lo_specified = pp_boundary.queryarr("lo", bc_str_2d[0]);
    bool bc_hi_specified = pp_boundary.queryarr("hi", bc_str_2d[1]);

#ifdef PRINT_LOW
    amrex::Print() << "\nboundary conditions strings, bc_str: \n";
    for (auto& i: bc_str_2d)
    {
        for (auto& j: i) 
        {
            amrex::Print() << j << "  ";
        }
        amrex::Print() << "\n";
    }
    amrex::Print() << "\n";
#endif

        
    for (std::size_t i = 0; i < 2; ++i) 
    {
        SortBoundaryTypeArrayString(bc_str_2d[i], bcType_2d[i], bcAny_2d[i], map_bcAny_2d[i]);
    }

    bc_str_2d.clear();

    /* WARNING: assert both sides to be periodic */
    for (std::size_t i = 0; i < 1; ++i) 
    {
        for (std::size_t j = 0; j < AMREX_SPACEDIM; ++j) 
        {
            if(map_bcAny_2d[i][j] == "periodic" or map_bcAny_2d[i+1][j] == "periodic") {
               if(map_bcAny_2d[i][j] != "periodic" or map_bcAny_2d[i+1][j] !="periodic") { //raise a warning

                    std::string dim;
                    if(j==0) { dim= 'x';} else if(j==1) {dim='y';} else if(j==2) {dim = 'z';}

                    std::stringstream warnMsg;
                    warnMsg << "Note that only one of the ends of a domain along a direction cannot be periodic; both must be. \n"
                    << "Input specification for direction '" << dim << "' violates that, and therefore, both sides are assumed to be periodic. \n";
                    c_Code::GetInstance().RecordWarning("Boundary Conditions", warnMsg.str());

                    map_bcAny_2d[i][j] = "periodic";    map_bcAny_2d[i+1][j] = "periodic";
                    bcType_2d[i][j] = "per";            bcType_2d[i+1][j] = "per";
                    bcAny_2d[i][j] = 0x0;               bcAny_2d[i+1][j] = 0x0;
               }
            }
        }
    }

#ifdef PRINT_LOW
    for (std::size_t i = 0; i < 2; ++i) 
    {
        for (std::size_t j = 0; j < AMREX_SPACEDIM; ++j) 
        {
            
            amrex::Print() << "\nboundary type, bcType_2d: " << bcType_2d[i][j] << "\n";
            amrex::Print() << "boundary map, map_bcAny_2d: " << map_bcAny_2d[i][j] << "\n";
            amrex::Print() << "(bcAny_2d) bracket ";

            process_std_any(bcAny_2d[i][j]);

            amrex::Print() << "\n";
        }
    }
    amrex::Print() << "\n";
#endif


    /* loop over map_bcAny_2d and fill in the map of function parser names and set number of function parser names. */
    int c=0;
    for (std::size_t i = 0; i < 2; ++i) 
    {
        for (std::size_t j = 0; j < AMREX_SPACEDIM; ++j) 
        {
            if(map_bcAny_2d[i][j] == "inhomogeneous_function") 
            {
                if( bcType_2d[i][j] == "rob") 
                {
                    std::string main_str = std::any_cast<std::string>(bcAny_2d[i][j]);

                    for(auto& str: robin_substrings) 
                    {
                        std::string str_with_subscript = main_str + str;
                        map_function_parser_name[ str_with_subscript ] = c; 
                        ++c;
                    }
                } 
                else 
                {
                    map_function_parser_name[ std::any_cast<std::string>(bcAny_2d[i][j]) ] = c; 
                    ++c;
                }
            }
        }
    }
    num_function_parsers = map_function_parser_name.size();

#ifdef PRINT_LOW
    amrex::Print() << "\nmap_function_parser_name: \n";
    for (auto& i: map_function_parser_name)
    {
        amrex::Print() << i.first << "  " << i.second << "\n";
    }
    amrex::Print() << "\nnum_function_parsers: " << num_function_parsers << "\n";
    amrex::Print() << "\n\n";
#endif    

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_BoundaryConditions::ReadBoundaryConditionsType()************************\n\n";
#endif    
}


void 
c_BoundaryConditions::DefineMacroVariableVectorSizes()
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_BoundaryConditions::DefineMacroVariableVectorSizes()************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif


    m_macro_str_function.resize(num_function_parsers);
    m_p_macro_parser.resize(num_function_parsers);
    m_p_mf.resize(num_function_parsers);



#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_BoundaryConditions::DefineMacroVariableVectorSizes()************************\n";
#endif
}


void 
c_BoundaryConditions::ReadBoundaryConditionsParser(std::string macro_str, int macro_num)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_BoundaryConditions::ReadBoundaryConditionsParser(*)************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif


    ParmParse pp_boundary("boundary");

    std::string macro_functionXYZ = macro_str + "_function(x,y,z)";
    bool specified = false;

    if (pp_boundary.query( macro_functionXYZ.c_str(), m_macro_str_function[macro_num]) ) {
        specified = true;
    }
    if(!specified) {
        std::string warnMsg = "Boundary Conditions: function parser '" + macro_functionXYZ + "' is not specified in the input file.\n";
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(specified==true, warnMsg);
    } 
    else 
    {
        Store_parserString(pp_boundary, macro_functionXYZ.c_str(),  m_macro_str_function[macro_num]);

        m_p_macro_parser[macro_num] = std::make_unique<amrex::Parser>(
                                           makeParser( m_macro_str_function[macro_num], m_parser_varname_vector));
    }

#ifdef PRINT_LOW
    if(specified) {
       amrex::Print() << "\nReading macro: " << macro_str << " with number: " << macro_num << "\n";
       amrex::Print() << "function_parser_name: " << macro_functionXYZ << "\n";
       amrex::Print() << "function_parser: " <<  m_macro_str_function[macro_num] << "\n\n";
    }
#endif

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_BoundaryConditions::ReadBoundaryConditionsParser(*)************************\n\n";
#endif    
}



void 
c_BoundaryConditions::InitData()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t{************************c_BoundaryConditions::InitData()************************\n";
    amrex::Print() << "\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    const int Ncomp1=1;
    const int Nghost1=1;

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;
    auto& geom = rGprop.geom;

    for (auto it: map_function_parser_name)
    {
        auto macro_str = it.first;
        auto macro_num = it.second;
        int parser_varname_size = 3;

        m_p_mf[macro_num] = std::make_unique<amrex::MultiFab>(ba, dm, Ncomp1, Nghost1); //cell-centered multifab

        Multifab_Manipulation::InitializeMacroMultiFabUsingParser_3vars(m_p_mf[macro_num].get(), m_p_macro_parser[macro_num]->compile<3>(), geom);
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t}************************c_BoundaryConditions::InitData()************************\n\n";
#endif    
}
