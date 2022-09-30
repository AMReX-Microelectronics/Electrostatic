#include "BoundaryConditions.H"

#include "../Utils/WarpXUtil.H"
#include "Utils/WarpXConst.H"

#include "Code.H"
#include "GeometryProperties.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <ctype.h>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

using namespace amrex;

c_BoundaryConditions::c_BoundaryConditions ()
{
    DefineBoundaryTypeMap();
//    DefineDefaultValueMap();
    ReadData();
} 


c_BoundaryConditions::~c_BoundaryConditions ()
{
//    for (auto& elem : m_p_mf) {
//        elem.release();
//    }
} 
//
//
void 
c_BoundaryConditions::ReadData()
{ 

    ReadBoundaryConditionsType();

    DefineMacroVariableVectorSizes();

//    ReadBoundaryConditionsParser();
    

//    std::map<std::string,amrex::Real>::iterator it_default;
//
//    for (auto it: map_param_all)
//    {
//        amrex::Real default_val;
//
//        it_default = map_default_value.find(it.first);
//
//        if (it_default == map_default_value.end()) {
//            default_val = 0.0;
//        }
//        else {
//            default_val = map_default_value[it.first];
//        }
//        ReadMacroparam(it.first, default_val);
//    }
}
//
//

//bcAny_lo_cast(int c, std::any bcAny_lo) 
//{   
//    switch(map_bcAny_lo[c])
//    {   
//        case inhomogeneous_constant: 
//        {
//            return std::any_cast<std::double>(bcAny_lo[c]);
//        }
//        case inhomogeneous_function:
//        {
//            return std::any_cast<std::string>(bcAny_lo[c]);
//        }
//        case homogeneous:
//        {
//            return std::any_cast<std::double>(bcAny_lo[c]);
//        }
//    }
//}
//
template<class T, class F>
inline std::pair<const std::type_index, std::function<void(std::any const&)>>
    to_any_visitor(F const &f)
{
    return {
        std::type_index(typeid(T)),
        [g = f](std::any const &a)
        {
            if constexpr (std::is_void_v<T>)
                g();
            else
                g(std::any_cast<T const&>(a));
        }
    };
}
 
static std::unordered_map<
    std::type_index, std::function<void(std::any const&)>>
    any_visitor {
        to_any_visitor<float>([](float x){ std::cout << "contains a value: " << x; }),
        to_any_visitor<double>([](double x){ std::cout << "contains a value: " << x; }),
        to_any_visitor<std::string>([](std::string s)
            { std::cout << "contains a function parser name: "<< s; })
        // ... add more handlers for your types ...
    };

inline void process(const std::any& a)
{
    if (const auto it = any_visitor.find(std::type_index(a.type()));
        it != any_visitor.cend()) {
        it->second(a);
    } 
    else {
//        std::cout << "unregistered type! "<< std::quoted(a.type().name());
        std::cout << "contents do not matter! ";
    }
}


void
c_BoundaryConditions::SortBoundaryTypeArrayString(const amrex::Vector<std::string>& bc_str, std::array< std::string, AMREX_SPACEDIM >& bcType, std::array< std::any, AMREX_SPACEDIM >& bcAny,  std::map<int,std::string>& map_bcAny)
{ 
    int c=0;
    for (auto str: bc_str)
    { 

        std::string first_three_letters = str.substr(0,3);
        bcType[c] = first_three_letters;
//        amrex::Print() << "string: " << str  << " first_three_letters: " << first_three_letters<< "\n";

        if(bcType[c] == "dir" or bcType[c] == "neu" or bcType[c] == "rob")
        { 
            std::string fourth_char = str.substr(3,1);
            std::string last_char = str.substr(str.length()-1);
//            amrex::Print() << "fourth_char: " << fourth_char << "\n";
//            amrex::Print() << "last_char: " << last_char << "\n";

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
               
//                amrex::Print() << "bracketed_str: " << bracketed_str << "\n";

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
//                    amrex::Print() << "inhomo constant: " << bracketed_str << "\n";
//                    amrex::Print() << "bcAny_lo: " << std::any_cast<std::double>(bcAny_lo[c]) << "\n";
                }
                else 
                {
                    map_bcAny[c] = "inhomogeneous_function"; 
                    bcAny[c] = bracketed_str;
//                    amrex::Print() << "inhomo function with parameter name: " << bracketed_str << "\n";
//                    amrex::Print() << "bcAny_lo: " << std::any_cast<std::string>(bcAny_lo[c]) << "\n";
                }
            }
            else if(fourth_char != "(" and last_char != ")") 
            {
                map_bcAny[c] = "homogeneous"; 
                bcAny[c] = 0.0;
//                amrex::Print() << "homo constant 0.0 " << "\n";
//                amrex::Print() << "bcAny_lo: " << std::any_cast<std::double>(bcAny_lo[c]) << "\n";
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
}

void 
c_BoundaryConditions::ReadBoundaryConditionsType()
{ 

    amrex::Vector<amrex::Vector<std::string>> bc_str_2d(2);

    amrex::ParmParse pp_boundary("boundary");
    bool bc_lo_specified = pp_boundary.queryarr("lo", bc_str_2d[0]);
    bool bc_hi_specified = pp_boundary.queryarr("hi", bc_str_2d[1]);

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
        
    for (std::size_t i = 0; i < 2; ++i) 
    {
        SortBoundaryTypeArrayString(bc_str_2d[i], bcType_2d[i], bcAny_2d[i], map_bcAny_2d[i]);
    }

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

    /* loop over map_bcAny_2d and determine num_function_parser */
    amrex::Print() << "printing map_bcAny_2d: \n";
    for (auto& i: map_bcAny_2d)
    {
        for (auto& j: i) 
        {
            amrex::Print() << j.first << "  " << j.second << "\n";
        }
        amrex::Print() << "\n";
    }
    amrex::Print() << "\n";

    /* Printing */
    for (std::size_t i = 0; i < 2; ++i) 
    {
        for (std::size_t j = 0; j < AMREX_SPACEDIM; ++j) 
        {
            
            amrex::Print() << "\nboundary type: " << bcType_2d[i][j] << "\n";
            amrex::Print() << "boundary map: " << map_bcAny_2d[i][j] << "\n";
            amrex::Print() << "bracket ";

            process(bcAny_2d[i][j]);

            amrex::Print() << "\n";
        }
    }
    amrex::Print() << "\n\n:";
    
}


void 
c_BoundaryConditions::DefineMacroVariableVectorSizes()
{ 
    //m_macro_type.resize(num_function_parsers, "constant"); //initialized to constant
    //m_macro_value.resize(num_function_parsers);
    //m_macro_str_function.resize(num_function_parsers);
    //m_p_macro_parser.resize(num_function_parsers);
    //m_p_mf.resize(num_function_parsers);
}


//template < class T >
//void 
//c_BoundaryConditions::ReadBoundaryConditionsParser(std::string macro_str, 
//                                        T default_value)
//{
//
//    auto macro_num = map_param_all[macro_str];
//    m_macro_value[macro_num] = default_value;
//
//    ParmParse pp_macroscopic("macroscopic");
//
//    bool specified = false; /** epsilon is the permittivity */
//    std::string macro_functionXYZ = macro_str+"_function(x,y,z)";
//    if (queryWithParser(pp_macroscopic, macro_str.c_str() , m_macro_value[macro_num]) ) {
//        m_macro_type[macro_num] = "constant";
//        specified = true;
//    }
//    if (pp_macroscopic.query( macro_functionXYZ.c_str(), m_macro_str_function[macro_num]) ) {
//        m_macro_type[macro_num] = "parse_" + macro_str + "_function";
//        specified = true;
//    }
//    if (!specified) {
//        std::stringstream warnMsg;
//        warnMsg << "Macroscopic parameter '" << macro_str << "' is not specified in the input file. The default value of " 
//                <<  m_macro_value[macro_num]
//                << " is used.";
//        c_Code::GetInstance().RecordWarning("Macroscopic properties", warnMsg.str());
//    }
//
//    if (m_macro_type[macro_num] == "parse_" + macro_str + "_function") /** initialization of permittivity with a parser */
//    { 
//        Store_parserString(pp_macroscopic, macro_functionXYZ.c_str(),  m_macro_str_function[macro_num]);
//
//        m_p_macro_parser[macro_num] = std::make_unique<amrex::Parser>(
//                                           makeParser( m_macro_str_function[macro_num], {"x","y","z"}));
//    }
//
//}



//void 
//c_BoundaryConditions::InitData()
//{
//
//    const int Ncomp1=1;
//    const int Nghost1=1;
//    const int Nghost0=0;
//
//    auto& rCode = c_Code::GetInstance();
//    auto& rGprop = rCode.get_GeometryProperties();
//    auto& ba = rGprop.ba;
//    auto& dm = rGprop.dm;
//    auto& geom = rGprop.geom;
//
//    for (auto it: map_param_all)
//    {
//        auto str = it.first;
//        DefineAndInitializeMacroparam(str, ba, dm, geom, Ncomp1, map_num_ghostcell[str]);
//    }
//
//    //auto& eps = get_mf("epsilon");
//    //const auto& eps_arr = eps[0].array();
//    //amrex::Print() << "eps_0,0,0 :  " << eps_arr(0,0,0) << "\n";
//    //amrex::Print() << "eps_15,49,49:  " << eps_arr(15,49,49) << "\n";
//    //amrex::Print() << "eps_24,49,49:  " << eps_arr(24,49,49) << "\n";
//    //amrex::Print() << "eps_25,49,49:  " << eps_arr(25,49,49) << "\n";
//    //amrex::Print() << "eps_49,49,49:  " << eps_arr(49,49,49) << "\n";
//    //amrex::Print() << "eps_74,49,49:  " << eps_arr(74,49,49) << "\n";
//    //amrex::Print() << "eps_75,49,49:  " << eps_arr(75,49,49) << "\n";
//    //amrex::Print() << "eps_85,49,49:  " << eps_arr(85,49,49) << "\n";
//    //auto& rho = get_mf("charge_density");
//    //const auto& rho_arr = rho[0].array();
//    //amrex::Print() << "rho_0,0,0 :  " << rho_arr(0,0,0) << "\n";
//    //amrex::Print() << "rho_15,49,49:  " << rho_arr(15,49,49) << "\n";
//    //amrex::Print() << "rho_24,49,49:  " << rho_arr(24,49,49) << "\n";
//    //amrex::Print() << "rho_25,49,49:  " << rho_arr(25,49,49) << "\n";
//    //amrex::Print() << "rho_49,49,49:  " << rho_arr(49,49,49) << "\n";
//    //amrex::Print() << "rho_74,49,49:  " << rho_arr(74,49,49) << "\n";
//    //amrex::Print() << "rho_75,49,49:  " << rho_arr(75,49,49) << "\n";
//    //amrex::Print() << "rho_85,49,49:  " << rho_arr(85,49,49) << "\n";
//}
//
//
//void 
//c_BoundaryConditions::DefineAndInitializeMacroparam(std::string macro_str, 
//                                                       amrex::BoxArray& ba, 
//                                                       amrex::DistributionMapping& dm, 
//                                                       amrex::Geometry& geom, 
//                                                       int Ncomp, 
//                                                       int Nghost)
//{
//
//    auto macro_num = map_param_all[macro_str];
//    //amrex::Print()  << " Initializing macro_str: " << macro_str << " macro_num: " << macro_num << " macro_type: " << m_macro_type[macro_num] << "\n";
//
//    m_p_mf[macro_num] = std::make_unique<amrex::MultiFab>(ba, dm, Ncomp, Nghost); //cell-centered multifab
//
//    if (m_macro_type[macro_num] == "constant") {
//    //    amrex::Print()  << macro_num << " parse function is constant and set to : " << m_macro_value[macro_num] << "\n";
//        m_p_mf[macro_num] -> setVal(m_macro_value[macro_num]);
//
//    } else if (m_macro_type[macro_num] == "parse_" + macro_str + "_function") {
//    //    amrex::Print() << macro_num << " parse function is used with name: " << m_macro_type[macro_num] << "\n";
//
//        InitializeMacroMultiFabUsingParser(m_p_mf[macro_num].get(), m_p_macro_parser[macro_num]->compile<3>(), geom);
//
//    }
//
//}
//
//
//void 
//c_BoundaryConditions::InitializeMacroMultiFabUsingParser (
//                       amrex::MultiFab *macro_mf,
//                       amrex::ParserExecutor<3> const& macro_parser,
//                       amrex::Geometry& geom)
//{
//
//    auto dx = geom.CellSizeArray();
//    amrex::Print() << "dx: " << dx[0] << " " << dx[1] << " " << dx[2] << "\n";
//    auto& real_box = geom.ProbDomain();
//
//    auto iv = macro_mf->ixType().toIntVect();
//
//    amrex::Print() << "iv: " << iv[0] << " " << iv[1] << " " << iv[2] << "\n";
//    amrex::Print() << "real_box_lo: " << real_box.lo(0) << " " << real_box.lo(1) << " " << real_box.lo(2) << "\n";
//
//    for ( amrex::MFIter mfi(*macro_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
//
//        const auto& tb = mfi.tilebox( iv, macro_mf->nGrowVect() ); /** initialize ghost cells in addition to valid cells.
//                                                                       auto = amrex::Box
//                                                                    */
//        auto const& mf_array =  macro_mf->array(mfi); //auto = amrex::Array4<amrex::Real>
//
//        amrex::ParallelFor (tb,
//            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
//
//                amrex::Real fac_x = (1._rt - iv[0]) * dx[0] * 0.5_rt;
//                amrex::Real x = i * dx[0] + real_box.lo(0) + fac_x;
//
//                amrex::Real fac_y = (1._rt - iv[1]) * dx[1] * 0.5_rt;
//                amrex::Real y = j * dx[1] + real_box.lo(1) + fac_y;
//
//                amrex::Real fac_z = (1._rt - iv[2]) * dx[2] * 0.5_rt;
//                amrex::Real z = k * dx[2] + real_box.lo(2) + fac_z;
//
//                mf_array(i,j,k) = macro_parser(x,y,z);
//        });
//    }
//
//}
