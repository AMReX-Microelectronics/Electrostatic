#include "MacroscopicProperties.H"

#include "../Utils/WarpXUtil.H"
#include "Utils/WarpXConst.H"

#include "Code.H"
#include "GeometryProperties.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

using namespace amrex;

c_MacroscopicProperties::c_MacroscopicProperties ()
{
    DefineParameterNameMap();
    DefineDefaultValueMap();
    ReadData();
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

    num_params = ReadParameterMap();
    
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
}


int 
c_MacroscopicProperties::ReadParameterMap()
{ 

    amrex::Vector< std::string > m_fields_to_define;

    amrex::ParmParse pp_macro("macroscopic");

    bool varnames_specified = pp_macro.queryarr("fields_to_define", m_fields_to_define);

    std::map<std::string,int>::iterator it_map_param_all;

    int c=0;
    for (auto it: m_fields_to_define)
    {
    //    amrex::Print() << "reading field to define " << it  << "\n";

        it_map_param_all = map_param_all.find(it);

        if (it_map_param_all == map_param_all.end()) {
                map_param_all[it] = c;
                ++c;
        }
    }
    m_fields_to_define.clear();

    //amrex::Print() <<  " map_param_all: \n";
    //for (auto it: map_param_all)
    //{
    //    amrex::Print() <<  it.first << "   " << it.second << "\n";
    //}
    //amrex::Print() << "total parameters to define (final): " << map_param_all.size() << "\n\n";

    return map_param_all.size();

}


void 
c_MacroscopicProperties::DefineMacroVariableVectorSizes()
{ 
    m_macro_type.resize(num_params, "constant"); //initialized to constant
    m_macro_value.resize(num_params);
    m_macro_str_function.resize(num_params);
    m_p_macro_parser.resize(num_params);
    m_p_mf.resize(num_params);
}


void 
c_MacroscopicProperties::InitData()
{

    const int Ncomp1=1;
    const int Nghost1=1;
    const int Nghost0=0;

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;
    auto& geom = rGprop.geom;

    DefineAndInitializeMacroparam("epsilon", ba, dm, geom, Ncomp1, Nghost1);
    DefineAndInitializeMacroparam("charge_density", ba, dm, geom, Ncomp1, Nghost0);
    DefineAndInitializeMacroparam("phi", ba, dm, geom, Ncomp1, Nghost1);

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
}


template < class T >
void 
c_MacroscopicProperties::ReadMacroparam(std::string macro_str, 
                                        T default_value)
{

    auto macro_num = map_param_all[macro_str];
    m_macro_value[macro_num] = default_value;

    ParmParse pp_macroscopic("macroscopic");

    bool specified = false; /** epsilon is the permittivity */
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

    if (m_macro_type[macro_num] == "parse_" + macro_str + "_function") /** initialization of permittivity with a parser */
    { 
        Store_parserString(pp_macroscopic, macro_functionXYZ.c_str(),  m_macro_str_function[macro_num]);

        m_p_macro_parser[macro_num] = std::make_unique<amrex::Parser>(
                                           makeParser( m_macro_str_function[macro_num], {"x","y","z"}));
    }

}


void 
c_MacroscopicProperties::DefineAndInitializeMacroparam(std::string macro_str, 
                                                       amrex::BoxArray& ba, 
                                                       amrex::DistributionMapping& dm, 
                                                       amrex::Geometry& geom, 
                                                       int Ncomp, 
                                                       int Nghost)
{

    auto macro_num = map_param_all[macro_str];
    //amrex::Print()  << " Initializing macro_str: " << macro_str << " macro_num: " << macro_num << " macro_type: " << m_macro_type[macro_num] << "\n";

    m_p_mf[macro_num] = std::make_unique<amrex::MultiFab>(ba, dm, Ncomp, Nghost); //cell-centered multifab

    if (m_macro_type[macro_num] == "constant") {
    //    amrex::Print()  << macro_num << " parse function is constant and set to : " << m_macro_value[macro_num] << "\n";
        m_p_mf[macro_num] -> setVal(m_macro_value[macro_num]);

    } else if (m_macro_type[macro_num] == "parse_" + macro_str + "_function") {
    //    amrex::Print() << macro_num << " parse function is used with name: " << m_macro_type[macro_num] << "\n";

        InitializeMacroMultiFabUsingParser(m_p_mf[macro_num].get(), m_p_macro_parser[macro_num]->compile<3>(), geom);

    }

}


void 
c_MacroscopicProperties::InitializeMacroMultiFabUsingParser (
                       amrex::MultiFab *macro_mf,
                       amrex::ParserExecutor<3> const& macro_parser,
                       amrex::Geometry& geom)
{

    auto dx = geom.CellSizeArray();
    amrex::Print() << "dx: " << dx[0] << " " << dx[1] << " " << dx[2] << "\n";
    auto& real_box = geom.ProbDomain();

    auto iv = macro_mf->ixType().toIntVect();

    amrex::Print() << "iv: " << iv[0] << " " << iv[1] << " " << iv[2] << "\n";
    amrex::Print() << "real_box_lo: " << real_box.lo(0) << " " << real_box.lo(1) << " " << real_box.lo(2) << "\n";

    for ( amrex::MFIter mfi(*macro_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const auto& tb = mfi.tilebox( iv, macro_mf->nGrowVect() ); /** initialize ghost cells in addition to valid cells.
                                                                       auto = amrex::Box
                                                                    */
        auto const& mf_array =  macro_mf->array(mfi); //auto = amrex::Array4<amrex::Real>

        amrex::ParallelFor (tb,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                amrex::Real fac_x = (1._rt - iv[0]) * dx[0] * 0.5_rt;
                amrex::Real x = i * dx[0] + real_box.lo(0) + fac_x;

                amrex::Real fac_y = (1._rt - iv[1]) * dx[1] * 0.5_rt;
                amrex::Real y = j * dx[1] + real_box.lo(1) + fac_y;

                amrex::Real fac_z = (1._rt - iv[2]) * dx[2] * 0.5_rt;
                amrex::Real z = k * dx[2] + real_box.lo(2) + fac_z;

                mf_array(i,j,k) = macro_parser(x,y,z);
        });
    }

}
