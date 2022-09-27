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
    ReadData();
} 


c_MacroscopicProperties::~c_MacroscopicProperties ()
{
//    for (auto& elem : m_p_mf) {
//        elem.release();
//    }
} 


int 
c_MacroscopicProperties::ReadParameterMap()
{ 

   param_map["epsilon"] = 0;
   param_map["charge_density"] = 1;
   param_map["phi"] = 2;

   return param_map.size();

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
c_MacroscopicProperties::ReadData()
{ 

    num_params = ReadParameterMap();
    
    DefineMacroVariableVectorSizes();

    ReadMacroparam("epsilon", PhysConst::ep0);
    ReadMacroparam("charge_density", 0.0);
    ReadMacroparam("phi", 0.0);

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

}


template < class T >
void 
c_MacroscopicProperties::ReadMacroparam(std::string macro_str, 
                                        T default_value)
{

    auto macro_num = param_map[macro_str];
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

    auto macro_num = param_map[macro_str];

    m_p_mf[macro_num] = std::make_unique<amrex::MultiFab>(ba, dm, Ncomp, Nghost); //cell-centered multifab

    if (m_macro_type[macro_num] == "constant") {
        m_p_mf[macro_num] -> setVal(m_macro_value[macro_num]);

    } else if (m_macro_type[macro_num] == "parse_epsilon_function") {

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
    auto& real_box = geom.ProbDomain();

    auto iv = macro_mf->ixType().toIntVect();

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
