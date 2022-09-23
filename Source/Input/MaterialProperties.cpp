#include "MaterialProperties.H"

#include "../Utils/WarpXUtil.H"
#include "Code.H"
#include "GeometryProperties.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

using namespace amrex;

c_MaterialProperties::c_MaterialProperties ()
{
    ReadData();
} 


c_MaterialProperties::~c_MaterialProperties ()
{
    m_p_epsilon_mf.release();
} 


void 
c_MaterialProperties::ReadData()
{    
     ReadPermittivity();
}


void 
c_MaterialProperties::InitData()
{
    InitializePermittivity();
}


void 
c_MaterialProperties::ReadPermittivity()
{

    ParmParse pp_material("material");
    /** permittivity (epsilon) */
    bool epsilon_specified = false;
    if (queryWithParser(pp_material, "epsilon", m_epsilon)) {
        m_epsilon_s = "constant";
        epsilon_specified = true;
    }
    if (pp_material.query("epsilon_function(x,y,z)", m_str_epsilon_function) ) {
        m_epsilon_s = "parse_epsilon_function";
        epsilon_specified = true;
    }
    if (!epsilon_specified) {
        std::stringstream warnMsg;
        Print() << "Oops " << m_epsilon << "is used\n";
        warnMsg << "Material permittivity is not specified. Using default vacuum value of " 
                <<  m_epsilon 
                << " in the simulation.";
        c_Code::GetInstance().RecordWarning("Material properties", warnMsg.str());
    }

    /** initialization of permittivity with a parser */
    if (m_epsilon_s == "parse_epsilon_function") {
        Store_parserString(pp_material, "epsilon_function(x,y,z)", m_str_epsilon_function);
        m_p_epsilon_parser = std::make_unique<amrex::Parser>(
                                             makeParser(m_str_epsilon_function,{"x","y","z"}));
    }

}


void 
c_MaterialProperties::InitializePermittivity()
{

    auto& pCode = c_Code::GetInstance();
    auto& pGeom = pCode.get_GeometryProperties();
    auto ba = pGeom.ba;
    auto dm = pGeom.dm;

    m_p_epsilon_mf = std::make_unique<amrex::MultiFab>(ba, dm, Ncomp1, Nghost1); //cell-centered multifab

    if (m_epsilon_s == "constant") {

        m_p_epsilon_mf->setVal(m_epsilon);

    } else if (m_epsilon_s == "parse_epsilon_function") {

        InitializeMacroMultiFabUsingParser(m_p_epsilon_mf.get(), m_p_epsilon_parser->compile<3>(), level0);

    }
    m_p_epsilon_parser.release();

}


void 
c_MaterialProperties::InitializeMacroMultiFabUsingParser (
                       amrex::MultiFab *macro_mf,
                       amrex::ParserExecutor<3> const& macro_parser,
                       const int lev)
{
//    WarpX& warpx = WarpX::GetInstance();
//    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx_lev = warpx.Geom(lev).CellSizeArray();
//    const amrex::RealBox& real_box = warpx.Geom(lev).ProbDomain();
//    amrex::IntVect iv = macro_mf->ixType().toIntVect();
//    for ( amrex::MFIter mfi(*macro_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
//        // Initialize ghost cells in addition to valid cells
//
//        const amrex::Box& tb = mfi.tilebox( iv, macro_mf->nGrowVect());
//        amrex::Array4<amrex::Real> const& macro_fab =  macro_mf->array(mfi);
//        amrex::ParallelFor (tb,
//            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
//                // Shift x, y, z position based on index type
//                amrex::Real fac_x = (1._rt - iv[0]) * dx_lev[0] * 0.5_rt;
//                amrex::Real x = i * dx_lev[0] + real_box.lo(0) + fac_x;
//#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
//                amrex::Real y = 0._rt;
//                amrex::Real fac_z = (1._rt - iv[1]) * dx_lev[1] * 0.5_rt;
//                amrex::Real z = j * dx_lev[1] + real_box.lo(1) + fac_z;
//#else
//                amrex::Real fac_y = (1._rt - iv[1]) * dx_lev[1] * 0.5_rt;
//                amrex::Real y = j * dx_lev[1] + real_box.lo(1) + fac_y;
//                amrex::Real fac_z = (1._rt - iv[2]) * dx_lev[2] * 0.5_rt;
//                amrex::Real z = k * dx_lev[2] + real_box.lo(2) + fac_z;
//#endif
//                // initialize the macroparameter
//                macro_fab(i,j,k) = macro_parser(x,y,z);
//        });
//
//    }
}
