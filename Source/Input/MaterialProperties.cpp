#include "MacroscopicProperties.H"

#include "../Utils/WarpXUtil.H"
#include "Code.H"
#include "GeometryProperties.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

using namespace amrex;

c_MacroscopicProperties::c_MacroscopicProperties ()
{
    ReadData();
} 


c_MacroscopicProperties::~c_MacroscopicProperties ()
{
    m_p_epsilon_mf.release();
} 


void 
c_MacroscopicProperties::ReadData()
{    
     ReadPermittivity();
}


void 
c_MacroscopicProperties::InitData()
{
    InitializePermittivity();
}


void 
c_MacroscopicProperties::ReadPermittivity()
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
        warnMsg << "Permittivity (epsilon) is not specified. Using default vacuum value of " 
                <<  m_epsilon 
                << " in the simulation.";
        c_Code::GetInstance().RecordWarning("Macroscopic properties", warnMsg.str());
    }

    /** initialization of permittivity with a parser */
    if (m_epsilon_s == "parse_epsilon_function") {
        Store_parserString(pp_material, "epsilon_function(x,y,z)", m_str_epsilon_function);
        m_p_epsilon_parser = std::make_unique<amrex::Parser>(
                                             makeParser(m_str_epsilon_function,{"x","y","z"}));
    }

}


void 
c_MacroscopicProperties::InitializePermittivity()
{

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;

    m_p_epsilon_mf = std::make_unique<amrex::MultiFab>(ba, dm, Ncomp1, Nghost1); //cell-centered multifab

    if (m_epsilon_s == "constant") {

        m_p_epsilon_mf->setVal(m_epsilon);

    } else if (m_epsilon_s == "parse_epsilon_function") {

        InitializeMacroMultiFabUsingParser(m_p_epsilon_mf.get(), m_p_epsilon_parser->compile<3>());

    }

}


void 
c_MacroscopicProperties::InitializeMacroMultiFabUsingParser (
                       amrex::MultiFab *macro_mf,
                       amrex::ParserExecutor<3> const& macro_parser)
{

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto dx = rGprop.geom.CellSizeArray();
    auto& real_box = rGprop.geom.ProbDomain();

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
