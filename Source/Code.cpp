#include "Code.H"
#include "Utils/WarpXUtil.H"
#include "Input/MaterialProperties.H"


c_Code::c_Code ()
{

    ReadInputParameters();

}


c_Code::~c_Code ()
{

 //

}


void 
c_Code::ReadInputParameters ()
{

     
     ParseGeometryInput();
     m_p_MaterialProperties = std::make_unique<c_MaterialProperties>();

}


void 
c_Code::InitializeData ()
{
 
    m_p_MaterialProperties->InitData();

}

void 
c_Code::ParseGeometryInput() 
{

    amrex::Vector<int> n_cell;
    amrex::Vector<amrex::Real> prob_lo(AMREX_SPACEDIM);
    amrex::Vector<amrex::Real> prob_hi(AMREX_SPACEDIM);
    amrex::Vector<int> mg;
    amrex::Vector<int> bf;

    amrex::ParmParse pp_domain("domain");

    getArrWithParser(pp_domain, "prob_lo", prob_lo, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_lo.size() == AMREX_SPACEDIM);

    getArrWithParser(pp_domain, "prob_hi", prob_hi, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_hi.size() == AMREX_SPACEDIM);

    pp_domain.getarr("n_cell", n_cell, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(n_cell.size() == AMREX_SPACEDIM);

    pp_domain.queryarr("max_grid_size", mg);

    pp_domain.queryarr("blocking_factor", bf);
    bf.resize(std::max(static_cast<int>(bf.size()),1),8);

    pp_domain.addarr("m_cell", n_cell);
    pp_domain.addarr("prob_lo", prob_lo);
    pp_domain.addarr("prob_hi", prob_hi);
    pp_domain.addarr("blocking_factor", bf);
    pp_domain.addarr("blocking_factor", mg);

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        domain.n_cell[i] = n_cell[i];
        domain.prob_lo[i] = prob_lo[i];
        domain.prob_hi[i] = prob_hi[i];
    }
    domain.max_grid_size = mg;
    domain.blocking_factor = bf;

}
