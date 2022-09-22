#include "Input.H"

void ReadBasicInput(GpuArray<int, AMREX_SPACEDIM>& n_cell, 
                    int& max_grid_size, 
                    GpuArray<Real, AMREX_SPACEDIM>& prob_lo, 
                    GpuArray<Real, AMREX_SPACEDIM>& prob_hi, 
                    int& nstep, 
                    int& plot_int,
                    int& plot_elems)
{
    // ParmParse is way of reading inputs from the inputs file
    // pp.get means we require the inputs file to have it
    // pp.query means we optionally need the inputs file to have it - but we must supply a default here
    ParmParse pp;

    amrex::Vector<int> temp_int(AMREX_SPACEDIM);
    if (pp.queryarr("n_cell",temp_int)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            n_cell[i] = temp_int[i];
        }
    }

    pp.get("max_grid_size",max_grid_size);

    amrex::Vector<Real> temp(AMREX_SPACEDIM);
    if (pp.queryarr("prob_lo",temp)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            prob_lo[i] = temp[i];
        }
    }

    if (pp.queryarr("prob_hi",temp)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            prob_hi[i] = temp[i];
        }
    }

    nstep = 10;
    pp.query("nstep",nstep);

    plot_int = -1;
    pp.query("plot_int",plot_int);

    pp.get("plot_elems",plot_elems);
}

void ReadMaterialInput(Real& eps_0,
                       Real& eps_r)
{
    ParmParse pp;
    std::unique_ptr<amrex::MultiFab> m_eps_mf;
    std::string m_epsilon_s = "constant";
    std::string m_str_eps_function;

    if (pp_macroscopic.query("epsilon_function(x,y,z)", m_str_epsilon_function) ) {
        m_epsilon_s = "parse_epsilon_function";
        m_epsilon_specified = true;
    }
    if (m_epsilon_s == "parse_epsilon_function") {
        Store_parserString(pp_macroscopic, "epsilon_function(x,y,z)", m_str_epsilon_function);
        m_epsilon_parser = std::make_unique<amrex::Parser>(
                                 makeParser(m_str_epsilon_function,{"x","y","z"}));
    }
    // Initialize epsilon
    if (m_epsilon_s == "constant") {

        m_eps_mf->setVal(m_epsilon);

    } else if (m_epsilon_s == "parse_epsilon_function") {

        InitializeMacroMultiFabUsingParser(m_eps_mf.get(), m_epsilon_parser->compile<3>(), lev);

    }
    // Material Properties
    pp.get("eps_0",eps_0); // permittivity of vaccum
    pp.get("eps_r",eps_r); // relative permittivity of the material
}

void ReadMLMGInput(Real& ascalar,
                   Real& bscalar,
                   int& set_verbose,
                   int& max_order)
{
    ParmParse pp;
    pp.get("mlmg_ascalar",ascalar);
    pp.get("mlmg_bscalar",bscalar);
    set_verbose = 0;
    pp.query("mlmg_set_verbose",set_verbose);
    max_order = 2;
    pp.query("mlmg_max_order",max_order);
}

void
MacroscopicProperties::InitializeMacroMultiFabUsingParser (
                       amrex::MultiFab *macro_mf,
                       amrex::ParserExecutor<3> const& macro_parser,
                       const int lev)
{
    WarpX& warpx = WarpX::GetInstance();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx_lev = warpx.Geom(lev).CellSizeArray();
    const amrex::RealBox& real_box = warpx.Geom(lev).ProbDomain();
    amrex::IntVect iv = macro_mf->ixType().toIntVect();
    for ( amrex::MFIter mfi(*macro_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        // Initialize ghost cells in addition to valid cells

        const amrex::Box& tb = mfi.tilebox( iv, macro_mf->nGrowVect());
        amrex::Array4<amrex::Real> const& macro_fab =  macro_mf->array(mfi);
        amrex::ParallelFor (tb,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Shift x, y, z position based on index type
                amrex::Real fac_x = (1._rt - iv[0]) * dx_lev[0] * 0.5_rt;
                amrex::Real x = i * dx_lev[0] + real_box.lo(0) + fac_x;
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                amrex::Real y = 0._rt;
                amrex::Real fac_z = (1._rt - iv[1]) * dx_lev[1] * 0.5_rt;
                amrex::Real z = j * dx_lev[1] + real_box.lo(1) + fac_z;
#else
                amrex::Real fac_y = (1._rt - iv[1]) * dx_lev[1] * 0.5_rt;
                amrex::Real y = j * dx_lev[1] + real_box.lo(1) + fac_y;
                amrex::Real fac_z = (1._rt - iv[2]) * dx_lev[2] * 0.5_rt;
                amrex::Real z = k * dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // initialize the macroparameter
                macro_fab(i,j,k) = macro_parser(x,y,z);
        });

    }
}

