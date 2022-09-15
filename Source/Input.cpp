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
