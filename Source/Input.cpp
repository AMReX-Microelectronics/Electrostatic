#include "Input.H"

void ReadBasicInput(GpuArray<int, AMREX_SPACEDIM>& n_cell, 
                    int& max_grid_size, 
                    GpuArray<Real, AMREX_SPACEDIM>& prob_lo, 
                    GpuArray<Real, AMREX_SPACEDIM>& prob_hi, 
                    int& nstep, 
                    int& plot_int, 
                    Real& dt) 
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
    }

    nsteps = 10;
    pp.query("nsteps",nsteps);

    plot_int = -1;
    pp.query("plot_int",plot_int);

    pp.get("dt",dt);
}

void ReadMaterialInput(Real& eps_0,
                       Real& eps_r)
{
    ParmParse pp;
    // Material Properties
    pp.get("eps_0",eps_0); // permittivity of vaccum
    pp.get("eps_r",eps_r); // relative permittivity of the material
}

void ReadMLMGInput(Real& alpha,
                   Real& beta,
                   int& set_verbose_param)
{
    ParmParse pp;
    pp.get("alpha",alpha);
    pp.get("beta",beta);
    set_verbose_param = 1;
    pp.query("set_verbose_param",set_verbose_param);
}


void DefineGlobalVariables() 
{
    const Real q = 1.602e-19;
    const Real pi = 3.14159265;
    const Real kb = 1.380649e-23;
}
