#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H> 
#include <AMReX_MultiFab.H> 
#include <AMReX_VisMF.H>
#include "Global.H"
#include "myfunc.H"
#include "Input.H"
#include "Initialize.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
    return 0;
}

void main_main ()
{

    Real total_step_strt_time = ParallelDescriptor::second();

    //*** SIMULATION INPUT START

    //Basic simulation input
    amrex::GpuArray<int, AMREX_SPACEDIM> n_cell;     // number of cells in each dimension
    int max_grid_size; // size of maximum box dimension used to divide the domain
    amrex::GpuArray<Real, AMREX_SPACEDIM> prob_lo; // physical lo coordinate
    amrex::GpuArray<Real, AMREX_SPACEDIM> prob_hi; // physical hi coordinate
    int nsteps; // total steps in simulation
    int plot_int; // plot file writing interval
    int plot_elems;
    ReadBasicInput(n_cell, max_grid_size, prob_lo, prob_hi, nsteps, plot_int, plot_elems);
    Print() << "prob_lo: " << prob_lo[0] << " " << prob_lo[1] << " " << prob_lo[2] << "\n";
    Print() << "prob_hi: " << prob_hi[0] << " " << prob_hi[1] << " " << prob_hi[2] << "\n";
    
    //Material input
    Real eps_0, eps_r;
    ReadMaterialInput(eps_0, eps_r); 

    //Material input
    Real mlmg_ascalar, mlmg_bscalar;
    int mlmg_set_verbose;
    int mlmg_max_order;
    ReadMLMGInput(mlmg_ascalar, mlmg_bscalar, mlmg_set_verbose, mlmg_max_order); 

    //*** SIMULATION INPUT END


    //*** DOMAIN SETUP START

    // create ba: BoxArray and geom: Geometry
    BoxArray ba; // a list of boxes that cover the domain
    Geometry geom; // contains info such as the physical domain size, number of points in the domain, and periodicity

    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0)); // domain low 
    IntVect dom_hi(AMREX_D_DECL(n_cell[0]-1, n_cell[1]-1, n_cell[2]-1)); // domain high

    Box domain(dom_lo, dom_hi); // Make a single box that is the entire domain

    ba.define(domain); // initialize the boxarray 'ba' from the single box 'domain'

    ba.maxSize(max_grid_size); // break up ba into chunks no larger than 'max_grid_size' along a direction

    RealBox real_box({AMREX_D_DECL( prob_lo[0], prob_lo[1], prob_lo[2])}, 
                     {AMREX_D_DECL( prob_hi[0], prob_hi[1], prob_hi[2])});  //physical domain 

    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)}; // 0: not periodic, 1: periodic

    geom.define(domain, real_box, CoordSys::cartesian, is_periodic); //define the geom object

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray(); // obtain cell size array dx from the geom object
    Print() << "dx: " << dx[0] << " " << dx[1] << " " << dx[2] << "\n";


    DistributionMapping dm(ba); //distribute boxes over MPI processors

    //*** DOMAIN SETUP END


    //*** MULTIFAB ARRAYS START

    int Nghost1 = 1; //number of ghost cells for each array
    int Nghost0 = 0; //number of ghost cells for each array
    int Ncomp1 = 1; //number of components for each array

    MultiFab Permittivity(ba, dm, Ncomp1, Nghost1);
    MultiFab PoissonPhi(ba, dm, Ncomp1, Nghost1);
    MultiFab rho(ba, dm, Ncomp1, Nghost0); //charge density
    MultiFab eps_cc(ba, dm, Ncomp1, Nghost1); //permittivity at cell-centers
    MultiFab Plt(ba, dm, plot_elems, Nghost0); //4: number of elements to print

    //*** MULTIFAB ARRAYS END


    //*** POISSON SOLVER SETUP SETUP AND RUN START

    // Multi-Level Multi-Grid solver is used which solves for phi in (A*alpha_cc - B * div (beta_face grad)) phi = rho
    // For Poisson solver, A and alpha_cc are set to zero.
    // see AMReX_MLLinOp.H for the defaults, accessors, and mutators
    LPInfo info;
    MLABecLaplacian mlabec({geom}, {ba}, {dm}, info); // Implicit solve using MLABecLaplacian class

    // Force singular system to be solvable
    mlabec.setEnforceSingularSolvable(false); 

    // set order of stencil
    mlabec.setMaxOrder(mlmg_max_order);  

    // build array of boundary conditions needed by MLABecLaplacian
    // see Src/Boundary/AMReX_LO_BCTYPES.H for supported types
    std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_bc_lo;
    std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_bc_hi; 

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          mlmg_bc_lo[idim] = mlmg_bc_hi[idim] = LinOpBCType::Dirichlet;
    } 

    // tell the solver what the domain boundary conditions are
    mlabec.setDomainBC(mlmg_bc_lo,mlmg_bc_hi);

    // initialize phi to zero including the ghost cells
    PoissonPhi.setVal(0.); 
    // set Dirichlet BC by reading in the ghost cell values
    int amrlev = 0; //refers to the setcoarsest level of the solve
    mlabec.setLevelBC(amrlev, &PoissonPhi);

    // set scaling factors 
    mlabec.setScalars(mlmg_ascalar, mlmg_bscalar); 

    // set alpha_cc, and beta_face coefficients
    MultiFab alpha_cc(ba, dm, 1, 0); // cell-centered 
    alpha_cc.setVal(0.); // fill in alpha_cc multifab to the value of 0.0
    mlabec.setACoeffs(amrlev, alpha_cc); 

    std::array< MultiFab, AMREX_SPACEDIM > beta_face; // beta_face is the permittivity in Poisson's equation
    // beta_face lives on faces so we make an array of face-centered MultiFabs
    AMREX_D_TERM(beta_face[0].define(convert(ba,IntVect(AMREX_D_DECL(1,0,0))), dm, 1, 0);,
                 beta_face[1].define(convert(ba,IntVect(AMREX_D_DECL(0,1,0))), dm, 1, 0);,
                 beta_face[2].define(convert(ba,IntVect(AMREX_D_DECL(0,0,1))), dm, 1, 0););
    
    InitializePermittivity(eps_cc, prob_lo, prob_hi, geom, eps_0, eps_r); //initialize permittivity to cell centers
    AveragePermittivityToCellFaces(eps_cc, beta_face); //converts from cell-centered permittivity to face-center and store in beta_face

    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(beta_face));

    //initialize charge density rho (cell-centered). This is RHS of Poisson solver.
    InitializeCharge(rho, prob_lo, prob_hi, geom); 

    //build an MLMG solver
    MLMG mlmg(mlabec); 
    mlmg.setVerbose(mlmg_set_verbose); 

    Real time = 0.0; //starting time in the simulation

    // call to Poisson solve
    // relative and absolute tolerances for linear solve
    const Real tol_rel = 1.e-10;
    const Real tol_abs = 0.0;
    mlmg.solve({&PoissonPhi}, {&rho}, tol_rel, tol_abs); 

    //*** POISSON SOLVER SETUP AND RUN END


    //*** PRINT PLOT FILES START

    amrex::Print() << "\n ========= Printing plot files ========== \n"<< std::endl;

    // Write a plotfile of the initial data if plot_int > 0
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,4);
        MultiFab::Copy(Plt, eps_cc, 0, 0, Ncomp1, Nghost0);
        MultiFab::Copy(Plt, PoissonPhi, 0, 1, Ncomp1, Nghost0);
        MultiFab::Copy(Plt, rho, 0, 2, Ncomp1, Nghost0);
        WriteSingleLevelPlotfile(pltfile, Plt, {"eps_cc", "phi", "rho"}, geom, time, 0);
    }

    //*** PRINT PLOT FILES END


    //*** RUN DIAGNOSTICS START

    // MultiFab memory usage
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    amrex::Long min_fab_megabytes  = amrex::TotalBytesAllocatedInFabsHWM()/1048576;
    amrex::Long max_fab_megabytes  = min_fab_megabytes;

    ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
    ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

    amrex::Print() << "High-water FAB megabyte spread across MPI nodes: ["
                   << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";

    min_fab_megabytes  = amrex::TotalBytesAllocatedInFabs()/1048576;
    max_fab_megabytes  = min_fab_megabytes;

    ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
    ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

    amrex::Print() << "Curent     FAB megabyte spread across MPI nodes: ["
                   << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";

    
    Real total_step_stop_time = ParallelDescriptor::second() - total_step_strt_time;
    ParallelDescriptor::ReduceRealMax(total_step_stop_time);

    amrex::Print() << "Total run time " << total_step_stop_time << " seconds\n";

    //*** RUN DIAGNOSTICS END
    

}
