//#include <AMReX_PlotFileUtil.H>
//#include <AMReX_ParmParse.H>
//#include <AMReX_MLABecLaplacian.H>
//#include <AMReX_MLMG.H> 
//#include <AMReX_MultiFab.H> 
//#include <AMReX_VisMF.H>
#include "Code.H"
//#include "Global.H"
//#include "myfunc.H"
//#include "Input.H"
//#include "Initialize.H"
//#include "Poisson_MLMG.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    c_Code p_Code;
    p_Code.InitializeData();

    auto& n_cell = p_Code.domain.n_cell;
    auto& prob_lo = p_Code.domain.prob_lo;
				auto& prob_hi = p_Code.domain.prob_hi;
    auto& max_grid_size = p_Code.domain.max_grid_size;
    Print() << "n_cell: " << n_cell[0] << " " << n_cell[1] << " " << n_cell[2] << "\n";
    Print() << "prob_lo: " << prob_lo[0] << " " << prob_lo[1] << " " << prob_lo[2] << "\n";
    Print() << "prob_hi: " << prob_hi[0] << " " << prob_hi[1] << " " << prob_hi[2] << "\n";
    Print() << "max_grid_size: " << max_grid_size[0]  << "\n";


//    main_main();

    amrex::Finalize();
    return 0;
}

//void main_main ()
//{
//
//    Real total_step_strt_time = ParallelDescriptor::second();
//
//    //*** SIMULATION INPUT START
//
//    //Basic simulation input
//    amrex::GpuArray<int, AMREX_SPACEDIM> n_cell;     // number of cells in each dimension
//    int max_grid_size; // size of maximum box dimension used to divide the domain
//    amrex::GpuArray<Real, AMREX_SPACEDIM> prob_lo; // physical lo coordinate
//    amrex::GpuArray<Real, AMREX_SPACEDIM> prob_hi; // physical hi coordinate
//    int plot_int; // plot file writing interval
//    int plot_elems;
//    ReadBasicInput(n_cell, max_grid_size, prob_lo, prob_hi, nsteps, plot_int, plot_elems);
//    Print() << "prob_lo: " << prob_lo[0] << " " << prob_lo[1] << " " << prob_lo[2] << "\n";
//    Print() << "prob_hi: " << prob_hi[0] << " " << prob_hi[1] << " " << prob_hi[2] << "\n";
//    
//    //Material input
//    Real mlmg_ascalar, mlmg_bscalar;
//    int mlmg_set_verbose;
//    int mlmg_max_order;
//    ReadMLMGInput(mlmg_ascalar, mlmg_bscalar, mlmg_set_verbose, mlmg_max_order); 
//
//    //*** SIMULATION INPUT END
//
//
//    //*** DOMAIN SETUP START
//
//    // create ba: BoxArray and geom: Geometry
//    BoxArray ba; // a list of boxes that cover the domain
//    Geometry geom; // contains info such as the physical domain size, number of points in the domain, and periodicity
//    InitializeBoxArrayAndGeom(ba, geom, prob_lo, prob_hi, n_cell, max_grid_size);
//    DistributionMapping dm(ba); //distribute boxes over MPI processors
//
//    //*** DOMAIN SETUP END
//
//
//    //*** MULTIFAB ARRAYS DEFINE START
//
//    int Nghost1 = 1; //number of ghost cells for each array
//    int Nghost0 = 0; //number of ghost cells for each array
//    int Ncomp1 = 1; //number of components for each array
//
//    MultiFab eps_cc(ba, dm, Ncomp1, Nghost1); //permittivity at cell-centers
//    MultiFab rho(ba, dm, Ncomp1, Nghost0); //charge density
//    MultiFab PoissonPhi(ba, dm, Ncomp1, Nghost1);
//    MultiFab Plt(ba, dm, plot_elems, Nghost0); //multifabs for plotting
////TD<decltype(eps_cc[0])>eps_cc_type;
//    InitializePermittivity(eps_cc, prob_lo, prob_hi, geom, eps_0, eps_r); //initialize permittivity (cell-centered)
//    InitializeCharge(rho, prob_lo, prob_hi, geom); //initialize charge density rho (cell-centered). 
//
//    //*** MULTIFAB ARRAYS DEFINE END
//
//
//    //*** POISSON SOLVER SETUP AND RUN START
//
//    LPInfo info;
//    MLABecLaplacian mlabec({geom}, {ba}, {dm}, info); // Implicit solve using MLABecLaplacian class
//
//    setupMLABecForPoissonProblem(mlabec, mlmg_ascalar, mlmg_bscalar, mlmg_max_order,
//                                 ba, dm,
//                                 PoissonPhi, eps_cc, rho);
//    MLMG mlmg(mlabec); //build an MLMG solver
//    mlmg.setVerbose(mlmg_set_verbose);
//    Real time = 0.0; //starting time in the simulation
//    const Real tol_rel = 1.e-10; //relative tolerance
//    const Real tol_abs = 0.0; //absolute tolerance, 0 means not known.
//    mlmg.solve({&PoissonPhi}, {&rho}, tol_rel, tol_abs); //call to Poisson solve
// 
//    //*** POISSON SOLVER SETUP AND RUN END
//    
//
//    //*** POSTPROCESS START
//    Array< MultiFab, AMREX_SPACEDIM > E; 
//    GetEFieldUsingMLMG(mlmg,E,ba,dm);
//
//    amrex::Vector<Array<MultiFab, AMREX_SPACEDIM>> minus_epsGradPhi;
//    GetFluxesUsingMLMG(mlmg, minus_epsGradPhi, ba,dm);
//
//    ComputeSumOfFluxes(minus_epsGradPhi);
//
//    //*** POSTPROCESS END
//
//    //*** PRINT PLOT FILES START
//
//    amrex::Print() << "\n ========= Printing plot files ========== \n"<< std::endl;
//
//    // Write a plotfile of the initial data if plot_int > 0
//    if (plot_int > 0)
//    {
//        int step = 0;
//        const std::string& pltfile = amrex::Concatenate("plt",step,4);
//        MultiFab::Copy(Plt, eps_cc, 0, 0, Ncomp1, Nghost0);
//        MultiFab::Copy(Plt, PoissonPhi, 0, 1, Ncomp1, Nghost0);
//        MultiFab::Copy(Plt, rho, 0, 2, Ncomp1, Nghost0);
//        MultiFab::Copy(Plt, E[0], 0, 3, Ncomp1, Nghost0);
//        MultiFab::Copy(Plt, E[1], 0, 4, Ncomp1, Nghost0);
//        MultiFab::Copy(Plt, E[2], 0, 5, Ncomp1, Nghost0);
//        MultiFab::Copy(Plt, minus_epsGradPhi[0][0], 0, 6, Ncomp1, Nghost0);
//        MultiFab::Copy(Plt, minus_epsGradPhi[0][1], 0, 7, Ncomp1, Nghost0);
//        MultiFab::Copy(Plt, minus_epsGradPhi[0][2], 0, 8, Ncomp1, Nghost0);
//        WriteSingleLevelPlotfile(pltfile, Plt, {"eps_cc", "phi", "rho", "Ex", "Ey", "Ez", "-epsGradPhi_x", "-epsGradPhi_y", "-epsGradPhi_z"}, geom, time, 0);
//    }
//
//    //*** PRINT PLOT FILES END
//
//
//    //*** RUN DIAGNOSTICS START
//
//    // MultiFab memory usage
//    const int IOProc = ParallelDescriptor::IOProcessorNumber();
//
//    amrex::Long min_fab_megabytes  = amrex::TotalBytesAllocatedInFabsHWM()/1048576;
//    amrex::Long max_fab_megabytes  = min_fab_megabytes;
//
//    ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
//    ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);
//
//    amrex::Print() << "High-water FAB megabyte spread across MPI nodes: ["
//                   << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
//
//    min_fab_megabytes  = amrex::TotalBytesAllocatedInFabs()/1048576;
//    max_fab_megabytes  = min_fab_megabytes;
//
//    ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
//    ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);
//
//    amrex::Print() << "Curent     FAB megabyte spread across MPI nodes: ["
//                   << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
//
//    
//    Real total_step_stop_time = ParallelDescriptor::second() - total_step_strt_time;
//    ParallelDescriptor::ReduceRealMax(total_step_stop_time);
//
//    amrex::Print() << "Total run time " << total_step_stop_time << " seconds\n";
//
//    //*** RUN DIAGNOSTICS END
//    
//
//}
