//#include <AMReX_PlotFileUtil.H>
//#include <AMReX_ParmParse.H>
//#include <AMReX_MLABecLaplacian.H>
//#include <AMReX_MLMG.H> 
//#include <AMReX_MultiFab.H> 
#include <AMReX_VisMF.H>

#include "Code.H"

//#include "Global.H"
//#include "myfunc.H"
//#include "Input.H"
//#include "Initialize.H"
//#include "Poisson_MLMG.H"

using namespace amrex;

template<typename T>
class TD;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    amrex::Real total_step_strt_time = ParallelDescriptor::second();

    c_Code pCode; 

    pCode.InitData();

    pCode.Solve();

    pCode.PostProcess();

    pCode.PrintGlobalWarnings("THE END"); //Print warning messages at the end of the simulation

//    //*** PRINT PLOT FILES START
//
//    MultiFab Plt(ba, dm, plot_elems, Nghost0); //multifabs for plotting
//    amrex::Print() << "\n ========= Printing plot files ========== \n"<< std::endl;
//
//    // Write a plotfile of the initial data if plot_int > 0
//    if (plot_int > 0)
//    {
//        int step = 0;
//        amrex::Real time = 0.0; //starting time in the simulation
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

    amrex::Finalize();
    return 0;
}

//void main_main ()
//{
//
//    Array< MultiFab, AMREX_SPACEDIM > E; 
//    amrex::Vector<Array<MultiFab, AMREX_SPACEDIM>> minus_epsGradPhi;
//
//    
//
//}
