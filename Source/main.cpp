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


//#if defined PRINT_HIGH
//
//    #define PRINT_MEDIUM = true
//    #define PRINT_LOW = true
//
//#elif defined PRINT_MEDIUM
//
//    #define PRINT_LOW = true
//
//#endif

using namespace amrex;

template<typename T>
class TD;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

#ifdef PRINT_NAME
    amrex::Print() << "\n\nSTART************************main(*)************************\n";
    amrex::Print() << "in file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif


    amrex::Real total_run_time = ParallelDescriptor::second();

    c_Code pCode; 

    pCode.InitData();

    pCode.PrintGlobalWarnings("the initialization step"); //Print warning at this stage

    pCode.Solve();

    pCode.PostProcess();

    pCode.Output();


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

    
    total_run_time = ParallelDescriptor::second() - total_run_time;
    ParallelDescriptor::ReduceRealMax(total_run_time);

    amrex::Print() << "Total run time " << total_run_time << " seconds\n";

    //*** RUN DIAGNOSTICS END

#ifdef PRINT_NAME
    amrex::Print() << "END************************main(*)************************\n";
#endif
    amrex::Finalize();
    return 0;
}
