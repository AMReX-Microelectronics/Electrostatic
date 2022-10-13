
#include "Code.H"
#include "CodeUtil.H"
#include "Utils/SelectWarpXUtils/WarpXUtil.H"
#include "Utils/SelectWarpXUtils/WarpXProfilerWrapper.H"


#include <AMReX_TinyProfiler.H>


using namespace amrex;

int main (int argc, char* argv[])
{

    amrex::Initialize(argc,argv);

#ifdef PRINT_NAME
    amrex::Print() << "\n\n************************************************************\n";
    amrex::Print() << "{**************************main(*)**************************\n";
    amrex::Print() << "in file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    amrex::Print() << "************************************************************\n";
#endif

    amrex::Real initial_time = ParallelDescriptor::second();
    
    {
        WARPX_PROFILE_VAR("main()", pmain);

        c_Code pCode; 

        pCode.InitData();

        pCode.PrintGlobalWarnings("the initialization step"); //Print warning at this stage

        //pCode.Output();

        pCode.Solve();

        pCode.PostProcess();

        pCode.Output();
 
        WARPX_PROFILE_VAR_STOP(pmain);

    } //destructor of c_Code is called here and all objects it created are destructed in reverse order.

    PrintRunDiagnostics(initial_time);

#ifdef PRINT_NAME
    amrex::Print() << "\n\n**********************************************************\n";
    amrex::Print() << "}**************************main(*)**************************\n";
    amrex::Print() << "**********************************************************\n\n";
#endif

    amrex::Finalize();

}
