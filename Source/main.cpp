//#include <AMReX_PlotFileUtil.H>
//#include <AMReX_ParmParse.H>
//#include <AMReX_MLABecLaplacian.H>
//#include <AMReX_MLMG.H> 
//#include <AMReX_MultiFab.H> 
#include <AMReX_VisMF.H>

#include "Code.H"
#include "CodeUtil.H"

using namespace amrex;

template<typename T>
class TD;

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

    c_Code pCode; 

    pCode.InitData();

    pCode.PrintGlobalWarnings("the initialization step"); //Print warning at this stage

//    pCode.Output();

    pCode.Solve();

    pCode.PostProcess();

    pCode.Output();

    PrintRunDiagnostics(initial_time);

#ifdef PRINT_NAME
    amrex::Print() << "\n\n**********************************************************\n";
    amrex::Print() << "}**************************main(*)**************************\n";
    amrex::Print() << "**********************************************************\n\n";
#endif

    amrex::Finalize();
    return 0;
}
