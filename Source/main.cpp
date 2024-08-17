
#include <AMReX_TinyProfiler.H>

#include "Code.H"
#include "CodeUtil.H"
#include "Utils/SelectWarpXUtils/WarpXProfilerWrapper.H"
#include "Utils/SelectWarpXUtils/WarpXUtil.H"

using namespace amrex;

int main(int argc, char *argv[])
{
    amrex::Initialize(argc, argv);

    amrex::Real initial_time = ParallelDescriptor::second();

    {
        WARPX_PROFILE_VAR("main()", pmain);

        c_Code pCode;
        amrex::ParmParse pp;

        pCode.InitData();

        // pCode.EstimateOfRequiredMemory();

        pCode.PrintGlobalWarnings("the initialization step");

        pCode.Solve_PostProcess_Output();

        pCode.Cleanup();

        WARPX_PROFILE_VAR_STOP(pmain);
    }

    PrintRunDiagnostics(initial_time);

    amrex::Finalize();
}
