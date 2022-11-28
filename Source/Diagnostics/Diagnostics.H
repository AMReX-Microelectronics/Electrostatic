#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include "Diagnostics_fwd.H"
#include "Diagnostics_Using_EB.H"

//#include <AMReX_REAL.H>
//#include <AMReX_GpuContainers.H>
//#include <AMReX_Vector.H>
//#include <AMReX_Array.H>
//#include <AMReX_IntVect.H>
//#include <AMReX_BoxArray.H>
//#include <AMReX_Geometry.H>
//#include <AMReX_MultiFab.H>
//
//#include <AMReX.H>
//
//#include <string>
//#include <any>

class
c_Diagnostics: private c_Diagnostics_Using_EB
{
public:
    c_Diagnostics();
    ~c_Diagnostics();

    int specify_using_eb;

    void InitData();
    void ComputeAndWriteDiagnostics(int step, amrex::Real time);

private:

    void ReadDiagnostics();

    std::string m_foldername_str;
    std::string m_diagnostics_str = "_diag";

};
#endif
