#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include "Diagnostics_fwd.H"

#ifdef AMREX_USE_EB
#include "Diagnostics_Using_EB.H"
#else
#include <AMReX_REAL.H>

#include <string>
#endif

class c_Diagnostics
#ifdef AMREX_USE_EB
    : private c_Diagnostics_Using_EB
#endif
{
   public:
    c_Diagnostics();
    ~c_Diagnostics();

    int specify_using_eb;

    void InitData();
    void ComputeAndWriteDiagnostics(int step, amrex::Real time);
    void ComputeFieldDiagnostics(int step, amrex::Real time);

   private:
    void ReadDiagnostics();

    std::string m_foldername_str;
    std::string m_diagnostics_str = "/diag";
};
#endif
