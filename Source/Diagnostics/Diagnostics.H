#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include "Diagnostics_fwd.H"
#include "Diagnostics_Using_EB.H"

class
c_Diagnostics: private c_Diagnostics_Using_EB
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
    std::string m_diagnostics_str = "_diag";

};
#endif
