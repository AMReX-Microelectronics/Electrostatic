#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>

class MaterialProperties
{
public: 
    MaterialProperties();
    void ReadData();
    void InitData();
    
    amrex::MultiFab& getpermittivity() {return (*m_eps_mf);}

    void InitializeMacroMultiFabUsingParser (amrex::MultiFab *macro_mf,
                                             amrex::ParserExecutor<3> const& macro_parser,
                                             const int lev);
    std::string m_epsilon_s = "constant";
    amrex::Real m_epsilon = PhysConst::ep0;
    std::unique_ptr<amrex::Parser> m_epsilon_parser;
private:

    /** Pointer to multifab of permittivity, eps */
    std::unique_ptr<amrex::MultiFab> m_epsilon_mf;
    /** String for storing parser function */
    std::string m_str_epsilon_function;
}
