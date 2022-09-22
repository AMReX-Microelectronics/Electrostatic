#ifndef CODE_H_
#define CODE_H_

#include "Input/GeometryProperties_fwd.H"
#include "Input/MaterialProperties_fwd.H"

#include <AMReX.H>
//#include <AMReX_AmrCore.H>
//#include <AMReX_Array.H>
//#include <AMReX_Config.H>
#include <AMReX_GpuContainers.H>
//#include <AMReX_IntVect.H>
//#include <AMReX_LayoutData.H>
//#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
//#include <AMReX_RealBox.H>
//#include <AMReX_RealVect.H>
//#include <AMReX_Vector.H>
//#include <AMReX_VisMF.H>

//#include <array>
//#include <iostream>
//#include <limits>
//#include <memory>
//#include <optional>
//#include <string>
//#include <vector>


//#include <AMReX_AmrCoreFwd.H>

//struct s_Domain
//{
//
//    amrex::GpuArray<int, AMREX_SPACEDIM> n_cell;     // number of cells in each dimension
//    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo; // physical lo coordinate
//    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi; // physical hi coordinate
//    amrex::Vector<int> max_grid_size;
//    amrex::Vector<int> blocking_factor;
//
//};

class c_Code
{

public:

    static c_Code& GetInstance ();
    static void ResetInstance ();

    c_Code ();
    ~c_Code ();
    
//    s_Domain domain;
    void ReadData ();
    void InitData ();

    c_GeometryProperties& GetGeometryPropertiesPointer () { return *m_pGeometryProperties;}
    c_MaterialProperties& GetMaterialPropertiesPointer () { return *m_pMaterialProperties;}

private:
    static c_Code* m_instance;

    std::unique_ptr<c_GeometryProperties> m_pGeometryProperties;
    std::unique_ptr<c_MaterialProperties> m_pMaterialProperties;

};
#endif
