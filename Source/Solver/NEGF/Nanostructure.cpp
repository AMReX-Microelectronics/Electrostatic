#include "Nanostructure.H"

//#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
//#include "../../Utils/SelectWarpXUtils/TextMsg.H"
//#include "../../Utils/CodeUtils/CodeUtil.H"
//#include "Code.H"
//#include "Input/GeometryProperties/GeometryProperties.H"
//#include "Input/MacroscopicProperties/MacroscopicProperties.H"
//#include "Input/BoundaryConditions/BoundaryConditions.H"
//
//#include <AMReX.H>
//#include <AMReX_ParmParse.H>
//#include <AMReX_Parser.H>
//#include <AMReX_RealBox.H>
//#include <AMReX_MultiFab.H>
////#include <AMReX_MultiFabUtil.H>


using namespace amrex;

c_Nanostructure::c_Nanostructure
                               (const amrex::Geometry            & geom,
                                const amrex::DistributionMapping & dmap,
                                const amrex::BoxArray            & ba)
    : amrex::ParticleContainer<realPD::NUM, intPD::NUM,
                               realPA::NUM, intPA::NUM> (geom, dmap, ba)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Nanostructure_Atom_Container() Constructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    ReadNanostructureProperties();
  
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_Nanostructure_Atom_Container() Constructor************************\n";
#endif
}

void 
c_Nanostructure::ReadNanostructureProperties() 
{
   ReadAtomLocations(); 
}

void 
c_Nanostructure::ReadAtomLocations() 
{

}
//c_Nanostructure_Atom_Container::~c_Nanostructure_Atom_Container()
//    : amrex::ParticleContainer<realPD::NUM, intPD::NUM,
//                               realPA::NUM, intPA::NUM> ~()
//{
//#ifdef PRINT_NAME
//    amrex::Print() << "\n\n\t\t\t{************************c_Nanostructure_Atom_Container() Constructor************************\n";
//    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
//#endif
//
//#ifdef PRINT_NAME
//    amrex::Print() << "\t\t\t}************************c_Nanostructure_Atom_Container() Constructor************************\n";
//#endif
//}
