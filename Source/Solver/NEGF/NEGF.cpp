#include "NEGF.H"

//#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
//#include "../../Utils/SelectWarpXUtils/TextMsg.H"
//#include "../../Utils/CodeUtils/CodeUtil.H"
#include "../../Code.H"
#include "../../Input/GeometryProperties/GeometryProperties.H"
//#include "Input/MacroscopicProperties/MacroscopicProperties.H"
//#include "Input/BoundaryConditions/BoundaryConditions.H"
//
#include <AMReX.H>
//#include <AMReX_ParmParse.H>
//#include <AMReX_Parser.H>
//#include <AMReX_RealBox.H>
//#include <AMReX_MultiFab.H>
////#include <AMReX_MultiFabUtil.H>


using namespace amrex;

c_NEGFSolver::c_NEGFSolver()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_NEGFSolver() Constructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    ReadData();
  
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_NEGFSolver() Constructor************************\n";
#endif
}


c_NEGFSolver::~c_NEGFSolver()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_NEGFSolver() Destructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    vp_CNT.clear();
    vp_Graphene.clear();
    vp_Silicon.clear();
  
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_NEGFSolver() Destructor************************\n";
#endif
}

void 
c_NEGFSolver::ReadData() 
{

    amrex::Print() << "\n##### NEGF Solver #####\n\n";

    amrex::ParmParse pp_negf("negf");
    num_NS = 0;

    amrex::Vector<std::string> temp_vec;
    bool NS_pecified = pp_negf.queryarr("NS_names", temp_vec);
    int c=0;
    for (auto it: temp_vec)
    {
        if (std::find(vec_NS_names.begin(), vec_NS_names.end(), it) == vec_NS_names.end()) 
        {
           vec_NS_names.push_back(it);
           ++c;
        }
    }
    temp_vec.clear();

    amrex::Print() << "##### negf.NS_names: ";
    for (auto it: vec_NS_names) amrex::Print() << it << "  ";
    amrex::Print() << "\n";
    
}

void 
c_NEGFSolver::InitData() 
{
    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();

    _geom = &rGprop.geom;
    _dm   = &rGprop.dm;
    _ba   = &rGprop.ba;
    //_n_cell = &rGprop.n_cell;

    std::string type;
    int c=0;
    for (auto name: vec_NS_names)
    {
        
        amrex::Print() << "\n##### NANOSTRUCTURE PROPERTIES #####\n\n";
     
        amrex::ParmParse pp_ns(name);
        pp_ns.get("type", type);
        amrex::Print() << "##### name: " << name << "\n";
        amrex::Print() << "##### type: " << type << "\n";

        switch (map_NSType_enum[type])
        {
            case s_NS::Type::CNT:
            {
                using T = c_CNT;
                vp_CNT.push_back(std::make_unique<c_Nanostructure<T>>(*_geom, *_dm, *_ba, name));
                break;
            }
            case s_NS::Type::Graphene:
            {
                amrex::Abort("NS_type " + type + " is not yet defined.");
            }
            case s_NS::Type::Silicon:
            {
                amrex::Abort("NS_type " + type + " is not yet defined.");
            }
            default:
            {
                amrex::Abort("NS_type " + type + " is not supported.");
            }
        }
        ++c;
    }

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
