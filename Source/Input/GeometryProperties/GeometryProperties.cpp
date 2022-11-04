#include "GeometryProperties.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
//#include "Code.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_RealBox.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBSupport.H>


using namespace amrex;

c_GeometryProperties::c_GeometryProperties ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_GeometryProperties Constructor()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    ReadData();

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_GeometryProperties Constructor()************************\n";
#endif
} 


c_GeometryProperties::~c_GeometryProperties ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_GeometryProperties Destructor()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_GeometryProperties Destructor()************************\n";
#endif
} 

void 
c_GeometryProperties::ReadData()
{    
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_GeometryProperties::ReadData()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

     ParseBasicDomainInput();

    if(embedded_boundary_flag) eb.ReadGeometry();

#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t}************************c_GeometryProperties::ReadData()************************\n";
#endif
}


void 
c_GeometryProperties::InitData()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t{************************c_GeometryProperties::InitData()************************\n";
    amrex::Print() << "\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    InitializeBoxArrayAndDistributionMap();

    if(embedded_boundary_flag) CreateEmbeddedBoundary();

#ifdef PRINT_NAME
    amrex::Print() << "\t\t}************************c_GeometryProperties::InitData()************************\n";
#endif
}


void
c_GeometryProperties::ParseBasicDomainInput()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_GeometryProperties::ParseBasicDomainInput()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif
    std::string prt = "\t\t\t\t";

    amrex::Vector<int> num_cell;
    amrex::Vector<amrex::Real> prob_min(AMREX_SPACEDIM);
    amrex::Vector<amrex::Real> prob_max(AMREX_SPACEDIM);
    amrex::Vector<amrex::Real> mg{AMREX_D_DECL(128,128,128)}; //default values
    amrex::Vector<amrex::Real> bf{AMREX_D_DECL(8,8,8)};
    amrex::Vector<amrex::Real> periodicity{AMREX_D_DECL(0,0,0)};
    std::string coord_sys_str = "cartesian";
    coord_sys =  amrex::CoordSys::cartesian; //default
    embedded_boundary_flag = 0;

    amrex::ParmParse pp_domain("domain");

    getArrWithParser(pp_domain, "prob_lo", prob_min, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_lo.size() == AMREX_SPACEDIM);

    getArrWithParser(pp_domain, "prob_hi", prob_max, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_hi.size() == AMREX_SPACEDIM);

    pp_domain.getarr("n_cell", num_cell, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(n_cell.size() == AMREX_SPACEDIM);

    pp_domain.queryarr("max_grid_size", mg);

    pp_domain.queryarr("blocking_factor", bf);

    pp_domain.queryarr("is_periodic", periodicity);

    pp_domain.query("coord_sys", coord_sys_str);

    pp_domain.query("embedded_boundary", embedded_boundary_flag);

    //pp_domain.addarr("n_cell", num_cell);
    //pp_domain.addarr("prob_lo", prob_min);
    //pp_domain.addarr("prob_hi", prob_max);
    //pp_domain.addarr("max_grid_size", mg);
    //pp_domain.addarr("blocking_factor", bf);

    for (int i=0; i<AMREX_SPACEDIM; ++i) 
    {
        n_cell[i] = num_cell[i];
        prob_lo[i] = prob_min[i]; //Converting vector to GpuArray
        prob_hi[i] = prob_max[i]; 
        max_grid_size[i] = mg[i];  //Converting Vector to IntVect
        blocking_factor[i] = bf[i]; 
        is_periodic[i] = periodicity[i]; 
    }
    if(coord_sys_str == "cartesian") 
    {
        coord_sys =  amrex::CoordSys::cartesian;
    }
    else if(coord_sys_str == "radial") 
    {
        coord_sys = amrex::CoordSys::RZ;
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i) 
    {
        amrex::Print() << prt << "\n";
        amrex::Print() << prt << "direction: " << i << "\n";
        amrex::Print() << prt << "prob_lo: " << prob_lo[i] << "\n";
        amrex::Print() << prt << "prob_hi: " << prob_hi[i] << "\n";
        amrex::Print() << prt << "max_grid_size: " << max_grid_size[i] << "\n";
        amrex::Print() << prt << "blocking_factor: " << blocking_factor[i] << "\n";
        amrex::Print() << prt << "is_periodic: " << is_periodic[i] << "\n";
    }
    amrex::Print() << prt << "\n";
    amrex::Print() << prt << "coord_sys: " << coord_sys << "\n";
    amrex::Print() << prt << "embedded_boundary_flag: " << embedded_boundary_flag << "\n";

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_GeometryProperties::ParseBasicDomainInput()************************\n";
#endif
}


void 
c_GeometryProperties::InitializeBoxArrayAndDistributionMap()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_GeometryProperties::InitializeBoxArrayAndDistributionMap()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    amrex::IntVect dom_lo(AMREX_D_DECL(       0,        0,        0)); // domain low
    amrex::IntVect dom_hi(AMREX_D_DECL(n_cell[0]-1, n_cell[1]-1, n_cell[2]-1)); // domain high

    amrex::Box domain(dom_lo, dom_hi); // Make a single box that is the entire domain

    ba.define(domain); // initialize the boxarray 'ba' from the single box 'domain'

    ba.maxSize(max_grid_size); // break up ba into chunks no larger than 'max_grid_size' along a direction

    amrex::RealBox real_box({AMREX_D_DECL( prob_lo[0], prob_lo[1], prob_lo[2])},
                    {AMREX_D_DECL( prob_hi[0], prob_hi[1], prob_hi[2])});  //physical domain

    geom.define(domain, real_box, coord_sys, is_periodic); //define the geom object

    dm.define(ba);

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_GeometryProperties::InitializeBoxArrayAndDistributionMap()************************\n";
#endif
}


void 
c_GeometryProperties::CreateEmbeddedBoundary()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_GeometryProperties::CreateEmbeddedBoundary()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    if(eb.specify_input_using_eb2) 
    {
        EB2::Build(geom, eb.required_coarsening_level, eb.max_coarsening_level);
    }
    else {
        eb.ConstructFinalObject("Hey",geom);
    }


    const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
    const EB2::Level& eb_level = eb_is.getLevel(geom);

    // number of ghost cells for each of the 3 EBSupport types
    Vector<int> ng_ebs = {2,2,2};

    // This object provides access to the EB database in the format of basic AMReX objects
    // such as BaseFab, FArrayBox, FabArray, and MultiFab
    eb.pFactory = amrex::makeEBFabFactory(&eb_level, ba, dm, ng_ebs, eb.support);
    
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_GeometryProperties::CreateEmbeddedBoundary()************************\n";
#endif
}
