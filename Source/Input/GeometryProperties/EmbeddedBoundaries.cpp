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

c_EmbeddedBoundaries::c_EmbeddedBoundaries()
{

}

void
c_EmbeddedBoundaries::ReadGeometry()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_GeometryProperties::ReadGeometry()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif
    std::string prt = "\t\t\t\t";

    //setting default values
    required_coarsening_level = 0; // max amr level (at present 0)
    max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible
    support = EBSupport::full;
    std::string eb_support_str = "full";
    specify_input_using_eb2 = 0;

    amrex::ParmParse pp_ebgeom("ebgeom");
    pp_ebgeom.query("required_coarsening_level", required_coarsening_level);
    pp_ebgeom.query("max_coarsening_level", max_coarsening_level);
    pp_ebgeom.query("support", eb_support_str);
    support = map_eb_support[eb_support_str];
    pp_ebgeom.query("specify_input_using_eb2", specify_input_using_eb2);
    amrex::Print() << "\n";
    amrex::Print() << prt << "ebgeom.required_coarsening_level: " << required_coarsening_level << "\n";
    amrex::Print() << prt << "ebgeom.max_coarsening_level: " << max_coarsening_level << "\n";
    amrex::Print() << prt << "ebgeom.support: " << eb_support_str << "\n";
    amrex::Print() << prt << "ebgeom.specify_input_using_eb2: " << specify_input_using_eb2 << "\n";

    if(!specify_input_using_eb2) 
    {
        num_basic_objects = 0;
        amrex::Vector< std::string > vec_basic_objects;
        bool basic_objects_specified = pp_ebgeom.queryarr("basic_objects", vec_basic_objects);
//        int c=0;
//        for (auto it: vec_basic_objects)
//        {
//            if (map_basic_objects_rank.find(it) == map_basic_objects_rank.end()) {
//                    map_basic_objects_rank[it] = c;
//                    ++c;
//            }
//        }
        int c=0;
        for (auto it: vec_basic_objects)
        {
            if (map_basic_objects_type.find(it) == map_basic_objects_type.end()) {
                amrex::ParmParse pp_object(it);

                pp_object.get("geom_type", map_basic_objects_type[it]);

                ReadObjectInfo(it, map_basic_objects_type[it], pp_object);
                ++c;
            }
        }
        vec_basic_objects.clear();
        num_basic_objects = c;
        amrex::Print() << "\ntotal number of basic objects: " << num_basic_objects << "\n";
        amrex::Print() << "\nObject Info:\n";
        for (auto it: map_basic_objects_info)
        {
            std::string name = it.first; 
            std::string geom_type = map_basic_objects_type[name];
            PrintObjectInfo(name, geom_type, it.second);
        }

        amrex::Vector< std::string > vec_final_objects;
        bool final_objects_specified = pp_ebgeom.queryarr("final_objects", vec_final_objects);
        c=0;
        for (auto it: vec_final_objects)
        {
            if (map_final_objects_rank.find(it) == map_final_objects_rank.end()) {
                map_final_objects_rank[it] = c;
                amrex::ParmParse pp_final_object(it);
                std::string construct_main;
                pp_final_object.get("construct", construct_main);
                
                //ConstructFinalObject(construct_main);
                  
                ++c;
            }
        }
        vec_final_objects.clear();
        num_final_objects = map_final_objects_rank.size();

        amrex::Print() << "\ntotal number of final objects: " << num_final_objects << "\n";
        for (auto it: map_final_objects_rank)
        {
            amrex::Print() << "final object no: " << it.second  << " name: " << it.first << "\n";
        }
        

//
//        EB2::Build(gshop, geom, required_coarsening_level,
//                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening,
//                   a_extend_domain_face, a_num_coarsen_opt);
//
    }
    //pp_domain.query("embedded_boundary", embedded_boundary_flag);
    //amrex::Vector<int> num_cell;
    //amrex::Vector<amrex::Real> prob_min(AMREX_SPACEDIM);
    //amrex::Vector<amrex::Real> prob_max(AMREX_SPACEDIM);
    //amrex::Vector<amrex::Real> mg{AMREX_D_DECL(128,128,128)}; //default values
    //amrex::Vector<amrex::Real> bf{AMREX_D_DECL(8,8,8)};
    //amrex::Vector<amrex::Real> periodicity{AMREX_D_DECL(0,0,0)};
    //std::string coord_sys_str = "cartesian";
    //coord_sys =  amrex::CoordSys::cartesian; //default
    //embedded_boundary_flag = 0;


    //getArrWithParser(pp_domain, "prob_lo", prob_min, 0, AMREX_SPACEDIM);
    //AMREX_ALWAYS_ASSERT(prob_lo.size() == AMREX_SPACEDIM);

    //getArrWithParser(pp_domain, "prob_hi", prob_max, 0, AMREX_SPACEDIM);
    //AMREX_ALWAYS_ASSERT(prob_hi.size() == AMREX_SPACEDIM);

    //pp_domain.getarr("n_cell", num_cell, 0, AMREX_SPACEDIM);
    //AMREX_ALWAYS_ASSERT(n_cell.size() == AMREX_SPACEDIM);

    //pp_domain.queryarr("max_grid_size", mg);

    //pp_domain.queryarr("blocking_factor", bf);

    //pp_domain.queryarr("is_periodic", periodicity);

    //pp_domain.query("coord_sys", coord_sys_str);

    //pp_domain.query("embedded_boundary", embedded_boundary_flag);

    //for (int i=0; i<AMREX_SPACEDIM; ++i) 
    //{
    //    n_cell[i] = num_cell[i];
    //    prob_lo[i] = prob_min[i]; //Converting vector to GpuArray
    //    prob_hi[i] = prob_max[i]; 
    //    max_grid_size[i] = mg[i];  //Converting Vector to IntVect
    //    blocking_factor[i] = bf[i]; 
    //    is_periodic[i] = periodicity[i]; 
    //}


//#ifdef PRINT_LOW
//    for (int i=0; i<AMREX_SPACEDIM; ++i) 
//    {
//        amrex::Print() << prt << "\n";
//        amrex::Print() << prt << "direction: " << i << "\n";
//        amrex::Print() << prt << "prob_lo: " << prob_lo[i] << "\n";
//        amrex::Print() << prt << "prob_hi: " << prob_hi[i] << "\n";
//        amrex::Print() << prt << "max_grid_size: " << max_grid_size[i] << "\n";
//        amrex::Print() << prt << "blocking_factor: " << blocking_factor[i] << "\n";
//        amrex::Print() << prt << "is_periodic: " << is_periodic[i] << "\n";
//    }
//    amrex::Print() << prt << "\n";
//    amrex::Print() << prt << "coord_sys: " << coord_sys << "\n";
//    amrex::Print() << prt << "embedded_boundary_flag: " << embedded_boundary_flag << "\n";
//#endif

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_GeometryProperties::ReadGeometry()************************\n";
#endif
}

void
c_EmbeddedBoundaries::ReadObjectInfo(std::string object_name, std::string object_type, amrex::ParmParse pp_object)
{

    switch (map_object_type_enum[object_type]) 
    {  
        case s_ObjectType::object::sphere:
        {
            RealArray center;
            pp_object.get("sphere_center", center);
            Real radius;
            pp_object.get("sphere_radius", radius);
            bool has_fluid_inside;
            pp_object.get("sphere_has_fluid_inside", has_fluid_inside);

//            s_Sphere temp = {center, radius, has_fluid_inside};

            EB2::SphereIF sphere(radius, center, has_fluid_inside);

            map_basic_objects_info[object_name] = sphere;  
            break;
        }
        case s_ObjectType::object::box:
        {
            RealArray lo;
            pp_object.get("box_lo", lo);
        
            RealArray hi;
            pp_object.get("box_hi", hi);
        
            bool has_fluid_inside;
            pp_object.get("box_has_fluid_inside", has_fluid_inside);
//          s_Box temp = {lo, hi, has_fluid_inside};

            EB2::BoxIF box(lo, hi, has_fluid_inside);

            map_basic_objects_info[object_name] = box;  
            break;
        }
    }
}

void
c_EmbeddedBoundaries::PrintObjectInfo(std::string object_name, std::string object_type, std::any object_info)
{
//
//    switch (map_object_type_enum[object_type]) 
//    {  
//        case s_ObjectType::object::sphere:
//        {
//            amrex::Print() << "\nObject name: " << object_name << "\n";
//            amrex::Print() << "Object type: sphere \n";
//            auto sphere = std::any_cast<EB2::SphereIF>(object_info); 
//            amrex::Print() << "Sphere center: " << sphere.m_center[0] << "  " << sphere.m_center[1] << "  " << sphere.m_center[2] << "\n";
//            amrex::Print() << "Sphere radius: " << sphere.m_radius << "\n";
//            amrex::Print() << "Sphere has_fluid_inside: " << sphere.m_sign << "\n";
//            break;
//        }
//        case s_ObjectType::object::box:
//        {
//        //    amrex::Print() << "\nObject name: " << object_name << "\n";
//        //    amrex::Print() << "Object type: box \n";
//
//        //    auto box = std::any_cast<EB2::BoxIF>(object_info); 
//        //    amrex::Print() << "Box lo: " << box.lo[0] << "  " << box.lo[1] << "  " << box.lo[2] << "\n";
//        //    amrex::Print() << "Box hi: " << box.hi[0] << "  " << box.hi[1] << "  " << box.hi[2] << "\n";
//        //    amrex::Print() << "Box has_fluid_inside: " << box.has_fluid_inside << "\n";
//            break;
//        }
//    }
}


//void
//c_EmbeddedBoundaries::BuildEBObjects(std::string object_name, std::string object_type, std::any object_info)
//{
//
//    switch (map_object_type_enum[object_type]) 
//    {  
//        case s_ObjectType::object::sphere:
//        {
//            amrex::Print() << "\nObject name: " << object_name << "\n";
//            amrex::Print() << "Object type: sphere \n";
//            s_Sphere sphere_info = std::any_cast<s_Sphere>(object_info); 
//            amrex::Print() << "Sphere center: " << sphere_info.center[0] << "  " << sphere_info.center[1] << "  " << sphere_info.center[2] << "\n";
//            amrex::Print() << "Sphere radius: " << sphere_info.radius << "\n";
//            amrex::Print() << "Sphere has_fluid_inside: " << sphere_info.has_fluid_inside << "\n";
//            break;
//        }
//        case s_ObjectType::object::box:
//        {
//            amrex::Print() << "\nObject name: " << object_name << "\n";
//            amrex::Print() << "Object type: box \n";
//
//            s_Box box_info = std::any_cast<s_Box>(object_info); 
//            amrex::Print() << "Box lo: " << box_info.lo[0] << "  " << box_info.lo[1] << "  " << box_info.lo[2] << "\n";
//            amrex::Print() << "Box hi: " << box_info.hi[0] << "  " << box_info.hi[1] << "  " << box_info.hi[2] << "\n";
//            amrex::Print() << "Box has_fluid_inside: " << box_info.has_fluid_inside << "\n";
//            break;
//        }
//    }
//}

template<typename T>
class TD;

void
c_EmbeddedBoundaries::ConstructFinalObject(std::string construct_main, amrex::Geometry geom,amrex::BoxArray ba, amrex::DistributionMapping dm)
{

//    std::string str = construct_main;
//    
//    std::map<std::string,std::string>::iterator it;
//
//    it = map_basic_objects_type.find(str);
//    if(it != map_basic_objects_type.end())
//    {
//    }
//    auto object_name = it->first;
//    auto geom_type = it->second;

//Geometry1
    auto object_name1 = "Sph1";
    auto object_name2 = "Box1";

    auto ob1 = std::any_cast<amrex::EB2::SphereIF>(map_basic_objects_info[object_name1]);
    auto ob2 = std::any_cast<amrex::EB2::BoxIF>(map_basic_objects_info[object_name2]);

    auto cubesphere1 = amrex::EB2::makeIntersection(ob1, ob2);
    auto gshop1 = amrex::EB2::makeShop(cubesphere1);

    amrex::EB2::Build(gshop1, geom, required_coarsening_level, max_coarsening_level);

    const EB2::IndexSpace& eb_is1 = EB2::IndexSpace::top();
    const EB2::Level& eb_level1 = eb_is1.getLevel(geom);

    // number of ghost cells for each of the 3 EBSupport types
    Vector<int> ng_ebs = {2,2,2};

    // This object provides access to the EB database in the format of basic AMReX objects
    // such as BaseFab, FArrayBox, FabArray, and MultiFab
    pFactory1 = amrex::makeEBFabFactory(&eb_level1, ba, dm, ng_ebs, support);

//Geometry2
    auto object_name3 = "Sph2";
    auto object_name4 = "Box2";

    auto ob3 = std::any_cast<amrex::EB2::SphereIF>(map_basic_objects_info[object_name3]);
    auto ob4 = std::any_cast<amrex::EB2::BoxIF>(map_basic_objects_info[object_name4]);

    auto cubesphere2 = amrex::EB2::makeIntersection(ob3, ob4);
    auto gshop2 = amrex::EB2::makeShop(cubesphere2);

    amrex::EB2::Build(gshop2, geom, required_coarsening_level, max_coarsening_level);

    const EB2::IndexSpace& eb_is2 = EB2::IndexSpace::top();
    const EB2::Level& eb_level2 = eb_is2.getLevel(geom);

    pFactory2 = amrex::makeEBFabFactory(&eb_level2, ba, dm, ng_ebs, support);

//Combined Geometry

    
    auto unionobject = amrex::EB2::makeUnion(cubesphere1, cubesphere2);
    auto gshop = amrex::EB2::makeShop(unionobject);

    amrex::EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);

    const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
    const EB2::Level& eb_level = eb_is.getLevel(geom);
    pFactory = amrex::makeEBFabFactory(&eb_level, ba, dm, ng_ebs, support);
}
