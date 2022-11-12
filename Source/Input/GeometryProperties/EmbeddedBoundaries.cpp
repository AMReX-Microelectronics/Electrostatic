#include "GeometryProperties.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/CodeUtils/CodeUtil.H"

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

template<typename T>
class TD;

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
    specify_separate_surf_beta = 0;
    specify_inhomogeneous_dirichlet = 0;

    amrex::ParmParse pp_ebgeom("ebgeom");
    pp_ebgeom.query("required_coarsening_level", required_coarsening_level);
    pp_ebgeom.query("max_coarsening_level", max_coarsening_level);
    pp_ebgeom.query("support", eb_support_str);
    support = map_eb_support[eb_support_str];
    pp_ebgeom.query("specify_input_using_eb2", specify_input_using_eb2);
    queryWithParser(pp_ebgeom,"specify_separate_surf_beta", specify_separate_surf_beta);
    pp_ebgeom.query("specify_inhomo_dir", specify_inhomogeneous_dirichlet);

    amrex::Print() << "\n";
    amrex::Print() << prt << "ebgeom.required_coarsening_level: " << required_coarsening_level << "\n";
    amrex::Print() << prt << "ebgeom.max_coarsening_level: " << max_coarsening_level << "\n";
    amrex::Print() << prt << "ebgeom.support: " << eb_support_str << "\n";
    amrex::Print() << prt << "ebgeom.specify_input_using_eb2: " << specify_input_using_eb2 << "\n";
    amrex::Print() << prt << "ebgeom.specify_inhomo_dir: " << specify_inhomogeneous_dirichlet << "\n";
    amrex::Print() << prt << "ebgeom.specify_separate_surf_beta: " << specify_separate_surf_beta << "\n";
    if(specify_separate_surf_beta == 0) 
    {
       getWithParser(pp_ebgeom,"surf_beta", surf_beta);
       amrex::Print() << prt << "ebgeom.surf_beta: " << surf_beta << "\n";
    }

    if(!specify_input_using_eb2) 
    {
        num_objects = 0;
        bool basic_objects_specified = pp_ebgeom.queryarr("objects", vec_object_names);
        int c=0;
        for (auto it: vec_object_names)
        {
            if (map_basic_objects_type.find(it) == map_basic_objects_type.end()) {
                amrex::ParmParse pp_object(it);

                pp_object.get("geom_type", map_basic_objects_type[it]);
                ReadObjectInfo(it, map_basic_objects_type[it], pp_object);

                if(specify_inhomogeneous_dirichlet == 1) 
                {
                    getWithParser(pp_object,"surf_soln", map_basic_objects_soln[it]);
                    amrex::Print()  << "surf_soln: " << map_basic_objects_soln[it] << "\n";
                } 

                if(specify_separate_surf_beta == 1) 
                {
                    getWithParser(pp_object,"surf_beta", map_basic_objects_beta[it]);
                    amrex::Print()  << "surf_beta: " << map_basic_objects_beta[it] << "\n";
                }

                ++c;
            }
        }
        num_objects = c;
        amrex::Print() << "\ntotal number of objects: " << num_objects << "\n";

        //amrex::Print() << "\nObject Info:\n";
//        for (auto it: map_basic_objects_info)
//        {
//            std::string name = it.first; 
//            std::string geom_type = map_basic_objects_type[name];
//            PrintObjectInfo(name, geom_type, it.second);
//        }

    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_GeometryProperties::ReadGeometry()************************\n";
#endif
}


void
c_EmbeddedBoundaries::ReadObjectInfo(std::string object_name, std::string object_type, amrex::ParmParse pp_object)
{

    amrex::Print() << "\nObject name: " << object_name << "\n";
    amrex::Print() << "Object type: " << object_type << "\n";

    switch (map_object_type_enum[object_type]) 
    {  
        case s_ObjectType::object::sphere:
        {
            amrex::Vector<amrex::Real> center;
            getArrWithParser(pp_object, "sphere_center", center, 0, AMREX_SPACEDIM);

            amrex::Print() << "Sphere center: ";
            for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << center[i] << "  ";
            amrex::Print() << "\n";

         
            amrex::Real radius;
            getWithParser(pp_object, "sphere_radius", radius);

            amrex::Print() << "Sphere radius: " << radius << "\n";


            bool has_fluid_inside;
            pp_object.get("sphere_has_fluid_inside", has_fluid_inside);

            amrex::Print() << "Sphere has_fluid_inside?: " << has_fluid_inside << "\n";


            EB2::SphereIF sphere(radius, vecToArr(center), has_fluid_inside);

            map_basic_objects_info[object_name] = sphere;  
            break;
        }
        case s_ObjectType::object::box:
        {
            amrex::Vector<amrex::Real> lo;
            getArrWithParser(pp_object,"box_lo", lo,0,AMREX_SPACEDIM);
        
            amrex::Print() << "box_lo: ";
            for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << lo[i] << "  ";
            amrex::Print() << "\n";


            amrex::Vector<amrex::Real> hi;
            getArrWithParser(pp_object,"box_hi", hi,0,AMREX_SPACEDIM);
        
            amrex::Print() << "box_hi: ";
            for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << hi[i] << "  ";
            amrex::Print() << "\n";


            bool has_fluid_inside;
            pp_object.get("has_fluid_inside", has_fluid_inside);

            amrex::Print() << "Box has_fluid_inside?: " << has_fluid_inside << "\n";


            EB2::BoxIF box(vecToArr(lo), vecToArr(hi), has_fluid_inside);

            map_basic_objects_info[object_name] = box;  
            break;
        }
        case s_ObjectType::object::cylinder:
        {
            amrex::Vector<amrex::Real> center;
            getArrWithParser(pp_object, "center", center, 0, AMREX_SPACEDIM);

            amrex::Print() << "cylinder center: ";
            for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << center[i] << "  ";
            amrex::Print() << "\n";

            amrex::Real radius;
            getWithParser(pp_object,"radius", radius);

            amrex::Print() << "cylinder radius: " << radius << "\n";

            int direction;
            pp_object.get("direction", direction);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(direction >=0 && direction < 3,
                                             "cyl_cavity_direction is invalid");

            amrex::Print() << "cylinder direction: " << direction << "\n";

            amrex::Real height=-1;
            queryWithParser(pp_object,"height", height);

            amrex::Print() << "cylinder height: " << height << "\n";
    
            bool cyl_has_fluid_inside;
            pp_object.get("has_fluid_inside", cyl_has_fluid_inside);

            EB2::CylinderIF cyl(radius, height, direction, vecToArr(center), cyl_has_fluid_inside);

            map_basic_objects_info[object_name] = cyl;  
            break;
        }
        case s_ObjectType::object::cntfet_contact:
        {
            amrex::Vector<amrex::Real> lo;
            getArrWithParser(pp_object,"box_lo", lo,0,AMREX_SPACEDIM);

            amrex::Print() << "box_lo: ";
            for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << lo[i] << "  ";
            amrex::Print() << "\n";


            amrex::Vector<amrex::Real> hi;
            getArrWithParser(pp_object,"box_hi", hi,0,AMREX_SPACEDIM);

            amrex::Print() << "box_hi: ";
            for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << hi[i] << "  ";
            amrex::Print() << "\n";


            amrex::Vector<amrex::Real> center;
            center.resize(AMREX_SPACEDIM);
            for(int idim=0; idim < AMREX_SPACEDIM; ++idim) 
            {
               center[idim] = lo[idim] + (hi[idim] - lo[idim])/2.;  
            }
            queryArrWithParser(pp_object,"cyl_cavity_center", center,0,AMREX_SPACEDIM);

            amrex::Print() << "cyl_cavity_center: ";
            for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << center[i] << "  ";
            amrex::Print() << "\n";


            amrex::Real radius;
            getWithParser(pp_object,"cyl_cavity_radius", radius);

            amrex::Print() << "cyl_cavity_radius: " << radius << "\n";


            int direction;
            pp_object.get("cyl_cavity_direction", direction);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(direction >=0 && direction < 3,
                                             "cyl_cavity_direction is invalid");

            amrex::Print() << "cyl_cavity_direction: " << direction << "\n";


            amrex::Real height = hi[direction] - lo[direction];
            queryWithParser(pp_object,"cyl_cavity_height", height);

            amrex::Print() << "cyl_cavity_height: " << height << "\n";
    
            bool box_has_fluid_inside=0;
            bool cyl_has_fluid_inside=1;

            EB2::CylinderIF cyl(radius, height, direction, vecToArr(center), cyl_has_fluid_inside);

            EB2::BoxIF box(vecToArr(lo), vecToArr(hi), box_has_fluid_inside);

            auto cntfet_contact = amrex::EB2::makeIntersection(cyl, box);
            //TD<decltype(cntfet_contact)> cntfet_contact_type;

            map_basic_objects_info[object_name] = cntfet_contact;  
            break;
        }
    }
}

void
c_EmbeddedBoundaries::BuildObjects(amrex::Geometry geom,amrex::BoxArray ba, amrex::DistributionMapping dm)
{

    Vector<int> ng_ebs = {2,2,2};
    
    pFactory.resize(num_objects);
    if(specify_inhomogeneous_dirichlet == 1)  m_p_soln_mf.resize(num_objects);
    if(specify_separate_surf_beta == 1)  m_p_beta_mf.resize(num_objects);

    int c=0;
    for (auto it: map_basic_objects_type)
    {
        std::string name = it.first; 
        std::string geom_type = it.second;
        amrex::Print() << "name: " << name << ", geom_type: " << geom_type << "\n";

        switch (map_object_type_enum[geom_type]) 
        {  
            case s_ObjectType::object::sphere:
            {
                using ObjectType = amrex::EB2::SphereIF;
                BuildSingleObject<ObjectType>(name, c, geom, ba, dm);
                break;
            }
            case s_ObjectType::object::box:
            {
                using ObjectType = amrex::EB2::BoxIF;
                BuildSingleObject<ObjectType>(name, c, geom, ba, dm);
                break;
            }
            case s_ObjectType::object::cylinder:
            {
                using ObjectType = amrex::EB2::CylinderIF;
                BuildSingleObject<ObjectType>(name, c, geom, ba, dm);
                break;
            }
            case s_ObjectType::object::cntfet_contact:
            {
                using ObjectType = cntfet_contact_type;
                BuildSingleObject<ObjectType>(name, c, geom, ba, dm);
                break;
            }
        }

        if(specify_separate_surf_beta == 1) 
        {
            m_p_beta_mf[c] = std::make_unique<amrex::MultiFab>(ba, dm, 1, 0, MFInfo(), *pFactory[c]); 
            Multifab_Manipulation::SpecifyValueOnlyOnCutcells(*m_p_beta_mf[c], map_basic_objects_beta[name]);
        }
        if(specify_inhomogeneous_dirichlet == 1) 
        {
            m_p_soln_mf[c] = std::make_unique<amrex::MultiFab>(ba, dm, 1, 0, MFInfo(), *pFactory[c]); 
            //(*m_p_soln_mf[c]).setVal(-1); 
            Multifab_Manipulation::SpecifyValueOnlyOnCutcells(*m_p_soln_mf[c], map_basic_objects_soln[name]);
        }
        ++c;
    }

    if(num_objects == 2) 
    {  
        auto name1 = vec_object_names[0];  
        auto geom_type1 = map_basic_objects_type[name1];  
        auto name2 = vec_object_names[1];
        auto geom_type2 = map_basic_objects_type[name2];  

        if( (map_object_type_enum[geom_type1] == s_ObjectType::object::cntfet_contact) && 
            (map_object_type_enum[geom_type2] == s_ObjectType::object::cntfet_contact) ) 
        {  
            using ObjectType1 = cntfet_contact_type;
            using ObjectType2 = cntfet_contact_type;

            BuildUnionObject<ObjectType1, ObjectType2>(name1, name2, geom, ba, dm);
        }
        else if ( (map_object_type_enum[geom_type1] == s_ObjectType::object::sphere) &&
                  (map_object_type_enum[geom_type2] == s_ObjectType::object::sphere) )
        {
            using ObjectType1 = amrex::EB2::SphereIF;
            using ObjectType2 = amrex::EB2::SphereIF;

            BuildUnionObject<ObjectType1, ObjectType2>(name1, name2, geom, ba, dm);
        }
        else if ( (map_object_type_enum[geom_type1] == s_ObjectType::object::box) &&
                  (map_object_type_enum[geom_type2] == s_ObjectType::object::box) )
        {
            using ObjectType1 = amrex::EB2::BoxIF;
            using ObjectType2 = amrex::EB2::BoxIF;

            BuildUnionObject<ObjectType1, ObjectType2>(name1, name2, geom, ba, dm);
        }
        else if ( (map_object_type_enum[geom_type1] == s_ObjectType::object::cylinder) &&
                  (map_object_type_enum[geom_type2] == s_ObjectType::object::cylinder) )
        {
            using ObjectType1 = amrex::EB2::CylinderIF;
            using ObjectType2 = amrex::EB2::CylinderIF;

            BuildUnionObject<ObjectType1, ObjectType2>(name1, name2, geom, ba, dm);
        }
        else 
        {
            amrex::Abort("Error: 1) For more than 1 objects, one must code the operation such as union, intersection, etc.\
                          2) At present, union operation is performed only when two geometries are of type cntfet_contact,\
                             amrex::EB2::BoxIF, amrex::EB2::SphereIF");
        }   
    }
}


template<typename ObjectType>
void
c_EmbeddedBoundaries::BuildSingleObject(std::string name, int c, amrex::Geometry geom,amrex::BoxArray ba, amrex::DistributionMapping dm)
{
    auto object = std::any_cast<ObjectType>(map_basic_objects_info[name]);
    auto gshop = amrex::EB2::makeShop(object);
    amrex::EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
    const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
    const EB2::Level& eb_level = eb_is.getLevel(geom);
    Vector<int> ng_ebs = {2,2,2};
    //pFactory.push_back(amrex::makeEBFabFactory(&eb_level, ba, dm, ng_ebs, support));
    pFactory[c] = amrex::makeEBFabFactory(&eb_level, ba, dm, ng_ebs, support);
}


template<typename ObjectType1, typename ObjectType2>
void
c_EmbeddedBoundaries::BuildUnionObject(std::string name1, std::string name2, amrex::Geometry geom,amrex::BoxArray ba, amrex::DistributionMapping dm)
{
    auto object1 = std::any_cast<ObjectType1>(map_basic_objects_info[name1]);
    auto object2 = std::any_cast<ObjectType2>(map_basic_objects_info[name2]);
    auto union_object = amrex::EB2::makeUnion(object1, object2);

    auto gshop = amrex::EB2::makeShop(union_object);
    amrex::EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
    const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
    const EB2::Level& eb_level = eb_is.getLevel(geom);
    Vector<int> ng_ebs = {2,2,2};

    pFactory.push_back(amrex::makeEBFabFactory(&eb_level, ba, dm, ng_ebs, support));
}
