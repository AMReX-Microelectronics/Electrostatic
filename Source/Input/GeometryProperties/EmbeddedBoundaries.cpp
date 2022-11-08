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
//        for (auto it: map_basic_objects_info)
//        {
//            std::string name = it.first; 
//            std::string geom_type = map_basic_objects_type[name];
//            PrintObjectInfo(name, geom_type, it.second);
//        }

        amrex::Vector< std::string > vec_final_objects;
        bool final_objects_specified = pp_ebgeom.queryarr("final_objects", vec_final_objects);
        c=0;
        for (auto it: vec_final_objects)
        {
            if (map_final_objects_rank.find(it) == map_final_objects_rank.end()) {
                map_final_objects_rank[it] = c;
                amrex::ParmParse pp_final_object(it);
                std::string construct_instruction_top;
                pp_final_object.get("construct",construct_instruction_top);
                //amrex::Print() << "\ntop construct instruction: " << construct_instruction_top << "\n";                
                DecodeConstructInstructionTree(construct_instruction_top);
                  
                ++c;
            }
        }
        vec_final_objects.clear();
        num_final_objects = map_final_objects_rank.size();

        amrex::Print() << "\ntotal number of final objects: " << num_final_objects << "\n";
        for (auto it: map_final_objects_rank)
        {
            amrex::Print() << "final object: " << it.second  << " name: " << it.first << "\n";
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

    amrex::Print() << "\nObject name: " << object_name << "\n";
    amrex::Print() << "Object type: " << object_type << "\n";

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

            amrex::Print() << "Sphere center: " << center[0] << "  " << center[1] << "  " << center[2] << "\n";
            amrex::Print() << "Sphere radius: " << radius << "\n";
            amrex::Print() << "Sphere has_fluid_inside?: " << has_fluid_inside << "\n";

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

            amrex::Print() << "Box lo: " << lo[0] << "  " << lo[1] << "  " << lo[2] << "\n";
            amrex::Print() << "Box hi: " << hi[0] << "  " << hi[1] << "  " << hi[2] << "\n";
            amrex::Print() << "Box has_fluid_inside?: " << has_fluid_inside << "\n";

            EB2::BoxIF box(lo, hi, has_fluid_inside);

            map_basic_objects_info[object_name] = box;  
            break;
        }
    }
}

//void
//c_EmbeddedBoundaries::PrintObjectInfo(std::string object_name, std::string object_type, std::any object_info)
//{
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
//}


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

//void TraverseString(string &str, int N)
//{
//    for (auto &ch : str) {
//        cout<< ch<< " ";
//    }
//}

// Function to reverse a string
std::string reverseStr(std::string str)
{
    int n = str.length();
    for (int i = 0; i < n / 2; i++)  std::swap(str[i], str[n - i - 1]);
    return str;
}

int findLocationVector(std::string str, std::string target, amrex::Vector<int>& loc_vec)
{
    //amrex::Print() << "string: " << str << ", string_to_find: " << target << "\n"; 
    int occurrences = 0;
    std::string::size_type pos = 0;

    while ((pos = str.find(target, pos )) != std::string::npos) 
    {
        loc_vec.push_back(pos);  
        ++ occurrences;
        pos += target.length();
    }
    return occurrences;
}


int findNumberOfOccurrencesBelowMaxIndex(std::string str, std::string target, int MaxIndex)
{
    int occurrences = 0;
    std::string::size_type pos = 1;
    while ((pos = str.find(target, pos )) != std::string::npos && pos < MaxIndex) 
    {
        ++ occurrences;
        pos += target.length();
    }
    return occurrences;
}

int findNumberOfOccurrencesAboveMaxIndex(std::string str, std::string target, int MaxIndex)
{
    int occurrences = 0;
    std::string::size_type pos = MaxIndex;
    while ((pos = str.find(target, pos )) != std::string::npos && pos < str.length()) 
    {
        ++ occurrences;
        pos += target.length();
    }
    return occurrences;
}

c_InstructionTreeNode* 
c_EmbeddedBoundaries::RecursivelyDecodeInstructions(std::string str, int level)
{
    if(str == "*") return NULL;

    c_InstructionTreeNode* root = new c_InstructionTreeNode(str);

    root->instruction = str;
    root->tree_level = level;

    std::string left_str = "*";
    std::string right_str = "*";

    root->map_basic_objects_type = &map_basic_objects_type;
    root->map_basic_objects_info = &map_basic_objects_info;
    root->map_object_type_enum = &map_object_type_enum;

    //amrex::Print() << "\ninstruction: " << str << "\n";

    amrex::Vector<int> loc_Comma;
    int total_commas = 0;
    total_commas = findLocationVector(str, ",", loc_Comma);
    //amrex::Print() << "total_commas: " << total_commas << "\n";
    //for (auto &loc : loc_Comma) amrex::Print() << " commas at: " << loc << "\n";

    for (auto &icomma: loc_Comma) 
    {
        //amrex::Print() << "icomma: "<< icomma << "\n"; 

        int num_OBs_left = findNumberOfOccurrencesBelowMaxIndex(str,"(",icomma);
        int num_CBs_left = findNumberOfOccurrencesBelowMaxIndex(str,")",icomma);
        int num_OBs_right = findNumberOfOccurrencesAboveMaxIndex(str,"(",icomma);
        int num_CBs_right = findNumberOfOccurrencesAboveMaxIndex(str,")",icomma);

        //amrex::Print() << "num_OBs_left: " << num_OBs_left << ", num_CBs_left: " << num_CBs_left << "\n";
        //amrex::Print() << "num_OBs_right: " << num_OBs_right << ", num_CBs_right: " << num_CBs_right << "\n";

        if(total_commas > 1) {
            if((num_OBs_left-1 == num_CBs_left) && (num_OBs_right == num_CBs_right-1)) 
            {
                //amrex::Print() << "num_OBs_left==num_CBs_left && num_OBs_right == num_CBs_right\n";
                left_str = str.substr(2, (icomma-1) - 2 + 1);  
                right_str = str.substr(icomma+1, str.length() - icomma -2);  
                //amrex::Print() << "left_str: " << left_str << "\n";
                //amrex::Print() << "right_str: " << right_str << "\n";
                root->operation = str.substr(0,1);
            }
        }
        else 
        {
           if((num_OBs_left==num_CBs_right) && (num_CBs_left == 0) && (num_OBs_right == 0)) 
           {
               if(num_OBs_left > 1) {
                   left_str = str.substr(2, str.length()-2-1);  
                   right_str = "*";
                   root->operation = str.substr(0,1);
               } 
               else if(num_OBs_left == 1) {
                   left_str = str.substr(2, (icomma-1) - 2 + 1);  
                   right_str = str.substr(icomma+1, str.length() - icomma -2);  
                   root->operation = str.substr(0,1);
               }
               else if (num_OBs_left == 0) {
                   left_str = str.substr(2, (icomma-1) - 2 + 1);  
                   right_str = str.substr(icomma+1, str.length() - icomma -2);  
                   root->operation = "*";
               }
           }
        }
    }

    root->left = RecursivelyDecodeInstructions(left_str, root->tree_level + 1);
    root->right = RecursivelyDecodeInstructions(right_str, root->tree_level + 1);

    return root;
}

void 
c_EmbeddedBoundaries::PrintInstructionTree(c_InstructionTreeNode* root) 
{
    if(root==NULL) return;
    amrex::Print() << "root->tree_level: " << root->tree_level << ", instruction: " << root->instruction << ", operation: " << root->operation << "\n";
    PrintInstructionTree(root->left);
    PrintInstructionTree(root->right);
}


template<class BasicObject>
BasicObject* c_InstructionTreeNode::GetABasicObject()
{
    auto op = this->operation;

    BasicObject object;
     
    auto object_name = this->instruction;
    auto object_type = (*map_basic_objects_type)[object_name];

    switch ((*map_object_type_enum)[object_type]) 
    {  
        case s_ObjectType::object::sphere:
        {
            object = std::any_cast<amrex::EB2::SphereIF>((*map_basic_objects_info)[object_name]);
            break;
        }
        case s_ObjectType::object::box:
        {
            object = std::any_cast<amrex::EB2::BoxIF>((*map_basic_objects_info)[object_name]);
            break;
        }
    }
    return &object;
}


template<typename CompoundObject, typename ObjectType1, typename ObjectType2> 
CompoundObject*
c_EmbeddedBoundaries::RecursivelyPerformObjectConstruction(c_InstructionTreeNode* root) 
{
    CompoundObject compound_object;
    ObjectType1* object1;
    ObjectType2* object2;

    if(root!=NULL) {
    
        if(root->left->operation == "*")
        {  
            object1 = root->left->GetABasicObject<ObjectType1>();
        }
        else 
        {
            object1 = RecursivelyPerformObjectConstruction<ObjectType1>(root->left);
        }
    
        if(root->right->operation == "*")
        {
            object2 = root->right->GetABasicObject<ObjectType2>();
        } 
        else 
        {
            object2 = RecursivelyPerformObjectConstruction<ObjectType2>(root->right);
        }
    
        if(root->operation == "I") 
        {
            compound_object = amrex::EB2::makeIntersection(*object1, *object2); 
        }
    }
    return &compound_object;
}


void
c_EmbeddedBoundaries::DecodeConstructInstructionTree(std::string top)
{

    amrex::Vector<int> loc_OBO;
    int TOBO=0;
    TOBO = findLocationVector(top, "(", loc_OBO);
    //amrex::Print() << "Total open bracket occurances: " << TOBO << "\n";
    //for (auto &loc : loc_OBO) amrex::Print() << " occurrence at: " << loc << "\n";
 

    amrex::Vector<int> loc_CBO;
    int TCBO=0;
    TCBO = findLocationVector(top, ")", loc_CBO);
    //amrex::Print() << "Total close bracket occurances: " << TCBO << "\n";
    //for (auto &loc : loc_CBO) amrex::Print() << " occurrence at: " << loc << "\n";
    
    if( TOBO == TCBO )
    {
       // amrex::Print() << "Assert equal number of open and close brackets. \n";
    }  
 
   int root_level = 0;
   itree = RecursivelyDecodeInstructions(top,root_level);

   amrex::Print() << "\nPrinting Instruction Tree: \n";
   PrintInstructionTree(itree);

//   auto* final_object = RecursivelyPerformObjectConstruction<amrex::EB2::IntersectionIF<{amrex::EB2::SphereIF, amrex::EB2::BoxIF},amrex::EB2::SphereIF,amrex::EB2::BoxIF>(itree);

   //auto gshop1 = amrex::EB2::makeShop(final_object);
   //amrex::EB2::Build(gshop1, geom, required_coarsening_level, max_coarsening_level);
   //const EB2::IndexSpace& eb_is1 = EB2::IndexSpace::top();
   //const EB2::Level& eb_level1 = eb_is1.getLevel(geom);
   //Vector<int> ng_ebs = {2,2,2};
   //pFactory1 = amrex::makeEBFabFactory(&eb_level1, ba, dm, ng_ebs, support);

}


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

//    TD<decltype(unionobject)> type;
//    amrex::EB2::UnionIF<amrex::EB2::IntersectionIF<amrex::EB2::SphereIF, amrex::EB2::BoxIF>, amrex::EB2::IntersectionIF<amrex::EB2::SphereIF, amrex::EB2::BoxIF>>
    amrex::EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);

    const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
    const EB2::Level& eb_level = eb_is.getLevel(geom);
    pFactory = amrex::makeEBFabFactory(&eb_level, ba, dm, ng_ebs, support);
}
