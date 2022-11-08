#ifndef EMBEDDED_BOUNDARIES_H_
#define EMBEDDED_BOUNDARIES_H_

#include "EmbeddedBoundaries_fwd.H"

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>

#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBSupport.H>


#include <string>
#include <any>


//struct
//s_Sphere
//{
//    amrex::RealArray center;
//    amrex::Real radius;
//    bool has_fluid_inside;
//};
//
//struct
//s_Box
//{
//    amrex::RealArray lo;
//    amrex::RealArray hi;
//    bool has_fluid_inside;
//};
//
struct s_ObjectType {
    enum object { sphere, box, cylinder, torus, parser, stl };
};

template <typename O>
class c_InstructionTreeNode 
{
    public:
    std::string instruction;
    std::string operation;
    int tree_level;
    O object;

    c_InstructionTreeNode* left;
    c_InstructionTreeNode* right;

    c_InstructionTreeNode(std::string str) 
    {
        instruction = str;
        operation = "*";
        tree_level=0;
        left = NULL;
        right = NULL;
    }
};


class
c_EmbeddedBoundaries
{
public:
    c_EmbeddedBoundaries();

    amrex::EBSupport support;
std::unique_ptr<amrex::EBFArrayBoxFactory> pFactory;
std::unique_ptr<amrex::EBFArrayBoxFactory> pFactory1;
std::unique_ptr<amrex::EBFArrayBoxFactory> pFactory2;
//amrex::Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > pFactory;

    void ReadGeometry();

    int required_coarsening_level;
    int max_coarsening_level;
    int specify_input_using_eb2;

    void ConstructFinalObject(std::string construct_main, amrex::Geometry geom,amrex::BoxArray ba, amrex::DistributionMapping dm);
private:
    void ReadObjectInfo(std::string object_name, std::string object_type, amrex::ParmParse pp_object);
//    void PrintObjectInfo(std::string object_name, std::string object_type, std::any object_info);

    int num_basic_objects;
    std::map<std::string,std::string> map_basic_objects_type;
    std::map<std::string,std::any> map_basic_objects_info;

    void DecodeConstructInstructionTree(std::string top);


    std::map<std::string, s_ObjectType::object>
    map_object_type_enum = {
                               {"sphere", s_ObjectType::object::sphere},
                               {"box", s_ObjectType::object::box},
                               {"cylinder", s_ObjectType::object::cylinder},
                               {"torus", s_ObjectType::object::torus},
                               {"parser", s_ObjectType::object::parser},
                               {"stl", s_ObjectType::object::stl},
                           };

    int num_final_objects;
    std::map<std::string,int> map_final_objects_rank;

    std::map<int,int> map_basic_objects;

    template <class O>
    c_InstructionTreeNode* RecursivelyDecodeInstructions(std::string str, int level);

    c_InstructionTreeNode* itree;

    std::map<std::string,amrex::EBSupport>
    map_eb_support = {
                         {"full", amrex::EBSupport::full},
             		     {"basic", amrex::EBSupport::basic},
        	    	     {"volume",amrex::EBSupport::volume},
                         {"Full", amrex::EBSupport::full},
        		         {"Basic", amrex::EBSupport::basic},
            		     {"Volume",amrex::EBSupport::volume},
                     };
};


#endif
