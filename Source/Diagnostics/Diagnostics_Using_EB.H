#ifndef DIAGNOSTICS_USING_EB_H_
#define DIAGNOSTICS_USING_EB_H_

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBSupport.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <any>
#include <string>

#include "../Code_Definitions.H"
#include "../Input/GeometryProperties/EmbeddedBoundaries.H"
#include "../Output/Field_Plot_Essentials.H"
#include "Diagnostics_Using_EB_fwd.H"

struct s_Cylinder
{
    amrex::RealArray center;
    amrex::Real radius;
    int axial_direction;
    int theta_reference_direction;
    amrex::Real height = -1;
    int has_fluid_inside = 0;
    int num_params_plot_single_level;
    std::map<std::string, int> map_param_all;
};

struct s_Plane
{
    int direction;
    amrex::Real location;
    int num_params_plot_single_level;
    std::map<std::string, int> map_param_all;
};

class c_Diagnostics_Using_EB : private c_Field_Plot_Essentials
{
    void Define_SliceGeometry(amrex::Geometry &sg, const int dir);

   protected:
    const amrex::Geometry *geom;
    const amrex::BoxArray *ba;
    const amrex::DistributionMapping *dm;
    void SetGeometry(const amrex::Geometry *GEOM, const amrex::BoxArray *BA,
                     const amrex::DistributionMapping *DM,
                     int specify_using_eb);
    amrex::EBSupport support;
    int required_coarsening_level;
    int max_coarsening_level;
    int num_objects;
    amrex::Vector<std::string> vec_object_names;

    std::map<std::string, std::string> map_object_type;
    std::map<std::string, std::any> map_object_info;
    std::map<std::string, std::unique_ptr<amrex::EBFArrayBoxFactory>>
        map_object_pfactory;

    std::map<std::string, s_ObjectType::object> map_object_type_enum = {
        //{"sphere"         , s_ObjectType::object::sphere},
        //{"box"            , s_ObjectType::object::box},
        {"cylinder", s_ObjectType::object::cylinder},
        {"plane", s_ObjectType::object::plane},
        //{"torus"          , s_ObjectType::object::torus},
        //{"parser"         , s_ObjectType::object::parser},
        //{"stl"            , s_ObjectType::object::stl},
        //{"cntfet_contact" , s_ObjectType::object::cntfet_contact},

        //{"Sphere"         , s_ObjectType::object::sphere},
        //{"Box"            , s_ObjectType::object::box},
        {"Cylinder", s_ObjectType::object::cylinder},
        {"Plane", s_ObjectType::object::plane}
        //{"Torus"          , s_ObjectType::object::torus},
        //{"Parser"         , s_ObjectType::object::parser},
        //{"Stl"            , s_ObjectType::object::stl},
        //{"Cntfet_contact" , s_ObjectType::object::cntfet_contact},

        //{"STL"            , s_ObjectType::object::stl},
        //{"CNTFET_contact" , s_ObjectType::object::cntfet_contact},
    };

    std::map<std::string, amrex::EBSupport> map_eb_support = {
        {"full", amrex::EBSupport::full},
        {"basic", amrex::EBSupport::basic},
        {"volume", amrex::EBSupport::volume},
        {"Full", amrex::EBSupport::full},
        {"Basic", amrex::EBSupport::basic},
        {"Volume", amrex::EBSupport::volume},
    };

    void ReadEBDiagnostics();
    void ReadEBObjectInfo(std::string object_name, std::string object_type,
                          amrex::ParmParse pp_object);

    template <typename IFType>
    void ObtainSingleObjectFactory(std::string name, IFType object_IF);

    void CreateFactory();

    void ComputeAndWriteEBDiagnostics(int step, amrex::Real time);
};
#endif
