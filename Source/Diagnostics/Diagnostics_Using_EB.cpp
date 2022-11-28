#include "Diagnostics_Using_EB.H"

#include "../Code.H"
#include "../Input/GeometryProperties/GeometryProperties.H"
#include "../Input/MacroscopicProperties/MacroscopicProperties.H"
#include "../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../Utils/CodeUtils/CodeUtil.H"
#include "../Output/Output.H"

#include <AMReX_PlotFileUtil.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_RealBox.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBSupport.H>
#include <AMReX_EB2_IndexSpace_STL.H>


using namespace amrex;
//template<typename T>
//class TD;

void
c_Diagnostics_Using_EB::ReadEBDiagnostics()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t\t{************************c_Diagnostics_Using_EB()::ReadEBDiagnostics()************************\n";
    amrex::Print() << "\t\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif


    amrex::ParmParse pp_diag("diag");
    //setting default values
    support = EBSupport::full;
    required_coarsening_level = 0; // max amr level (at present 0)
    max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible
    std::string eb_support_str = "full";

    pp_diag.query("required_coarsening_level", required_coarsening_level);
    pp_diag.query("max_coarsening_level", max_coarsening_level);
    pp_diag.query("support", eb_support_str);

    amrex::Print() << "\n##### DIAGNOSTIC EB PROPERTIES #####\n\n";
    amrex::Print() << "##### diag.required_coarsening_level: " << required_coarsening_level << "\n";
    amrex::Print() << "##### diag.max_coarsening_level: " << max_coarsening_level << "\n";
    amrex::Print() << "##### diag.support: " << eb_support_str << "\n";

    num_objects = 0;
    bool basic_objects_specified = pp_diag.queryarr("objects", vec_object_names);
    int c=0;
    for (auto it: vec_object_names)
    {
        if (map_object_type.find(it) == map_object_type.end()) {
            amrex::ParmParse pp_object(it);

            pp_object.get("geom_type", map_object_type[it]);
            ReadEBObjectInfo(it, map_object_type[it], pp_object);

            ++c;
        }
    }
    num_objects = c;
    amrex::Print() << "\n##### total number of diagnostic objects: " << num_objects << "\n";

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t\t}************************c_Diagnostics_Using_EB()::ReadEBDiagnostics()************************\n";
#endif
}

void
c_Diagnostics_Using_EB::ReadEBObjectInfo(std::string object_name, std::string object_type, amrex::ParmParse pp_object)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t\t{************************c_Diagnostics_Using_EB()::ReadEBObjectInfo()************************\n";
    amrex::Print() << "\t\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    amrex::Print() << "\n##### Object name: " << object_name << "\n";
    amrex::Print() << "##### Object type: " << object_type << "\n";

    switch (map_object_type_enum[object_type]) 
    {  
        case s_ObjectType::object::cylinder:
        {
            s_Cylinder cyl;

            amrex::Vector<amrex::Real> center;
            getArrWithParser(pp_object, "center", center, 0, AMREX_SPACEDIM);
            cyl.center = vecToArr(center);

            amrex::Print() << "##### cylinder center: ";
            for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << cyl.center[i] << "  ";
            amrex::Print() << "\n";

            getWithParser(pp_object,"radius", cyl.radius);

            amrex::Print() << "##### cylinder radius: " << cyl.radius << "\n";

            pp_object.get("axial_direction", cyl.axial_direction);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cyl.axial_direction >=0 && cyl.axial_direction < 3,
                                             "cyl.axial_direction is invalid");

            pp_object.get("theta_reference_direction", cyl.theta_reference_direction);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cyl.theta_reference_direction >=0 && cyl.theta_reference_direction < 3,
                                             "cyl.theta_reference_direction is invalid");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cyl.axial_direction !=cyl.theta_reference_direction,
                                             "cyl.theta_reference_direction is invalid");

            amrex::Print() << "##### cylinder axial_direction: " << cyl.axial_direction << "\n";
            amrex::Print() << "##### cylinder theta_reference_direction: " << cyl.theta_reference_direction << "\n";

            queryWithParser(pp_object,"height", cyl.height);

            amrex::Print() << "##### cylinder height: " << cyl.height << "\n";

            pp_object.get("has_fluid_inside", cyl.has_fluid_inside);
            amrex::Print() << "##### cylinder has_fluid_inside: " << cyl.has_fluid_inside << "\n";

            bool varnames_specified = pp_object.queryarr("fields_to_plot", cyl.fields_to_plot);
            amrex::Print() << "##### cylinder fields_to_plot: ";
            for (auto field: cyl.fields_to_plot) amrex::Print() << field << "  ";
            amrex::Print() << "\n";
           
            map_object_info[object_name] = cyl;  
            break;

        }
        default:
        {
            amrex::Abort("geom_type " + object_type + " not supported in diagnostics.");
            break;
        }
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t\t}************************c_Diagnostics_Using_EB()::ReadEBObjectInfo()************************\n";
#endif
}


void
c_Diagnostics_Using_EB::SetGeometry(const amrex::Geometry* GEOM, 
                                    const amrex::BoxArray* BA, 
                                    const amrex::DistributionMapping* DM)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t\t{************************c_Diagnostics_Using_EB()::SetGeometry()************************\n";
    amrex::Print() << "\t\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    geom = GEOM;
    ba = BA;
    dm = DM;

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t\t}************************c_Diagnostics_Using_EB()::SetGeometry()************************\n";
#endif
}


void
c_Diagnostics_Using_EB::CreateFactory()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t\t{************************c_Diagnostics_Using_EB()::CreateFactory()************************\n";
    amrex::Print() << "\t\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    amrex::Print() <<  "Before index space size : " << amrex::EB2::IndexSpace::size() << "\n";
    for (auto it: map_object_info)
    {
        std::string name = it.first; 
        std::string object_type = map_object_type[name];

        switch (map_object_type_enum[object_type]) 
        {  
            case s_ObjectType::object::cylinder:
            {
                auto ob = std::any_cast<s_Cylinder>(map_object_info[name]); 

                using IFType = amrex::EB2::CylinderIF;
                 
                IFType object_IF(ob.radius, ob.height, ob.axial_direction, 
                                 ob.center, ob.has_fluid_inside);
                
                ObtainSingleObjectFactory<IFType>(name, object_IF);

                //amrex::EB2::IndexSpace::clear();
                break;
            }
        }
    }
    amrex::Print() << "After index space size : " << amrex::EB2::IndexSpace::size() << "\n";


#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t\t}************************c_Diagnostics_Using_EB()::CreateFactory()************************\n";
#endif
}


template<typename IFType>
void
c_Diagnostics_Using_EB::ObtainSingleObjectFactory(std::string name, IFType object_IF)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_Diagnostics_Using_EB()::ObtainSingleObjectFactory()***********************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    auto gshop = amrex::EB2::makeShop(object_IF);
    amrex::EB2::Build(gshop, *geom, required_coarsening_level, max_coarsening_level);

    const auto& eb_is = EB2::IndexSpace::top();
    const auto& eb_level = eb_is.getLevel(*geom);
    Vector<int> ng_ebs = {2,2,2};
     
    map_object_pfactory[name] = amrex::makeEBFabFactory(&eb_level, *ba, *dm, ng_ebs, support);
    //amrex::EB2::IndexSpace::erase(&EB2::IndexSpace::top());

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_Diagnostics_Using_EB()::ObtainSingleObjectFactory()************************\n";
#endif
}


void
c_Diagnostics_Using_EB::ComputeAndWriteEBDiagnostics(int step, amrex::Real time, std::string folder_diag_str)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t\t{************************c_Diagnostics_Using_EB()::ComputeAndWriteEBDiagnostics()************************\n";
    amrex::Print() << "\t\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

       auto& rCode = c_Code::GetInstance();
       auto& rMprop = rCode.get_MacroscopicProperties();
       auto& rGprop = rCode.get_GeometryProperties();
       auto& n_cell = rGprop.n_cell;

       for(auto name: vec_object_names)
       {
          std::string filename_prefix_str = folder_diag_str + name + m_plt_str;

          std::string plot_file_name = amrex::Concatenate(filename_prefix_str, step, m_plt_name_digits);

          auto object_type = map_object_type[name];

          switch (map_object_type_enum[object_type])
          {
              case s_ObjectType::object::cylinder:
              {
                  auto info = std::any_cast<s_Cylinder>(map_object_info[name]);
                  auto p_factory = map_object_pfactory[name].get();

                  std::unique_ptr<amrex::MultiFab> pFactoryMF 
                           = std::make_unique<amrex::MultiFab>(*ba, *dm, 1, 0, MFInfo(), *p_factory);

                  for(auto field_name: info.fields_to_plot)
                  {
                      pFactoryMF->setVal(0.);

                      amrex::MultiFab* pFieldMF = NULL; 
                      if ( Evaluate_TypeOf_MacroStr(field_name) == 0 )/*field belongs to MacroscopicProperties*/
                      {
                          pFieldMF = rMprop.get_p_mf(field_name);
                      }

                      Multifab_Manipulation::CopyValuesIntoAMultiFabOnCutcells(*pFactoryMF, *pFieldMF);        

                      //int total_cutcells = Multifab_Manipulation::GetTotalNumberOfCutcells(*pFactoryMF);
                      //amrex::Print() << "total_cutcells: " << total_cutcells << "\n";

                      //const int cyl_ncell_axis = n_cell[info.axial_direction];
                      //const int cyl_ncell_ring = int (total_cutcells/(cyl_ncell_axis));
                      //amrex::Print() << "cyl_ncell_axis "<< cyl_ncell_axis << ", cyl_ncell_ring: " << cyl_ncell_ring << "\n";

                      //Multifab_Manipulation::
                      // Copy_3DCartesian_To_2DAzimuthalLongitudinal(*geom,
                      //                                             *pFactoryMF, 
                      //                                             &cyl_surf_grid, 
                      //                                             cyl_ncell_ring, 
                      //                                             cyl_ncell_axis,
                      //                                             info.center,
                      //                                             info.axial_direction,
                      //                                             info.theta_reference_direction);
                      //
                      ////cyl_surf_grid(1,1)=1;
                      //amrex::Print() << "cyl_surf_grid value: (1,1) : " << cyl_surf_grid(1,1) << " (1+cyl_ncell_ring, 1): " << cyl_surf_grid(1+cyl_ncell_ring, 1) << "\n";
                      amrex::EB_WriteSingleLevelPlotfile( plot_file_name,
                                                      *pFactoryMF, {field_name},
                                                       *geom,
                                                       time, step,
                                                       "HyperCLaw-V1.1", m_default_level_prefix, "Cell");
                  }
                  break;
              }
          }
       }
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t\t}************************c_Diagnostics_Using_EB()::ComputeAndWriteEBDiagnostics()************************\n";
#endif
}
