#include "Common_Properties.H"
//#include <AMReX_Particles.H>
//
#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"

amrex::Vector<int> vec_avg_indices;

void
c_Common_Properties:: ReadNanostructureProperties ()
{

    amrex::ParmParse pp_ns(name);

    getWithParser(pp_ns,"num_unitcells", num_unitcells);
    amrex::Print() << "##### num_unitcells: " << num_unitcells << "\n";

    amrex::Vector<amrex::Real> vec_offset;
    getArrWithParser(pp_ns, "offset", vec_offset, 0, AMREX_SPACEDIM);
    offset = vecToArr(vec_offset);

    amrex::Print() << "##### offset: ";
    for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << offset[i] << "  ";
    amrex::Print() << "\n";

    std::string avg_type_str = "all";
    pp_ns.query("field_averaging_type", avg_type_str);
    if(avg_type_str == "all" or avg_type_str == "All" or avg_type_str == "ALL")
    {
        avg_type = s_AVG_TYPE::ALL;
    }
    if(avg_type_str == "specific" or avg_type_str == "Specific" or avg_type_str == "SPECIFIC")
    {
        avg_type = s_AVG_TYPE::SPECIFIC;
    }
    amrex::Print() << "##### field_averaging_type: " << avg_type_str << " enum id: " << avg_type << "\n";
   
    if(avg_type == s_AVG_TYPE::SPECIFIC) {
       pp_ns.getarr("atom_indices_for_averaging", vec_avg_indices);
       amrex::Print() << "##### atom_indices_for_averaging: \n";
       for (int i=0; i<vec_avg_indices.size(); ++i) 
       {
           amrex::Print() << vec_avg_indices[i] << "  ";
       }
    }
    
}
