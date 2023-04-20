#include "Common_Properties.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"

#include <cmath>
#include <math.h>
#include<stdlib.h>

/*Explicit specializations*/
template class c_Common_Properties<ComplexType[NUM_MODES]>; //of c_CNT
template class c_Common_Properties<ComplexType[NUM_MODES][NUM_MODES]>; //c_Graphene
//template struct MatrixBlock<ComplexType[NUM_MODES]>;

template<typename T>
void
c_Common_Properties<T>:: ReadNanostructureProperties ()
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


template<typename T>
void 
c_Common_Properties<T>::DefineMatrixPartition(int num_proc) {}


template<typename T>
void 
c_Common_Properties<T>::AllocateArrays () 
{
    h_Ha_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());

    h_Sigma_glo_data.resize({0,0},{NUM_CONTACTS,NUM_ENERGY_PTS_REAL},The_Pinned_Arena());

    #if AMREX_USE_GPU
    d_GR_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    d_A_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    #else
    h_GR_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    h_A_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    #endif
}


template<typename T>
void 
c_Common_Properties<T>:: ConstructHamiltonian () {}


template<typename T>
void 
c_Common_Properties<T>:: AddPotentialToHamiltonian () 
{
}

//template<typename T>
//void 
//c_Common_Properties<T>:: DefineIntegrationPaths ()
//{
//
//
//}
//
//
//template<typename T>
//void 
//c_Common_Properties<T>:: DefineSelfEnergy ()
//{
//
//
//}
//
//
//template<typename T>
//void 
//c_Common_Properties<T>::DeallocateArrays () 
//{
//
//}
//
//
//template<typename T>
//void 
//c_Common_Properties<T>::ComputeChargeDensity () 
//{
//
//}


AMREX_GPU_HOST_DEVICE
ComplexType get_Gamma(ComplexType Sigma)
{
   ComplexType val(-2.*Sigma.imag(), 0.);
   return val;
}


AMREX_GPU_HOST_DEVICE
ComplexType conjugate(ComplexType a)
{
   ComplexType a_conj(a.real(), -1.*a.imag());
   return a_conj;
}


//template<typename T>
//MatrixBlock<T>
//MatrixBlock<T>::operator*(const amrex::Real m)
//{
//   MatrixBlock<T> result;
//   for(int i=0; i < NUM_MODES; ++i)
//   {
//       result.block[i] = this->block[i]*m;
//   }
//   return result;
//   /*if S = A*1; In this case &A is this, 1 is m, and S is the returned value*/
//}
//
//template<typename T>
//MatrixBlock<T>
//MatrixBlock<T>::operator*(MatrixBlock<T>& rhs)
//{
//   MatrixBlock<T> result;
//   return result;
//}
