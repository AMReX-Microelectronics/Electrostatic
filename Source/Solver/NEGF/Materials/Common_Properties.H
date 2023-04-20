#ifndef COMMON_PROPERTIES_H_
#define COMMON_PROPERTIES_H_

#include "../../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../../Code_Definitions.H"
#include "../../../Utils/SelectWarpXUtils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_REAL.H>
#include <AMReX.H>
#include<AMReX_TableData.H>
#include <AMReX_GpuComplex.H>

using namespace amrex;

using ComplexType = amrex::GpuComplex<amrex::Real>;

using ComplexTable1D = TableData<ComplexType, 1>;
using ComplexTable2D = TableData<ComplexType, 2>;

template<typename T>
struct MatrixBlock
{
   public:
   T block;

   MatrixBlock<T> operator=(const amrex::Real c);
   MatrixBlock<T> operator=(const ComplexType c);

   MatrixBlock<T> operator*(const amrex::Real c);
   MatrixBlock<T> operator*(const MatrixBlock<T>& rhs);
   template<typename U>
   friend MatrixBlock<U> operator*(const amrex::Real c, const MatrixBlock<U>& rhs);

   MatrixBlock<T> operator+(const ComplexType c);
   MatrixBlock<T> operator+(const MatrixBlock<T>& rhs);
   template<typename U>
   friend MatrixBlock<U> operator+(const ComplexType c, const MatrixBlock<U>& rhs);

   MatrixBlock<T> operator-(const amrex::Real c);
   MatrixBlock<T> operator-(const ComplexType c);
   MatrixBlock<T> operator-(const MatrixBlock<T>& rhs);
   template<typename U>
   friend MatrixBlock<U> operator-(const ComplexType c, const MatrixBlock<U>& rhs);

   MatrixBlock<T> operator/(const ComplexType c);
   MatrixBlock<T> operator/(const MatrixBlock<T>& rhs);

   template<typename U>
   friend std::ostream& operator<<(std::ostream& stream, const MatrixBlock<U>& rhs);
};
/*BlkTable1D<T> 1DTable of blocks*/
template<typename T>
using BlkTable1D = TableData<MatrixBlock<T>, 1>;

/*BlkTable1D<T> 2DTable of blocks*/
template<typename T>
using BlkTable2D = TableData<MatrixBlock<T>, 2>;

struct s_AVG_TYPE 
{
    enum { 
          ALL,
          SPECIFIC
    };
};


template<typename T>
class c_Common_Properties
{

public:

    std::string name;
    int num_atoms;
    int num_unitcells;
    amrex::RealArray offset;
    int avg_type;
    amrex::Vector<int> vec_avg_indices;

    const int spin_degen = 2;

    amrex::Vector<amrex::Real> avg_gatherField;

    int Hsize_glo=0; /*global size of the Hamiltonian (each element can be a number or a block of numbers)*/

    int max_blkCol_perProc=0; /*maximum block columns per MPI/GPU process*/
    int blkCol_size_loc=0; /*size of the local block of columns*/

    int num_proc_with_blkCol=0;
    amrex::Vector<int> vec_cumu_blkCol_size; /*vec of cumulative number of blk column sizes*/
    amrex::Vector<int> vec_blkCol_gids; /*vec of block column global ids*/

    virtual void ReadNanostructureProperties();
    BlkTable1D<T> h_Ha_loc_data;
    BlkTable1D<T> h_Hb_loc_data;
    BlkTable1D<T> h_Hc_loc_data;
    BlkTable2D<T> h_Sigma_glo_data;
    #ifdef AMREX_USE_GPU
    BlkTable2D<T> d_GR_loc_data;
    BlkTable2D<T> d_A_loc_data;
    #else
    BlkTable2D<T> h_GR_loc_data;
    BlkTable2D<T> h_A_loc_data;
    #endif
    amrex::GpuArray<int,NUM_CONTACTS> global_contact_index = {0, Hsize_glo-1};
    amrex::GpuArray<int,NUM_CONTACTS> contact_transmission_index = {Hsize_glo-1, 0};

    virtual void DefineMatrixPartition(int num_proc);
    virtual void AllocateArrays();
    virtual void ConstructHamiltonian ();
    virtual void AddPotentialToHamiltonian ();

    //virtual void DefineIntegrationPaths ();
    //virtual void DefineSelfEnergy();
    //virtual void DeallocateArrays();
    //virtual void ComputeChargeDensity();

    AMREX_GPU_HOST_DEVICE
    ComplexType get_Gamma(ComplexType Sigma);
    
    AMREX_GPU_HOST_DEVICE     
    ComplexType conjugate(ComplexType a);

    //AMREX_GPU_HOST_DEVICE
    //ComplexType Compute_SurfaceGreensFunction();

    //AMREX_GPU_HOST_DEVICE
    //ComplexType get_Sigma();
};

#endif
