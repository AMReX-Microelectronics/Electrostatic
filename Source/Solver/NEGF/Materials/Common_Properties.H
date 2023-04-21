#ifndef COMMON_PROPERTIES_H_
#define COMMON_PROPERTIES_H_

#include "Matrix_Block.H"
#include "../../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../../Utils/SelectWarpXUtils/WarpXConst.H"

#include<AMReX_TableData.H>

struct s_AVG_TYPE 
{
    enum { 
          ALL,
          SPECIFIC
    };
};


template<typename T>
class c_Common_Properties: public MatrixBlock<T>
{

public:

    using ComplexTable1D = TableData<ComplexType, 1>;
    using ComplexTable2D = TableData<ComplexType, 2>;

    //template<typename U>
    using BlkTable1D = TableData<MatrixBlock<T>, 1>;
    //template<typename U>
    using BlkTable2D = TableData<MatrixBlock<T>, 2>;

    std::string name;
    int num_atoms;
    int num_unitcells;
    amrex::RealArray offset;
    int avg_type;
    amrex::Vector<int> vec_avg_indices;


    int num_proc;
    int my_rank;
    int Hsize_glo=0; 
    /*global size of the Hamiltonian (each element can be a number or a block of numbers)*/
    int offDiag_repeatBlkSize=0;
    /*repeated block size of offdiagonal elements of Hamiltonian*/
 
    int max_blkCol_perProc=0; /*maximum block columns per MPI/GPU process*/
    int blkCol_size_loc=0; /*size of the local block of columns*/

    int num_proc_with_blkCol=0;
    amrex::Vector<int> vec_cumu_blkCol_size; /*vec of cumulative number of blk column sizes*/
    amrex::Vector<int> vec_blkCol_gids; /*vec of block column global ids*/

    const int spin_degen = 2;
    amrex::Vector<amrex::Real> avg_gatherField;

    virtual void ReadNanostructureProperties();
    BlkTable1D h_Ha_loc_data;
    BlkTable1D h_Hb_loc_data;
    BlkTable1D h_Hc_loc_data;

    BlkTable1D h_tau_glo_data;
    BlkTable2D h_Sigma_glo_data;
    #ifdef AMREX_USE_GPU
    BlkTable2D d_GR_loc_data;
    BlkTable2D d_A_loc_data;
    #else
    BlkTable2D h_GR_loc_data;
    BlkTable2D h_A_loc_data;
    #endif
    amrex::GpuArray<int,NUM_CONTACTS> global_contact_index = {0, Hsize_glo-1};
    amrex::GpuArray<int,NUM_CONTACTS> contact_transmission_index = {Hsize_glo-1, 0};
    amrex::Real U_contact[NUM_CONTACTS] = {0.,0.};

    virtual int get_Hsize() {return 0;}
    virtual int get_offDiag_repeatBlkSize() {return 0;}
    virtual void DefineMatrixPartition();
    virtual void AllocateArrays();
    virtual void ConstructHamiltonian ();
    void AddPotentialToHamiltonian ();
    void Update_ContactPotential ();
    //virtual void DefineIntegrationPaths ();
    virtual void Define_tau();
    virtual void Compute_SurfaceGreensFunction (MatrixBlock<T> gr, const ComplexType EmU);
    virtual void Define_SelfEnergy();
    //virtual void DeallocateArrays();
    //virtual void ComputeChargeDensity();

    AMREX_GPU_HOST_DEVICE
    ComplexType get_Gamma(ComplexType Sigma);
    
    AMREX_GPU_HOST_DEVICE     
    ComplexType conjugate(ComplexType a);

    ComplexTable1D h_E_RealPath_data;

    void Write_BlkTable1D_asaf_E(const BlkTable1D& Arr_data,
                                 std::string filename,
                                 std::string header);

    void Write_BlkTable2D_asaf_E(const BlkTable2D& Arr_data,
                                 std::string filename,
                                 std::string header);

    void Print_BlkTable1D_loc(const BlkTable1D& Tab1D_data);
    void Print_BlkTable2D_loc(const BlkTable2D& Tab2D_data);
    void SetZero_BlkTable1D(BlkTable1D& Tab1D_data);
    void SetZero_BlkTable2D(BlkTable2D& Tab2D_data);
    //AMREX_GPU_HOST_DEVICE
    //ComplexType Compute_SurfaceGreensFunction();

    //AMREX_GPU_HOST_DEVICE
    //ComplexType get_Sigma();
};

#endif
