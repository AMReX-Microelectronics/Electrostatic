#ifndef NEGF_COMMON_H_
#define NEGF_COMMON_H_

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
class c_NEGF_Common: public MatrixBlock<T>
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

    amrex::Vector<int> MPI_recv_count;
    amrex::Vector<int> MPI_disp;

    amrex::Vector<int> vec_cumu_blkCol_size; /*vec of cumulative number of blk column sizes*/
    amrex::Vector<int> vec_blkCol_gids; /*vec of block column global ids*/

    const int spin_degen = 2;
    amrex::Vector<amrex::Real> avg_gatherField;


    /*Tables required to store Hamiltonian */
    BlkTable1D h_Ha_loc_data;
    BlkTable1D h_Hb_loc_data;
    BlkTable1D h_Hc_loc_data;

    BlkTable1D h_tau_glo_data;

    /*Tables required to store essential quantities and observables */
    BlkTable2D h_Sigma_glo_data;
    #ifdef AMREX_USE_GPU
    BlkTable2D d_GR_loc_data;
    BlkTable2D d_A_loc_data;
    #else
    BlkTable2D h_GR_loc_data;
    BlkTable2D h_A_loc_data;
    #endif
    BlkTable1D h_Rho0_loc_data;

    /*Tables required to compute Green and Spectral Functions */
    BlkTable1D h_Alpha_loc_data;
    BlkTable1D h_Alpha_glo_data;
    BlkTable1D h_Xtil_glo_data;
    BlkTable1D h_Ytil_glo_data;
    BlkTable1D h_X_glo_data;
    BlkTable1D h_Y_glo_data;
    BlkTable1D h_X_loc_data;
    BlkTable1D h_Y_loc_data;

    BlkTable1D h_Alpha_contact_data;
    BlkTable1D h_X_contact_data;
    BlkTable1D h_Y_contact_data;
    //BlkTable1D h_Adiag_loc_data;

    #ifdef AMREX_USE_GPU
    BlkTable1D d_Alpha_loc_data;
    BlkTable1D d_X_loc_data;
    BlkTable1D d_Y_loc_data;
    BlkTable1D d_Xtil_glo_data;
    BlkTable1D d_Ytil_glo_data;

    BlkTable1D d_Alpha_contact_data;
    BlkTable1D d_X_contact_data;
    BlkTable1D d_Y_contact_data;
    //BlkTable1D d_Adiag_loc_data;
    #endif





    amrex::GpuArray<int,NUM_CONTACTS> global_contact_index; 
    amrex::GpuArray<int,NUM_CONTACTS> contact_transmission_index;
    amrex::Real U_contact[NUM_CONTACTS] = {0.,0.};

    amrex::Real E_U_contact[NUM_CONTACTS] = {0.,0.};

    bool flag_noneq_exists = false;
    ComplexType E_zPlus = 0.;
    amrex::Real E_valence_min = -14;
    amrex::Real E_pole_max = 10.;
    ComplexType E_contour_left = 0.;
    ComplexType E_contour_right = 0.;
    ComplexType E_rightmost = 0.;
    ComplexType E_zeta = 0.;
    ComplexType E_eta = 0.;
    amrex::Real E_f = 0.;
    amrex::Vector<ComplexType> E_f_vec;
    int num_enclosed_poles = 0.;
    amrex::Vector<ComplexType> E_poles_vec;
    amrex::Real Fermi_tail_factor = 14.;

    amrex::Real kT_contact[NUM_CONTACTS] = {PhysConst::kb_eVperK*298.,PhysConst::kb_eVperK*298.};
    amrex::Real kT_min=PhysConst::kb_eVperK*298.;
    amrex::Real kT_max=PhysConst::kb_eVperK*298.;

    amrex::Real mu_contact[NUM_CONTACTS] = {0.,0.};
    amrex::Real mu_min=0.;
    amrex::Real mu_max=0.;

    void GaussLegendre(int degree, amrex::Vector<amrex::Real>& x, 
                                   amrex::Vector<amrex::Real>& weight);

    class c_IntegrationPath
    {
        private:
        ComplexType E_min  = 0.;
        ComplexType E_max  = 0.;
        int type_id     = 0; /*0=line, 1=circle*/
        //ComplexType E_center  = 0.; 
        //amrex::Real E_radius  = 0.;
        
        public:
        int num_pts = 0;
        amrex::Vector<ComplexType> weight_vec;
        amrex::Vector<ComplexType> mul_factor_vec;
        amrex::Vector<ComplexType> E_vec;

        //void Define_GaussLegendrePoints();
        void Define_RegularPoints(ComplexType min, ComplexType max, int pts) 
        { 
            E_min = min; 
            E_max = max;
            num_pts = pts;

            ComplexType dE = (E_max-E_min)/static_cast<amrex::Real>(num_pts-1);

            weight_vec.resize(num_pts);
            mul_factor_vec.resize(num_pts);   
            E_vec.resize(num_pts);

            for (int i=0; i<num_pts; ++i) 
            {
                E_vec[i] = E_min + static_cast<amrex::Real>(i)*dE;  
                weight_vec[i] = dE;
                mul_factor_vec[i] = 1.;
            }
        }
    };



    amrex::GpuArray<c_IntegrationPath,3> ContourPath_Rho0;  
    amrex::GpuArray<c_IntegrationPath,3> ContourPath_RhoEq;  
    c_IntegrationPath ContourPath_RhoEq_Direct;
    c_IntegrationPath ContourPath_RhoNonEq;

    virtual void ReadNanostructureProperties();
    virtual int get_Hsize() {return 0;}
    virtual int get_offDiag_repeatBlkSize() {return 0;}
    void DefineMatrixPartition();

    MPI_Datatype MPI_BlkType;
    virtual void Define_MPI_BlkType ();
    virtual void AllocateArrays();
    virtual void ConstructHamiltonian ();
    void AddPotentialToHamiltonian ();
    void Update_ContactPotential ();
    virtual void Define_EnergyLimits ();
    virtual void Define_IntegrationPaths ();
    virtual void Define_ContactInfo();
    virtual void Compute_SurfaceGreensFunction (MatrixBlock<T>& gr, const ComplexType EmU);
    virtual void Define_SelfEnergy();
    void Compute_DensityOfStates();
    //virtual void DeallocateArrays();
    //virtual void ComputeChargeDensity();

    void Allocate_TemporaryArraysForGFComputation();
    void Deallocate_TemporaryArraysForGFComputation();

    AMREX_GPU_HOST_DEVICE
    void get_Sigma_at_contacts(amrex::GpuArray<MatrixBlock<T>,NUM_CONTACTS>& Sigma, 
                               ComplexType E);

    AMREX_GPU_HOST_DEVICE
    MatrixBlock<T> get_Gamma(const MatrixBlock<T>& Sigma);
    
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
