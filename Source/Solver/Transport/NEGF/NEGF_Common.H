#ifndef NEGF_COMMON_H_
#define NEGF_COMMON_H_

#include "Matrix_Block.H"
#include "IntegrationPath.H"
#include "../../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../../Utils/SelectWarpXUtils/WarpXConst.H"

#include<AMReX_TableData.H>

#include <fstream>

struct s_AVG_TYPE 
{
    enum { 
          ALL,
          SPECIFIC
    };
};


struct s_Potential_Profile
{
    enum Type
    { 
          CONSTANT,
          LINEAR,
	  POINT_CHARGE
    };
};

struct s_Position3D
{
   amrex::Array<amrex::Real,AMREX_SPACEDIM> dir = {AMREX_D_DECL(0.,0.,0.)};
};

template<typename T>
class c_NEGF_Common: public MatrixBlock<T>, 
                     public c_IntegrationPath
{

public:
    
    using ComplexTable1D = TableData<ComplexType, 1>;
    using ComplexTable2D = TableData<ComplexType, 2>;
    using RealTable1D    = TableData<amrex::Real, 1>;
    using RealTable2D    = TableData<amrex::Real, 2>;
    using BlkTable1D     = TableData<MatrixBlock<T>, 1>;
    using BlkTable2D     = TableData<MatrixBlock<T>, 2>;

    std::string name;
    int num_atoms = 0;
    int num_atoms_per_unitcell = 0;
    int num_unitcells = 0;
    int num_field_sites = 0;
    int average_field_flag = 1;
    int num_atoms_per_field_site = 1;
    int num_atoms_to_avg_over = 1;
    int primary_transport_dir = 1; /*Y*/
    amrex::Real cell_volume = 1.;

    amrex::RealArray offset;
    int avg_type;
    amrex::Vector<int> vec_avg_indices;
    amrex::Gpu::DeviceVector<int> gpuvec_avg_indices;
   
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
    amrex::Vector<int> block_degen_vec;
    amrex::Gpu::DeviceVector<int> block_degen_gpuvec;

    amrex::Vector<amrex::Real> Potential;
    amrex::Vector<amrex::Real> PTD; /*PTD: primary transport dir*/

    /*Tables required to store Hamiltonian */
    BlkTable1D h_minusHa_loc_data;
    BlkTable1D h_Hb_loc_data;
    BlkTable1D h_Hc_loc_data;

    BlkTable1D h_tau_glo_data;

    /*Tables required to store essential quantities and observables */
    BlkTable2D h_GR_loc_data;
    BlkTable2D h_A_loc_data;

    RealTable1D h_LDOS_loc_data;
    RealTable1D h_Transmission_loc_data;
    BlkTable1D  h_Rho0_loc_data;
    BlkTable1D  h_GR_atPoles_loc_data;
    BlkTable1D  h_RhoEq_loc_data;
    BlkTable1D  h_RhoNonEq_loc_data;
    RealTable1D h_Current_loc_data;


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
    BlkTable1D h_Sigma_contact_data;
    BlkTable1D h_Fermi_contact_data;

    const int num_traces = 2;
    amrex::Gpu::HostVector<amrex::Real> h_Trace_r;
    amrex::Gpu::HostVector<amrex::Real> h_Trace_i;

    #ifdef AMREX_USE_GPU
    BlkTable2D d_GR_loc_data;
    BlkTable2D d_A_loc_data;
    BlkTable1D d_Rho0_loc_data;
    BlkTable1D d_GR_atPoles_loc_data;
    BlkTable1D d_RhoEq_loc_data;
    BlkTable1D d_RhoNonEq_loc_data;

    BlkTable1D d_Alpha_loc_data;
    BlkTable1D d_Xtil_glo_data;
    BlkTable1D d_Ytil_glo_data;
    BlkTable1D d_X_loc_data;
    BlkTable1D d_Y_loc_data;

    BlkTable1D d_Alpha_contact_data;
    BlkTable1D d_X_contact_data;
    BlkTable1D d_Y_contact_data;
    BlkTable1D d_Sigma_contact_data;
    BlkTable1D d_Fermi_contact_data;

    amrex::Gpu::DeviceVector<amrex::Real> d_Trace_r;
    amrex::Gpu::DeviceVector<amrex::Real> d_Trace_i;
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
    amrex::Vector<int> eq_integration_pts;
    amrex::Vector<int> noneq_integration_pts;
    amrex::Real Fermi_tail_factor = 20.;

    amrex::Real Contact_Temperature[NUM_CONTACTS] = {298., 298.};
    amrex::Real Contact_Electrochemical_Potential[NUM_CONTACTS] = {0., 0.};

    int flag_impose_potential = 0;
    std::string potential_profile_type_str = "constant";
    std::map<std::string, s_Potential_Profile::Type>
    map_PotentialProfile =
        {
           {"constant"        , s_Potential_Profile::Type::CONSTANT},
           {"linear"          , s_Potential_Profile::Type::LINEAR},
           {"point_charge"    , s_Potential_Profile::Type::POINT_CHARGE}
        };

    amrex::Real kT_contact[NUM_CONTACTS] = {PhysConst::kb_eVperK*298.,PhysConst::kb_eVperK*298.};
    amrex::Real kT_min=PhysConst::kb_eVperK*298.;
    amrex::Real kT_max=PhysConst::kb_eVperK*298.;

    amrex::Real mu_contact[NUM_CONTACTS] = {0.,0.};
    amrex::Real mu_min=0.;
    amrex::Real mu_max=0.;

    amrex::Vector<c_IntegrationPath> ContourPath_Rho0;  
    amrex::Vector<c_IntegrationPath> ContourPath_RhoEq;  
    amrex::Vector<c_IntegrationPath> ContourPath_RhoNonEq;
    c_IntegrationPath ContourPath_DOS;


    /*For Broyden*/
    int Broyden_Step = 1;
    amrex::Real Broyden_fraction = 0.1;
    amrex::Real Broyden_Original_Fraction = 0.1;
    amrex::Real Broyden_Norm = 1.0;
    int Broyden_Correction_Step = 0;
    int Broyden_NormSumIsIncreasing_Step = 0;
    int Broyden_Reset_Step = 0;
    
    int Broyden_Threshold_MaxStep = 400;
    int Broyden_Threshold_NormSumIncreaseStep = 4;
    int Broyden_Threshold_CorrectionStep = 2;
    amrex::Real Broyden_Threshold_MinFraction = 5.e-4;
    amrex::Real Broyden_Fraction_Decrease_Factor = 5.;
    amrex::Real Broyden_Scalar_Decrease_Factor   = 2.;

    amrex::Real Broyden_Scalar = 1.0;
    amrex::Real Broyden_NormSum_Curr = 1.e10;
    amrex::Real Broyden_NormSum_Prev = 1.e100;
    RealTable1D h_n_curr_in_data;
    RealTable1D d_n_curr_in_data;

    RealTable1D n_start_in_data;
    RealTable1D n_curr_out_data;
    RealTable1D n_prev_in_data;
    RealTable1D F_curr_data;
    RealTable1D Norm_data;

    RealTable1D delta_F_curr_data;
    RealTable1D delta_n_curr_data;

    RealTable2D Jinv_curr_data;

    amrex::Vector<RealTable1D*> W_Broyden;
    amrex::Vector<RealTable1D*> V_Broyden;


    void Set_Broyden();
    void Reset_Broyden();
    void Restart_Broyden();

    /*For Output*/
    const int negf_plt_name_digits = 4;
    std::string step_foldername_str;
    std::string step_filename_prefix_str;
    std::string current_filename_str; /*current here means charge current, I */
    std::ofstream outfile_I;

    std::string step_filename_str;
    std::string iter_foldername_str;
    std::string iter_filename_prefix_str;
    std::string iter_filename_str;
    int write_at_iter = 0;

    void Set_StepFilenameString(const int step);
    void Set_IterationFilenameString();


    void Initialize_ChargeAtFieldSites(const amrex::Real value);
    virtual void ReadNanostructureProperties();
    virtual int get_Hsize() {return 0;}
    virtual int get_offDiag_repeatBlkSize() {return 0;}
    virtual void Set_Material_Specific_Parameters();
    virtual void Generate_AtomLocations(amrex::Vector<s_Position3D>& pos);
    virtual void Define_PotentialProfile();

    void DefineMatrixPartition();

    MPI_Datatype MPI_BlkType;
    virtual void Define_MPI_BlkType ();
    virtual void AllocateArrays();
    virtual void ConstructHamiltonian ();
    void AddPotentialToHamiltonian ();
    void Update_ContactPotential ();
    void Update_ContactElectrochemicalPotential ();
    virtual void Define_EnergyLimits ();
    virtual void Define_IntegrationPaths ();
    void Update_IntegrationPaths ();
    virtual void Define_ContactInfo();
    virtual 
    AMREX_GPU_HOST_DEVICE
    void Compute_SurfaceGreensFunction (MatrixBlock<T>& gr, const ComplexType EmU);
    void Compute_DensityOfStates();
    void Compute_Rho0();
    void Compute_GR_atPoles();
    void Compute_RhoEq();
    void Compute_RhoNonEq();
    void Compute_InducedCharge();
    void Compute_Current();
    //void Compute_ChargeDensity();
    void GuessNewCharge_SimpleMixingAlg();
    void GuessNewCharge_ModifiedBroydenSecondAlg();
    void GuessNewCharge_ModifiedBroydenSecondAlg_WithCorrection();
    void GuessNewCharge_Broyden_FirstAlg();
    

    void Set_Arrays_OfSize_NumFieldSites();

    //virtual void ComputeChargeDensity();
    //virtual void DeallocateArrays();

    AMREX_GPU_HOST_DEVICE
    ComplexType FermiFunction(ComplexType E_minus_Mu, const amrex::Real kT);

    void Allocate_TemporaryArraysForGFComputation();
    void Deallocate_TemporaryArraysForGFComputation();

    void get_Sigma_at_contacts(BlkTable1D& h_Sigma_contact_data, 
                               ComplexType E);
   
    //AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
    //MatrixBlock<T> get_Gamma(const MatrixBlock<T>& Sigma);

    ComplexTable1D h_E_RealPath_data;

    template<typename VectorType, typename TableType> 
    void Write_Table1D(const amrex::Vector<VectorType>& Vec,
                       const TableType& Arr_data,
                       std::string filename,
                       std::string header);

    void Write_BlkTable1D_asaf_E(const amrex::Vector<ComplexType>& E_vec,
                                 const BlkTable1D& Arr_data,
                                 std::string filename,
                                 std::string header);

    void Write_BlkTable2D_asaf_E(const BlkTable2D& Arr_data,
                                 std::string filename,
                                 std::string header);

    void Write_InducedCharge (const std::string filename_prefix);
    void Write_ChargeNorm(const std::string filename_prefix);
    void Write_Current (const int step);

    template<typename U> 
    void Print_Table1D_loc(const U& Tab1D_data);
    template<typename U> 
    void Print_Table2D_loc(const U& Tab2D_data);

    template<typename U, typename V> 
    void SetVal_Table1D (U& Tab1D_data, V val);

    template<typename U, typename V> 
    void SetVal_Table2D(U& Tab2D_data, V val);

    virtual void Deallocate();

};


template<typename T> 
template<typename U, typename V> 
void 
c_NEGF_Common<T>:: SetVal_Table1D (U& Tab1D_data, V val)
{
    auto tlo = Tab1D_data.lo();
    auto thi = Tab1D_data.hi();
    
    auto const& Tab1D = Tab1D_data.table();
    
    for (int i = tlo[0]; i < thi[0]; ++i)
    {
        Tab1D(i) = val;   
    }
}


template<typename T> 
template<typename U, typename V> 
void 
c_NEGF_Common<T>:: SetVal_Table2D (U& Tab2D_data, V val)
{
    auto tlo = Tab2D_data.lo();
    auto thi = Tab2D_data.hi();

    auto const& Tab2D = Tab2D_data.table();

    for (int i = tlo[0]; i < thi[0]; ++i)
    {
        for (int j = tlo[1]; j < thi[1]; ++j) //slow moving index. printing slow
        {
            Tab2D(i,j) = val;
        }
    }
}

#endif
