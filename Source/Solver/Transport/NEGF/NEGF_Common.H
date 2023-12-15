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

struct s_Norm {
    enum Type { Absolute, Relative };
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

    virtual ~c_NEGF_Common();
    /* This is the name of the nanostructure specified by user in the input file. */
    std::string name;

    /* rank number of current process */
    int my_rank;

    /* total number of processors*/
    int num_proc;

    /* Total number of atoms. 
     * This number can be different from num_field_sites.
     * E.g. in case of mode-space approximation for a carbon nanotube,
     * we can treat each ring of a nanotube as a field site.
     * So the num_atoms are m times larger than num_field_sites, 
     * where m is the number of atoms in a ring.
     */
    int num_atoms = 0;

    /* The following two quantities are specific to material structure,
     * and are useful in defining atom locations. 
     */
    int num_atoms_per_unitcell = 0;
    int num_unitcells = 0;
    
    /*number of all sites*/
    int num_field_sites = 0; 

    /* number of field sites this processor gets
     * Note: a field site can be shared by multiple processors,
     * e.g., a ring of a carbon nanotube can be split across, 
     * say, 4 processors, yet it can be represented as 1 field site
     * in mode-space approximation.
     * In this case, each of the 4 processors contain the field site 
     * information.
     */
    int num_local_field_sites = 0; 

    /* This is the offset of the current processors local field site id 
     * with respect to global field site id,
     * e.g., if a processor contains field sites [4, 5, 6, 7]
     * their local site ids are [0, 1, 2, 3], while the site id offset is 4.
     */
    int site_id_offset = 0;

    /* avg_type can be s_AVG_TYPE::ALL or s_AVG_TYPE::SPECIFIC */
    int avg_type;

    /* The following host and device vectors represent indices of atoms to average over, while 
     * averaging data represented on atoms of a material to corresponding field site representation, 
     * e.g. data represented on  atoms of a ring of a nanotube to average data per ring.
     */
    amrex::Vector<int> vec_avg_indices;
    amrex::Gpu::DeviceVector<int> gpuvec_avg_indices;

    /* whether to average over atoms to represent field sites */
    int average_field_flag = 1;

    /* This number can be > 1, e.g. a carbon ring treated as a field site in mode-space approximation, 
     * having more than 1 atoms in a ring
     */
    int num_atoms_per_field_site = 1;

    /* This number is > 1, if field site representation is constructed by averaging over more than 1 atoms. */
    int num_atoms_to_avg_over = 1; 

    /* 1 indicates Y*/
    int primary_transport_dir = 1; 

    /* computational cell volume */
    amrex::Real cell_volume = 1.;
    amrex::RealArray offset;
   
    /* global size of the Hamiltonian (each element can be a number or a block of numbers)*/
    int Hsize_glo=0; 

    /* repeated block size of offdiagonal elements of Hamiltonian*/
    int offDiag_repeatBlkSize=0;
 
    /* maximum block columns per MPI/GPU process*/
    int max_blkCol_perProc=0; 

    /* size of the local block of columns*/
    int blkCol_size_loc=0; 

    /* Number of processors that hold nonzero number of block columns*/
    int num_proc_with_blkCol=0;

    /* vec of cumulative number of block column sizes*/
    amrex::Vector<int> vec_cumu_blkCol_size; 

    /* vec of block column global ids*/
    amrex::Vector<int> vec_blkCol_gids; 

    /* spin degeneracy */
    const int spin_degen = 2;

    /* Host and device vectors of block degeneracy vectors. 
     * This is relevant in mode-space approximation, when each block in 
     * the Hamiltonian is not a number but an array representing distinct modes
     * and elements of an array have degeneracy associated with it.
     */
    amrex::Vector<int> block_degen_vec;
    amrex::Gpu::DeviceVector<int> block_degen_gpuvec;

    /* MPI_recv_count stores blkCol_size_loc of each processor.
     * MPI_recv_disp stores the offset id of block columns for each processor.
     * These vectors are defined for each processor and are useful during
     * MPI_Gatherv of of induced charge computed by each processor.
     */
    amrex::Vector<int> MPI_recv_count;
    amrex::Vector<int> MPI_recv_disp;

    /* MPI_send_count stores num_local_field_sites of each processor.
     * MPI_send_disp stores site_id_offset of each processor.
     * These vectors are defined only for IOProcessor.
     */
    amrex::Vector<int> MPI_send_count;
    amrex::Vector<int> MPI_send_disp;

    int NS_data_offset=0; 
    /* offset in a global 1D table data of multiple structures
     * e.g. if there are 2 nanotubes, then for the first and second nanotubes
     * it is 0 and Hsize_glo of the first nanotube.
     */


    /*Tables required to store Hamiltonian */
    BlkTable1D h_minusHa_loc_data;
    BlkTable1D h_Hb_loc_data;
    BlkTable1D h_Hc_loc_data;

    BlkTable1D h_tau_glo_data;

    /*Tables required to store essential quantities and observables */

    #ifdef COMPUTE_GREENS_FUNCTION_OFFDIAG_ELEMS
    BlkTable2D h_GR_loc_data;
    #else
    BlkTable1D h_GR_loc_data;
    #endif
    #ifdef COMPUTE_SPECTRAL_FUNCTION_OFFDIAG_ELEMS
    BlkTable2D h_A_loc_data;
    #else
    BlkTable1D h_A_loc_data;
    #endif

    RealTable1D h_DOS_loc_data;
    RealTable1D h_Transmission_loc_data;
    BlkTable1D  h_Rho0_loc_data;
    BlkTable1D  h_GR_atPoles_loc_data;
    BlkTable1D  h_RhoEq_loc_data;
    BlkTable1D  h_RhoNonEq_loc_data;
    RealTable1D h_RhoInduced_loc_data;
    RealTable1D h_U_loc_data;
    amrex::Vector<amrex::Real> h_PTD_glo_vec; /*PTD: primary transport direction*/

    RealTable1D h_n_curr_in_glo_data; //temporarily used for initialization
    RealTable1D h_n_curr_in_loc_data;
    RealTable1D d_n_curr_in_loc_data;
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

    RealTable1D h_NonEq_Integrand_data;
    RealTable1D h_NonEq_Integrand_Source_data;
    RealTable1D h_NonEq_Integrand_Drain_data;
    amrex::Real E_at_max_noneq_integrand = 0;

    #ifdef AMREX_USE_GPU
    #ifdef COMPUTE_GREENS_FUNCTION_OFFDIAG_ELEMS
    BlkTable2D d_GR_loc_data;
    #else
    BlkTable1D d_GR_loc_data;
    #endif
    #ifdef COMPUTE_SPECTRAL_FUNCTION_OFFDIAG_ELEMS
    BlkTable2D d_A_loc_data;
    #else
    BlkTable1D d_A_loc_data;
    #endif
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

    RealTable1D d_NonEq_Integrand_data;
    RealTable1D d_NonEq_Integrand_Source_data;
    RealTable1D d_NonEq_Integrand_Drain_data;
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
    int total_noneq_integration_pts=0;
    amrex::Vector<amrex::Real> noneq_percent_intercuts;

    int num_noneq_paths=1;
    bool flag_compute_flatband_dos = false;
    int flatband_dos_integration_pts = 100;
    amrex::Vector<amrex::Real> flatband_dos_integration_limits = {-1., 1.};
    amrex::Real Fermi_tail_factor_lower = 14.;
    amrex::Real Fermi_tail_factor_upper = 14.;

    amrex::Real Contact_Temperature[NUM_CONTACTS] = {298., 298.};

    int flag_contact_mu_specified = 1;
    bool flag_EC_potential_updated = false;
    amrex::Real Contact_Electrochemical_Potential[NUM_CONTACTS] = {0., 0.};
    std::string Contact_Parser_String[NUM_CONTACTS];
    std::string Gate_String = "Gate";

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
    amrex::Vector<c_IntegrationPath> ContourPath_DOS;

    /*For Restart*/
    int flag_initialize_charge_distribution=0;
    std::string charge_distribution_filename;

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
    bool flag_adaptive_integration_limits = false;
    bool flag_correct_integration_limits = false;
    bool flag_write_integrand_main = false;
    bool flag_write_integrand_iter = false;
    int integrand_correction_interval = 500;
    int write_integrand_interval = 500;
    amrex::Vector<amrex::Real> kT_window_around_singularity = {1., 1.};
    

    void Set_StepFilenameString(const int step);
    void Set_IterationFilenameString(const int iter);

    amrex::Real initial_charge = 1e-3;
    int num_recursive_parts = 1;
    int Hsize_recur_part=-1; 

    void Initialize_ChargeAtFieldSites();
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
    void Generate_NonEq_Paths ();
    void Find_NonEq_Percent_Intercuts_Adaptively();
    virtual void Define_ContactInfo();
    virtual 
    AMREX_GPU_HOST_DEVICE
    void Compute_SurfaceGreensFunction (MatrixBlock<T>& gr, const ComplexType EmU);
    void Compute_DensityOfStates(std::string DOS_foldername, bool flag_write_LDOS);
    void Write_FermiFunction (const amrex::Vector<ComplexType> E_vec, std::string filename);
    void Compute_Rho0();
    void Compute_GR_atPoles();
    void Compute_RhoEq();
    void Compute_RhoNonEq();
    void Compute_RhoNonEqOld();
    void Compute_InducedCharge(RealTable1D& n_curr_out_data);
    void Compute_Current();
    //void GuessNewCharge_SimpleMixingAlg();
    
    void Gatherv_NEGFComputed_LocalCharge(RealTable1D& n_curr_out_glo_data);
    void Fetch_InputLocalCharge_FromNanostructure(RealTable1D& container_data, 
                                                  const int disp, 
                                                  const int data_size);
    void Scatterv_BroydenComputed_GlobalCharge(RealTable1D& n_curr_in_glo_data);

    void Set_Arrays_OfSize_NumFieldSites();

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

    template<typename TableType>
    void
    Write_ChargeComponents(std::string filename,
                           const TableType& h_RhoEq_data,
                           const TableType& h_RhoNonEq_data,
                           const TableType& h_Rho0_data,
                           const TableType& h_RhoInduced_data);

    template<typename TableType> 
    void Write_Integrand(const amrex::Vector<ComplexType>& Vec,
                         const TableType& Arr_Channel_data,
                         const TableType& Arr_Source_data,
                         const TableType& Arr_Drain_data,
                         std::string filename);

    void Write_BlkTable1D_asaf_E(const amrex::Vector<ComplexType>& E_vec,
                                 const BlkTable1D& Arr_data,
                                 std::string filename,
                                 std::string header);

    void Write_BlkTable2D_asaf_E(const BlkTable2D& Arr_data,
                                 std::string filename,
                                 std::string header);

    void Write_PotentialAtSites(const std::string filename_prefix);
    void Write_InputInducedCharge (const std::string filename_prefix, const RealTable1D& n_curr_in_data);
    void Write_InducedCharge (const std::string filename_prefix, 
                              const RealTable1D& n_curr_out_data);

    void Write_ChargeNorm(const std::string filename_prefix, 
                          const RealTable1D& Norm_data);

    void Write_Current (const int step, 
                        const amrex::Real Vds,
                        const amrex::Real Vgs,
                        const int Broyden_Step,
                        const int max_iter,
                        const amrex::Real Broyden_fraction,
                        const int Broyden_Scalar);

    template<typename U> 
    void Print_Table1D_loc(const U& Tab1D_data);
    template<typename U> 
    void Print_Table2D_loc(const U& Tab2D_data);

    template<typename TableType>
    void Read_Table1D(int assert_size,
                      TableType& Tab1D_data,
                      std::string filename);

    template<typename U, typename V> 
    void SetVal_Table1D (U& Tab1D_data, V val);

    template<typename U, typename V> 
    void SetVal_Table2D(U& Tab2D_data, V val);

    //virtual void Deallocate();

    void Define_MPISendCountAndDisp();
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
