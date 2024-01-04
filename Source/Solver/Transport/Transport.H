#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "Transport_fwd.H"
#include "Nanostructure.H"

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>

#include<AMReX_TableData.H>

#include <string>

struct s_NS {
    enum Type { CNT, Graphene, Silicon };
};

struct s_Algorithm { 
    enum Type { broyden_first, broyden_second, simple_mixing };
};

class 
c_TransportSolver
{

public:
    amrex::Vector<std::unique_ptr<c_Nanostructure<c_CNT>>> vp_CNT;
    amrex::Vector<std::unique_ptr<c_Nanostructure<c_Graphene>>> vp_Graphene;
    //amrex::Vector<std::unique_ptr<c_Nanostructure<c_Silicon>>> vp_Silicon;

    c_TransportSolver();
    ~c_TransportSolver();

    void InitData();
    void Solve(const int step, const amrex::Real time);
    
    #ifdef BROYDEN_PARALLEL
    void Set_Broyden_Parallel();
    void Reset_Broyden_Parallel();
    void Execute_Broyden_Modified_Second_Algorithm_Parallel();
    void Execute_Broyden_Modified_Second_Algorithm_Parallel_SkipGPU();
    #endif

private:
    using RealTable1D    = TableData<amrex::Real, 1>;
    using RealTable2D    = TableData<amrex::Real, 2>;

    int num_NS=0;
    int m_step=0;
    amrex::Geometry  dummy_geom;
    amrex::BoxArray  dummy_ba;
    amrex::DistributionMapping  dummy_dm;
    amrex::Real total_max_time_across_all_steps[6] = {0., 0., 0., 0., 0., 0.};
    int total_iter = 0;

    amrex::Geometry const * _geom = &dummy_geom;
    amrex::BoxArray const * _ba = &dummy_ba;
    amrex::DistributionMapping const * _dm = &dummy_dm;
    //const amrex::GpuArray<int, AMREX_SPACEDIM> * _n_cell;
    amrex::Vector< std::string > vec_NS_names;
   
    amrex::Real Vds = 0;
    amrex::Real Vgs = 0;
    int max_iter = 1;

    void ReadData();

    std::map<std::string, s_NS::Type>
    map_NSType_enum = 
        {
           {"carbon nanotube" , s_NS::Type::CNT},
           {"Carbon Nanotube" , s_NS::Type::CNT},
           {"CNT"             , s_NS::Type::CNT},
           {"graphene"        , s_NS::Type::Graphene},
           {"Graphene"        , s_NS::Type::Graphene},
           {"silicon"         , s_NS::Type::Silicon},
           {"Silicon"         , s_NS::Type::Silicon}
        };

    std::string NS_gather_field_str;
    std::string NS_deposit_field_str;
    amrex::Real NS_initial_deposit_value;

    int use_selfconsistent_potential;
    int use_negf;
    const int negf_plt_name_digits = 4;

    std::string negf_foldername_str = "output/negf";
    std::string common_foldername_str = "output/negf/transport_common";
    std::string common_step_folder_str;
    std::string inverse_jacobian_filename;

    void Set_CommonStepFolder(const int step);

    int flag_reset_with_previous_charge_distribution = 0;
    int flag_initialize_inverse_jacobian = 0;
    bool flag_write_LDOS = false;
    bool flag_write_LDOS_iter = false;
    int  write_LDOS_iter_period = 1e6;
    /*For Broyden*/
    std::string Algorithm_Type = "broyden_second";
    std::map<std::string, s_Algorithm::Type>
    map_AlgorithmType =
        {
           {"broyden_first"        , s_Algorithm::Type::broyden_first},
           {"Broyden_first"        , s_Algorithm::Type::broyden_first},
           {"broyden_second"        , s_Algorithm::Type::broyden_second},
           {"Broyden_second"        , s_Algorithm::Type::broyden_second},
           {"simple_mixing"         , s_Algorithm::Type::simple_mixing},
        };
    
    std::string Broyden_Norm_Type = "relative";
    std::map<std::string, s_Norm::Type>
    map_NormType =
        {
           {"absolute"        , s_Norm::Type::Absolute},
           {"relative"        , s_Norm::Type::Relative},
           {"Absolute"        , s_Norm::Type::Absolute},
           {"Relative"        , s_Norm::Type::Relative}
        };

    int         num_field_sites_all_NS = 0;
    int         Broyden_Step = 1;
    amrex::Real Broyden_fraction = 0.1;
    amrex::Real Broyden_max_norm = 1.e-5;
    amrex::Real Broyden_Original_Fraction = 0.1;
    amrex::Real Broyden_Scalar = 1.0;
    amrex::Real Broyden_NormSum_Prev = 1.e100;
    int         Broyden_Threshold_MaxStep = 200;
    amrex::Real Broyden_Norm = 0.;
    amrex::Real Broyden_NormSum_Curr = 0.;

    int total_proc;
    int my_rank;
    int num_procs_with_sites;
    int site_size_loc_all_NS = 0;
    amrex::Vector<int> site_size_loc_cumulative;
    //std::map<std::string, amrex::Vector<int>> MPI_recv_count;
    //std::map<std::string, amrex::Vector<int>> MPI_recv_disp;
    //amrex::Vector<int> MPI_recv_count;
    //amrex::Vector<int> MPI_recv_disp;

    RealTable1D h_n_curr_in_data; 
    RealTable1D h_n_curr_out_data; 
    RealTable1D h_n_prev_in_data;
    RealTable1D h_F_curr_data;
    RealTable1D h_delta_F_curr_data;
    RealTable1D h_Norm_data;

    #ifdef BROYDEN_PARALLEL
    RealTable1D n_curr_in_glo_data;
    RealTable1D h_intermed_vector_data;
    RealTable1D n_curr_out_glo_data;
    RealTable1D Norm_glo_data;
    amrex::Gpu::HostVector<amrex::Real> h_Intermed_values_vec = {0., 0., 0.}; 
    /*stores Broyden_Norm, Broyden_NormSum_Curr,  Broyden_Denom*/
        #ifdef BROYDEN_SKIP_GPU_OPTIMIZATION
        RealTable1D h_sum_vector_data;
        RealTable2D h_Wmat_data;
        RealTable2D h_VmatTran_data;
        #else
        RealTable1D d_n_curr_in_data;
        RealTable1D d_n_curr_out_data;
        RealTable1D d_n_prev_in_data;
        RealTable1D d_F_curr_data;
        RealTable1D d_delta_F_curr_data;
        RealTable1D d_Norm_data;
        RealTable1D d_sum_vector_data;
        RealTable1D d_intermed_vector_data;

        RealTable2D d_Wmat_data;
        RealTable2D d_VmatTran_data;
        amrex::Gpu::DeviceVector<amrex::Real> d_Intermed_values_vec = {0., 0., 0.};
        /*stores Broyden_Norm, Broyden_NormSum_Curr,  Broyden_Denom*/
        #endif
    #else
    RealTable2D h_Jinv_curr_data;
    RealTable1D h_n_start_in_data;
    RealTable1D h_delta_n_curr_data;
    amrex::Vector<RealTable1D*> W_Broyden;
    amrex::Vector<RealTable1D*> V_Broyden;
    #endif

    template<typename NSType>
    void Set_TerminalBiasesAndContactPotential(NSType const& NS);

    template<typename NSType>
    void CopyFromNS_ChargeComputedFromNEGF(NSType const& NS);

    template<typename NSType>
    void CopyToNS_ChargeComputedUsingSelfConsistencyAlgorithm(NSType const& NS);

    template<typename NSType>
    void Write_MoreDataAndComputeCurrent(NSType const& NS, 
                                         std::string const& write_filename,
                                         bool const compute_current_flag);

    void Perform_SelfConsistencyAlgorithm();

    void Reset_ForNextBiasStep();

    void Obtain_maximum_time(amrex::Real const* total_time_counter_diff);


    #ifdef BROYDEN_PARALLEL
    void Deallocate_Broyden_Parallel();
    void Define_Broyden_Partition();
    void Define_MPI_Vector_Type_and_MPI_Vector_Sum();
    void Free_MPI_Vector_Type_and_MPI_Vector_Sum();
    MPI_Datatype MPI_Vector_Type;
    MPI_Op Vector_Add;

    static void Vector_Add_Func_wrapper(double* A, double* B, int* veclen, MPI_Datatype* dtype) 
    {
        reinterpret_cast<c_TransportSolver*>(A,B,veclen,dtype)->Vector_Add_Func(A,B,veclen,dtype);
    }

    void Vector_Add_Func(double* A, double* B, int* veclen, MPI_Datatype* dtype)
    {
        /* Note: https://www.mpich.org/static/docs/v3.2/www3/MPI_Op_create.html
         * The calling list for the user function type is
         * typedef void (MPI_User_function) ( void * A, void * B, int * veclen, MPI_Datatype * ); 
         * where the operation is B[i] = A[i] op B[i], for i=0,...,veclen-1. 
         * A pointer to the datatype given to the MPI collective computation routine (i.e., MPI_Reduce, MPI_Allreduce, 
         * MPI_Scan, or MPI_Reduce_scatter) is also passed to the user-specified routine. 
         */
        //double* Avec = A;
        //double* Bvec = B;
        for (int i=0; i < *veclen; ++i)
        {
            B[i] += A[i];
        }
    }
    template<typename NSType>
    void Create_Global_Output_Data(NSType const& NS);
    #else
    void Set_Broyden();
    void Reset_Broyden();
    void Execute_Broyden_First_Algorithm();
    void Execute_Broyden_Modified_Second_Algorithm();
    void Execute_Simple_Mixing_Algorithm();
    void Deallocate_Broyden_Serial();
    #endif

    void SetVal_RealTable1D (RealTable1D& Tab1D_data, amrex::Real val);    
    void SetVal_RealTable2D(RealTable2D& Tab2D_data, amrex::Real val);

    template<typename TableType>
    void Read_Table1D(int assert_size,
                      TableType& Tab1D_data,
                      std::string filename);
    
    template<typename TableType>
    void Read_Table2D(int assert_size,
                      TableType& Tab2D_data,
                      std::string filename);

    template<typename VectorType, typename TableType>
    void Write_Table1D(const amrex::Vector<VectorType>& Vec,
                       const TableType& Arr_data,
                       std::string filename,
                       std::string header);

    template<typename TableType>
    void Write_Table2D(const TableData<TableType, 2>& Arr_data,
                       std::string filename,
                       std::string header);
};

#endif
