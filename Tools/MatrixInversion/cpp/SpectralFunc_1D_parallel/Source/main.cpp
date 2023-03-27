
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_MFParallelForC.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuComplex.H>
#include<AMReX_Print.H>
#include<AMReX_TableData.H>
#include <AMReX_TinyProfiler.H>
#include <AMReX_GpuUtility.H>
#include <AMReX_Config.H>
#include <AMReX_Loop.H>

#include <cmath>
#include <iomanip>
#include <math.h> 
#include<stdlib.h>
#define NUM_CONTACTS 2
using namespace amrex;
using MatrixDType = amrex::GpuComplex<amrex::Real>;
using Matrix1D = TableData<MatrixDType, 1>;
using Matrix2D = TableData<MatrixDType, 2>;

template<typename T>
class TD;

namespace MathConst
{
    static constexpr amrex::Real pi = static_cast<amrex::Real>(3.14159265358979323846);
}

namespace PhysConst
{
    static constexpr auto q_e   = static_cast<amrex::Real>( 1.602176634e-19 );
    static constexpr auto ep0   = static_cast<amrex::Real>( 8.8541878128e-12 );
//    static constexpr auto hbar  = static_cast<amrex::Real>( 1.054571817e-34 );
}

AMREX_GPU_HOST_DEVICE 
MatrixDType get_beta(amrex::Real gamma, int M, int J) 
{
   amrex::GpuComplex arg(0., -MathConst::pi*J/M);
   return 2 * gamma * cos(-1*arg.imag());// * exp(arg); 
}

AMREX_GPU_HOST_DEVICE 
MatrixDType get_Gamma(MatrixDType Sigma)
{
   amrex::GpuComplex val(-2.*Sigma.imag(), 0.);
   return val;
}


AMREX_GPU_HOST_DEVICE 
MatrixDType Compute_SurfaceGreensFunction(MatrixDType E, 
		                          amrex::Real U, 
					  MatrixDType beta, 
					  amrex::Real gamma) 
{

     MatrixDType EmU = E-U;

     MatrixDType EmU_sq = pow(EmU,2.);

     amrex::Real gamma_sq = pow(gamma,2.);

     MatrixDType Factor = EmU_sq + gamma_sq - pow(beta,2.);

     MatrixDType Sqrt = sqrt(pow(Factor,2.) - 4. * EmU_sq * gamma_sq);

     MatrixDType Denom = 2. * gamma_sq * EmU; 

     //amrex::Print() << "EmU: " << EmU << "\n";
     //amrex::Print() << "Factor: " << Factor << "\n";
     //amrex::Print() << "Sqrt: "  << Sqrt << "\n";
     //amrex::Print() << "Denom: " << Denom << "\n";
     //amrex::Print() << "Numerator: " << Factor+Sqrt << "\n";
     //amrex::Print() << "Value: " << (Factor+Sqrt)/Denom << "\n";

     return (Factor + Sqrt) / Denom;

}


AMREX_GPU_HOST_DEVICE 
MatrixDType get_Sigma(MatrixDType E,
                      amrex::Real U,
                      MatrixDType beta,
                      amrex::Real gamma) 
{

   MatrixDType gr = Compute_SurfaceGreensFunction(E, U, beta, gamma);

   MatrixDType Sigma = pow(gamma,2.)*gr;

   amrex::GpuComplex val(640.5, -3.1977);
   return val;

}



AMREX_GPU_HOST_DEVICE 
MatrixDType conjugate(MatrixDType a) 
{
   amrex::GpuComplex a_conj(a.real(), -1.*a.imag());
   return a_conj;
}


void PrintTable_loc(Matrix2D& G_blk)
{

    auto tlo = G_blk.lo();
    auto thi = G_blk.hi();
	
    auto const& G = G_blk.table();

	std::cout << "tlo: " << tlo[0] << " " << tlo[1] 
	          << "thi: " << thi[0] << " " << thi[1] << "\n";

    for (int i = tlo[0]; i < thi[0]; ++i) 
    { 
        for (int j = tlo[1]; j < thi[1]; ++j) //slow moving index. printing slow
	{ 
            std::cout << std::setw(12) << std::setprecision(2) << G(i,j);
        }
        std::cout << "\n";
    }
}


void PrintTable(Matrix2D& h_G_loc_data, 
		int N_total, 
		const amrex::Vector<int>& cumu_blk_size,
                const amrex::Vector<int>& vec_col_gids,
                const int num_proc_with_blk)
{

    /*Gatherv G_loc into G_glo*/
    amrex::Real G_locglo_gatherv_time_beg = amrex::second();

    int num_proc = ParallelDescriptor::NProcs(); 
    int my_rank = ParallelDescriptor::MyProc();
    int num_cols_loc = vec_col_gids.size();

    Matrix2D h_G_glo_data;
    Table2D<MatrixDType> h_G_glo;
    int* recv_count;
    int* disp;

    if(my_rank == ParallelDescriptor::IOProcessorNumber()) 
    {
        h_G_glo_data.resize({0,0},{N_total, N_total}, The_Pinned_Arena());
        h_G_glo = h_G_glo_data.table();
        recv_count = new int[num_proc];
        disp = new int [num_proc];
        for(int p=0; p < num_proc; ++p) {

           if(p < num_proc_with_blk) {
              recv_count[p] = cumu_blk_size[p+1] - cumu_blk_size[p];
              disp[p] = cumu_blk_size[p];
           }
           else {
              recv_count[p] = 0;
              disp[p] = 0;
           }
           //amrex::Print() << "p/recv_count/disp: " << p << " " << recv_count[p] << " " << disp[p] << "\n";
        }
    }

    MPI_Datatype column_type;

    MPI_Type_vector(1, N_total+1, N_total+1, MPI_DOUBLE_COMPLEX, &column_type);

    MPI_Type_commit(&column_type);

    auto h_G_loc = h_G_loc_data.table();

    MPI_Gatherv(&h_G_loc(0,0),
                num_cols_loc,
                column_type,
                &h_G_glo(0,0),
                recv_count,
                disp,
                column_type,
                ParallelDescriptor::IOProcessorNumber(),
                ParallelDescriptor::Communicator());
    ParallelDescriptor::Barrier();  

    if(my_rank == ParallelDescriptor::IOProcessorNumber()) 
    { 
        delete [] recv_count;
        delete [] disp;
    } 
    MPI_Type_free(&column_type);

    amrex::Real G_locglo_gatherv_time = amrex::second() - G_locglo_gatherv_time_beg;
    amrex::Print() << "\nG_locglo_gatherv_time: " << std::setw(15) << G_locglo_gatherv_time<< "\n";



    if(my_rank == ParallelDescriptor::IOProcessorNumber()) 
    {
        auto tlo = h_G_glo_data.lo();
        auto thi = h_G_glo_data.hi();
            
        auto const& G = h_G_glo_data.table();

            amrex::Print() << "tlo: " << tlo[0] << " " << tlo[1] << "\n"
                           << "thi: " << thi[0] << " " << thi[1] << "\n";

        for (int i = tlo[0]; i < thi[0]; ++i) 
        { 
            for (int j = tlo[1]; j < thi[1]; ++j) //slow moving index. printing slow
           	{ 
                amrex::Print() << std::setw(12) << std::setprecision(2) << G(i,j);
            }
            amrex::Print() << "\n";
        }

	h_G_glo_data.clear();
    }

}


//template<std::size_t N_total>
void Obtain_GreensAndSpectralFunctions(int N_total,
		             const Matrix1D& h_Alpha_loc_data, 
                             const Matrix1D& h_B_data, 
                             const Matrix1D& h_C_data,
                             const Matrix1D& h_E_data,
                             Matrix2D& h_Sigma_glo_data, 
                             Matrix2D& hd_G_loc_data, 
                             Matrix2D& hd_A_loc_data, 
                             const amrex::Vector<int>& cumu_blk_size,
                             const amrex::Vector<int>& vec_col_gids,
                             const int num_proc_with_blk,
			     const amrex::GpuArray<int, NUM_CONTACTS>& global_contact_index)
{
/* Six step procedure:
 * 1) MPI_Allgatherv local Hamiltonian diagonal elements in Alpha_loc_data into Alpha_glo_data, both allocated on the PinnedArena.
 * 2) Compute on all processes, 1DTable arrays Xtil, Ytil, X, Y recursively, on the PinnedArena and then delete Alpha_glo_data. 
 * These arrays have a size equal to the column (or row) size of the full (or global) matrix to be inverted.
 * 3) All processes copy the respective portion of X and Y needed later into arrays allocated on managed memory. 
 * These portions have a size equal to the number of columns in the local matrix block. 
 * 4) Copy needed arrays on to the device. These include Xtil, Ytil, X_local, Y_local, Alpha_loc. 
 * 5) Launce parallelFor, assign each gpu-thread to compute a column of the inverted matrix. 
 * Each GPU is responsible for computing it's own block.
 * 6) Do a stream synchronise and copy the inverted matrix block on to the CPU.
 */ 
    int num_proc = ParallelDescriptor::NProcs(); 
    int my_rank = ParallelDescriptor::MyProc();
    int num_cols_loc = vec_col_gids.size();
    int ioproc = ParallelDescriptor::IOProcessorNumber();  
    int num_traces=1;


    auto thi_Sigma_glo = h_Sigma_glo_data.hi();
    int num_contacts = thi_Sigma_glo[0];
    int num_EnPts    = thi_Sigma_glo[1];
    amrex::Print() << "num_contacts/num_EnPts: " << num_contacts << " " << num_EnPts << "\n";

    Matrix1D h_Alpha_glo_data({0},{N_total}, The_Pinned_Arena());
    int recv_count[num_proc];
    int disp[num_proc];

    for(int p=0; p < num_proc; ++p) {

       if(p < num_proc_with_blk) {
          recv_count[p] = cumu_blk_size[p+1] - cumu_blk_size[p];
          disp[p] = cumu_blk_size[p];
       }
       else {
          recv_count[p] = 0;
          disp[p] = 0;
       }
    }
    Matrix1D h_Xtil_glo_data({0},{N_total},The_Pinned_Arena());
    Matrix1D h_Ytil_glo_data({0},{N_total},The_Pinned_Arena());
    Matrix1D h_X_glo_data({0},{N_total},The_Pinned_Arena());
    Matrix1D h_Y_glo_data({0},{N_total},The_Pinned_Arena());
    Matrix1D h_Alpha_contact_data({0},{num_contacts},The_Pinned_Arena()); 

    Matrix1D h_X_loc_data({0},{num_cols_loc},The_Pinned_Arena()); 
    Matrix1D h_Y_loc_data({0},{num_cols_loc},The_Pinned_Arena()); 
    Matrix1D h_X_contact_data({0},{num_contacts},The_Pinned_Arena()); 
    Matrix1D h_Y_contact_data({0},{num_contacts},The_Pinned_Arena()); 

    #if (defined AMREX_USE_GPU && GPU_MATRIX_OP)
    Matrix1D d_Alpha_loc_data({0},{num_cols_loc}, The_Arena()); 
    Matrix1D  d_Xtil_glo_data({0},{N_total}, The_Arena()); 
    Matrix1D  d_Ytil_glo_data({0},{N_total}, The_Arena()); 
    Matrix1D d_X_loc_data({0},{num_cols_loc}, The_Arena()); 
    Matrix1D d_Y_loc_data({0},{num_cols_loc}, The_Arena()); 
    Matrix1D d_Alpha_contact_data({0},{num_contacts}, The_Arena()); 
    Matrix1D d_X_contact_data({0},{num_contacts}, The_Arena()); 
    Matrix1D d_Y_contact_data({0},{num_contacts}, The_Arena()); 
    Matrix2D d_Sigma_glo_data({0,0},{num_contacts,num_EnPts}, The_Arena()); 
    Matrix1D d_Adiag_loc_data({0}, {num_cols_loc}, The_Arena());
    Matrix1D d_T_loc_data({0}, {num_EnPts}, The_Arena());
    amrex::Gpu::DeviceVector<amrex::Real> d_Trace_r(num_traces, amrex::Real(0.0));
    amrex::Gpu::DeviceVector<amrex::Real> d_Trace_i(num_traces, amrex::Real(0.0));
    #endif
    amrex::Gpu::HostVector<amrex::Real> h_Trace_r(num_traces);
    amrex::Gpu::HostVector<amrex::Real> h_Trace_i(num_traces);
    Matrix1D h_T_loc_data({0}, {num_EnPts}, The_Arena());


    Table1D<MatrixDType> h_Alpha_glo = h_Alpha_glo_data.table();
    auto const& h_Alpha_loc = h_Alpha_loc_data.table();
    auto const& h_B = h_B_data.const_table();
    auto const& h_C = h_C_data.const_table();
    auto const& h_E = h_E_data.const_table();
    auto const& h_Ytil_glo = h_Ytil_glo_data.table();
    auto const& h_Xtil_glo = h_Xtil_glo_data.table();
    auto const& h_Y_glo = h_Y_glo_data.table();
    auto const& h_X_glo = h_X_glo_data.table();
    auto const& h_Alpha_contact = h_Alpha_contact_data.table();
    auto const& h_Y_loc = h_Y_loc_data.table();
    auto const& h_X_loc = h_X_loc_data.table();
    auto const& h_Y_contact = h_Y_contact_data.table();
    auto const& h_X_contact = h_X_contact_data.table();

    #if(defined AMREX_USE_GPU && GPU_MATRIX_OP)
    auto const& Alpha         = d_Alpha_loc_data.const_table();
    auto const& Xtil_glo      = d_Xtil_glo_data.const_table();
    auto const& Ytil_glo      = d_Ytil_glo_data.const_table();
    auto const& X             = d_X_loc_data.const_table();
    auto const& Y             = d_Y_loc_data.const_table();
    auto const& Alpha_contact = d_Alpha_contact_data.const_table();
    auto const& X_contact     = d_X_contact_data.const_table();
    auto const& Y_contact     = d_Y_contact_data.const_table();
    auto const& Sigma         = d_Sigma_glo_data.const_table();
    auto const& Adiag_loc     = d_Adiag_loc_data.table(); //Spectral function
    auto const& T_loc         = d_T_loc_data.table(); //Spectral function

    amrex::Real* trace_r      = d_Trace_r.dataPtr();
    amrex::Real* trace_i      = d_Trace_i.dataPtr();

    //auto const& G_contact     = d_G_contact_data.table();
    #else 
    auto const& Alpha         = h_Alpha_loc_data.const_table();
    auto const& Xtil_glo      = h_Xtil_glo_data.const_table();
    auto const& Ytil_glo      = h_Ytil_glo_data.const_table();
    auto const& X             = h_X_loc_data.const_table();
    auto const& Y             = h_Y_loc_data.const_table();
    auto const& Alpha_contact = h_Alpha_contact_data.const_table();
    auto const& X_contact     = h_X_contact_data.const_table();
    auto const& Y_contact     = h_Y_contact_data.const_table();
    auto const& Sigma         = h_Sigma_glo_data.const_table();
    auto const& T_loc         = h_T_loc_data.const_table();

    amrex::Real* trace_r      = h_Trace_r.dataPtr();
    amrex::Real* trace_i      = h_Trace_i.dataPtr();
    #endif
    auto const& G_loc    =    hd_G_loc_data.table(); //Green's function
    auto const& A_loc    =    hd_A_loc_data.table(); //Spectral function


    for(int e=0; e<num_EnPts; ++e) {

        //BL_PROFILE_VAR("steps_1_2_3", steps_1_2_3);

        /*Step 1*/	

        MPI_Allgatherv(&h_Alpha_loc(0),
                        num_cols_loc,
                        MPI_DOUBLE_COMPLEX,
                       &h_Alpha_glo(0),
                        recv_count,
                        disp,
                        MPI_DOUBLE_COMPLEX,
                        ParallelDescriptor::Communicator());


        /*Step 2*/	

        h_Y_glo(0) = 0;
        for (int n = 1; n < N_total; ++n)
        {
        int p = (n-1)%2;
            h_Ytil_glo(n) = h_C(p) / ( h_Alpha_glo(n-1) - h_Y_glo(n-1) );
            h_Y_glo(n) = h_B(p) * h_Ytil_glo(n);
        }

        h_X_glo(N_total-1) = 0;
        for (int n = N_total-2; n > -1; n--)
        {
        int p = n%2;
            h_Xtil_glo(n) = h_B(p)/(h_Alpha_glo(n+1) - h_X_glo(n+1));
            h_X_glo(n) = h_C(p)*h_Xtil_glo(n);
        }

        for (int c = 0; c < num_contacts; ++c) 
        {
            int n = global_contact_index[c];
            h_Alpha_contact(c) = h_Alpha_glo(n);
        } 

        /*Step 3*/	

        for (int c = 0; c < num_cols_loc; ++c) 
        {
            int n = c + cumu_blk_size[my_rank];
            h_Y_loc(c) = h_Y_glo(n);
            h_X_loc(c) = h_X_glo(n);
        } 

        for (int c = 0; c < num_contacts; ++c) 
        {
            int n = global_contact_index[c];
            h_Y_contact(c) = h_Y_glo(n);
            h_X_contact(c) = h_X_glo(n);
        } 

        //BL_PROFILE_VAR_STOP(steps_1_2_3);

        /*Step 4*/	
        #if (defined AMREX_USE_GPU && GPU_MATRIX_OP)
        //BL_PROFILE_VAR("step4_CpuToGpuCopy_time", step4_CpuToGpuCopy_time);
	
        d_Alpha_loc_data.copy(h_Alpha_loc_data); 
        d_Xtil_glo_data.copy(h_Xtil_glo_data); 
        d_Ytil_glo_data.copy(h_Ytil_glo_data); 
        d_X_loc_data.copy(h_X_loc_data); 
        d_Y_loc_data.copy(h_Y_loc_data); 
        d_Alpha_contact_data.copy(h_Alpha_contact_data); 
        d_X_contact_data.copy(h_X_contact_data); 
        d_Y_contact_data.copy(h_Y_contact_data); 
        d_Sigma_glo_data.copy(h_Sigma_glo_data); 

        amrex::Gpu::streamSynchronize();

        //BL_PROFILE_VAR_STOP(step4_CpuToGpuCopy_time);
        #endif
 
        /*Step 5*/	
        //BL_PROFILE_VAR("step5_pFor_time", step5_pFor_time);

        int cumulative_columns = cumu_blk_size[my_rank]; 
	int e_id = e;
	MatrixDType E = h_E(e_id);
        amrex::ParallelFor(num_cols_loc, [=] AMREX_GPU_DEVICE (int n) noexcept
        {
            int n_glo = n + cumulative_columns; /*global column number*/

            G_loc(n_glo,n) =  1./(Alpha(n) - X(n) - Y(n)); 

            for (int m = n_glo; m > 0; m--)
            {   
                G_loc(m-1,n) =  -Ytil_glo(m)*G_loc(m,n);
            }
            for (int m = n_glo; m < N_total-1; ++m)
            {   
                G_loc(m+1,n) = -Xtil_glo(m)*G_loc(m,n);
            }

            //MatrixDType G_contact_kk contains G_glo(k_glo,k_glo) element
            //MatrixDType G_contact_nk contains G_glo(n_glo,k_glo) element
            MatrixDType A_kn_arr[NUM_CONTACTS];

            for (int k=0; k < num_contacts; ++k) 
            {
                int k_glo = global_contact_index[k];
                MatrixDType G_contact_kk =  1./(Alpha_contact(k) - X_contact(k) - Y_contact(k)); 
                
                MatrixDType temp = G_contact_kk;
                for (int m = k_glo; m < n_glo; ++m)
                {   
                    temp = -Xtil_glo(m)*temp;
                }
                for (int m = k_glo; m > n_glo; m--)
                {   
                    temp = -Ytil_glo(m)*temp;
                }
                MatrixDType G_contact_nk = temp;

                MatrixDType A_kn  = G_contact_kk * get_Gamma(Sigma(k,e_id)) * 
                                    conjugate( G_contact_nk );

                A_loc(k_glo, n) += A_kn;
                for (int m = k_glo+1; m < N_total; ++m)
                {   
                    A_kn = -Xtil_glo(m-1)*A_kn;	
                    A_loc(m,n) +=  A_kn;
                }
                for (int m = k_glo-1; m >= 0; m--)
                {   
                    A_kn = -Ytil_glo(m+1)*A_kn;	
                    A_loc(m,n) +=  A_kn;
                }
                A_kn_arr[k] = A_kn;
            }  
            
            T_loc(0) = get_Gamma(Sigma(0,e_id)) * A_kn_arr[1];

            amrex::HostDevice::Atomic::Add(&(trace_r[0]), Adiag_loc(n).real());
            amrex::HostDevice::Atomic::Add(&(trace_i[0]), Adiag_loc(n).imag());
        });

        //BL_PROFILE_VAR_STOP(step5_pFor_time);

        #ifdef AMREX_USE_GPU 
        amrex::Gpu::streamSynchronize();
        h_T_loc_data.copy(d_T_loc_data); 

        amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_Trace_r.begin(), d_Trace_r.end(), h_Trace_r.begin());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_Trace_i.begin(), d_Trace_i.end(), h_Trace_i.begin());
        #endif 
   
        for (int t=0; t< num_traces; ++t)
        {
            amrex::ParallelDescriptor::ReduceRealSum(h_Trace_r[t]);
            amrex::ParallelDescriptor::ReduceRealSum(h_Trace_i[t]);
            amrex::Print() << "Trace: "<< t << " " << h_Trace_r[t] << " "<< h_Trace_i[t]<< "\n";
        }

    }

    /*deallocate memory*/
    #if(defined AMREX_USE_GPU && GPU_MATRIX_OP)
    d_Alpha_loc_data.clear();
    d_Xtil_glo_data.clear();
    d_Ytil_glo_data.clear();
    d_X_loc_data.clear();
    d_Y_loc_data.clear();
    d_Alpha_contact_data.clear();
    d_X_contact_data.clear();
    d_Y_contact_data.clear();
    d_Adiag_loc_data.clear();
    d_T_loc_data.clear();
    d_Trace_r.clear();
    d_Trace_i.clear();
    #endif 
    h_Alpha_glo_data.clear();
    h_Xtil_glo_data.clear();
    h_Ytil_glo_data.clear();
    h_X_glo_data.clear();
    h_Y_glo_data.clear();
    h_X_loc_data.clear();
    h_Y_loc_data.clear();
    h_Alpha_contact_data.clear();
    h_X_contact_data.clear();
    h_Y_contact_data.clear();
    h_T_loc_data.clear();
    h_Trace_r.clear();
    h_Trace_i.clear();

}

int main (int argc, char* argv[])
{

    amrex::Initialize(argc,argv);
    int num_proc = ParallelDescriptor::NProcs(); 
    int my_rank = ParallelDescriptor::MyProc();
    amrex::Print() << "total number of procs: " << num_proc << "\n";

    //Construct Tridiagonal Dummy Hamiltonian
    std::array<amrex::Real,3> point_charge_loc {0., 0., 1e-9};

    int N_total = 4; /*matrix size*/
    amrex::ParmParse pp;
    pp.query("N_total", N_total);

    //ParallelDescriptor::Bcast(&N_total, 1, ParallelDescriptor::IOProcessorNumber());
    if(my_rank == 1) {
	    std::cout << "Rank 1 reporting, total matrix size: " << N_total << "\n";
    }

    amrex::Print() << "total matrix size: " << N_total << "\n";
    const int MAX_THRESHOLD_BLK_SIZE = 40000; /*matrix size*/

    bool flag_fixed_blk_size = false;

    int max_blk_size = 2;   /*max block size, may change depending on flag_fixed_blk_size*/

    if(flag_fixed_blk_size) amrex::Print() << "max_blk_size is fixed by user\n";
    else {
         amrex::Print() << "max_blk_size is computed at run-time\n";
         max_blk_size = ceil(static_cast<amrex::Real>(N_total)/num_proc);
         if(max_blk_size > MAX_THRESHOLD_BLK_SIZE) {
            max_blk_size = MAX_THRESHOLD_BLK_SIZE;
            /*assert that use larger number of procs*/
         }
    }
    amrex::Print() << "maximum block size (in number of columns): " << max_blk_size << "\n";

    const int num_blk_glo = ceil(static_cast<amrex::Real>(N_total)/max_blk_size);
    amrex::Print() << "global number of blocks: " << num_blk_glo << "\n";
     
    const int num_proc_with_blk = num_blk_glo;
    /*if num_proc_with_blk >= num_proc, assert.*/

    amrex::Vector<int> cumu_blk_size(num_proc_with_blk + 1);

    cumu_blk_size[0] = 0;
    for(int p=1; p < num_proc_with_blk; ++p)
    {
        cumu_blk_size[p] = cumu_blk_size[p-1] + max_blk_size;
        /*note that all proc except the last proc in "all procs with blocks" have N_blk number of columns. The last proc in "all procs with blocks" may have number of columns smaller than N_blk.
        */
    }
    cumu_blk_size[num_proc_with_blk] = cumu_blk_size[num_proc_with_blk-1] 
                                     + N_total - (num_proc_with_blk-1)*max_blk_size;

    //amrex::Print() << "cumulative blk size: \n";
    //for(int p=0; p < num_proc_with_blk+1; ++p)
    //{
    //    amrex::Print() << cumu_blk_size[p] << "\n";
    //}


    amrex::Vector<int>  vec_col_gids;
   
    if(my_rank < num_proc_with_blk) 
    {

        int blk_gid = my_rank;
        int blk_size = cumu_blk_size[blk_gid+1] - cumu_blk_size[blk_gid];

        for (int c=0; c < blk_size; ++c) 
        {
            int col_gid = cumu_blk_size[blk_gid] + c;
            vec_col_gids.push_back(col_gid);
        } 

        //if(my_rank == 1) {
        //    std::cout      << "(proc) My rank: " << my_rank << "\n";
        //    std::cout << " blk_size: " << vec_col_gids.size() << "\n";
        //    std::cout << " vec_col_gids: ";
        //    for(auto& col_gid: vec_col_gids) {
        //       std::cout << std::setw(5) << col_gid << " ";
        //    } 
        //    std::cout <<  "\n ";
        //}
    }


    /*setting A, the diagonal 1D vector, which is partitioned in the form of some blocks*/
    /*in general A can be 2D*/
    int num_cols_loc = vec_col_gids.size();
    Matrix1D h_Alpha_loc_data({0},{num_cols_loc}, The_Pinned_Arena());

    /*R this is the size of repeated B and C vectors in Hamiltonian in the mode-space approximation for CNT*/
    const int R = 2; 

    /*In general, B and C are 2D blocks of different sizes 2D blocks*/
    Matrix1D h_B_data({0},{R},The_Pinned_Arena());
    Matrix1D h_C_data({0},{R},The_Pinned_Arena());
    auto const& h_B = h_B_data.table();
    auto const& h_C = h_C_data.table();

    /*specifying A block vector for CNT with point charge potential*/

    auto blk_size = vec_col_gids.size();
    int c=0;
    auto const& h_Alpha_loc = h_Alpha_loc_data.table();

    for(auto& col_gid: vec_col_gids) 
    {
        amrex::Real layer_loc_y = -10e-9 
		+ (static_cast<amrex::Real>(col_gid)/(N_total-1))*20e-9;        
        
        amrex::Real r = pow( (  pow((0.          - point_charge_loc[0]),2)
                              + pow((layer_loc_y - point_charge_loc[1]),2)
                              + pow((0.          - point_charge_loc[2]),2) ), 0.5);

        h_Alpha_loc(c) = -(PhysConst::q_e)/(4. * MathConst::pi * PhysConst::ep0 * r);

        //amrex::Print() << "(proc 0) column, A: " << col_gid << " "<< A_loc(c) << "\n";
        ++c;  
    }

    /*below are some parameters for (17,0) carbon nanotube*/
    amrex::Real gamma = 2.5; //eV
    int M=17;
    MatrixDType beta = get_beta(gamma,M,6);
    amrex::Print() << "beta: " << beta << "\n";
    amrex::Real U_contact[NUM_CONTACTS] = {0.,0.};

    /*specifying B and C vectors for the special case of mode-space Hamiltonian of CNT*/
    for (std::size_t i = 0; i < R; ++i)
    {
       if(i%2 == 0) {
          h_B(i) = conjugate(beta);
          h_C(i) = beta;
       }
       else {
          h_B(i) = gamma;
          h_C(i) = gamma;
       } 
    }
    
    /*define energy grid and sigma*/
    int num_contacts = NUM_CONTACTS;
    amrex::GpuArray<int,NUM_CONTACTS> global_contact_index = {0, N_total-1};
    int num_EnPts = 10;
    amrex::Real EnRange[2] = {-1., 1.}; //eV
    Matrix1D h_E_data({0},{num_EnPts},The_Pinned_Arena());
    auto const& h_E = h_E_data.table();
    amrex::Real deltaE = (EnRange[1]-EnRange[0])/(num_EnPts-1);
    for (std::size_t e = 0; e < num_EnPts; ++e)
    {
        amrex::GpuComplex E(EnRange[0] + e*deltaE, 1e-14);
	h_E(e) = E;
    }


    Matrix2D h_Sigma_glo_data({0,0},{num_contacts,num_EnPts},The_Pinned_Arena());
    auto const& h_Sigma = h_Sigma_glo_data.table();

    for (std::size_t c = 0; c < num_contacts; ++c)
    {
        amrex::Print() << "contact: " << c << "\n";    
        for (std::size_t e = 0; e < num_EnPts; ++e)
        {
            h_Sigma(c,e) = get_Sigma(h_E(e), U_contact[c], beta, gamma);
	    //amrex::Print() << "e: " << e <<  " " << h_E(e) << "  "<< h_Sigma(c,e) << "\n";
        }
    }

    //PrintTable_loc(h_Sigma_glo_data);

    /*allocate local G*/
    #if(defined AMREX_USE_GPU && GPU_MATRIX_OP)
    Matrix2D d_G_loc_data({0,0},{N_total, num_cols_loc}, The_Arena());
    Matrix2D d_A_loc_data({0,0},{N_total, num_cols_loc}, The_Arena());
    
    Obtain_GreensAndSpectralFunctions(N_total, h_Alpha_loc_data, h_B_data, h_C_data, h_E_data,h_Sigma_glo_data,
		                               d_G_loc_data, d_A_loc_data,   
					       cumu_blk_size, vec_col_gids, num_proc_with_blk, global_contact_index);   
    #else 
    Matrix2D h_G_loc_data({0,0},{N_total, num_cols_loc}, The_Pinned_Arena());
    Matrix2D h_A_loc_data({0,0},{N_total, num_cols_loc}, The_Pinned_Arena());

    Obtain_GreensAndSpectralFunctions(N_total, h_Alpha_loc_data, h_B_data, h_C_data, h_E_data,h_Sigma_glo_data,
		            h_G_loc_data, h_A_loc_data,
                            cumu_blk_size, vec_col_gids, num_proc_with_blk, global_contact_index);   
    #endif



    int print_matrix_flag = false;
    pp.query("print_matrix", print_matrix_flag);

    if(print_matrix_flag) {
        /**copy G_loc from device to host*/
        
        #if (defined AMREX_USE_GPU && GPU_MATRIX_OP)
        BL_PROFILE_VAR("step6_CopyGpuToCpu", step6_copyGpuToCpu);
        
	amrex::Gpu::streamSynchronize();
        Matrix2D h_G_loc_data({0,0},{N_total, num_cols_loc}, The_Pinned_Arena());
        Matrix2D h_A_loc_data({0,0},{N_total, num_cols_loc}, The_Pinned_Arena());
        h_G_loc_data.copy(d_G_loc_data); //copy from gpu to cpu
        h_A_loc_data.copy(d_A_loc_data); //copy from gpu to cpu

        BL_PROFILE_VAR_STOP(step6_copyGpuToCpu);
        #endif
        
	amrex::Print() << "G_glo: \n";
        PrintTable(h_G_loc_data,N_total,cumu_blk_size, vec_col_gids, num_proc_with_blk);
	amrex::Print() << "A_glo: \n";
        PrintTable(h_A_loc_data,N_total,cumu_blk_size, vec_col_gids, num_proc_with_blk);

        #if (defined AMREX_USE_GPU && GPU_MATRIX_OP)
        h_G_loc_data.clear();
        h_A_loc_data.clear();
        #endif
    }
    

    /*deallocate memory*/
    h_Alpha_loc_data.clear();
    h_B_data.clear();
    h_C_data.clear();
    #if (defined AMREX_USE_GPU && GPU_MATRIX_OP)
    d_G_loc_data.clear();
    d_A_loc_data.clear();
    #else 
    h_G_loc_data.clear();
    h_A_loc_data.clear();
    #endif
    h_Sigma_glo_data.clear();
    h_E_data.clear();


    amrex::Finalize();

}
