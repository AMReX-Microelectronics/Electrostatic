
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

#include <cmath>
#include <iomanip>
#include <math.h> 
#include<stdlib.h>
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

MatrixDType get_beta(amrex::Real gamma, int M, int J) 
{
   amrex::GpuComplex arg(0., -MathConst::pi*J/M);
   return 2 * gamma * cos(-1*arg.imag()) * exp(arg); 
}

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

void PrintTable_glo(Matrix2D& G_glo_data)
{

    auto tlo = G_glo_data.lo();
    auto thi = G_glo_data.hi();
	
    auto const& G = G_glo_data.table();

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
}

//template<std::size_t N_total>
void MatInv_BlockTriDiagonal(int N_total,
		             const Matrix1D& h_A_loc_data, 
                             const Matrix1D& h_B_data, 
                             const Matrix1D& h_C_data,
                             Matrix2D& h_G_loc_data, 
                             const amrex::Vector<int>& cumu_blk_size,
                             const amrex::Vector<int>& vec_col_gids,
                             const int num_proc_with_blk)
{
/* Gatherv local A 1D table, A_loc_data, from all procs into global 1D table, A_glo_data, at the root.
 * Root processes Xtil, Ytil, X, Y 1Dtables.
 * Root deletes A_glo_data. 
 * Root broadcasts Xtil and Ytil to all procs.
 * Root Scatterv X and Y to other procs.
 * Each processor launches a GPU kernel where each thread works on each local column of G (the inverse we compute).
 * stream or device synchronise.
 */ 

    amrex::Real mat_inv_beg_time_cpu = amrex::second();


    int num_proc = ParallelDescriptor::NProcs(); 
    int my_rank = ParallelDescriptor::MyProc();
    int num_cols_loc = vec_col_gids.size();

    auto const& h_A_loc = h_A_loc_data.table();

    int ioproc = ParallelDescriptor::IOProcessorNumber();  
    Matrix1D h_A_glo_data({0},{N_total}, The_Pinned_Arena());
    Table1D<MatrixDType> h_A_glo = h_A_glo_data.table();
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

    MPI_Allgatherv(&h_A_loc(0),
                    num_cols_loc,
                    MPI_DOUBLE_COMPLEX,
                   &h_A_glo(0),
                    recv_count,
                    disp,
                    MPI_DOUBLE_COMPLEX,
                    ParallelDescriptor::Communicator());

    //if(my_rank == 0) {
    //   for (int i=0; i <N_total; ++i) {
    //       amrex::Print() << i << "  " <<  A_glo(i) << "\n";
    //   }
    //}
     
    auto const& h_B = h_B_data.const_table();
    auto const& h_C = h_C_data.const_table();


    Matrix1D h_Xtil_glo_data({0},{N_total},The_Pinned_Arena());
    Matrix1D h_Ytil_glo_data({0},{N_total},The_Pinned_Arena());
    Matrix1D h_X_glo_data({0},{N_total},The_Pinned_Arena());
    Matrix1D h_Y_glo_data({0},{N_total},The_Pinned_Arena());
    auto const& h_Ytil_glo = h_Ytil_glo_data.table();
    auto const& h_Xtil_glo = h_Xtil_glo_data.table();
    auto const& h_Y_glo = h_Y_glo_data.table();
    auto const& h_X_glo = h_X_glo_data.table();

    //if(my_rank == ParallelDescriptor::IOProcessorNumber()) 
    //{ 
        amrex::Real recursion_time_beg = amrex::second();

        h_Y_glo(0) = 0;
        for (int n = 1; n < N_total; ++n)
        {
        int p = (n-1)%2;
            h_Ytil_glo(n) = h_C(p) / ( h_A_glo(n-1) - h_Y_glo(n-1) );
            h_Y_glo(n) = h_B(p) * h_Ytil_glo(n);
        }

        h_X_glo(N_total-1) = 0;
        for (int n = N_total-2; n > -1; n--)
        {
        int p = n%2;
            h_Xtil_glo(n) = h_B(p)/(h_A_glo(n+1) - h_X_glo(n+1));
            h_X_glo(n) = h_C(p)*h_Xtil_glo(n);
        }
        //amrex::Print() << "\nYtil & Y: \n";
        //for (std::size_t n = 0; n < N_total; ++n)
        //{   
        //    amrex::Print() << std::setw(25)<< Ytil_glo(n) << std::setw(25) << Y_glo(n)<< "\n";
        //}
        //amrex::Print() << "\nXtil & X: \n";
        //for (int n = N_total-1; n > -1; n--)
        //{   
        //    amrex::Print() << std::setw(25)<< Xtil_glo(n) << std::setw(25) << X_glo(n)<< "\n";
        //}

        amrex::Real recursion_time = amrex::second() - recursion_time_beg;

        h_A_glo_data.clear();


        Matrix1D h_X_loc_data({0},{num_cols_loc}); //located on GPU 
        Matrix1D h_Y_loc_data({0},{num_cols_loc}); //located on GPU 
        auto const& h_Y_loc = h_Y_loc_data.table();
        auto const& h_X_loc = h_X_loc_data.table();

        for (int e = 0; e < num_cols_loc; ++e) /*loop over elements of a block*/
        {
            int n = e + cumu_blk_size[my_rank]; 
            h_Y_loc(e) = h_Y_glo(n);
            h_X_loc(e) = h_X_glo(n);
        } 

        h_X_glo_data.clear();
        h_Y_glo_data.clear();

    //if(my_rank == 1) 
    //{
    //    std::cout << "num_cols_loc: " << num_cols_loc << "\n";
    //    std::cout << "\n(proc 1) Y: \n";
    //    for (int e = 0; e < num_cols_loc; ++e) /*loop over elements of a block*/
    //    {
    //        int n = e + cumu_blk_size[my_rank]; 
    //        std::cout  << "e, n, Yloc: " << e << " " << n << std::setw(25) << Y_loc(e)<< "\n";
    //    } 
    //    std::cout << "\n(proc 1) X: \n";
    //    for (int e = num_cols_loc-1; e > -1; --e) /*loop over elements of a block*/
    //    {
    //        int n = e + cumu_blk_size[my_rank]; 
    //        std::cout << "e, n, Yloc: " << e << " " << n << std::setw(25) << X_loc(e)<< "\n";
    //    } 
    //}
    amrex::Real mat_inv_time_cpu = amrex::second() - mat_inv_beg_time_cpu;

    amrex::Print() << "\nRecursion time (cpu): " << std::setw(15) << recursion_time << "\n";
    amrex::Print() << "Matrix inversion overall cpu time: " << std::setw(15) << mat_inv_time_cpu << "\n";
    amrex::Print() << "Fraction of recursion/overall_cpu times: " << recursion_time/mat_inv_time_cpu << "\n";


    amrex::Real gpu_overall_beg_time = amrex::second();
    amrex::Real gpu_parallelFor_time = 0;

    if(num_cols_loc > 0) {

        Matrix2D d_G_loc_data({0,0},{N_total, num_cols_loc}, The_Device_Arena());
        Matrix1D d_A_loc_data({0},{num_cols_loc}, The_Device_Arena()); 
        d_A_loc_data.copy(h_A_loc_data); 

	Matrix1D  d_Xtil_glo_data({0},{N_total}, The_Device_Arena()); 
        Matrix1D  d_Ytil_glo_data({0},{N_total}, The_Device_Arena()); 

	d_Xtil_glo_data.copy(h_Xtil_glo_data); 
	d_Ytil_glo_data.copy(h_Ytil_glo_data); 

	Matrix1D d_X_loc_data({0},{num_cols_loc}, The_Device_Arena()); 
        Matrix1D d_Y_loc_data({0},{num_cols_loc}, The_Device_Arena()); 
	d_X_loc_data.copy(h_X_loc_data); 
	d_Y_loc_data.copy(h_Y_loc_data); 
        amrex::Gpu::streamSynchronize();

        auto const& A        =    d_A_loc_data.const_table();
        auto const& Xtil_glo =    d_Xtil_glo_data.const_table();
        auto const& Ytil_glo =    d_Ytil_glo_data.const_table();
        auto const& X        =    d_X_loc_data.const_table();
        auto const& Y        =    d_Y_loc_data.const_table();

        auto const& G_loc = d_G_loc_data.table();

        int cumulative_columns = cumu_blk_size[my_rank]; 

        amrex::Real gpu_parallelFor_beg_time = amrex::second();

        amrex::ParallelFor(num_cols_loc, [=] AMREX_GPU_DEVICE (int n) noexcept
        {
            int n_glo = n + cumulative_columns; /*global column number*/

            G_loc(n_glo,n) =  1./(A(n) - X(n) - Y(n)); 

            for (int m = n_glo; m > 0; m--)
            {   
                G_loc(m-1,n) =  -Ytil_glo(m)*G_loc(m,n);
            }
            for (int m = n_glo; m < N_total-1; ++m)
            {   
                G_loc(m+1,n) = -Xtil_glo(m)*G_loc(m,n);
            }
        });

        gpu_parallelFor_time = amrex::second() - gpu_parallelFor_beg_time;

        amrex::Gpu::streamSynchronize();
        h_G_loc_data.copy(d_G_loc_data); //copy from gpu to cpu
    }

    amrex::Real gpu_overall_time = amrex::second() - gpu_overall_beg_time;

    amrex::Print() << "\nMatrix inversion gpu parallelFor time: " 
		       << std::setw(15) << gpu_parallelFor_time << "\n";
    amrex::Print() << "Matrix inversion overall gpu time: " << std::setw(15) << gpu_overall_time << "\n";
    amrex::Print() << "Fraction of gpu_pFor/overall_gpu times: " 
	           << gpu_parallelFor_time/gpu_overall_time << "\n";

    amrex::Real total_cpu_gpu_time = mat_inv_time_cpu + gpu_overall_time;
    amrex::Print() << "\nTotal cpu gpu time: " 
	           << total_cpu_gpu_time << "\n";

    amrex::Print() << "\nRecursion/total, (overall_cpu-recursion)/total, gpu_pFor/total, (overall_gpu-gpu_pFor)/total time: " << std::setw(20) 
	           << recursion_time/total_cpu_gpu_time << std::setw(20) 
	           << (mat_inv_time_cpu-recursion_time)/total_cpu_gpu_time << std::setw(20) 
	           << gpu_parallelFor_time/total_cpu_gpu_time << std::setw(20) 
		   << (gpu_overall_time - gpu_parallelFor_time)/total_cpu_gpu_time << "\n";
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
    Matrix1D h_A_blk_loc({0},{num_cols_loc}, The_Pinned_Arena());

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
    auto const& h_A_loc = h_A_blk_loc.table();

    for(auto& col_gid: vec_col_gids) 
    {
        amrex::Real layer_loc_y = -10e-9 
		+ (static_cast<amrex::Real>(col_gid)/(N_total-1))*20e-9;        
        
        amrex::Real r = pow( (  pow((0.          - point_charge_loc[0]),2)
                              + pow((layer_loc_y - point_charge_loc[1]),2)
                              + pow((0.          - point_charge_loc[2]),2) ), 0.5);

        h_A_loc(c) = -(PhysConst::q_e)/(4. * MathConst::pi * PhysConst::ep0 * r);

        //amrex::Print() << "(proc 0) column, A: " << col_gid << " "<< A_loc(c) << "\n";
        ++c;  
    }

    /*below are some parameters for (17,0) carbon nanotube*/
    amrex::Real gamma = 2.5; //eV
    int M=17;
    MatrixDType beta = get_beta(gamma,M,6);

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
  
    /*allocate local G*/
    Matrix2D h_G_blk_loc({0,0},{N_total, num_cols_loc}, The_Pinned_Arena());

    MatInv_BlockTriDiagonal(N_total, h_A_blk_loc, h_B_data, h_C_data, h_G_blk_loc, 
                            cumu_blk_size, vec_col_gids, num_proc_with_blk);   





    ///*Gatherv G_blk_loc into G_blk_glo*/
    //amrex::Real G_locglo_gatherv_time_beg = amrex::second();

    //Matrix2D h_G_glo_data;
    //Table2D<MatrixDType> h_G_glo;
    //int* recv_count;
    //int* disp;

    //if(my_rank == ParallelDescriptor::IOProcessorNumber()) 
    //{
    //    h_G_glo_data.resize({0,0},{N_total, N_total}, The_Pinned_Arena());
    //    h_G_glo = h_G_glo_data.table();
    //    recv_count = new int[num_proc];
    //    disp = new int [num_proc];
    //    for(int p=0; p < num_proc; ++p) {

    //       if(p < num_proc_with_blk) {
    //          recv_count[p] = cumu_blk_size[p+1] - cumu_blk_size[p];
    //          disp[p] = cumu_blk_size[p];
    //       }
    //       else {
    //          recv_count[p] = 0;
    //          disp[p] = 0;
    //       }
    //       //amrex::Print() << "p/recv_count/disp: " << p << " " << recv_count[p] << " " << disp[p] << "\n";
    //    }
    //}

    //MPI_Datatype column_type;

    //MPI_Type_vector(1, N_total+1, N_total+1, MPI_DOUBLE_COMPLEX, &column_type);

    //MPI_Type_commit(&column_type);

    //auto h_G_loc = h_G_blk_loc.table();

    //MPI_Gatherv(&h_G_loc(0,0),
    //            num_cols_loc,
    //            column_type,
    //            &h_G_glo(0,0),
    //            recv_count,
    //            disp,
    //            column_type,
    //            ParallelDescriptor::IOProcessorNumber(),
    //            ParallelDescriptor::Communicator());
    //ParallelDescriptor::Barrier();  

    //if(my_rank == ParallelDescriptor::IOProcessorNumber()) 
    //{ 
    //    delete [] recv_count;
    //    delete [] disp;
    //} 
    //MPI_Type_free(&column_type);

    //amrex::Real G_locglo_gatherv_time = amrex::second() - G_locglo_gatherv_time_beg;
    //amrex::Print() << "\nG_locglo_gatherv_time: " << std::setw(15) << G_locglo_gatherv_time<< "\n";

    //if(my_rank == ParallelDescriptor::IOProcessorNumber()) 
    //{
    //    amrex::Print() << "G_glo:\n";
    //    PrintTable_glo(h_G_glo_data);
    //}
    ////if(my_rank == 1) 
    ////{
    ////    amrex::Print() << "(rank 1) G_loc:\n";
    ////    PrintTable_loc(h_G_blk_loc);
    ////}
 
    amrex::Finalize();

}
