
#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Config.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_GpuUtility.H>
#include <AMReX_Loop.H>
#include <AMReX_MFParallelForC.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_TableData.H>
#include <AMReX_TinyProfiler.H>
#include <AMReX_Vector.H>
#include <math.h>
#include <stdlib.h>

#include <cmath>
#include <iomanip>
#define NUM_CONTACTS 2
#define ZPLUS 1e-8
using namespace amrex;
using MatrixDType = amrex::GpuComplex<amrex::Real>;
using Matrix1D = TableData<MatrixDType, 1>;
using Matrix2D = TableData<MatrixDType, 2>;

template <typename T>
class TD;

struct CondensedHamiltonianElem
{
    MatrixDType Xi_s;
    MatrixDType Xi;
    MatrixDType Pi;
};

namespace MathConst
{
static constexpr amrex::Real pi =
    static_cast<amrex::Real>(3.14159265358979323846);
}

namespace PhysConst
{
static constexpr auto q_e = static_cast<amrex::Real>(1.602176634e-19);
static constexpr auto ep0 = static_cast<amrex::Real>(8.8541878128e-12);
//    static constexpr auto hbar  = static_cast<amrex::Real>( 1.054571817e-34 );
}  // namespace PhysConst

AMREX_GPU_HOST_DEVICE
MatrixDType get_beta(amrex::Real gamma, int M, int J)
{
    amrex::GpuComplex arg(0., -MathConst::pi * J / M);
    return 2 * gamma * cos(-1 * arg.imag());  // * exp(arg);
}

AMREX_GPU_HOST_DEVICE
MatrixDType conjugate(MatrixDType a)
{
    amrex::GpuComplex a_conj(a.real(), -1. * a.imag());
    return a_conj;
}

AMREX_GPU_HOST_DEVICE
MatrixDType Compute_SurfaceGreensFunction_Quadratic(MatrixDType E,
                                                    amrex::Real U,
                                                    MatrixDType beta,
                                                    amrex::Real gamma)
{
    MatrixDType EmU = E - U;

    MatrixDType EmU_sq = pow(EmU, 2.);

    amrex::Real gamma_sq = pow(gamma, 2.);

    MatrixDType Factor = EmU_sq + gamma_sq - pow(beta, 2.);

    MatrixDType Sqrt = sqrt(pow(Factor, 2.) - 4. * EmU_sq * gamma_sq);

    MatrixDType Denom = 2. * gamma_sq * EmU;

    amrex::Print() << "EmU: " << EmU << "\n";
    // amrex::Print() << "Factor: " << Factor << "\n";
    // amrex::Print() << "Sqrt: "  << Sqrt << "\n";
    // amrex::Print() << "Denom: " << Denom << "\n";
    // amrex::Print() << "Numerator: " << Factor-Sqrt << "\n";
    amrex::Print() << "Quadratic Value: " << (Factor - Sqrt) / Denom << "\n";

    return (Factor - Sqrt) / Denom;
}

AMREX_GPU_HOST_DEVICE
void CondenseHamiltonian(CondensedHamiltonianElem& cond, MatrixDType E,
                         amrex::Real U, MatrixDType beta, amrex::Real gamma)
{
    MatrixDType H_tilde_11 = pow(E - U, -1.);

    cond.Xi_s = U + beta * conjugate(beta) * H_tilde_11;

    cond.Xi = cond.Xi_s + gamma * conjugate(gamma) * H_tilde_11;

    cond.Pi = gamma * beta * H_tilde_11;

    amrex::Print() << "E: "    << E << "\n";
    amrex::Print() << "U: "    << U << "\n";
    amrex::Print() << "EmU: "  << E-U << "\n";
    amrex::Print() << "Xi_s: " << cond.Xi_s << "\n";
    amrex::Print() << "Xi: "   << cond.Xi << "\n";
    amrex::Print() << "Pi: "   << cond.Pi << "\n";
}

AMREX_GPU_HOST_DEVICE
MatrixDType Compute_SurfaceGreensFunction_EigenfunctionTechnique(
    MatrixDType E, CondensedHamiltonianElem& cond)
{
    MatrixDType EmXi = E - cond.Xi;

    MatrixDType Factor = EmXi * pow(cond.Pi, -1.);

    MatrixDType Eigval = 0.5 * Factor - 0.5 * sqrt(pow(Factor, 2.) - 4.);

    MatrixDType Sigma_dash = cond.Pi * Eigval;

    MatrixDType Denomenator = E - cond.Xi_s - Sigma_dash;

    amrex::Print() << "E: " << E << "\n";
    // amrex::Print() << "Factor Eig: " << Factor << "\n";
    // amrex::Print() << "Eigval: " << Eigval << "\n";
    amrex::Print() << "Eigenfunction Value: " << pow(Denomenator, -1.) << "\n";

    return pow(Denomenator, -1.);
}

AMREX_GPU_HOST_DEVICE
MatrixDType Compute_SurfaceGreensFunction_DecimationTechnique(
    MatrixDType E, const CondensedHamiltonianElem& cond)
{
    MatrixDType mu_0 = E - cond.Xi_s;

    MatrixDType nu_0 = E - cond.Xi;

    MatrixDType gamma_0 = cond.Pi;

    MatrixDType mu_old = mu_0 - gamma_0 * conjugate(gamma_0) / nu_0;

    MatrixDType nu_old = nu_0 - gamma_0 * conjugate(gamma_0) / nu_0;

    MatrixDType gamma_old = pow(gamma_0, 2.) / nu_0;

    MatrixDType zeta_old = pow(conjugate(gamma_0), 2.) / nu_0;

    MatrixDType mu_new = mu_old;

    MatrixDType nu_new = nu_old;

    MatrixDType gamma_new = gamma_old;

    MatrixDType zeta_new = zeta_old;

    for (int i = 1; i < 10; ++i)
    {
        mu_new = mu_old - gamma_old * zeta_old / nu_old;

        nu_new = nu_old - 2. * gamma_old * zeta_old / nu_old;

        gamma_new = pow(gamma_old, 2.) / nu_old;

        zeta_new = pow(zeta_old, 2.) / nu_old;

        mu_old = mu_new;

        nu_old = nu_new;

        gamma_old = gamma_new;

        zeta_old = zeta_new;
    }

    amrex::Print() << "E: " << E << "\n";
    amrex::Print() << "Decimation Value: " << pow(mu_new, -1.) << "\n";

    return pow(mu_new, -1.);
}

AMREX_GPU_HOST_DEVICE
MatrixDType get_Sigma(MatrixDType E, amrex::Real U, MatrixDType beta,
                      amrex::Real gamma)
{
    CondensedHamiltonianElem cond;
    CondenseHamiltonian(cond, E, U, beta, gamma);
    Compute_SurfaceGreensFunction_DecimationTechnique(E, cond);
    Compute_SurfaceGreensFunction_EigenfunctionTechnique(E, cond);

    MatrixDType gr = Compute_SurfaceGreensFunction_Quadratic(E, U, beta, gamma);

    MatrixDType Sigma = pow(gamma, 2.) * gr;

    // amrex::GpuComplex val(640.5, -3.1977);
    return Sigma;
}

AMREX_GPU_HOST_DEVICE
MatrixDType get_Gamma(MatrixDType Sigma)
{
    amrex::GpuComplex val(-2. * Sigma.imag(), 0.);
    return val;
}

void Write1DArrayVsE(const Matrix1D& E_data, Matrix1D& Arr_data,
                     std::string filename, std::string header)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\n Root Writing " << filename << "\n";
        std::ofstream outfile;
        outfile.open(filename.c_str());

        auto const& E = E_data.const_table();
        auto thi = E_data.hi();
        auto const& Arr = Arr_data.const_table();

        outfile << header << "\n";
        for (int e = 0; e < thi[0]; ++e)
        {
            outfile << std::setw(10) << E(e).real() << std::setw(15)
                    << Arr(e).real() << std::setw(15) << Arr(e).imag() << "\n";
        }

        outfile.close();
    }
}

void Write2DArrayVsE(const Matrix1D& E_data, Matrix2D& Arr_data,
                     std::string filename, std::string header)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\n Root Writing " << filename << "\n";
        std::ofstream outfile;
        outfile.open(filename.c_str());

        auto const& E = E_data.const_table();
        auto thi = E_data.hi();
        auto const& Arr = Arr_data.const_table();

        outfile << header << "\n";
        for (int e = 0; e < thi[0]; ++e)
        {
            outfile << std::setw(10) << E(e).real() << std::setw(15)
                    << Arr(0, e).real() << std::setw(15) << Arr(0, e).imag()
                    << std::setw(15) << Arr(1, e).real() << std::setw(15)
                    << Arr(1, e).imag() << "\n";
        }

        outfile.close();
    }
}

void PrintTable_loc(Matrix2D& G_blk)
{
    auto tlo = G_blk.lo();
    auto thi = G_blk.hi();

    auto const& G = G_blk.table();

    std::cout << "tlo: " << tlo[0] << " " << tlo[1] << "thi: " << thi[0] << " "
              << thi[1] << "\n";

    for (int i = tlo[0]; i < thi[0]; ++i)
    {
        for (int j = tlo[1]; j < thi[1];
             ++j)  // slow moving index. printing slow
        {
            std::cout << std::setw(12) << std::setprecision(2) << G(i, j);
        }
        std::cout << "\n";
    }
}

void PrintTable(Matrix2D& h_G_loc_data, int N_total,
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

    if (my_rank == ParallelDescriptor::IOProcessorNumber())
    {
        h_G_glo_data.resize({0, 0}, {N_total, N_total}, The_Pinned_Arena());
        h_G_glo = h_G_glo_data.table();
        recv_count = new int[num_proc];
        disp = new int[num_proc];
        for (int p = 0; p < num_proc; ++p)
        {
            if (p < num_proc_with_blk)
            {
                recv_count[p] = cumu_blk_size[p + 1] - cumu_blk_size[p];
                disp[p] = cumu_blk_size[p];
            }
            else
            {
                recv_count[p] = 0;
                disp[p] = 0;
            }
            // amrex::Print() << "p/recv_count/disp: " << p << " " <<
            // recv_count[p] << " " << disp[p] << "\n";
        }
    }

    MPI_Datatype column_type;

    MPI_Type_vector(1, N_total + 1, N_total + 1, MPI_DOUBLE_COMPLEX,
                    &column_type);

    MPI_Type_commit(&column_type);

    auto h_G_loc = h_G_loc_data.table();

    MPI_Gatherv(&h_G_loc(0, 0), num_cols_loc, column_type, &h_G_glo(0, 0),
                recv_count, disp, column_type,
                ParallelDescriptor::IOProcessorNumber(),
                ParallelDescriptor::Communicator());
    ParallelDescriptor::Barrier();

    if (my_rank == ParallelDescriptor::IOProcessorNumber())
    {
        delete[] recv_count;
        delete[] disp;
    }
    MPI_Type_free(&column_type);

    amrex::Real G_locglo_gatherv_time =
        amrex::second() - G_locglo_gatherv_time_beg;
    amrex::Print() << "\nG_locglo_gatherv_time: " << std::setw(15)
                   << G_locglo_gatherv_time << "\n";

    if (my_rank == ParallelDescriptor::IOProcessorNumber())
    {
        auto tlo = h_G_glo_data.lo();
        auto thi = h_G_glo_data.hi();

        auto const& G = h_G_glo_data.table();

        amrex::Print() << "tlo: " << tlo[0] << " " << tlo[1] << "\n"
                       << "thi: " << thi[0] << " " << thi[1] << "\n";

        for (int i = tlo[0]; i < thi[0]; ++i)
        {
            for (int j = tlo[1]; j < thi[1];
                 ++j)  // slow moving index. printing slow
            {
                amrex::Print()
                    << std::setw(12) << std::setprecision(2) << G(i, j);
            }
            amrex::Print() << "\n";
        }

        h_G_glo_data.clear();
    }
}

void Print_GreensAndSpectralFunction(Matrix2D& hd_G_loc_data,
                                     Matrix2D& hd_A_loc_data, int N_total,
                                     const amrex::Vector<int>& cumu_blk_size,
                                     const amrex::Vector<int>& vec_col_gids,
                                     const int num_proc_with_blk)
{
#ifdef AMREX_USE_GPU
    int num_cols_loc = vec_col_gids.size();
    /**copy G_loc from device to host*/
    Matrix2D h_G_loc_data({0, 0}, {N_total, num_cols_loc}, The_Pinned_Arena());
    Matrix2D h_A_loc_data({0, 0}, {N_total, num_cols_loc}, The_Pinned_Arena());
    h_G_loc_data.copy(hd_G_loc_data);  // copy from gpu to cpu
    h_A_loc_data.copy(hd_A_loc_data);  // copy from gpu to cpu
    amrex::Print() << "G_glo: \n";
    PrintTable(h_G_loc_data, N_total, cumu_blk_size, vec_col_gids,
               num_proc_with_blk);
    amrex::Print() << "A_glo: \n";
    PrintTable(h_A_loc_data, N_total, cumu_blk_size, vec_col_gids,
               num_proc_with_blk);
    h_G_loc_data.clear();
    h_A_loc_data.clear();
#else
    amrex::Print() << "G_glo: \n";
    PrintTable(hd_G_loc_data, N_total, cumu_blk_size, vec_col_gids,
               num_proc_with_blk);
    amrex::Print() << "A_glo: \n";
    PrintTable(hd_A_loc_data, N_total, cumu_blk_size, vec_col_gids,
               num_proc_with_blk);
#endif
}

void Obtain_GreensAndSpectralFunctions(
    int N_total, const Matrix1D& h_U_loc_data, const Matrix1D& h_B_data,
    const Matrix1D& h_C_data, const Matrix1D& h_E_glo_data,
    Matrix2D& h_Sigma_glo_data, Matrix2D& hd_G_loc_data,
    Matrix2D& hd_A_loc_data, const amrex::Vector<int>& cumu_blk_size,
    const amrex::Vector<int>& vec_col_gids, const int num_proc_with_blk,
    const amrex::GpuArray<int, NUM_CONTACTS>& global_contact_index,
    const amrex::GpuArray<int, NUM_CONTACTS>& contact_transmission_index,
    int print_matrix_flag)
{
    /* Six step procedure:
     * 1) MPI_Allgatherv local Hamiltonian diagonal elements in Alpha_loc_data
     * into Alpha_glo_data, both allocated on the PinnedArena. 2) Compute on all
     * processes, 1DTable arrays Xtil, Ytil, X, Y recursively, on the
     * PinnedArena and then delete Alpha_glo_data. These arrays have a size
     * equal to the column (or row) size of the full (or global) matrix to be
     * inverted. 3) All processes copy the respective portion of X and Y needed
     * later into arrays allocated on managed memory. These portions have a size
     * equal to the number of columns in the local matrix block. 4) Copy needed
     * arrays on to the device. These include Xtil, Ytil, X_local, Y_local,
     * Alpha_loc. 5) Launce parallelFor, assign each gpu-thread to compute a
     * column of the inverted matrix. Each GPU is responsible for computing it's
     * own block. 6) Do a stream synchronise and copy the inverted matrix block
     * on to the CPU.
     */
    int num_proc = ParallelDescriptor::NProcs();
    int my_rank = ParallelDescriptor::MyProc();
    int num_cols_loc = vec_col_gids.size();
    int ioproc = ParallelDescriptor::IOProcessorNumber();
    int num_traces = 2;

    auto thi_Sigma_glo = h_Sigma_glo_data.hi();
    int num_contacts = thi_Sigma_glo[0];
    int num_EnPts = thi_Sigma_glo[1];
    amrex::Print() << "num_contacts/num_EnPts: " << num_contacts << " "
                   << num_EnPts << "\n";

    int recv_count[num_proc];
    int disp[num_proc];

    for (int p = 0; p < num_proc; ++p)
    {
        if (p < num_proc_with_blk)
        {
            recv_count[p] = cumu_blk_size[p + 1] - cumu_blk_size[p];
            disp[p] = cumu_blk_size[p];
        }
        else
        {
            recv_count[p] = 0;
            disp[p] = 0;
        }
    }
    /*allocating local arrays on the host*/
    Matrix1D h_Alpha_loc_data({0}, {num_cols_loc}, The_Pinned_Arena());
    Matrix1D h_Xtil_glo_data({0}, {N_total}, The_Pinned_Arena());
    Matrix1D h_Ytil_glo_data({0}, {N_total}, The_Pinned_Arena());
    Matrix1D h_X_glo_data({0}, {N_total}, The_Pinned_Arena());
    Matrix1D h_Y_glo_data({0}, {N_total}, The_Pinned_Arena());
    Matrix1D h_Alpha_contact_data({0}, {num_contacts}, The_Pinned_Arena());

    Matrix1D h_X_loc_data({0}, {num_cols_loc}, The_Pinned_Arena());
    Matrix1D h_Y_loc_data({0}, {num_cols_loc}, The_Pinned_Arena());
    Matrix1D h_X_contact_data({0}, {num_contacts}, The_Pinned_Arena());
    Matrix1D h_Y_contact_data({0}, {num_contacts}, The_Pinned_Arena());
    Matrix1D h_Adiag_loc_data({0}, {num_cols_loc}, The_Pinned_Arena());

#ifdef AMREX_USE_GPU
    /*allocating local arrays on the device*/
    Matrix1D d_Alpha_loc_data({0}, {num_cols_loc}, The_Arena());
    Matrix1D d_Xtil_glo_data({0}, {N_total}, The_Arena());
    Matrix1D d_Ytil_glo_data({0}, {N_total}, The_Arena());
    Matrix1D d_X_loc_data({0}, {num_cols_loc}, The_Arena());
    Matrix1D d_Y_loc_data({0}, {num_cols_loc}, The_Arena());
    Matrix1D d_Alpha_contact_data({0}, {num_contacts}, The_Arena());
    Matrix1D d_X_contact_data({0}, {num_contacts}, The_Arena());
    Matrix1D d_Y_contact_data({0}, {num_contacts}, The_Arena());
    Matrix2D d_Sigma_glo_data({0, 0}, {num_contacts, num_EnPts}, The_Arena());
    Matrix1D d_Adiag_loc_data({0}, {num_cols_loc}, The_Arena());
    amrex::Gpu::DeviceVector<amrex::Real> d_Trace_r(num_traces,
                                                    amrex::Real(0.0));
    amrex::Gpu::DeviceVector<amrex::Real> d_Trace_i(num_traces,
                                                    amrex::Real(0.0));
#endif
    amrex::Gpu::HostVector<amrex::Real> h_Trace_r(num_traces);
    amrex::Gpu::HostVector<amrex::Real> h_Trace_i(num_traces);

    Matrix1D h_T_loc_data({0}, {num_EnPts}, The_Pinned_Arena());
    Matrix1D h_DOS_loc_data({0}, {num_EnPts}, The_Pinned_Arena());

    /*get the reference of tables on the host*/
    Table1D<MatrixDType> h_Alpha_loc = h_Alpha_loc_data.table();

    auto const& h_U_loc = h_U_loc_data.const_table();
    auto const& h_B = h_B_data.const_table();
    auto const& h_C = h_C_data.const_table();
    auto const& h_E_glo = h_E_glo_data.const_table();
    auto const& h_Ytil_glo = h_Ytil_glo_data.table();
    auto const& h_Xtil_glo = h_Xtil_glo_data.table();
    auto const& h_Y_glo = h_Y_glo_data.table();
    auto const& h_X_glo = h_X_glo_data.table();
    auto const& h_Alpha_contact = h_Alpha_contact_data.table();
    auto const& h_Y_loc = h_Y_loc_data.table();
    auto const& h_X_loc = h_X_loc_data.table();
    auto const& h_Y_contact = h_Y_contact_data.table();
    auto const& h_X_contact = h_X_contact_data.table();
    auto const& h_T_loc = h_T_loc_data.table();
    auto const& h_DOS_loc = h_DOS_loc_data.table();
    auto const& h_Sigma_glo = h_Sigma_glo_data.const_table();

#ifdef AMREX_USE_GPU
    /*get the reference of tables on the host*/
    auto const& Alpha = d_Alpha_loc_data.const_table();
    auto const& Xtil_glo = d_Xtil_glo_data.const_table();
    auto const& Ytil_glo = d_Ytil_glo_data.const_table();
    auto const& X = d_X_loc_data.const_table();
    auto const& Y = d_Y_loc_data.const_table();
    auto const& Alpha_contact = d_Alpha_contact_data.const_table();
    auto const& X_contact = d_X_contact_data.const_table();
    auto const& Y_contact = d_Y_contact_data.const_table();
    auto const& Sigma_glo = d_Sigma_glo_data.const_table();
    auto const& Adiag_loc = d_Adiag_loc_data.table();  // Spectral function

    amrex::Real* trace_r = d_Trace_r.dataPtr();
    amrex::Real* trace_i = d_Trace_i.dataPtr();

// auto const& G_contact     = d_G_contact_data.table();
#else
    auto const& Alpha = h_Alpha_loc_data.const_table();
    auto const& Xtil_glo = h_Xtil_glo_data.const_table();
    auto const& Ytil_glo = h_Ytil_glo_data.const_table();
    auto const& X = h_X_loc_data.const_table();
    auto const& Y = h_Y_loc_data.const_table();
    auto const& Alpha_contact = h_Alpha_contact_data.const_table();
    auto const& X_contact = h_X_contact_data.const_table();
    auto const& Y_contact = h_Y_contact_data.const_table();
    auto const& Sigma_glo = h_Sigma_glo_data.const_table();
    auto const& Adiag_loc = h_Adiag_loc_data.table();  // Spectral function

    amrex::Real* trace_r = h_Trace_r.dataPtr();
    amrex::Real* trace_i = h_Trace_i.dataPtr();
#endif
    auto const& G_loc = hd_G_loc_data.table();  // Green's function
    auto const& A_loc = hd_A_loc_data.table();  // Spectral function

    for (int e = 0; e < num_EnPts; ++e)
    {
        // BL_PROFILE_VAR("steps_1_2_3", steps_1_2_3);

        /*Step 1*/

        for (int n = 0; n < num_cols_loc; ++n)
        {
            h_Alpha_loc(n) = h_E_glo(e);  // - h_U_loc(n);
        }
        for (int c = 0; c < num_contacts; ++c)
        {
            int n_glo = global_contact_index[c];
            if (n_glo >= cumu_blk_size[my_rank] &&
                n_glo < cumu_blk_size[my_rank + 1])
            {
                int n = n_glo - cumu_blk_size[my_rank];
                h_Alpha_loc(n) -= h_Sigma_glo(c, e);
            };
        }

        Matrix1D h_Alpha_glo_data({0}, {N_total}, The_Pinned_Arena());
        Table1D<MatrixDType> h_Alpha_glo = h_Alpha_glo_data.table();

        MPI_Allgatherv(&h_Alpha_loc(0), num_cols_loc, MPI_DOUBLE_COMPLEX,
                       &h_Alpha_glo(0), recv_count, disp, MPI_DOUBLE_COMPLEX,
                       ParallelDescriptor::Communicator());

        for (int c = 0; c < num_contacts; ++c)
        {
            int n_glo = global_contact_index[c];
            h_Alpha_contact(c) = h_Alpha_glo(n_glo);
        }
        // amrex::Print() << "n & Alpha_glo: " << "\n";
        // for(int n=0; n<N_total; ++n)
        //{
        //     amrex::Print() << n << " "<< h_Alpha_glo(n) << "\n";
        // }

        /*Step 2*/

        h_Y_glo(0) = 0;
        for (int n = 1; n < N_total; ++n)
        {
            int p = (n - 1) % 2;
            h_Ytil_glo(n) = h_C(p) / (h_Alpha_glo(n - 1) - h_Y_glo(n - 1));
            h_Y_glo(n) = h_B(p) * h_Ytil_glo(n);
        }

        h_X_glo(N_total - 1) = 0;
        for (int n = N_total - 2; n > -1; n--)
        {
            int p = n % 2;
            h_Xtil_glo(n) = h_B(p) / (h_Alpha_glo(n + 1) - h_X_glo(n + 1));
            h_X_glo(n) = h_C(p) * h_Xtil_glo(n);
        }

        h_Alpha_glo_data.clear();

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

        // BL_PROFILE_VAR_STOP(steps_1_2_3);

        for (int t = 0; t < num_traces; ++t)
        {
            h_Trace_r[t] = 0.;
            h_Trace_i[t] = 0.;
        }

/*Step 4*/
#ifdef AMREX_USE_GPU
        // BL_PROFILE_VAR("step4_CpuToGpuCopy_time", step4_CpuToGpuCopy_time);

        d_Alpha_loc_data.copy(h_Alpha_loc_data);
        d_Xtil_glo_data.copy(h_Xtil_glo_data);
        d_Ytil_glo_data.copy(h_Ytil_glo_data);
        d_X_loc_data.copy(h_X_loc_data);
        d_Y_loc_data.copy(h_Y_loc_data);
        d_Alpha_contact_data.copy(h_Alpha_contact_data);
        d_X_contact_data.copy(h_X_contact_data);
        d_Y_contact_data.copy(h_Y_contact_data);
        d_Sigma_glo_data.copy(h_Sigma_glo_data);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_Trace_r.begin(),
                         h_Trace_r.end(), d_Trace_r.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_Trace_i.begin(),
                         h_Trace_i.end(), d_Trace_i.begin());

        amrex::Gpu::streamSynchronize();

// BL_PROFILE_VAR_STOP(step4_CpuToGpuCopy_time);
#endif

        /*Step 5*/
        // BL_PROFILE_VAR("step5_pFor_time", step5_pFor_time);

        int cumulative_columns = cumu_blk_size[my_rank];
        int e_id = e;
        MatrixDType E = h_E_glo(e_id);

        amrex::ParallelFor(
            num_cols_loc,
            [=] AMREX_GPU_DEVICE(int n) noexcept
            {
                int n_glo = n + cumulative_columns; /*global column number*/

                G_loc(n_glo, n) = 1. / (Alpha(n) - X(n) - Y(n));

                for (int m = n_glo; m > 0; m--)
                {
                    G_loc(m - 1, n) = -Ytil_glo(m) * G_loc(m, n);
                }
                for (int m = n_glo; m < N_total - 1; ++m)
                {
                    G_loc(m + 1, n) = -Xtil_glo(m) * G_loc(m, n);
                }

                // MatrixDType G_contact_kk contains G_glo(k_glo,k_glo) element
                // MatrixDType G_contact_nk contains G_glo(n_glo,k_glo) element
                MatrixDType A_tk[NUM_CONTACTS];
                for (int m = 0; m < N_total; ++m)
                {
                    A_loc(m, n) = 0.;
                }
                for (int k = 0; k < num_contacts; ++k)
                {
                    int k_glo = global_contact_index[k];
                    MatrixDType G_contact_kk =
                        1. / (Alpha_contact(k) - X_contact(k) - Y_contact(k));

                    MatrixDType temp = G_contact_kk;
                    for (int m = k_glo; m < n_glo; ++m)
                    {
                        temp = -Xtil_glo(m) * temp;
                    }
                    for (int m = k_glo; m > n_glo; m--)
                    {
                        temp = -Ytil_glo(m) * temp;
                    }
                    MatrixDType G_contact_nk = temp;

                    MatrixDType A_kn = G_contact_kk *
                                       get_Gamma(Sigma_glo(k, e_id)) *
                                       conjugate(G_contact_nk);

                    A_loc(k_glo, n) += A_kn;
                    for (int m = k_glo + 1; m < N_total; ++m)
                    {
                        A_kn = -Xtil_glo(m - 1) * A_kn;
                        A_loc(m, n) += A_kn;
                    }
                    for (int m = k_glo - 1; m >= 0; m--)
                    {
                        A_kn = -Ytil_glo(m + 1) * A_kn;
                        A_loc(m, n) += A_kn;
                    }
                    A_tk[k] = 0.;
                    if (n_glo == contact_transmission_index[k])
                    {
                        A_tk[k] = A_kn;
                    }
                }

                MatrixDType T12 = get_Gamma(Sigma_glo(0, e_id)) * A_tk[1];

                Adiag_loc(n) = A_loc(n_glo, n);

                amrex::HostDevice::Atomic::Add(&(trace_r[0]),
                                               Adiag_loc(n).real());
                amrex::HostDevice::Atomic::Add(&(trace_i[0]),
                                               Adiag_loc(n).imag());

                amrex::HostDevice::Atomic::Add(&(trace_r[1]), T12.real());
                amrex::HostDevice::Atomic::Add(&(trace_i[1]), T12.imag());
            });

        // BL_PROFILE_VAR_STOP(step5_pFor_time);

#ifdef AMREX_USE_GPU
        amrex::Gpu::streamSynchronize();
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_Trace_r.begin(),
                         d_Trace_r.end(), h_Trace_r.begin());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_Trace_i.begin(),
                         d_Trace_i.end(), h_Trace_i.begin());
#endif

        for (int t = 0; t < num_traces; ++t)
        {
            amrex::ParallelDescriptor::ReduceRealSum(h_Trace_r[t]);
            amrex::ParallelDescriptor::ReduceRealSum(h_Trace_i[t]);
            // amrex::Print() << "Trace: "<< t << " " << h_Trace_r[t] << " "<<
            // h_Trace_i[t]<< "\n";
        }

        MatrixDType DOS_temp(h_Trace_r[0], h_Trace_i[0]);
        h_DOS_loc(e_id) = DOS_temp / (2. * MathConst::pi);

        MatrixDType T12_temp(h_Trace_r[1], h_Trace_i[1]);
        h_T_loc(e_id) = T12_temp;

        if (print_matrix_flag)
        {
            if (e == 0)
            {
                amrex::Print()
                    << "Printing G_R and A at E(e)=" << h_E_glo(e) << "\n";
                Print_GreensAndSpectralFunction(hd_G_loc_data, hd_A_loc_data,
                                                N_total, cumu_blk_size,
                                                vec_col_gids,
                                                num_proc_with_blk);
            }
        }
    }

    Write1DArrayVsE(h_E_glo_data, h_T_loc_data, "output/Transmission",
                    "E   T_r   T_i");
    Write1DArrayVsE(h_E_glo_data, h_DOS_loc_data, "output/DOS",
                    "E  DOS_r  DOS_i");

/*deallocate memory*/
#ifdef AMREX_USE_GPU
    d_Alpha_loc_data.clear();
    d_Xtil_glo_data.clear();
    d_Ytil_glo_data.clear();
    d_X_loc_data.clear();
    d_Y_loc_data.clear();
    d_Alpha_contact_data.clear();
    d_X_contact_data.clear();
    d_Y_contact_data.clear();
    d_Adiag_loc_data.clear();
    d_Trace_r.clear();
    d_Trace_i.clear();
#endif
    h_Alpha_loc_data.clear();
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
    h_DOS_loc_data.clear();
    h_Trace_r.clear();
    h_Trace_i.clear();
}

#ifdef ASSERT
void Check_SurfaceGreensFunctionAsserts()
{
    MatrixDType gr1 = Compute_SurfaceGreensFunction(0.1, 0, 0.01, 2.5);
    MatrixDType gr2 = Compute_SurfaceGreensFunction(0.1, 0, 0.015, 2.5);

    std::string assrt_msg = "Real parts are not equal: " + std::to_string(4.) +
                            " " + std::to_string(gr1.real()) + " " +
                            std::to_string(gr2.real()) + "\n";

    amrex::Print() << "gr1: " << gr1.real() << " " << gr1.imag() << "\n";
    amrex::Print() << "gr2: " << gr2.real() << " " << gr2.imag() << "\n";

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(gr1.real() == gr2.real(), assrt_msg);

    assrt_msg = "Imag parts are not equal: " + std::to_string(gr1.imag()) +
                " " + std::to_string(gr2.imag()) + "\n";

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(gr1.imag() == gr2.imag(), assrt_msg);
}
#endif

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

#ifdef ASSERT
    Check_SurfaceGreensFunctionAsserts();
#endif

    int num_proc = ParallelDescriptor::NProcs();
    int my_rank = ParallelDescriptor::MyProc();
    amrex::Print() << "total number of procs: " << num_proc << "\n";

    // Construct Tridiagonal Dummy Hamiltonian
    std::array<amrex::Real, 3> point_charge_loc{0., 0., 1e-9};

    int N_total = 4; /*matrix size*/
    amrex::ParmParse pp;
    pp.query("N_total", N_total);

    // ParallelDescriptor::Bcast(&N_total, 1,
    // ParallelDescriptor::IOProcessorNumber());
    if (my_rank == 1)
    {
        std::cout << "Rank 1 reporting, total matrix size: " << N_total << "\n";
    }

    amrex::Print() << "total matrix size: " << N_total << "\n";
    const int MAX_THRESHOLD_BLK_SIZE = 40000; /*matrix size*/

    bool flag_fixed_blk_size = false;

    int max_blk_size =
        2; /*max block size, may change depending on flag_fixed_blk_size*/

    if (flag_fixed_blk_size)
        amrex::Print() << "max_blk_size is fixed by user\n";
    else
    {
        amrex::Print() << "max_blk_size is computed at run-time\n";
        max_blk_size = ceil(static_cast<amrex::Real>(N_total) / num_proc);
        if (max_blk_size > MAX_THRESHOLD_BLK_SIZE)
        {
            max_blk_size = MAX_THRESHOLD_BLK_SIZE;
            /*assert that use larger number of procs*/
        }
    }
    amrex::Print() << "maximum block size (in number of columns): "
                   << max_blk_size << "\n";

    const int num_blk_glo =
        ceil(static_cast<amrex::Real>(N_total) / max_blk_size);
    amrex::Print() << "global number of blocks: " << num_blk_glo << "\n";

    const int num_proc_with_blk = num_blk_glo;
    /*if num_proc_with_blk >= num_proc, assert.*/

    amrex::Vector<int> cumu_blk_size(num_proc_with_blk + 1);

    cumu_blk_size[0] = 0;
    for (int p = 1; p < num_proc_with_blk; ++p)
    {
        cumu_blk_size[p] = cumu_blk_size[p - 1] + max_blk_size;
        /*note that all proc except the last proc in "all procs with blocks"
         * have N_blk number of columns. The last proc in "all procs with
         * blocks" may have number of columns smaller than N_blk.
         */
    }
    cumu_blk_size[num_proc_with_blk] = cumu_blk_size[num_proc_with_blk - 1] +
                                       N_total -
                                       (num_proc_with_blk - 1) * max_blk_size;

    // amrex::Print() << "cumulative blk size: \n";
    // for(int p=0; p < num_proc_with_blk+1; ++p)
    //{
    //     amrex::Print() << cumu_blk_size[p] << "\n";
    // }

    amrex::Vector<int> vec_col_gids;

    if (my_rank < num_proc_with_blk)
    {
        int blk_gid = my_rank;
        int blk_size = cumu_blk_size[blk_gid + 1] - cumu_blk_size[blk_gid];

        for (int c = 0; c < blk_size; ++c)
        {
            int col_gid = cumu_blk_size[blk_gid] + c;
            vec_col_gids.push_back(col_gid);
        }

        // if(my_rank == 1) {
        //     std::cout      << "(proc) My rank: " << my_rank << "\n";
        //     std::cout << " blk_size: " << vec_col_gids.size() << "\n";
        //     std::cout << " vec_col_gids: ";
        //     for(auto& col_gid: vec_col_gids) {
        //        std::cout << std::setw(5) << col_gid << " ";
        //     }
        //     std::cout <<  "\n ";
        // }
    }

    /*setting A, the diagonal 1D vector, which is partitioned in the form of
     * some blocks*/
    /*in general A can be 2D*/
    int num_cols_loc = vec_col_gids.size();
    Matrix1D h_U_loc_data({0}, {num_cols_loc}, The_Pinned_Arena());

    /*R this is the size of repeated B and C vectors in Hamiltonian in the
     * mode-space approximation for CNT*/
    const int R = 2;

    /*In general, B and C are 2D blocks of different sizes 2D blocks*/
    Matrix1D h_B_data({0}, {R}, The_Pinned_Arena());
    Matrix1D h_C_data({0}, {R}, The_Pinned_Arena());
    auto const& h_B = h_B_data.table();
    auto const& h_C = h_C_data.table();

    /*specifying A block vector for CNT with point charge potential*/

    auto blk_size = vec_col_gids.size();
    int c = 0;
    auto const& h_U_loc = h_U_loc_data.table();

    for (auto& col_gid : vec_col_gids)
    {
        amrex::Real layer_loc_y =
            -10e-9 +
            (static_cast<amrex::Real>(col_gid) / (N_total - 1)) * 20e-9;

        amrex::Real r = pow((pow((0. - point_charge_loc[0]), 2) +
                             pow((layer_loc_y - point_charge_loc[1]), 2) +
                             pow((0. - point_charge_loc[2]), 2)),
                            0.5);

        h_U_loc(c) =
            -(PhysConst::q_e) / (4. * MathConst::pi * PhysConst::ep0 * r);

        // amrex::Print() << "(proc 0) column, A: " << col_gid << " "<< A_loc(c)
        // << "\n";
        ++c;
    }

    /*below are some parameters for (17,0) carbon nanotube*/
    amrex::Real gamma = 2.5;  // eV
    int M = 17;
    MatrixDType beta = get_beta(gamma, M, 6);
    amrex::Print() << "beta: " << beta << "\n";
    amrex::Real U_contact[NUM_CONTACTS] = {0., 0.};

    /*specifying B and C vectors for the special case of mode-space Hamiltonian
     * of CNT*/
    for (std::size_t i = 0; i < R; ++i)
    {
        if (i % 2 == 0)
        {
            h_B(i) = -beta; /*negative sign because (E[I] - [H]) will have
                               negative B and C*/
            h_C(i) = -beta;
        }
        else
        {
            h_B(i) = -gamma;
            h_C(i) = -gamma;
        }
        // amrex::Print() << "i,B,C: " << i << " "<< h_B(i) << " " << h_C(i) <<
        // "\n";
    }

    /*define energy grid and sigma*/
    int num_contacts = NUM_CONTACTS;
    amrex::GpuArray<int, NUM_CONTACTS> global_contact_index = {0, N_total - 1};
    amrex::GpuArray<int, NUM_CONTACTS> contact_transmission_index = {N_total -
                                                                         1,
                                                                     0};
    int num_EnPts = 10;
    pp.query("num_EnPts", num_EnPts);

    amrex::Real EnRange[2] = {-6., 0.};  // eV
    Matrix1D h_E_glo_data({0}, {num_EnPts}, The_Pinned_Arena());
    auto const& h_E_glo = h_E_glo_data.table();

    amrex::Real deltaE = (EnRange[1] - EnRange[0]) / (num_EnPts - 1);

    for (std::size_t e = 0; e < num_EnPts; ++e)
    {
        amrex::GpuComplex E(EnRange[0] + e * deltaE, ZPLUS);
        h_E_glo(e) = E;
    }

    Matrix2D h_Sigma_glo_data({0, 0}, {num_contacts, num_EnPts},
                              The_Pinned_Arena());
    auto const& h_Sigma = h_Sigma_glo_data.table();

    for (std::size_t c = 0; c < num_contacts; ++c)
    {
        amrex::Print() << "contact: " << c << "\n";
        for (std::size_t e = 0; e < num_EnPts; ++e)
        {
            h_Sigma(c, e) = get_Sigma(h_E_glo(e), U_contact[c], beta, gamma);
        }
    }
    Write2DArrayVsE(h_E_glo_data, h_Sigma_glo_data, "output/Sigma",
                    "E Sigma_s_r Sigma_s_i Sigma_d_r Sigma_d_i");

    // PrintTable_loc(h_Sigma_glo_data);

    int print_matrix_flag = false;
    pp.query("print_matrix", print_matrix_flag);

/*allocate local G*/
#ifdef AMREX_USE_GPU
    Matrix2D d_G_loc_data({0, 0}, {N_total, num_cols_loc}, The_Arena());
    Matrix2D d_A_loc_data({0, 0}, {N_total, num_cols_loc}, The_Arena());

    Obtain_GreensAndSpectralFunctions(N_total, h_U_loc_data, h_B_data, h_C_data,
                                      h_E_glo_data, h_Sigma_glo_data,
                                      d_G_loc_data, d_A_loc_data, cumu_blk_size,
                                      vec_col_gids, num_proc_with_blk,
                                      global_contact_index,
                                      contact_transmission_index,
                                      print_matrix_flag);
#else
    Matrix2D h_G_loc_data({0, 0}, {N_total, num_cols_loc}, The_Pinned_Arena());
    Matrix2D h_A_loc_data({0, 0}, {N_total, num_cols_loc}, The_Pinned_Arena());

    Obtain_GreensAndSpectralFunctions(N_total, h_U_loc_data, h_B_data, h_C_data,
                                      h_E_glo_data, h_Sigma_glo_data,
                                      h_G_loc_data, h_A_loc_data, cumu_blk_size,
                                      vec_col_gids, num_proc_with_blk,
                                      global_contact_index,
                                      contact_transmission_index,
                                      print_matrix_flag);
#endif

    // if(print_matrix_flag) {
    //     /**copy G_loc from device to host*/
    //
    //     #ifdef AMREX_USE_GPU
    //     BL_PROFILE_VAR("step6_CopyGpuToCpu", step6_copyGpuToCpu);
    //
    //     Matrix2D h_G_loc_data({0,0},{N_total, num_cols_loc},
    //     The_Pinned_Arena()); Matrix2D h_A_loc_data({0,0},{N_total,
    //     num_cols_loc}, The_Pinned_Arena()); h_G_loc_data.copy(d_G_loc_data);
    //     //copy from gpu to cpu h_A_loc_data.copy(d_A_loc_data); //copy from
    //     gpu to cpu

    //    BL_PROFILE_VAR_STOP(step6_copyGpuToCpu);
    //    #endif
    //
    // amrex::Print() << "G_glo: \n";
    //    PrintTable(h_G_loc_data,N_total,cumu_blk_size, vec_col_gids,
    //    num_proc_with_blk);
    // amrex::Print() << "A_glo: \n";
    //    PrintTable(h_A_loc_data,N_total,cumu_blk_size, vec_col_gids,
    //    num_proc_with_blk);

    //    #ifdef AMREX_USE_GPU
    //    h_G_loc_data.clear();
    //    h_A_loc_data.clear();
    //    #endif
    //}

    /*deallocate memory*/
    h_U_loc_data.clear();
    h_B_data.clear();
    h_C_data.clear();
#ifdef AMREX_USE_GPU
    d_G_loc_data.clear();
    d_A_loc_data.clear();
#else
    h_G_loc_data.clear();
    h_A_loc_data.clear();
#endif
    h_Sigma_glo_data.clear();
    h_E_glo_data.clear();

    amrex::Finalize();
}
