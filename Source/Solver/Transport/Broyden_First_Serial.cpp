#include "Transport.H"

using namespace amrex;

#ifndef BROYDEN_PARALLEL
void c_TransportSolver::Execute_Broyden_First_Algorithm()
{
    if (ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\nBroydenStep: " << Broyden_Step
                       << ",  fraction: " << Broyden_fraction
                       << ",  scalar: " << Broyden_Scalar << "\n";

        auto const &n_curr_in = h_n_curr_in_data.table();
        auto const &n_curr_out = h_n_curr_out_data.table();
        auto const &n_prev_in = h_n_prev_in_data.table();
        auto const &F_curr = h_F_curr_data.table();
        auto const &delta_F_curr = h_delta_F_curr_data.table();
        auto const &delta_n_curr = h_delta_n_curr_data.table();
        auto const &Jinv_curr = h_Jinv_curr_data.table();
        auto const &Norm = h_Norm_data.table();

        RealTable1D h_sum_Fcurr_data({0}, {num_field_sites_all_NS},
                                     The_Pinned_Arena());
        RealTable1D h_sum_deltaFcurr_data({0}, {num_field_sites_all_NS},
                                          The_Pinned_Arena());
        RealTable1D h_delta_n_Jinv_data({0}, {num_field_sites_all_NS},
                                        The_Pinned_Arena());

        auto const &sum_Fcurr = h_sum_Fcurr_data.table();
        auto const &sum_deltaFcurr = h_sum_deltaFcurr_data.table();
        auto const &delta_n_Jinv = h_delta_n_Jinv_data.table();

        SetVal_RealTable1D(h_Norm_data, 0.);
        Broyden_NormSum_Curr = 0.;

        for (int l = 0; l < num_field_sites_all_NS; ++l)
        {
            sum_Fcurr(l) = 0;
            sum_deltaFcurr(l) = 0;
            delta_n_Jinv(l) = 0.;
        }

        switch (map_NormType[Broyden_Norm_Type])
        {
            case s_Norm::Type::Absolute:
            {
                for (int l = 0; l < num_field_sites_all_NS; ++l)
                {
                    amrex::Real Fcurr = n_curr_in(l) - n_curr_out(l);
                    Norm(l) = fabs(Fcurr);
                    Broyden_NormSum_Curr += pow(Fcurr, 2);
                }
                break;
            }
            case s_Norm::Type::Relative:
            {
                for (int l = 0; l < num_field_sites_all_NS; ++l)
                {
                    amrex::Real Fcurr = n_curr_in(l) - n_curr_out(l);
                    Norm(l) = fabs(Fcurr / (n_curr_in(l) + n_curr_out(l)));
                    Broyden_NormSum_Curr += pow(Norm(l), 2);
                }
                break;
            }
            default:
            {
                amrex::Abort("Norm Type " + Broyden_Norm_Type +
                             " is not yet defined.");
            }
        }
        Broyden_NormSum_Curr = sqrt(Broyden_NormSum_Curr);

        Broyden_Norm = Norm(0);
        int norm_index = 0;
        for (int l = 1; l < num_field_sites_all_NS; ++l)
        {
            if (Broyden_Norm < Norm(l))
            {
                Broyden_Norm = Norm(l);
                norm_index = l;
            }
        }
        amrex::Print() << "\nBroyden_NormSum_Curr: " << std::setw(20)
                       << Broyden_NormSum_Curr << "\n";
        amrex::Print() << "Broyden_NormSum_Prev: " << std::setw(20)
                       << Broyden_NormSum_Prev << ",   Difference: "
                       << (Broyden_NormSum_Curr - Broyden_NormSum_Prev) << "\n";
        amrex::Print() << "Broyden max norm: " << Broyden_Norm
                       << " at location: " << norm_index << "\n\n";

        Broyden_NormSum_Prev = Broyden_NormSum_Curr;

        for (int l = 0; l < num_field_sites_all_NS; ++l)
        {
            amrex::Real Fcurr = n_curr_in(l) - n_curr_out(l);
            delta_F_curr(l) = Fcurr - F_curr(l);
            F_curr(l) = Fcurr;
            delta_n_curr(l) = n_curr_in(l) - n_prev_in(l);
        }

        amrex::Print() << "n_curr_in, n_prev_in: " << n_curr_in(0) << " "
                       << n_prev_in(0) << "\n";

        int m = Broyden_Step - 1;
        if (m > 0)
        {
            for (int a = 0; a < num_field_sites_all_NS; ++a)
            {
                amrex::Real sum = 0.;
                for (int b = 0; b < num_field_sites_all_NS; ++b)
                {
                    sum += Jinv_curr(a, b) * delta_F_curr(b);
                }
                sum_deltaFcurr(a) = sum;
            }

            amrex::Real denom = 0.;
            for (int l = 0; l < num_field_sites_all_NS; ++l)
            {
                denom += delta_n_curr(l) * sum_deltaFcurr(l);
            }

            for (int b = 0; b < num_field_sites_all_NS; ++b)
            {
                amrex::Real sum = 0.;
                for (int a = 0; a < num_field_sites_all_NS; ++a)
                {
                    sum += delta_n_curr(a) * Jinv_curr(a, b);
                }
                delta_n_Jinv(b) = sum;
            }

            for (int a = 0; a < num_field_sites_all_NS; ++a)
            {
                for (int b = 0; b < num_field_sites_all_NS; ++b)
                {
                    Jinv_curr(a, b) += (delta_n_curr(a) - sum_deltaFcurr(a)) *
                                       delta_n_Jinv(b) / denom;
                }
            }
        }
        for (int a = 0; a < num_field_sites_all_NS; ++a)
        {
            amrex::Real sum = 0.;
            for (int b = 0; b < num_field_sites_all_NS; ++b)
            {
                sum += Jinv_curr(a, b) * F_curr(b);
            }
            sum_Fcurr(a) = sum;
        }

        for (int l = 0; l < num_field_sites_all_NS; ++l)
        {
            n_prev_in(l) = n_curr_in(l);
            n_curr_in(l) = n_prev_in(l) - sum_Fcurr(l);
        }
        amrex::Print() << "n_new_in: " << n_curr_in(0) << "\n";

        h_sum_Fcurr_data.clear();
        h_sum_deltaFcurr_data.clear();
        h_delta_n_Jinv_data.clear();

        Broyden_Step += 1;
    }

    MPI_Bcast(&Broyden_Norm, 1, MPI_DOUBLE,
              ParallelDescriptor::IOProcessorNumber(),
              ParallelDescriptor::Communicator());
}
#endif
