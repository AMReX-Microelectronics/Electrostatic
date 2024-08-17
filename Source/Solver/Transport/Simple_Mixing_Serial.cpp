#include "Transport.H"

using namespace amrex;

#ifndef BROYDEN_PARALLEL
void c_TransportSolver::Execute_Simple_Mixing_Algorithm()
{
    /*update h_RhoInduced_glo*/

    if (ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\nBroydenStep: " << Broyden_Step
                       << ",  fraction: " << Broyden_fraction
                       << ",  scalar: " << Broyden_Scalar << "\n";

        auto const &n_curr_in = h_n_curr_in_data.table();
        auto const &n_curr_out = h_n_curr_out_data.table();
        auto const &n_prev_in = h_n_prev_in_data.table();
        auto const &F_curr = h_F_curr_data.table();
        auto const &Norm = h_Norm_data.table();

        SetVal_RealTable1D(h_Norm_data, 0.);
        Broyden_NormSum_Curr = 0.;

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

        // if ((Broyden_NormSum_Curr - Broyden_NormSum_Prev) > 1.e-6)
        //{
        //     Broyden_NormSumIsIncreasing_Step +=1;
        //     amrex::Print() << "\nBroyden_NormSumIsIncreasing_Step: " <<
        //     Broyden_NormSumIsIncreasing_Step << "\n";
        // }
        // else
        //{
        //     Broyden_NormSumIsIncreasing_Step = 0;
        // }

        Broyden_NormSum_Prev = Broyden_NormSum_Curr;

        for (int l = 0; l < num_field_sites_all_NS; ++l)
        {
            F_curr(l) = n_curr_in(l) - n_curr_out(l);
        }

        for (int l = 0; l < num_field_sites_all_NS; ++l)
        {
            n_prev_in(l) = n_curr_in(l);
            n_curr_in(l) = n_prev_in(l) - Broyden_fraction * F_curr(l);
        }

        Broyden_Step += 1;
    }

    MPI_Bcast(&Broyden_Norm, 1, MPI_DOUBLE,
              ParallelDescriptor::IOProcessorNumber(),
              ParallelDescriptor::Communicator());
}
#endif
