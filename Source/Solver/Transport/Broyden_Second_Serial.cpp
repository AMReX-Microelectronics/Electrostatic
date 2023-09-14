#include "Transport.H"
#include "Broyden_Namespace.H"

using namespace amrex;
using namespace Broyden;

#ifndef BROYDEN_PARALLEL
void 
c_TransportSolver::Execute_Broyden_Modified_Second_Algorithm()
{

    if (ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\nBroydenStep: " << Broyden_Step  << ",  fraction: " << Broyden_fraction << ",  scalar: " << Broyden_Scalar<< "\n";

        auto const& n_curr_in  = h_n_curr_in_data.table();
        auto const& n_curr_out = h_n_curr_out_data.table();
        auto const& n_prev_in  = h_n_prev_in_data.table();
        auto const& F_curr     = h_F_curr_data.table();
        auto const& delta_F_curr   = h_delta_F_curr_data.table();
        auto const& delta_n_curr   = h_delta_n_curr_data.table();

        SetVal_RealTable1D(h_Norm_data, 0.);
        auto const& Norm       = h_Norm_data.table();

        RealTable1D sum_Fcurr_data({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        RealTable1D sum_deltaFcurr_data({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        RealTable1D W_curr_data({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        RealTable1D V_curr_data({0},{num_field_sites_all_NS}, The_Pinned_Arena());

        auto const& sum_Fcurr      = sum_Fcurr_data.table();
        auto const& sum_deltaFcurr = sum_deltaFcurr_data.table();
        auto const& W_curr   = W_curr_data.table();
        auto const& V_curr   = V_curr_data.table();

        amrex::Real denom = 0.;
        Broyden_NormSum_Curr = 0.;

        for(int l=0; l < num_field_sites_all_NS; ++l)
        {
            sum_deltaFcurr(l) = 0;
            sum_Fcurr(l) = 0;
            W_curr(l) = 0;
            V_curr(l) = 0;
        }

        switch(map_NormType[Broyden_Norm_Type])
        {
            case s_Norm::Type::Absolute:
            {
                for(int l=0; l < num_field_sites_all_NS; ++l)
                {
                    amrex::Real Fcurr = n_curr_in(l) - n_curr_out(l);
                    Norm(l) = fabs(Fcurr);
                    Broyden_NormSum_Curr += pow(Fcurr,2);
                }
                break;
            }
            case s_Norm::Type::Relative:
            {
                for(int l=0; l < num_field_sites_all_NS; ++l)
                {
                    amrex::Real Fcurr = n_curr_in(l) - n_curr_out(l);
                    Norm(l) = fabs(Fcurr/(n_curr_in(l) + n_curr_out(l)));
                    Broyden_NormSum_Curr += pow(Norm(l),2);
                }
                break;
            }
            default:
            {
                amrex::Abort("Norm Type " + Broyden_Norm_Type + " is not yet defined.");
            }
        }

        Broyden_NormSum_Curr = sqrt(Broyden_NormSum_Curr);

        Broyden_Norm = Norm(0);
        int norm_index = 0;
        for(int l=1; l < num_field_sites_all_NS; ++l)
        {
            if(Broyden_Norm < Norm(l))
            {
                Broyden_Norm = Norm(l);
                norm_index = l;
            }
        }
        amrex::Print() << "\nBroyden_NormSum_Curr: " << std::setw(20) << Broyden_NormSum_Curr << "\n";
        amrex::Print() <<   "Broyden_NormSum_Prev: " << std::setw(20) << Broyden_NormSum_Prev
                       << ",   Difference: " << (Broyden_NormSum_Curr - Broyden_NormSum_Prev) << "\n";
        amrex::Print() << "Broyden max norm: " << Broyden_Norm << " at location: " <<  norm_index << "\n\n";


        //if ((Broyden_NormSum_Curr - Broyden_NormSum_Prev) > 1.e-6)
        //{
        //    Broyden_NormSumIsIncreasing_Step +=1;
        //    amrex::Print() << "\nBroyden_NormSumIsIncreasing_Step: " << Broyden_NormSumIsIncreasing_Step << "\n";
        //}
        //else
        //{
        //    Broyden_NormSumIsIncreasing_Step = 0;
        //}
        //if (Broyden_Step > Broyden_Threshold_MaxStep)
        //{
        //      Restart_Broyden();
        //}
        //else if (Broyden_NormSumIsIncreasing_Step > Broyden_Threshold_NormSumIncreaseStep)
        //{
        //    for(int l=0; l < num_field_sites_all_NS; ++l)
        //    {
        //        n_curr_in(l) = n_prev_in(l);
        //        n_prev_in(l) = n_curr_in(l) - delta_n_curr(l);
        //        denom += pow(delta_F_curr(l),2.);
        //    }

        //    Broyden_Correction_Step += 1;

        //    if(Broyden_Correction_Step > Broyden_Threshold_CorrectionStep)
        //    {
        //      Restart_Broyden();
        //    }
        //    else
        //    {
        //        Broyden_Scalar = Broyden_Scalar/std::pow(Broyden_Scalar_Decrease_Factor,Broyden_Correction_Step);
        //        Broyden_Step -= 1;
        //        Broyden_NormSumIsIncreasing_Step = 0;
        //        amrex::Print() << "\n******Reducing Broyden scalar to: " << Broyden_Scalar << "\n";
        //    }
        //}
        //else
        //{
            Broyden_NormSum_Prev = Broyden_NormSum_Curr;

            for(int l=0; l < num_field_sites_all_NS; ++l)
            {
                amrex::Real Fcurr = n_curr_in(l) - n_curr_out(l);
                delta_F_curr(l) = Fcurr - F_curr(l);
                F_curr(l) = Fcurr;
                denom += pow(delta_F_curr(l),2.);
                delta_n_curr(l) = n_curr_in(l) - n_prev_in(l);
            }
            W_Broyden.push_back(new RealTable1D({0},{num_field_sites_all_NS}, The_Pinned_Arena()));
            V_Broyden.push_back(new RealTable1D({0},{num_field_sites_all_NS}, The_Pinned_Arena()));

        //}
        amrex::Print() << "W&V size: " << W_Broyden.size() << "\n";
        amrex::Print() << "n_curr_in, n_prev_in: " << n_curr_in(0) << " " << n_prev_in(0) << "\n";

        int m = Broyden_Step-1;
        if(m > 0)
        {
            for(int j=1; j <= m-1; ++j)
            {
                auto const& W_j = W_Broyden[j]->table();
                auto const& V_j = V_Broyden[j]->table();

                for(int a=0; a < num_field_sites_all_NS; ++a)
                {
                    amrex::Real sum = 0.;
                    for(int b=0; b < num_field_sites_all_NS; ++b)
                    {
                        sum += W_j(a)*V_j(b)*delta_F_curr(b);
                    }
                    sum_deltaFcurr(a) += sum;
                }
            }

            for(int l=0; l < num_field_sites_all_NS; ++l)
            {

                  amrex::Real delta_n = n_curr_in(l) - n_prev_in(l);
                  V_curr(l) = delta_F_curr(l)/denom;
                  W_curr(l) = -Broyden_fraction*delta_F_curr(l) + delta_n - sum_deltaFcurr(l);
            }

            W_Broyden[m]->copy(W_curr_data);
            V_Broyden[m]->copy(V_curr_data);

            auto const& W_m = W_Broyden[m]->table();
            auto const& V_m = V_Broyden[m]->table();

            for(int j=1; j <= m; ++j)
            {
                auto const& W_j = W_Broyden[j]->table();
                auto const& V_j = V_Broyden[j]->table();

                for(int a=0; a < num_field_sites_all_NS; ++a)
                {
                    amrex::Real sum = 0.;
                    for(int b=0; b < num_field_sites_all_NS; ++b)
                    {
                        sum += W_j(a)*V_j(b)*F_curr(b);
                    }
                    sum_Fcurr(a) += sum;
                }
            }
        }


        for(int l=0; l < num_field_sites_all_NS; ++l)
        {
            n_prev_in(l) = n_curr_in(l);
            n_curr_in(l) = n_prev_in(l) - Broyden_Scalar * Broyden_fraction*F_curr(l) - Broyden_Scalar * sum_Fcurr(l);
        }
        amrex::Print() << "n_new_in: " << n_curr_in(0) << "\n";

        sum_Fcurr_data.clear();
        sum_deltaFcurr_data.clear();
        W_curr_data.clear();
        V_curr_data.clear();

        Broyden_Step += 1;

    }

    MPI_Bcast(&Broyden_Norm,
               1,
               MPI_DOUBLE,
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

}
#endif
