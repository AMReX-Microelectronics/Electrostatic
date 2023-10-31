#include "Transport.H"
#include "Transport_Table_ReadWrite.H"

#include "../../Utils/SelectWarpXUtils/TextMsg.H"

using namespace amrex;

#ifndef BROYDEN_PARALLEL
void
c_TransportSolver:: Set_Broyden ()
{

    num_field_sites_all_NS = 0;
    for (int c=0; c < vp_CNT.size(); ++c)
    {
	    num_field_sites_all_NS += vp_CNT[c]->num_field_sites;
    }	
    amrex::Print() << "Number of field_sites at all nanostructures, num_field_sites_all_NS: " 
                   << num_field_sites_all_NS << "\n";

    if (ParallelDescriptor::IOProcessor())
    {
        Broyden_Step = 1;
        Broyden_Norm = 1.;
        Broyden_Scalar = 1.;
       // Broyden_Correction_Step = 0;
       // Broyden_NormSumIsIncreasing_Step = 0;
       // Broyden_Reset_Step   = 0;
        Broyden_NormSum_Curr = 1.e10;
        Broyden_NormSum_Prev = 1.e10;
        Broyden_fraction = Broyden_Original_Fraction;

        h_n_curr_in_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        h_n_curr_out_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        h_n_start_in_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        h_n_prev_in_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        h_F_curr_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        h_Norm_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        h_delta_F_curr_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        h_delta_n_curr_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());

        SetVal_RealTable1D(h_n_curr_in_data,0.);
        SetVal_RealTable1D(h_n_curr_out_data,0.);
        SetVal_RealTable1D(h_n_start_in_data,0.);
        SetVal_RealTable1D(h_n_prev_in_data,0.);
        SetVal_RealTable1D(h_F_curr_data,0.);
        SetVal_RealTable1D(h_Norm_data, 0.);
        SetVal_RealTable1D(h_delta_F_curr_data, 0.);
        SetVal_RealTable1D(h_delta_n_curr_data, 0.);

        /*Need generalization for multiple CNTs*/
        for (int c=0; c < vp_CNT.size(); ++c)
        {
            h_n_curr_in_data.copy(vp_CNT[c]->h_n_curr_in_data); 
        }
        h_n_start_in_data.copy(h_n_curr_in_data);

        switch(map_AlgorithmType[Algorithm_Type])
        {
            case s_Algorithm::Type::broyden_first:
            {
                h_Jinv_curr_data.resize({0,0},{num_field_sites_all_NS, num_field_sites_all_NS}, The_Pinned_Arena());
                SetVal_RealTable2D(h_Jinv_curr_data,0.);

                auto const& Jinv_curr    = h_Jinv_curr_data.table();


		        if(flag_initialize_inverse_jacobian) 
    	    	{
    		    int assert_size = std::pow(num_field_sites_all_NS, 2);

                    Read_Table2D(assert_size, h_Jinv_curr_data, inverse_jacobian_filename);
		        }
        		else 
		        {
                    for(int a=0; a < num_field_sites_all_NS; ++a) 
                    {
                        Jinv_curr(a,a) = Broyden_fraction;
                    }
        		}
                break;
            }
            case s_Algorithm::Type::broyden_second:
            {
                break;
            }
            case s_Algorithm::Type::simple_mixing:
            {
                break;
            }
            default:
            {
                amrex::Abort("In Set_Broyden: selfconsistency_algorithm, " + Algorithm_Type + ", is not yet defined.");
            }
        }
	
	
        amrex::Print() << "\nBroyden parameters are set to the following: \n";
        amrex::Print() << " Broyden_Step: "                          << Broyden_Step                          << "\n";
        amrex::Print() << " Broyden_Scalar: "                        << Broyden_Scalar                        << "\n";
        amrex::Print() << " Broyden_fraction: "                      << Broyden_fraction                      << "\n";
        //amrex::Print() << " Broyden_Correction_Step: "               << Broyden_Correction_Step               << "\n";
        //amrex::Print() << " Broyden_NormSumIsIncreasing_Step: "      << Broyden_NormSumIsIncreasing_Step      << "\n";
        //amrex::Print() << " Broyden_Reset_Step: "                    << Broyden_Reset_Step                    << "\n";
        amrex::Print() << " Broyden_NormSum_Curr: "                  << Broyden_NormSum_Curr                  << "\n";
        amrex::Print() << " Broyden_NormSum_Prev: "                  << Broyden_NormSum_Prev                  << "\n";
        amrex::Print() << " Broyden_Max_Norm: "                      << Broyden_Norm                          << "\n";
        amrex::Print() << " Broyden_Threshold_MaxStep: "             << Broyden_Threshold_MaxStep             << "\n";
        //amrex::Print() << " Broyden_Threshold_NormSumIncreaseStep: " << Broyden_Threshold_NormSumIncreaseStep << "\n";
        //amrex::Print() << " Broyden_Threshold_CorrectionStep: "      << Broyden_Threshold_CorrectionStep      << "\n";
        //amrex::Print() << " Broyden_Threshold_MinFraction: "         << Broyden_Threshold_MinFraction         << "\n";
        //amrex::Print() << " Broyden_Fraction_Decrease_Factor: "      << Broyden_Fraction_Decrease_Factor      << "\n";
        //amrex::Print() << " Broyden_Scalar_Decrease_Factor: "        << Broyden_Scalar_Decrease_Factor        << "\n\n";
    }

}


void
c_TransportSolver:: Reset_Broyden ()
{
    if(flag_reset_with_previous_charge_distribution == 0) 
    {
        SetVal_RealTable1D(h_n_curr_in_data, NS_initial_deposit_value);
        #if AMREX_USE_GPU
        d_n_curr_in_data.copy(h_n_curr_in_data);
        #endif
	amrex::Print() << "Input charge distribution is reset to: " << NS_initial_deposit_value << "\n";
    }

    if (ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\n\n\n\n**********************************Resetting Broyden**********************************\n";
        int size = W_Broyden.size();
        for(int j=0; j<size; ++j)
        {
            W_Broyden[j]->clear();
            V_Broyden[j]->clear();
        }
        W_Broyden.clear();
        V_Broyden.clear();

        SetVal_RealTable1D(h_n_curr_out_data, 0.);
        SetVal_RealTable1D(h_n_prev_in_data, 0.);
        SetVal_RealTable1D(h_F_curr_data, 0.);
        SetVal_RealTable1D(h_Norm_data, 0.);
        SetVal_RealTable1D(h_delta_F_curr_data, 0.);
        SetVal_RealTable1D(h_delta_n_curr_data, 0.);

        switch(map_AlgorithmType[Algorithm_Type])
        {
            case s_Algorithm::Type::broyden_first:
            {
                //SetVal_RealTable2D(Jinv_curr_data,0.);

                //auto const& Jinv_curr    = Jinv_curr_data.table();
                //for(int a=0; a < num_field_sites_all_NS; ++a) 
                //{
                //    Jinv_curr(a,a) = Broyden_fraction;
                //}
                break;
            }
            case s_Algorithm::Type::broyden_second:
            {
                break;
            }
            case s_Algorithm::Type::simple_mixing:
            {
                break;
            }
            default:
            {
                amrex::Abort("In Set_Broyden: selfconsistency_algorithm, " + Algorithm_Type + ", is not yet defined.");
            }
        }

        Broyden_Step = 1;
        Broyden_Norm = 1;
        Broyden_Scalar          = 1.;
        Broyden_NormSum_Curr    = 1.e10;
        Broyden_NormSum_Prev    = 1.e10;

        //Broyden_Correction_Step = 0;
        //Broyden_NormSumIsIncreasing_Step = 0;
        //Broyden_Reset_Step = 0;
        Broyden_fraction = Broyden_Original_Fraction;
        h_n_start_in_data.copy(h_n_curr_in_data);

        amrex::Print() << "\nBroyden parameters are reset to the following: \n";
        amrex::Print() << " Broyden_Step: "      << Broyden_Step << "\n";
        amrex::Print() << " Broyden_Scalar: "    << Broyden_Scalar << "\n";
        amrex::Print() << " Broyden_fraction: "  << Broyden_fraction << "\n";
        //amrex::Print() << " Broyden_Correction_Step: "  << Broyden_Correction_Step << "\n";
        //amrex::Print() << " Broyden_NormSumIsIncreasing_Step: "  << Broyden_NormSumIsIncreasing_Step << "\n";
        //amrex::Print() << " Broyden_Reset_Step: "       << Broyden_Reset_Step << "\n";
        amrex::Print() << " Broyden_NormSum_Curr: "     << Broyden_NormSum_Curr << "\n";
        amrex::Print() << " Broyden_NormSum_Prev: "     << Broyden_NormSum_Prev << "\n";
        amrex::Print() << " Broyden_Max_Norm: "         << Broyden_Norm << "\n\n";

    }
}


//void
//c_TransportSolver:: Restart_Broyden ()
//{
//
//    if (ParallelDescriptor::IOProcessor())
//    {
//        amrex::Print() << "\n\n\n\n**********************************Restarting Broyden**********************************\n";
//        int size = W_Broyden.size();
//        for(int j=0; j<size; ++j)
//        {
//            W_Broyden[j]->clear();
//            V_Broyden[j]->clear();
//        }
//        W_Broyden.clear();
//        V_Broyden.clear();
//        amrex::Print() << " W&V sizes before/after: " << size << ", "<< W_Broyden.size() << "\n";
//
//        SetVal_RealTable1D(h_n_curr_out_data, 0.);
//        SetVal_RealTable1D(h_n_prev_in_data, 0.);
//        SetVal_RealTable1D(h_F_curr_data, 0.);
//        SetVal_RealTable1D(h_Norm_data, 0.);
//        SetVal_RealTable1D(h_delta_F_curr_data, 0.);
//        SetVal_RealTable1D(h_delta_n_curr_data, 0.);
//
//
//        Broyden_Step = 0;
//        Broyden_Norm = 1;
//        Broyden_Scalar          = 1.;
//        Broyden_Correction_Step = 0;
//        Broyden_NormSumIsIncreasing_Step = 0;
//        Broyden_NormSum_Curr    = 1.e10;
//        Broyden_NormSum_Prev    = 1.e10;
//
//        SetVal_RealTable1D(h_n_curr_in_data, 0.);
//        h_n_curr_in_data.copy(h_n_start_in_data);
//
//        Broyden_Reset_Step += 1;
//        Broyden_fraction = Broyden_Original_Fraction/std::pow(Broyden_Fraction_Decrease_Factor,Broyden_Reset_Step);
//
//        if(Broyden_fraction < Broyden_Threshold_MinFraction)
//        {
//            amrex::Print() << "*Broyden_fraction is less than the threshold: " << Broyden_Threshold_MinFraction << "\n";
//            amrex::Print() << "*Charge density is reset to: " << NS_initial_deposit_value << "\n";
//
//            Broyden_fraction = Broyden_Original_Fraction;
//            Broyden_Reset_Step = 0;
//
//            SetVal_RealTable1D(h_n_curr_in_data, NS_initial_deposit_value);
//            SetVal_RealTable1D(h_n_start_in_data, NS_initial_deposit_value);
//
//        }
//
//        amrex::Print() << "\nBroyden parameters are reset to the following: \n";
//        amrex::Print() << " Broyden_Step: "      << Broyden_Step << "\n";
//        amrex::Print() << " Broyden_Scalar: "    << Broyden_Scalar << "\n";
//        amrex::Print() << " Broyden_fraction: "  << Broyden_fraction << "\n";
//        amrex::Print() << " Broyden_Correction_Step: "  << Broyden_Correction_Step << "\n";
//        amrex::Print() << " Broyden_NormSumIsIncreasing_Step: "  << Broyden_NormSumIsIncreasing_Step << "\n";
//        amrex::Print() << " Broyden_Reset_Step: "       << Broyden_Reset_Step << "\n";
//        amrex::Print() << " Broyden_NormSum_Curr: "     << Broyden_NormSum_Curr << "\n";
//        amrex::Print() << " Broyden_NormSum_Prev: "     << Broyden_NormSum_Prev << "\n";
//        amrex::Print() << " Broyden_Max_Norm: "         << Broyden_Norm << "\n\n";
//
//    }
//}


void
c_TransportSolver::Deallocate_Broyden_Serial ()
{
    if (ParallelDescriptor::IOProcessor())
    {
        h_n_curr_in_data.clear();
        h_n_prev_in_data.clear();
        h_F_curr_data.clear();
        h_delta_F_curr_data.clear();
        h_Norm_data.clear();
        h_delta_n_curr_data.clear();
    }
}
#endif
