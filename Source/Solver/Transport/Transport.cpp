#include "Transport.H"

#include "../../Utils/SelectWarpXUtils/WarpXConst.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"

#include "../../Utils/CodeUtils/CodeUtil.H"
#include "../../Code.H"
#include "../../Input/GeometryProperties/GeometryProperties.H"
#include "../../Input/MacroscopicProperties/MacroscopicProperties.H"
#include "../Electrostatics/MLMG.H"
#include "../Output/Output.H"
#include "../PostProcessor/PostProcessor.H"

#include <AMReX.H>
//#include <AMReX_ParmParse.H>
//#include <AMReX_Parser.H>
//#include <AMReX_RealBox.H>
//#include <AMReX_MultiFab.H>
////#include <AMReX_MultiFabUtil.H>


using namespace amrex;

c_TransportSolver::c_TransportSolver()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_TransportSolver() Constructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    ReadData();
  
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_TransportSolver() Constructor************************\n";
#endif
}


c_TransportSolver::~c_TransportSolver()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_TransportSolver() Destructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    vp_CNT.clear();
    vp_Graphene.clear();
    //vp_Silicon.clear();
    
    Deallocate_Broyden();
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_TransportSolver() Destructor************************\n";
#endif
}

void 
c_TransportSolver::ReadData() 
{

    amrex::Print() << "\n##### Transport Solver #####\n\n";

    amrex::ParmParse pp_transport("transport");
    num_NS = 0;

    amrex::Vector<std::string> temp_vec;
    bool NS_pecified = pp_transport.queryarr("NS_names", temp_vec);

    int c=0;
    for (auto it: temp_vec)
    {
        if (std::find(vec_NS_names.begin(), vec_NS_names.end(), it) == vec_NS_names.end()) 
        {
           vec_NS_names.push_back(it);
           ++c;
        }
    }
    temp_vec.clear();

    amrex::Print() << "##### transport.NS_names: ";
    for (auto it: vec_NS_names) amrex::Print() << it << "  ";
    amrex::Print() << "\n";

}


void 
c_TransportSolver::InitData() 
{
    amrex::Real negf_init_time = amrex::second();

    amrex::Print() << "\n##### TRANSPORT PROPERTIES #####\n\n";

    amrex::ParmParse pp_transport("transport");
    
    use_selfconsistent_potential = 0;
    use_negf = 0;
    pp_transport.query("use_selfconsistent_potential", use_selfconsistent_potential);
    pp_transport.query("use_negf", use_negf);
    amrex::Print() << "##### transport.use_selfconsistent_potential: " 
                   << use_selfconsistent_potential << "\n";
    amrex::Print() << "##### transport.use_negf: " 
                   << use_negf << "\n";

    NS_gather_field_str = "phi";
    NS_deposit_field_str = "charge_density";
    amrex::Real NS_initial_deposit_value = 0.;

    amrex::Real NS_Broyden_frac = 0.1;
    std::string NS_Broyden_norm_type = "relative";

    auto& rCode    = c_Code::GetInstance();

    amrex::ParmParse pp_plot("plot");
    std::string foldername_str = "output";
    pp_plot.query("folder_name", foldername_str);
    negf_foldername_str = foldername_str + "/negf";
    common_foldername_str = negf_foldername_str + "/transport_common";

    CreateDirectory(foldername_str);
    CreateDirectory(negf_foldername_str);
    CreateDirectory(common_foldername_str);
    

    if(rCode.use_electrostatic)
    {
        auto& rGprop = rCode.get_GeometryProperties();

        _geom = &rGprop.geom;
        _dm   = &rGprop.dm;
        _ba   = &rGprop.ba;

        pp_transport.get("NS_gather_field", NS_gather_field_str);
        amrex::Print() << "##### transport.NS_gather_field: " << NS_gather_field_str << "\n";
        if ( Evaluate_TypeOf_MacroStr(NS_gather_field_str) != 0 )
        {
            /*mf not defined in MacroscopicProperties*/
            amrex::Abort("NS_gather_field " + NS_gather_field_str + " not defined in Mprop.");
        } 

        pp_transport.get("NS_deposit_field", NS_deposit_field_str);
        amrex::Print() << "##### transport.NS_deposit_field: " << NS_deposit_field_str << "\n";
        if ( Evaluate_TypeOf_MacroStr(NS_deposit_field_str) != 0 )
        {
            /*mf not defined in MacroscopicProperties*/
            amrex::Abort("NS_deposit_field " + NS_deposit_field_str + " not defined in Mprop.");
        } 

        auto dxi = _geom->InvCellSizeArray();
        amrex::Real inv_vol = AMREX_D_TERM(dxi[0], *dxi[1], *dxi[2]);
        NS_initial_deposit_value = PhysConst::q_e*inv_vol;

        queryWithParser(pp_transport,"NS_initial_deposit_value", NS_initial_deposit_value);
        amrex::Print() << "##### transport.NS_initial_deposit_value: " << NS_initial_deposit_value << "\n";

        queryWithParser(pp_transport,"Broyden_fraction", NS_Broyden_frac);
        amrex::Print() << "##### Broyden_fraction: " << NS_Broyden_frac << "\n";

        queryWithParser(pp_transport,"Broyden_max_norm", Broyden_max_norm);
        amrex::Print() << "##### Broyden_max_norm: " << Broyden_max_norm << "\n";

        pp_transport.query("Broyden_norm_type", NS_Broyden_norm_type);
        amrex::Print() << "##### Broyden_norm_type: " << NS_Broyden_norm_type << "\n";

        pp_transport.query("Broyden_threshold_maxstep", Broyden_Threshold_MaxStep);
        amrex::Print() << "##### Broyden_Threshold_MaxStep: " << Broyden_Threshold_MaxStep << "\n";

        pp_transport.query("selfconsistency_algorithm", Algorithm_Type);
        amrex::Print() << "##### selfconsistency_algorithm: " << Algorithm_Type << "\n";

        pp_transport.query("reset_with_previous_charge_distribution", flag_reset_with_previous_charge_distribution);
        amrex::Print() << "##### reset_with_previous_charge_distribution: " << flag_reset_with_previous_charge_distribution  << "\n";

	flag_initialize_inverse_jacobian = 0;
        pp_transport.query("initialize_inverse_jacobian", flag_initialize_inverse_jacobian);
        amrex::Print() << "##### flag_initialize_inverse_jacobian: " << flag_initialize_inverse_jacobian  << "\n";

	if(flag_initialize_inverse_jacobian) 
	{
            amrex::ParmParse pp;
            int flag_restart = 0;
            pp.query("restart", flag_restart);

	    if(flag_restart) 
	    {		    
                int restart_step = 0;       
                getWithParser(pp,"restart_step", restart_step);

                std::string restart_folder_str  = amrex::Concatenate(common_foldername_str + "/step", restart_step-1, negf_plt_name_digits);
                /*eg. output/negf/transport_common/step0001 for step 1*/

                inverse_jacobian_filename  = restart_folder_str + "/Jinv.dat";
                pp_transport.query("inverse_jacobian_filename", inverse_jacobian_filename);
	    }	
	    else 
	    {
                pp_transport.get("inverse_jacobian_filename", inverse_jacobian_filename);
	    }
            amrex::Print() << "##### inverse_jacobian_filename: " << inverse_jacobian_filename << "\n";
        }

    }

    std::string type;

    for (auto name: vec_NS_names)
    {

        amrex::ParmParse pp_ns(name);
        pp_ns.get("type", type);
        amrex::Print() << "\n##### name: " << name << "\n";
        amrex::Print() << "##### type: " << type << "\n";

        switch (map_NSType_enum[type])
        {
            case s_NS::Type::CNT:
            {
                using T = c_CNT;
 
                vp_CNT.push_back(std::make_unique<c_Nanostructure<T>>(*_geom, *_dm, *_ba,
                                                                       name, 
                                                                       NS_gather_field_str, 
                                                                       NS_deposit_field_str, 
                                                                       NS_initial_deposit_value,
                                                                       NS_Broyden_frac,
                                                                       NS_Broyden_norm_type,
                                                                       use_negf,
								       negf_foldername_str
                                                                      ));
                break;
            }
            case s_NS::Type::Graphene:
            {
                using T = c_Graphene;
 
                vp_Graphene.push_back(std::make_unique<c_Nanostructure<T>>(*_geom, *_dm, *_ba,
                                                                            name, 
                                                                            NS_gather_field_str, 
                                                                            NS_deposit_field_str, 
                                                                            NS_initial_deposit_value,
                                                                            NS_Broyden_frac,
                                                                            NS_Broyden_norm_type,
                                                                            use_negf,
									    negf_foldername_str
                                                                           ));
                amrex::Abort("NS_type " + type + " is not yet defined.");
                break; 
            }
            case s_NS::Type::Silicon:
            {
                amrex::Abort("NS_type " + type + " is not yet defined.");
            }
            default:
            {
                amrex::Abort("NS_type " + type + " is not supported.");
            }
        }
    }


    num_field_sites_all_NS = 0;
    for (int c=0; c < vp_CNT.size(); ++c)
    {
	num_field_sites_all_NS += vp_CNT[c]->num_field_sites;
    }	
    /*Here add field sites from all other nanostructures*/
    amrex::Print() << "Number of field_sites at all nanostructures, num_field_sites_all_NS: " << num_field_sites_all_NS << "\n";

    Broyden_Original_Fraction = NS_Broyden_frac;
    Broyden_Norm_Type         = NS_Broyden_norm_type;
    amrex::Print() << "#####* Broyden_Original_Fraction: "  << Broyden_Original_Fraction  << "\n";
    amrex::Print() << "#####* Broyden_Norm_Type: "          << Broyden_Norm_Type          << "\n";
    

    Set_Broyden();

    if (ParallelDescriptor::IOProcessor())
    {
        for (int c=0; c < vp_CNT.size(); ++c)
        {
            n_curr_in_data.copy(vp_CNT[c]->h_n_curr_in_data); /*Need generalization for multiple CNTs*/
	}
        n_start_in_data.copy(n_curr_in_data);
    }
}


void
c_TransportSolver::Set_CommonStepFolder(const int step)
{

    common_step_folder_str  = amrex::Concatenate(common_foldername_str + "/step", step, negf_plt_name_digits);
    amrex::Print() << "common_step_folder_str: " << common_step_folder_str << "\n";
   
    CreateDirectory(common_step_folder_str);
    /*eg. output/negf/common/step0001 for step 1*/

}

void 
c_TransportSolver::Solve(const int step, const amrex::Real time) 
{

   auto& rCode    = c_Code::GetInstance();
   auto& rGprop = rCode.get_GeometryProperties();
   auto& rMprop = rCode.get_MacroscopicProperties();
   auto& rMLMG    = rCode.get_MLMGSolver();
   auto& rOutput  = rCode.get_Output();
   auto& rPostPro = rCode.get_PostProcessor();

   amrex::Real total_mlmg_solve_time = 0.;

   int max_iter = 1;
   int MAX_ITER_THRESHOLD = 800;


   for (int c=0; c < vp_CNT.size(); ++c)
   {
       vp_CNT[c]->Set_StepFilenameString(step);
       Set_CommonStepFolder(step);
   }

   amrex::Real Vds = 0.;
   amrex::Real Vgs = 0.;
   if(rCode.use_electrostatic) 
   {	
	   
       bool update_surface_soln_flag = true;	   

       do 
       {
	   if(max_iter > MAX_ITER_THRESHOLD) amrex::Abort("Iteration step is GREATER than the THRESHOLD" + MAX_ITER_THRESHOLD);

           amrex::Print() << "\n\nSelf-consistent iteration: " << max_iter << "\n";

           rMprop.ReInitializeMacroparam(NS_gather_field_str);
           rMLMG.UpdateBoundaryConditions(update_surface_soln_flag);

           auto mlmg_solve_time = rMLMG.Solve_PoissonEqn();
           total_mlmg_solve_time += mlmg_solve_time;
           amrex::Print() << "\nmlmg_solve_time: " << mlmg_solve_time << "\n";



           rPostPro.Compute();

           //rOutput.WriteOutput(max_iter+1000, time);

           for (int c=0; c < vp_CNT.size(); ++c)
           {

	       vp_CNT[c]->Set_IterationFilenameString(max_iter);

	       vp_CNT[c]->Write_InputInducedCharge(vp_CNT[c]->iter_filename_str, n_curr_in_data);  

	       vp_CNT[c]->Gather_MeshAttributeAtAtoms();  
               

	       if(update_surface_soln_flag && vp_CNT[c]->flag_contact_mu_specified == 0) 
	       {
                   amrex::Real V_contact[NUM_CONTACTS] = {0., 0.};

		   for(int k=0; k<NUM_CONTACTS; ++k) 
		   {    
                       V_contact[k] = rGprop.pEB->Read_SurfSoln(vp_CNT[c]->Contact_Parser_String[k]);
                       
		       amrex::Print() << "Updated terminal voltage: " << k << "  " << V_contact[k] << "\n";

		       vp_CNT[c]->Contact_Electrochemical_Potential[k] = vp_CNT[c]->E_f - V_contact[k];
		   }

		   Vds = V_contact[1] - V_contact[0];

                   Vgs = rGprop.pEB->Read_SurfSoln(vp_CNT[c]->Gate_String) - V_contact[0];

		   amrex::Print() << "Vds, Vgs: " << Vds << "  " << Vgs << "\n";
	       }

               vp_CNT[c]->Solve_NEGF();

	       vp_CNT[c]->Gather_NEGFComputed_Charge(n_curr_out_data); /*Need a strategy to gather data for multiple CNTs*/

           }

            if (Broyden_Step > Broyden_Threshold_MaxStep)
            {
                amrex::Abort("Broyden_Step has exceeded the Broyden_Threshold_MaxStep!");
            }

            switch(map_AlgorithmType[Algorithm_Type])
            {
                case s_Algorithm::Type::broyden_first:
                {
	            Execute_Broyden_First_Algorithm(); 
                    break;
                }
                case s_Algorithm::Type::broyden_second:
                {
	            Execute_Broyden_Modified_Second_Algorithm(); 
                    break;
                }
                case s_Algorithm::Type::simple_mixing:
                {
	            Execute_Simple_Mixing_Algorithm(); 
                    break;
                }
                default:
                {
                    amrex::Abort("selfconsistency_algorithm, " + Algorithm_Type + ", is not yet defined.");
                }
            }

           rMprop.ReInitializeMacroparam(NS_deposit_field_str);

           for (int c=0; c < vp_CNT.size(); ++c)
           {
               if (ParallelDescriptor::IOProcessor())
               {
                   vp_CNT[c]->h_n_curr_in_data.copy(n_curr_in_data); /*Need generalization for multiple CNTs*/
	       }

	       vp_CNT[c]->Broadcast_BroydenPredicted_Charge();

	       if(vp_CNT[c]->write_at_iter) 
	       {
                   vp_CNT[c]->Write_Data(vp_CNT[c]->iter_filename_str, n_curr_out_data, Norm_data);
	       }

               vp_CNT[c]->Deposit_AtomAttributeToMesh();
	   }
           update_surface_soln_flag = false;
           max_iter += 1;

       } while(Broyden_Norm > Broyden_max_norm);    


       for (int c=0; c < vp_CNT.size(); ++c)
       {

           vp_CNT[c]->Write_Data(vp_CNT[c]->step_filename_str, n_curr_out_data, Norm_data);

           if(map_AlgorithmType[Algorithm_Type] == s_Algorithm::Type::broyden_first) 
	   {
               Write_Table2D(Jinv_curr_data, common_step_folder_str + "/Jinv.dat", "Jinv");
	   }

           vp_CNT[c]->Compute_Current();

           vp_CNT[c]->Write_Current(step, Vds, Vgs, Broyden_Step, max_iter, Broyden_fraction, Broyden_Scalar);

           if(rCode.use_electrostatic)
           {
               Reset_Broyden();
           }

           //rMprop.ReInitializeMacroparam(NS_deposit_field_str);
       }

       amrex::Print() << "\nAverage mlmg time for self-consistency (s): " << total_mlmg_solve_time / max_iter << "\n";

   }
   else 
   {
       for (int c=0; c < vp_CNT.size(); ++c)
       {
           vp_CNT[c]->Solve_NEGF();

           vp_CNT[c]->Write_Data(vp_CNT[c]->step_filename_str, n_curr_out_data, Norm_data);

           vp_CNT[c]->Compute_Current();

           vp_CNT[c]->Write_Current(step, Vds, Vgs, Broyden_Step, max_iter, Broyden_fraction, Broyden_Scalar);
       }
   }

}


void
c_TransportSolver:: Set_Broyden ()
{

    if (ParallelDescriptor::IOProcessor())
    {

        Broyden_Step = 1;
        Broyden_Norm = 1.;
        Broyden_Scalar = 1.;
        Broyden_Correction_Step = 0;
        Broyden_NormSumIsIncreasing_Step = 0;
        Broyden_Reset_Step   = 0;
        Broyden_NormSum_Curr = 1.e10;
        Broyden_NormSum_Prev = 1.e10;
        Broyden_fraction = Broyden_Original_Fraction;

        n_curr_in_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        n_curr_out_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        n_start_in_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        n_prev_in_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        F_curr_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        Norm_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        delta_F_curr_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        delta_n_curr_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());

        SetVal_RealTable1D(n_curr_in_data,0.);
        SetVal_RealTable1D(n_curr_out_data,0.);
        SetVal_RealTable1D(n_start_in_data,0.);
        SetVal_RealTable1D(n_prev_in_data,0.);
        SetVal_RealTable1D(F_curr_data,0.);
        SetVal_RealTable1D(Norm_data, 0.);
        SetVal_RealTable1D(delta_F_curr_data, 0.);
        SetVal_RealTable1D(delta_n_curr_data, 0.);


        switch(map_AlgorithmType[Algorithm_Type])
        {
            case s_Algorithm::Type::broyden_first:
            {
                Jinv_curr_data.resize({0,0},{num_field_sites_all_NS, num_field_sites_all_NS}, The_Pinned_Arena());
                SetVal_RealTable2D(Jinv_curr_data,0.);

                auto const& Jinv_curr    = Jinv_curr_data.table();


		if(flag_initialize_inverse_jacobian) 
		{
		    int assert_size = std::pow(num_field_sites_all_NS, 2);

                    Read_Table2D(assert_size, Jinv_curr_data, inverse_jacobian_filename);
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
        amrex::Print() << " Broyden_Correction_Step: "               << Broyden_Correction_Step               << "\n";
        amrex::Print() << " Broyden_NormSumIsIncreasing_Step: "      << Broyden_NormSumIsIncreasing_Step      << "\n";
        amrex::Print() << " Broyden_Reset_Step: "                    << Broyden_Reset_Step                    << "\n";
        amrex::Print() << " Broyden_NormSum_Curr: "                  << Broyden_NormSum_Curr                  << "\n";
        amrex::Print() << " Broyden_NormSum_Prev: "                  << Broyden_NormSum_Prev                  << "\n";
        amrex::Print() << " Broyden_Max_Norm: "                      << Broyden_Norm                          << "\n";
        amrex::Print() << " Broyden_Threshold_MaxStep: "             << Broyden_Threshold_MaxStep             << "\n";
        amrex::Print() << " Broyden_Threshold_NormSumIncreaseStep: " << Broyden_Threshold_NormSumIncreaseStep << "\n";
        amrex::Print() << " Broyden_Threshold_CorrectionStep: "      << Broyden_Threshold_CorrectionStep      << "\n";
        amrex::Print() << " Broyden_Threshold_MinFraction: "         << Broyden_Threshold_MinFraction         << "\n";
        amrex::Print() << " Broyden_Fraction_Decrease_Factor: "      << Broyden_Fraction_Decrease_Factor      << "\n";
        amrex::Print() << " Broyden_Scalar_Decrease_Factor: "        << Broyden_Scalar_Decrease_Factor        << "\n\n";
    }

}


void
c_TransportSolver:: Reset_Broyden ()
{
    if(flag_reset_with_previous_charge_distribution == 0) 
    {
        SetVal_RealTable1D(n_curr_in_data, NS_initial_deposit_value);
        #if AMREX_USE_GPU
        d_n_curr_in_data.copy(n_curr_in_data);
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

        SetVal_RealTable1D(n_curr_out_data, 0.);
        SetVal_RealTable1D(n_prev_in_data, 0.);
        SetVal_RealTable1D(F_curr_data, 0.);
        SetVal_RealTable1D(Norm_data, 0.);
        SetVal_RealTable1D(delta_F_curr_data, 0.);
        SetVal_RealTable1D(delta_n_curr_data, 0.);

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
        Broyden_Correction_Step = 0;
        Broyden_NormSumIsIncreasing_Step = 0;
        Broyden_NormSum_Curr    = 1.e10;
        Broyden_NormSum_Prev    = 1.e10;

        Broyden_Reset_Step = 0;
        Broyden_fraction = Broyden_Original_Fraction;
        n_start_in_data.copy(n_curr_in_data);

        amrex::Print() << "\nBroyden parameters are reset to the following: \n";
        amrex::Print() << " Broyden_Step: "      << Broyden_Step << "\n";
        amrex::Print() << " Broyden_Scalar: "    << Broyden_Scalar << "\n";
        amrex::Print() << " Broyden_fraction: "  << Broyden_fraction << "\n";
        amrex::Print() << " Broyden_Correction_Step: "  << Broyden_Correction_Step << "\n";
        amrex::Print() << " Broyden_NormSumIsIncreasing_Step: "  << Broyden_NormSumIsIncreasing_Step << "\n";
        amrex::Print() << " Broyden_Reset_Step: "       << Broyden_Reset_Step << "\n";
        amrex::Print() << " Broyden_NormSum_Curr: "     << Broyden_NormSum_Curr << "\n";
        amrex::Print() << " Broyden_NormSum_Prev: "     << Broyden_NormSum_Prev << "\n";
        amrex::Print() << " Broyden_Max_Norm: "         << Broyden_Norm << "\n\n";

    }
}


void
c_TransportSolver:: Restart_Broyden ()
{

    if (ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\n\n\n\n**********************************Restarting Broyden**********************************\n";
        int size = W_Broyden.size();
        for(int j=0; j<size; ++j)
        {
            W_Broyden[j]->clear();
            V_Broyden[j]->clear();
        }
        W_Broyden.clear();
        V_Broyden.clear();
        amrex::Print() << " W&V sizes before/after: " << size << ", "<< W_Broyden.size() << "\n";

        SetVal_RealTable1D(n_curr_out_data, 0.);
        SetVal_RealTable1D(n_prev_in_data, 0.);
        SetVal_RealTable1D(F_curr_data, 0.);
        SetVal_RealTable1D(Norm_data, 0.);
        SetVal_RealTable1D(delta_F_curr_data, 0.);
        SetVal_RealTable1D(delta_n_curr_data, 0.);


        Broyden_Step = 0;
        Broyden_Norm = 1;
        Broyden_Scalar          = 1.;
        Broyden_Correction_Step = 0;
        Broyden_NormSumIsIncreasing_Step = 0;
        Broyden_NormSum_Curr    = 1.e10;
        Broyden_NormSum_Prev    = 1.e10;

        SetVal_RealTable1D(n_curr_in_data, 0.);
        n_curr_in_data.copy(n_start_in_data);

        Broyden_Reset_Step += 1;
        Broyden_fraction = Broyden_Original_Fraction/std::pow(Broyden_Fraction_Decrease_Factor,Broyden_Reset_Step);

        if(Broyden_fraction < Broyden_Threshold_MinFraction)
        {
            amrex::Print() << "*Broyden_fraction is less than the threshold: " << Broyden_Threshold_MinFraction << "\n";
            amrex::Print() << "*Charge density is reset to: " << NS_initial_deposit_value << "\n";

            Broyden_fraction = Broyden_Original_Fraction;
            Broyden_Reset_Step = 0;

            SetVal_RealTable1D(n_curr_in_data, NS_initial_deposit_value);
            SetVal_RealTable1D(n_start_in_data, NS_initial_deposit_value);

        }

        amrex::Print() << "\nBroyden parameters are reset to the following: \n";
        amrex::Print() << " Broyden_Step: "      << Broyden_Step << "\n";
        amrex::Print() << " Broyden_Scalar: "    << Broyden_Scalar << "\n";
        amrex::Print() << " Broyden_fraction: "  << Broyden_fraction << "\n";
        amrex::Print() << " Broyden_Correction_Step: "  << Broyden_Correction_Step << "\n";
        amrex::Print() << " Broyden_NormSumIsIncreasing_Step: "  << Broyden_NormSumIsIncreasing_Step << "\n";
        amrex::Print() << " Broyden_Reset_Step: "       << Broyden_Reset_Step << "\n";
        amrex::Print() << " Broyden_NormSum_Curr: "     << Broyden_NormSum_Curr << "\n";
        amrex::Print() << " Broyden_NormSum_Prev: "     << Broyden_NormSum_Prev << "\n";
        amrex::Print() << " Broyden_Max_Norm: "         << Broyden_Norm << "\n\n";

    }
}


void 
c_TransportSolver:: Execute_Broyden_First_Algorithm ()
{

    if (ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\nBroydenStep: " << Broyden_Step  << ",  fraction: " << Broyden_fraction << ",  scalar: " << Broyden_Scalar<< "\n";
	
        auto const& n_curr_in  = n_curr_in_data.table();
        auto const& n_curr_out = n_curr_out_data.table();
        auto const& n_prev_in  = n_prev_in_data.table();
        auto const& F_curr     = F_curr_data.table();
        auto const& delta_F_curr  = delta_F_curr_data.table();
        auto const& delta_n_curr   = delta_n_curr_data.table();
        auto const& Jinv_curr    = Jinv_curr_data.table();

        RealTable1D sum_Fcurr_data({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        RealTable1D sum_deltaFcurr_data({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        RealTable1D delta_n_Jinv_data({0},{num_field_sites_all_NS}, The_Pinned_Arena());

        auto const& sum_Fcurr      = sum_Fcurr_data.table();
        auto const& sum_deltaFcurr = sum_deltaFcurr_data.table();
        auto const& delta_n_Jinv   = delta_n_Jinv_data.table();

        SetVal_RealTable1D(Norm_data, 0.);
        auto const& Norm       = Norm_data.table();

        Broyden_NormSum_Curr = 0.;

        for(int l=0; l < num_field_sites_all_NS; ++l)
        {
            sum_Fcurr(l) = 0;
            sum_deltaFcurr(l) = 0;
            delta_n_Jinv(l) = 0.;
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


        Broyden_NormSum_Prev = Broyden_NormSum_Curr;

        for(int l=0; l < num_field_sites_all_NS; ++l)
        {
            amrex::Real Fcurr = n_curr_in(l) - n_curr_out(l);
            delta_F_curr(l) = Fcurr - F_curr(l);
            F_curr(l) = Fcurr;
            delta_n_curr(l) = n_curr_in(l) - n_prev_in(l);
        }

        amrex::Print() << "n_curr_in, n_prev_in: " << n_curr_in(0) << " " << n_prev_in(0) << "\n";

        int m = Broyden_Step-1;
        if(m > 0 ) 
	{
            for(int a=0; a < num_field_sites_all_NS; ++a)
            {
                amrex::Real sum = 0.;
                for(int b=0; b < num_field_sites_all_NS; ++b)
                {
                    sum += Jinv_curr(a,b)*delta_F_curr(b);
                }
                sum_deltaFcurr(a) = sum;
            }

	    amrex::Real denom = 0.;
            for(int l=0; l < num_field_sites_all_NS; ++l)
            {
                denom += delta_n_curr(l) * sum_deltaFcurr(l);
            }

            for(int b=0; b < num_field_sites_all_NS; ++b)
            {
                amrex::Real sum = 0.;
                for(int a=0; a < num_field_sites_all_NS; ++a)
                {
                    sum += delta_n_curr(a)*Jinv_curr(a,b);
                }
                delta_n_Jinv(b) = sum;
            }

            for(int a=0; a < num_field_sites_all_NS; ++a)
            {
                for(int b=0; b < num_field_sites_all_NS; ++b)
                {
                    Jinv_curr(a,b) += ( delta_n_curr(a) - sum_deltaFcurr(a) )  * delta_n_Jinv(b) / denom;
                }
            }
        }
        for(int a=0; a < num_field_sites_all_NS; ++a)
        {
            amrex::Real sum = 0.;
            for(int b=0; b < num_field_sites_all_NS; ++b)
            {
                sum += Jinv_curr(a,b)*F_curr(b);
            }
            sum_Fcurr(a) = sum;
        }

        for(int l=0; l < num_field_sites_all_NS; ++l)
        {
            n_prev_in(l) = n_curr_in(l);
            n_curr_in(l) = n_prev_in(l) - sum_Fcurr(l);
        }
        amrex::Print() << "n_new_in: " << n_curr_in(0) << "\n";

        sum_Fcurr_data.clear();
        sum_deltaFcurr_data.clear();
        delta_n_Jinv_data.clear();

        Broyden_Step += 1;
    }

    MPI_Bcast(&Broyden_Norm,
               1,
               MPI_DOUBLE,
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

}


void 
c_TransportSolver:: Execute_Simple_Mixing_Algorithm ()
{
    /*update h_RhoInduced_glo*/

    if (ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\nBroydenStep: " << Broyden_Step  << ",  fraction: " << Broyden_fraction << ",  scalar: " << Broyden_Scalar<< "\n";

        auto const& n_curr_in  = n_curr_in_data.table();
        auto const& n_curr_out = n_curr_out_data.table();
        auto const& n_prev_in  = n_prev_in_data.table();
        auto const& F_curr     = F_curr_data.table();

        SetVal_RealTable1D(Norm_data, 0.);
        auto const& Norm       = Norm_data.table();

        amrex::Real denom = 0.;
        Broyden_NormSum_Curr = 0.;

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

        if ((Broyden_NormSum_Curr - Broyden_NormSum_Prev) > 1.e-6)
        {
            Broyden_NormSumIsIncreasing_Step +=1;
            amrex::Print() << "\nBroyden_NormSumIsIncreasing_Step: " << Broyden_NormSumIsIncreasing_Step << "\n";
        }
        else
        {
            Broyden_NormSumIsIncreasing_Step = 0;
        }

        Broyden_NormSum_Prev = Broyden_NormSum_Curr;

        for(int l=0; l < num_field_sites_all_NS; ++l) 
        {
            F_curr(l) = n_curr_in(l) - n_curr_out(l);
        }

        for(int l=0; l < num_field_sites_all_NS; ++l) 
        {
            n_prev_in(l) = n_curr_in(l); 
            n_curr_in(l) = n_prev_in(l) - Broyden_fraction*F_curr(l);
        }

        Broyden_Step += 1;
    }

    MPI_Bcast(&Broyden_Norm,
               1,
               MPI_DOUBLE,
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

}



void 
c_TransportSolver::Execute_Broyden_Modified_Second_Algorithm()
{

    if (ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\nBroydenStep: " << Broyden_Step  << ",  fraction: " << Broyden_fraction << ",  scalar: " << Broyden_Scalar<< "\n";

        auto const& n_curr_in  = n_curr_in_data.table();
        auto const& n_curr_out = n_curr_out_data.table();
        auto const& n_prev_in  = n_prev_in_data.table();
        auto const& F_curr     = F_curr_data.table();
        auto const& delta_F_curr   = delta_F_curr_data.table();
        auto const& delta_n_curr   = delta_n_curr_data.table();

        SetVal_RealTable1D(Norm_data, 0.);
        auto const& Norm       = Norm_data.table();

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


        if ((Broyden_NormSum_Curr - Broyden_NormSum_Prev) > 1.e-6)
        {
            Broyden_NormSumIsIncreasing_Step +=1;
            amrex::Print() << "\nBroyden_NormSumIsIncreasing_Step: " << Broyden_NormSumIsIncreasing_Step << "\n";
        }
        else
        {
            Broyden_NormSumIsIncreasing_Step = 0;
        }

        if (Broyden_Step > Broyden_Threshold_MaxStep)
        {
              Restart_Broyden();
        }
        else if (Broyden_NormSumIsIncreasing_Step > Broyden_Threshold_NormSumIncreaseStep)
        {
            for(int l=0; l < num_field_sites_all_NS; ++l)
            {
                n_curr_in(l) = n_prev_in(l);
                n_prev_in(l) = n_curr_in(l) - delta_n_curr(l);
                denom += pow(delta_F_curr(l),2.);
            }

            Broyden_Correction_Step += 1;

            if(Broyden_Correction_Step > Broyden_Threshold_CorrectionStep)
            {
              Restart_Broyden();
            }
            else
            {
                Broyden_Scalar = Broyden_Scalar/std::pow(Broyden_Scalar_Decrease_Factor,Broyden_Correction_Step);
                Broyden_Step -= 1;
                Broyden_NormSumIsIncreasing_Step = 0;
                amrex::Print() << "\n******Reducing Broyden scalar to: " << Broyden_Scalar << "\n";
            }
        }
        else
        {
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

        }
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


void
c_TransportSolver::Deallocate_Broyden ()
{
    if (ParallelDescriptor::IOProcessor())
    {
        n_curr_in_data.clear();
        n_prev_in_data.clear();
        F_curr_data.clear();
        Norm_data.clear();

        delta_F_curr_data.clear();
        delta_n_curr_data.clear();
    }
}

void
c_TransportSolver:: SetVal_RealTable1D (RealTable1D& Tab1D_data, amrex::Real val)
{
    auto tlo = Tab1D_data.lo();
    auto thi = Tab1D_data.hi();

    auto const& Tab1D = Tab1D_data.table();

    for (int i = tlo[0]; i < thi[0]; ++i)
    {
        Tab1D(i) = val;
    }
}

void 
c_TransportSolver:: SetVal_RealTable2D (RealTable2D& Tab2D_data, amrex::Real val)
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


template<typename TableType>
void
c_TransportSolver::Read_Table1D(int assert_size,
                               TableType& Tab1D_data,
                               std::string filename)
{

    amrex::Print() << "Reading Table1D. filename: " << filename << "\n";

    std::ifstream infile;
    infile.open(filename.c_str());

    if(infile.fail())
    {
        amrex::Abort("Failed to read file " + filename);
    }
    else
    {
        auto const& Tab1D = Tab1D_data.table();
        auto thi = Tab1D_data.hi();
        auto tlo = Tab1D_data.lo();

        int filesize=0;
        std::string line;
        while(infile.peek()!=EOF)
        {
            std::getline(infile, line);
            filesize++;
        }
        amrex::Print() << "filesize: " << filesize << "\n";

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(filesize-1 == assert_size,
        "Assert size, " + std::to_string(assert_size)
         + ", is not equal to the filesize-1, " + std::to_string(filesize-1) + " !");

        infile.seekg(0, std::ios_base::beg);

        std::getline(infile, line);
        amrex::Print() << "file header: " << line << "\n";

        amrex::Real position, value;
        for(int l=tlo[0]; l < thi[0]; ++l)
        {
            infile >> position >> value;
            Tab1D(l) = value;
            amrex::Print() << "position/value: " << position << "    " << Tab1D(l) << "\n";
        }
        infile.close();
    }
}



template<typename TableType>
void
c_TransportSolver::Read_Table2D(int assert_size,
                               TableType& Tab2D_data,
                               std::string filename)
{

    amrex::Print() << "Reading Table2D. filename: " << filename << "\n";

    std::ifstream infile;
    infile.open(filename.c_str());

    if(infile.fail())
    {
        amrex::Abort("Failed to read file " + filename);
    }
    else
    {
        auto const& Tab2D = Tab2D_data.table();
        auto thi = Tab2D_data.hi();
        auto tlo = Tab2D_data.lo();

        int filesize=0;
        std::string line;
        while(infile.peek()!=EOF)
        {
            std::getline(infile, line);
            filesize++;
        }
        amrex::Print() << "filesize: " << filesize << "\n";

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(filesize-1 == assert_size,
        "Assert size, " + std::to_string(assert_size)
         + ", is not equal to the filesize-1, " + std::to_string(filesize-1) + " !");

        infile.seekg(0, std::ios_base::beg);

        std::getline(infile, line);
        amrex::Print() << "file header: " << line << "\n";

	amrex::Real value;
        for (int j = tlo[1]; j < thi[1]; ++j) //slow moving index.
        {
            for (int i = tlo[0]; i < thi[0]; ++i)
            {
                infile >> value;
		Tab2D(i,j) = value;
            }
        }
        infile.close();
    }
}


template<typename VectorType, typename TableType>
void
c_TransportSolver::Write_Table1D(const amrex::Vector<VectorType>& Vec,
                                const TableType& Arr_data,
                                std::string filename,
                                std::string header)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        //amrex::Print() << "\nRoot Writing " << filename << "\n";
        std::ofstream outfile;
        outfile.open(filename.c_str());

        auto const& Arr = Arr_data.const_table();
        auto thi = Arr_data.hi();

        outfile << header  << "\n";
        if(Vec.size() == thi[0]) {
            for (int e=0; e< thi[0]; ++e)
            {
                outfile << std::setprecision(15)
                        << std::setw(20) << Vec[e]
                        << std::setw(20) << Arr(e) << "\n";
            }
        }
        else {
            outfile << "Mismatch in the size of Vec and Table1D_data!"  << "\n";
        }
        outfile.close();
    }
}


template<typename TableType>
void
c_TransportSolver::Write_Table2D(const TableData<TableType, 2>& Tab_data,
                                 std::string filename,
                                 std::string header)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        //amrex::Print() << "\nRoot Writing " << filename << "\n";
        std::ofstream outfile;
        outfile.open(filename.c_str());

        auto const& Tab2D = Tab_data.const_table();
        auto thi = Tab_data.hi();
        auto tlo = Tab_data.lo();

        outfile << header  << "\n";
        for (int j = tlo[1]; j < thi[1]; ++j) //slow moving index.
        {
            for (int i = tlo[0]; i < thi[0]; ++i)
            {
                outfile  << std::setprecision(15) << Tab2D(i,j) << "\n";
            }
        }
        outfile.close();
    }
}
