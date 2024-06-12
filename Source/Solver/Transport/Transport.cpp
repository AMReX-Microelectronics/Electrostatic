#include "Transport.H"
#include "Transport_Table_ReadWrite.H"

#include "../../Utils/SelectWarpXUtils/WarpXConst.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"

#include "../../Utils/CodeUtils/CodeUtil.H"
#include "../../Code.H"
#include "../../Input/GeometryProperties/GeometryProperties.H"
#include "../../Input/MacroscopicProperties/MacroscopicProperties.H"
#include "../../Input/BoundaryConditions/BoundaryConditions.H"
#include "../Electrostatics/MLMG.H"
#include "../Output/Output.H"
#include "../PostProcessor/PostProcessor.H"

#include <AMReX.H>
#include <AMReX_GpuContainers.H>


using namespace amrex;

enum class s_NS_Type : int { CNT, Graphene, Silicon };
enum class s_Algorithm_Type : int { broyden_first, broyden_second, simple_mixing };
enum class Gate_Terminal_Type : int {EB, Boundary};

const std::map<std::string, s_NS_Type>
c_TransportSolver::map_NSType_enum =
{
   {"carbon nanotube" , s_NS_Type::CNT},
   {"Carbon Nanotube" , s_NS_Type::CNT},
   {"CNT"             , s_NS_Type::CNT},
   {"graphene"        , s_NS_Type::Graphene},
   {"Graphene"        , s_NS_Type::Graphene},
   {"silicon"         , s_NS_Type::Silicon},
   {"Silicon"         , s_NS_Type::Silicon}
};

const std::map<std::string, s_Algorithm_Type> 
c_TransportSolver::map_AlgorithmType =
{
   {"broyden_first"   , s_Algorithm_Type::broyden_first},
   {"Broyden_first"   , s_Algorithm_Type::broyden_first},
   {"broyden_second"  , s_Algorithm_Type::broyden_second},
   {"Broyden_second"  , s_Algorithm_Type::broyden_second},
   {"simple_mixing"   , s_Algorithm_Type::simple_mixing},
};

const std::map<std::string, Gate_Terminal_Type>
c_TransportSolver::map_S_GateTerminalType =
{
   {"eb"       , Gate_Terminal_Type::EB},
   {"EB"       , Gate_Terminal_Type::EB},
   {"boundary" , Gate_Terminal_Type::Boundary},
   {"Boundary" , Gate_Terminal_Type::Boundary}
};


c_TransportSolver::c_TransportSolver()
{
    ReadData();
}


void 
c_TransportSolver::Cleanup() 
{
    #ifdef BROYDEN_PARALLEL
    Free_MPIDerivedDataTypes();
    #endif
}


void 
c_TransportSolver::ReadData() 
{
    amrex::Print() << "\n##### Transport Solver #####\n\n";

    amrex::ParmParse pp_transport("transport");

    Read_NSNames(pp_transport);

    Read_NSTypes(pp_transport);

    Read_GatherAndDepositFields(pp_transport);

    Read_SelfConsistencyInput(pp_transport);

    Read_DOSInput(pp_transport);

    Read_GateTerminalType(pp_transport);

}

void 
c_TransportSolver::Read_NSNames(amrex::ParmParse& pp) 
{
    amrex::Vector<std::string> temp_vec;
    bool NS_names_specified = pp.queryarr("NS_names", temp_vec);

    if(NS_names_specified) 
    {
        int c=0;
        for (auto it: temp_vec)
        {
            /*handling duplicates*/
            if (std::find(vec_NS_names.begin(), vec_NS_names.end(), it) == vec_NS_names.end()) 
            {
                vec_NS_names.push_back(it);
                ++c;
            }
            else {
               amrex::Abort("NS_name " + it + " is duplicated!");
            }
        }
        temp_vec.clear();

        NS_num = vec_NS_names.size();
        amrex::Print() << "#####* Number of nanostructures, NS_num: " << NS_num << "\n";
    }
    else 
    {
        auto NS_num_isDefined = queryWithParser(pp,  "NS_num", NS_num);

        if(NS_num == 0) amrex::Abort("NS_num=0!, Either specify vector of NS_names or NS_num!");
        
        if(NS_num_isDefined) {
            for(int n=0; n < NS_num; ++n) 
            {
                vec_NS_names.push_back("NS_"+ std::to_string(n+1));
            }
        }
    }
    amrex::Print() << "#####* Number of nanostructures, NS_num: " << NS_num << "\n";
    amrex::Print() << "#####* Names: ";
    for (auto name: vec_NS_names) amrex::Print() << name << "\n";
    amrex::Print() << "\n";
}


void 
c_TransportSolver::Read_NSTypes(amrex::ParmParse& pp) 
{
    flag_isDefined_NS_type_default = pp.query("NS_type_default", NS_type_default);

    if(!flag_isDefined_NS_type_default) 
    {
        for (auto name: vec_NS_names)
        {
            amrex::ParmParse pp_ns(name);
            std::string type_str;
            pp_ns.get("type", type_str);
            map_NSNameToTypeStr[name] = type_str;
            amrex::Print() << "##### name & type: " << name << "  " << type_str << "\n";
        }
    }
}


void 
c_TransportSolver::InitData() 
{
    amrex::Real negf_init_time = amrex::second();

    amrex::Print() << "\n##### TRANSPORT PROPERTIES #####\n\n";

    amrex::ParmParse pp_transport("transport");
    Read_ControlFlags(pp_transport);

    Set_NEGFFolderDirectories();
    Create_NEGFFolderDirectories();

    auto& rCode    = c_Code::GetInstance();

    if(rCode.use_electrostatic)
    {
        auto& rGprop = rCode.get_GeometryProperties();
        _geom = &rGprop.geom;
        _dm   = &rGprop.dm;
        _ba   = &rGprop.ba;

        Assert_GatherAndDepositFields();

        if (!flag_isDefined_InitialDepositValue) Define_InitialDepositValue();

        Set_gate_terminal_type(gate_terminal_type_str);
    }
    num_field_sites_all_NS = Instantiate_Materials();

    if(rCode.use_electrostatic) Sum_ChargeDepositedByAllNS();

    Set_Broyden_Parallel();

}


void 
c_TransportSolver::Read_ControlFlags(amrex::ParmParse& pp) 
{
    pp.query("use_selfconsistent_potential", use_selfconsistent_potential);
    pp.query("use_negf", use_negf);
    amrex::Print() << "##### transport.use_selfconsistent_potential: " 
                   << use_selfconsistent_potential << "\n";
    amrex::Print() << "##### transport.use_negf: " 
                   << use_negf << "\n";
}


void 
c_TransportSolver::Set_NEGFFolderDirectories() 
{
    amrex::ParmParse pp_plot("plot");
    std::string foldername_str = "output";
    pp_plot.query("folder_name", foldername_str);
    CreateDirectory(foldername_str);

    negf_foldername_str = foldername_str + "/negf";
    common_foldername_str = negf_foldername_str + "/transport_common";
}


void 
c_TransportSolver::Create_NEGFFolderDirectories() 
{
    CreateDirectory(negf_foldername_str);
    CreateDirectory(common_foldername_str);
}


void 
c_TransportSolver::Read_GatherAndDepositFields(amrex::ParmParse& pp) 
{
    pp.get("NS_gather_field", NS_gather_field_str);
    amrex::Print() << "##### transport.NS_gather_field: " << NS_gather_field_str << "\n";

    pp.get("NS_deposit_field", NS_deposit_field_str);
    amrex::Print() << "##### transport.NS_deposit_field: " << NS_deposit_field_str << "\n";
}


void 
c_TransportSolver::Assert_GatherAndDepositFields() 
{
    if ( Evaluate_TypeOf_MacroStr(NS_gather_field_str) != 0 )
    {
        amrex::Abort("NS_gather_field " + NS_gather_field_str + " not defined in Mprop.");
    } 
    if ( Evaluate_TypeOf_MacroStr(NS_deposit_field_str) != 0 )
    {
        amrex::Abort("NS_deposit_field " + NS_deposit_field_str + " not defined in Mprop.");
    } 
}


void 
c_TransportSolver::Define_InitialDepositValue() 
{
    auto dxi = _geom->InvCellSizeArray();
    amrex::Real inv_vol = AMREX_D_TERM(dxi[0], *dxi[1], *dxi[2]);
    NS_initial_deposit_value = PhysConst::q_e*inv_vol;
}

void 
c_TransportSolver::Read_SelfConsistencyInput(amrex::ParmParse& pp) 
{
    flag_isDefined_InitialDepositValue = 
        queryWithParser(pp,"NS_initial_deposit_value", NS_initial_deposit_value);
    amrex::Print() << "##### transport.NS_initial_deposit_value: " << NS_initial_deposit_value << "\n";

    queryWithParser(pp,"Broyden_fraction", Broyden_Original_Fraction);
    amrex::Print() << "##### Broyden_fraction: " << Broyden_Original_Fraction << "\n";

    queryWithParser(pp,"Broyden_max_norm", Broyden_max_norm);
    amrex::Print() << "##### Broyden_max_norm: " << Broyden_max_norm << "\n";

    pp.query("Broyden_norm_type", Broyden_Norm_Type);
    amrex::Print() << "##### Broyden_norm_type: " << Broyden_Norm_Type << "\n";

    Broyden_Threshold_MaxStep = 200;
    pp.query("Broyden_threshold_maxstep", Broyden_Threshold_MaxStep);
    amrex::Print() << "##### Broyden_Threshold_MaxStep: " << Broyden_Threshold_MaxStep << "\n";

    pp.query("selfconsistency_algorithm", Algorithm_Type);
    amrex::Print() << "##### selfconsistency_algorithm: " << Algorithm_Type << "\n";

    pp.query("reset_with_previous_charge_distribution", flag_reset_with_previous_charge_distribution);
    amrex::Print() << "##### reset_with_previous_charge_distribution: " << flag_reset_with_previous_charge_distribution  << "\n";

    flag_initialize_inverse_jacobian = 0;
    pp.query("initialize_inverse_jacobian", flag_initialize_inverse_jacobian);
    amrex::Print() << "##### flag_initialize_inverse_jacobian: " << flag_initialize_inverse_jacobian  << "\n";

	if(flag_initialize_inverse_jacobian) Read_InverseJacobianFilename(pp);

}


void 
c_TransportSolver::Read_InverseJacobianFilename(amrex::ParmParse& pp) 
{
    amrex::ParmParse pp_default;
    int flag_restart = 0;
    pp_default.query("restart", flag_restart);

	if(flag_restart) 
	{		    
        int restart_step = 0;       
        getWithParser(pp_default,"restart_step", restart_step);

        std::string restart_folder_str  = 
            amrex::Concatenate(common_foldername_str + "/step", 
                               restart_step-1, negf_plt_name_digits);
        /*eg. output/negf/transport_common/step0001 for step 1*/

        inverse_jacobian_filename  = restart_folder_str + "/Jinv.dat";
        pp.query("inverse_jacobian_filename", inverse_jacobian_filename);
	}	
	else 
	{
        pp.get("inverse_jacobian_filename", inverse_jacobian_filename);
	}
    amrex::Print() << "##### inverse_jacobian_filename: " 
                   << inverse_jacobian_filename << "\n";
}


void 
c_TransportSolver::Read_DOSInput(amrex::ParmParse& pp) 
{
    pp.query("flag_compute_DOS", flag_compute_DOS);
    amrex::Print() << "##### flag_compute_DOS: " << flag_compute_DOS << "\n";

    pp.query("flag_write_LDOS", flag_write_LDOS);
    amrex::Print() << "##### flag_write_LDOS: " << flag_write_LDOS << "\n";

    pp.query("flag_write_LDOS_iter", flag_write_LDOS_iter);
    amrex::Print() << "##### flag_write_LDOS_iter: " << flag_write_LDOS_iter << "\n";

    pp.query("write_LDOS_iter_period", write_LDOS_iter_period);
    amrex::Print() << "##### write_LDOS_iter_period: " << write_LDOS_iter_period << "\n";

}


void 
c_TransportSolver::Read_GateTerminalType(amrex::ParmParse& pp) 
{
    pp.query("gate_terminal_type", gate_terminal_type_str);
}


void 
c_TransportSolver::Set_gate_terminal_type(const std::string gate_terminal_type_str) 
{
    if(doesKeyExist(gate_terminal_type_str))
    {
        amrex::Print() << "##### gate_terminal_type: " 
                       << gate_terminal_type_str   << "\n";
        gate_terminal_type = 
            c_TransportSolver::map_S_GateTerminalType.at(gate_terminal_type_str);
    }
    else {
        amrex::Abort("gate_terminal_type: " + Algorithm_Type 
                  + " is invalid. Valid type: EB or Boundary.");
    }
}


std::string 
c_TransportSolver::Get_NS_type_str(const std::string& name) 
{
    return flag_isDefined_NS_type_default ? NS_type_default : map_NSNameToTypeStr.at(name);
}


int
c_TransportSolver::Instantiate_Materials() 
{
    amrex::Vector<int> NS_field_sites_cumulative(1,0);
    int NS_id_counter = 0;

    for (auto name: vec_NS_names)
    {
        amrex::Print() << "##### Instantiating material: " << name << "\n";

        int field_sites = 0;
        int NS_field_sites_offset = NS_field_sites_cumulative.back();
        switch (c_TransportSolver::map_NSType_enum.at(Get_NS_type_str(name)))
        {
            case s_NS_Type::CNT:
            {
                using T = c_CNT;
                vp_CNT.push_back(std::make_unique<c_Nanostructure<T>>(*_geom, *_dm, *_ba,
                                                                       name, 
                                                                       NS_id_counter,
                                                                       NS_gather_field_str, 
                                                                       NS_deposit_field_str, 
                                                                       NS_initial_deposit_value,
                                                                       use_negf,
								                                       negf_foldername_str,
                                                                       NS_field_sites_offset
                                                                      ));
                field_sites = vp_CNT.back()->get_num_field_sites();
                break;
            }
            case s_NS_Type::Graphene:
            {
                using T = c_Graphene;
 
                vp_Graphene.push_back(std::make_unique<c_Nanostructure<T>>(*_geom, *_dm, *_ba,
                                                                            name, 
                                                                            NS_id_counter,
                                                                            NS_gather_field_str, 
                                                                            NS_deposit_field_str, 
                                                                            NS_initial_deposit_value,
                                                                            use_negf,
	                                    								    negf_foldername_str,
                                                                            NS_field_sites_offset
                                                                           ));
                field_sites = vp_Graphene.back()->get_num_field_sites();
                amrex::Abort("NS_type graphene is not yet defined.");
                break; 
            }
            case s_NS_Type::Silicon:
            {
                amrex::Abort("NS_type silicon is not yet defined.");
            }
            default:
            {
                amrex::Abort("NS_type " + Get_NS_type_str(name) + " is not supported.");
            }
        }
        int cumulative_sites = NS_field_sites_cumulative.back() + field_sites;
        NS_field_sites_cumulative.push_back(cumulative_sites);

        NS_id_counter++;
    }

    amrex::Print() << "NS_field_sites_cumulative: \n";
    for (auto offset: NS_field_sites_cumulative) 
    {
        amrex::Print() << " " << offset << "\n";
    }

    return NS_field_sites_cumulative.back();
}


void 
c_TransportSolver::Sum_ChargeDepositedByAllNS() 
{
    auto& rCode  = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& rMprop = rCode.get_MacroscopicProperties();

    amrex::MultiFab* p_mf_deposit = rMprop.get_p_mf(NS_deposit_field_str);
    p_mf_deposit->SumBoundary(rGprop.geom.periodicity());
}


void
c_TransportSolver::Set_CommonStepFolder(const int step)
{

    common_step_folder_str  = amrex::Concatenate(common_foldername_str + "/step", 
                                                 step, negf_plt_name_digits);

    amrex::Print() << "common_step_folder_str: " 
                   << common_step_folder_str << "\n";
   
    CreateDirectory(common_step_folder_str);
    /*eg. output/negf/common/step0001 for step 1*/

}


void 
c_TransportSolver::Solve(const int step, const amrex::Real time) 
{

    auto& rCode    = c_Code::GetInstance();
    auto& rMprop   = rCode.get_MacroscopicProperties();
    auto& rMLMG    = rCode.get_MLMGSolver();
    auto& rOutput  = rCode.get_Output();
    auto& rPostPro = rCode.get_PostProcessor();

    max_iter = 0;
    total_intg_pts_in_all_iter = 0;
    m_step = step;

    for (int c=0; c < vp_CNT.size(); ++c)
    {
       vp_CNT[c]->Set_StepFilenameString(step);
       Set_CommonStepFolder(step);
    }

    amrex::Real time_counter[7] = {0., 0., 0., 0., 0., 0., 0.};
    amrex::Real total_time_counter_diff[7] = {0., 0., 0., 0., 0., 0., 0.};

    if(rCode.use_electrostatic) 
    {	
        BL_PROFILE_VAR("Part1_to_6_sum", part1_to_6_sum_counter);

        bool update_surface_soln_flag = true;	   
        do 
        {
           amrex::Print() << "\n\n##### Self-Consistent Iteration: " << max_iter << " #####\n";
           if (Broyden_Step > Broyden_Threshold_MaxStep)
           {
               amrex::Abort("Broyden_Step has exceeded the Broyden_Threshold_MaxStep!");
           }
    
           //Part 1: Electrostatics
           time_counter[0] = amrex::second();
           

           rMprop.ReInitializeMacroparam(NS_gather_field_str);
           rMLMG.UpdateBoundaryConditions(update_surface_soln_flag);

           auto mlmg_solve_time = rMLMG.Solve_PoissonEqn();

           rPostPro.Compute();
           //rOutput.WriteOutput(max_iter+100, time);

           time_counter[1] = amrex::second();

           //Part 2: Gather
           for (int c=0; c < vp_CNT.size(); ++c)
           {
	           if( vp_CNT[c]->get_flag_write_at_iter() )
	           {
	                vp_CNT[c]->Set_IterationFilenameString(max_iter);
               }

	           vp_CNT[c]->Gather_MeshAttributeAtAtoms();  
            }
           time_counter[2] = amrex::second();

           //Part 3: NEGF
           int total_intg_pts_in_this_iter = 0;
           for (int c=0; c < vp_CNT.size(); ++c)
           {
	           if(update_surface_soln_flag)
	           {
                   Set_TerminalBiasesAndContactPotential(vp_CNT[c]);   
               }
               amrex::Print() << " Vds: " << Vds << " V, Vgs: " << Vgs << " V\n";

               #ifdef AMREX_USE_GPU
               vp_CNT[c]->Solve_NEGF(d_n_curr_out_data, max_iter);
               #else
               vp_CNT[c]->Solve_NEGF(h_n_curr_out_data, max_iter);
               #endif
               total_intg_pts_in_this_iter += vp_CNT[c]->get_Total_Integration_Pts();
               //CopyFromNS_ChargeComputedFromNEGF(vp_CNT[c]);
           }
           total_intg_pts_in_all_iter += total_intg_pts_in_this_iter;
           time_counter[3] = amrex::second();

           //Part 4: Self-consistency
           Perform_SelfConsistencyAlgorithm();

           time_counter[4] = amrex::second();

           //Part 5: Deposit
           rMprop.ReInitializeMacroparam(NS_deposit_field_str);

           for (int c=0; c < vp_CNT.size(); ++c)
           {
              CopyToNS_ChargeComputedUsingSelfConsistencyAlgorithm(vp_CNT[c]);

              vp_CNT[c]->Deposit_AtomAttributeToMesh();

	           if( vp_CNT[c]->get_flag_write_at_iter() ) 
               {
	               vp_CNT[c]->Write_InputInducedCharge(vp_CNT[c]->get_iter_filename(), 
                                                       n_curr_in_glo_data); 
               }
               if(ParallelDescriptor::IOProcessor())
               {
                   n_curr_in_glo_data.clear();
               }
	       }

           Sum_ChargeDepositedByAllNS();

           time_counter[5] = amrex::second();

           //Part 6: Write data
           for (int c=0; c < vp_CNT.size(); ++c)
           {
	          if( vp_CNT[c]->get_flag_write_at_iter() ) 
              {
                  bool compute_current = false;
                  Write_MoreDataAndComputeCurrent
                      (vp_CNT[c], 
                       vp_CNT[c]->get_iter_filename(), compute_current);
    	      }

              if(flag_write_LDOS_iter and (max_iter+1)%write_LDOS_iter_period == 0) 
              {
                  std::string iter_dos_foldername_str 
                  = amrex::Concatenate(vp_CNT[c]->get_iter_foldername() + "/DOS_iter", 
                          max_iter, negf_plt_name_digits);

                  CreateDirectory(iter_dos_foldername_str);
                  // e.g.: /negf/cnt/step0001_iter/LDOS_iter0002/
                  
                  vp_CNT[c]->Compute_DensityOfStates(iter_dos_foldername_str, 
                          flag_write_LDOS_iter);
              }
	       }
           time_counter[6] = amrex::second();

           update_surface_soln_flag = false;
           max_iter += 1;

           total_time_counter_diff[0] += time_counter[1] - time_counter[0];
           total_time_counter_diff[1] += time_counter[2] - time_counter[1];
           total_time_counter_diff[2] += time_counter[3] - time_counter[2];
           total_time_counter_diff[3] += time_counter[4] - time_counter[3];
           total_time_counter_diff[4] += time_counter[5] - time_counter[4];
           total_time_counter_diff[5] += time_counter[6] - time_counter[5];
           total_time_counter_diff[6] += (time_counter[3] - time_counter[2])/total_intg_pts_in_this_iter;

           amrex::Print() << " Times for: \n";
           amrex::Print() << " Electrostatics:   " << time_counter[1] - time_counter[0] << "\n";    
           amrex::Print() << " Gathering field:  " << time_counter[2] - time_counter[1] << "\n";    
           amrex::Print() << " NEGF:             " << time_counter[3] - time_counter[2] << "\n";    
           amrex::Print() << " NEGF (per intg pt.), intg_pts:" 
                                  << (time_counter[3] - time_counter[2])/total_intg_pts_in_this_iter << "  "
                                  << total_intg_pts_in_this_iter << "\n";    
           amrex::Print() << " Self-consistency: " << time_counter[4] - time_counter[3] << "\n";    
           amrex::Print() << " Deposit charge:   " << time_counter[5] - time_counter[4] << "\n";    
           amrex::Print() << " Writing at iter:  " << time_counter[6] - time_counter[5] << "\n";           
           amrex::Print() << " Total time (write excluded):   " << time_counter[5] - time_counter[0] << "\n";    

        } while(Broyden_Norm > Broyden_max_norm);    

        BL_PROFILE_VAR_STOP(part1_to_6_sum_counter);

        Obtain_maximum_time(total_time_counter_diff);


        /* LDOS computation is before current computation because
         * if total conductance is calculated during LDOS computation,
         * and it is written in the same file used for writing current
         * after current computation.*/
        if(flag_compute_DOS) 
        {
            for (int c=0; c < vp_CNT.size(); ++c)
            {
                std::string 
                dos_step_foldername_str 
                = amrex::Concatenate(vp_CNT[c]->get_step_foldername() + "/DOS_step", 
                        step, negf_plt_name_digits);

                CreateDirectory(dos_step_foldername_str);
                // e.g.: /negf/cnt/DOS_step0001/

                vp_CNT[c]->Compute_DensityOfStates(dos_step_foldername_str, flag_write_LDOS);
            }
        }

        //Part 7: current computation & writing data
        amrex::Real time_for_current = amrex::second();
        for (int c=0; c < vp_CNT.size(); ++c)
        {
            bool compute_current = true;
            Write_MoreDataAndComputeCurrent(vp_CNT[c],
                                            vp_CNT[c]->get_step_filename(), 
                                            compute_current);
        }

        amrex::Print() << "Time for current computation & writing data:   " 
                       << amrex::second() - time_for_current << "\n";    



        Reset_ForNextBiasStep();

   } //if use electrostatics
   else 
   {
       for (int c=0; c < vp_CNT.size(); ++c)
       {
           RealTable1D RhoInduced; /*this is not correct but added just so Solve_NEGF compiles*/
           vp_CNT[c]->Solve_NEGF(RhoInduced, 0);

           bool compute_current = true;

           Write_MoreDataAndComputeCurrent(vp_CNT[c], 
                                           vp_CNT[c]->get_step_filename(), 
                                           compute_current);
       }
   }

}


template<typename NSType>
void
c_TransportSolver:: Set_TerminalBiasesAndContactPotential(NSType const& NS)
{
    auto& rCode  = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& rBC    = rCode.get_BoundaryConditions();

    if(!NS->is_contact_mu_specified()) 
    {
         amrex::Real V_contact[NUM_CONTACTS] = {0., 0.};
        
         amrex::Vector<amrex::Real> ep(NUM_CONTACTS,0);

         const auto* ps = NS->get_Contact_Parser_String();
         for(int k=0; k<NUM_CONTACTS; ++k) 
         {    
             V_contact[k] = rGprop.pEB->Read_SurfSoln(ps[k]);
                
             amrex::Print() << "\n Updated terminal voltage: " << k << "  " << V_contact[k] << " V\n";

             ep[k] = NS->get_Fermi_level() - V_contact[k];
         }
         NS->set_Contact_Electrochemical_Potential(ep);
    
         Vds = V_contact[1] - V_contact[0];

         if(gate_terminal_type == Gate_Terminal_Type::EB) 
         {
            Vgs = rGprop.pEB->Read_SurfSoln(NS->get_Gate_String()) - V_contact[0];
         }
         else if (gate_terminal_type == Gate_Terminal_Type::Boundary)
         {
            Vgs = rBC.Read_BoundarySoln(NS->get_Gate_String()) - V_contact[0];
         }

    }
    else 
    {
         amrex::Real ep_s = NS->get_Source_Electrochemical_Potential();
         amrex::Real ep_d = NS->get_Drain_Electrochemical_Potential();

         Vds = ep_s - ep_d;
        
         amrex::Real GV = 0.;
         if(gate_terminal_type == Gate_Terminal_Type::EB) 
         {
            GV  = rGprop.pEB->Read_SurfSoln(NS->get_Gate_String());
         }
         else if (gate_terminal_type == Gate_Terminal_Type::Boundary) 
         {
            GV = rBC.Read_BoundarySoln(NS->get_Gate_String());
         }

         amrex::Print() << "\n Updated gate voltage: " << GV << " V\n";
         Vgs = GV - (NS->get_Fermi_level() - ep_s);
    }
}


template<typename NSType>
void
c_TransportSolver:: CopyFromNS_ChargeComputedFromNEGF(NSType const& NS)
{
   /* (?) Need a strategy to gather data for multiple CNTs*/
   /* (?) Future: if condition to check if proc was used for charge computation of this nanostructure*/
   /* (?) Future: for multiple nanotubes, then gather the data into a global array and then partition */
   #ifdef BROYDEN_SKIP_GPU_OPTIMIZATION
   h_n_curr_out_data.copy(NS->h_RhoInduced_loc_data); 
   #else
   d_n_curr_out_data.copy(NS->h_RhoInduced_loc_data); 
   #endif
}


template<typename NSType>
void
c_TransportSolver:: CopyToNS_ChargeComputedUsingSelfConsistencyAlgorithm(NSType const& NS)
{

    int begin = site_size_loc_cumulative[NS->get_NS_Id()]; 
    int end = site_size_loc_cumulative[NS->get_NS_Id() + 1];
    int site_size_loc = end - begin;

    auto const& h_n_curr_in      = h_n_curr_in_data.table();

    #ifdef AMREX_USE_GPU
    auto const& d_n_curr_in      = d_n_curr_in_data.table();

    h_n_curr_in_data.resize({0}, {site_size_loc}, The_Pinned_Arena());

    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, 
                          d_n_curr_in.p + begin , d_n_curr_in.p + end, h_n_curr_in.p);
    amrex::Gpu::streamSynchronize();
    #endif

    if(ParallelDescriptor::IOProcessor())
    {
        const int Hsize = NS->get_num_field_sites();
        n_curr_in_glo_data.resize({0},{Hsize}, The_Pinned_Arena());
    }
    auto const& n_curr_in_glo  = n_curr_in_glo_data.table();

    /*offset necessary for multiple NS*/
    MPI_Gatherv(&h_n_curr_in(0),
                 NS->MPI_recv_count[my_rank],
                 MPI_DOUBLE,
                &n_curr_in_glo(0),
                 NS->MPI_recv_count.data(),
                 NS->MPI_recv_disp.data(),
                 MPI_DOUBLE,
                 ParallelDescriptor::IOProcessorNumber(),
                 ParallelDescriptor::Communicator());

    //amrex::Print() << "n_curr_in_glo in CopyToNS: \n";
    //if (ParallelDescriptor::IOProcessor())
    //{
    //    for(int n=0; n < 5; ++n) {
    //       amrex::Print() << n << "  " <<  n_curr_in_glo(n) << "\n";
    //    }
    //}
    NS->Scatterv_BroydenComputed_GlobalCharge(n_curr_in_glo_data); 
}


template<typename NSType>
void
c_TransportSolver:: Write_MoreDataAndComputeCurrent(NSType const& NS, 
                                                    std::string const& write_filename,
                                                    bool const compute_current_flag)
{
    //Note: n_curr_out_glo was output from negf and input to broyden.
    //n_curr_in is the broyden predicted charge for next iteration.
    //NEGF->n_curr_out -> Broyden->n_curr_in -> Electrostatics -> NEGF.
	Create_Global_Output_Data(NS); 
    NS->Write_Data(write_filename, n_curr_out_glo_data, Norm_glo_data);

    if (ParallelDescriptor::IOProcessor())
    {
        n_curr_out_glo_data.clear();
        Norm_glo_data.clear();
    }

    if(compute_current_flag) 
    {
        NS->Compute_Current();
        NS->Write_Current(m_step, Vds, Vgs, total_intg_pts_in_all_iter, max_iter, Broyden_fraction, Broyden_Scalar);
    }
}


template<typename NSType>
void
c_TransportSolver:: Create_Global_Output_Data(NSType const& NS) 
{

    #ifdef BROYDEN_SKIP_GPU_OPTIMIZATION 
    auto const& h_n_curr_out = h_n_curr_out_data.table();
    auto const& h_Norm = h_Norm_data.table();
    #else 
    /*only select data need to be copied for multiple NS*/

    int begin = site_size_loc_cumulative[NS->get_NS_Id()]; 
    int end = site_size_loc_cumulative[NS->get_NS_Id() + 1];
    int site_size_loc = end - begin;

    amrex::Print() << "Create_Global: begin/end/site_size_loc: " << begin << " " << end << " " 
                                                                 << site_size_loc << "\n";

    h_n_curr_out_data.resize({0}, {site_size_loc}, The_Pinned_Arena());
    h_Norm_data.resize({0}, {site_size_loc}, The_Pinned_Arena());

    auto const& h_n_curr_out = h_n_curr_out_data.table();
    auto const& h_Norm = h_Norm_data.table();
    
    auto const& d_n_curr_out = d_n_curr_out_data.const_table();
    auto const& d_Norm = d_Norm_data.const_table();
    
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, 
                          d_n_curr_out.p + begin , d_n_curr_out.p + end, h_n_curr_out.p);

    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, 
                          d_Norm.p + begin , d_Norm.p + end, h_Norm.p);

    amrex::Gpu::streamSynchronize();

    amrex::Print() << "h_n_curr_out(0): " << h_n_curr_out(0) << "\n";
    amrex::Print() << "h_Norm(0): " << h_Norm(0) << "\n";
    #endif
 
    if (ParallelDescriptor::IOProcessor())
    {
        const int Hsize = NS->get_num_field_sites();
        n_curr_out_glo_data.resize({0},{Hsize}, The_Pinned_Arena());
        Norm_glo_data.resize({0},{Hsize}, The_Pinned_Arena());
    }
 
    auto const& n_curr_out_glo = n_curr_out_glo_data.table();
    auto const& Norm_glo       = Norm_glo_data.table();

    /*offset necessary for multiple NS*/
    MPI_Gatherv(&h_n_curr_out(0),
                NS->MPI_recv_count[my_rank],
                MPI_DOUBLE,
               &n_curr_out_glo(0),
                NS->MPI_recv_count.data(),
                NS->MPI_recv_disp.data(),
                MPI_DOUBLE,
                ParallelDescriptor::IOProcessorNumber(),
                ParallelDescriptor::Communicator());

    /*offset necessary for multiple NS*/
    MPI_Gatherv(&h_Norm(0),
                NS->MPI_recv_count[my_rank],
                MPI_DOUBLE,
               &Norm_glo(0),
                NS->MPI_recv_count.data(),
                NS->MPI_recv_disp.data(),
                MPI_DOUBLE,
                ParallelDescriptor::IOProcessorNumber(),
                ParallelDescriptor::Communicator());

    //if (ParallelDescriptor::IOProcessor())
    //{
    //    amrex::Print() << "\nPrinting n_curr_out_glo_data \n";
    //    for (int n=0; n <num_field_sites_all_NS; ++n)
    //    {
    //        amrex::Print()
    //        << n
    //        << std::setw(20) <<  n_curr_out_glo(n)
    //        << "\n";
    //    }
    //}
    
    #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
    h_n_curr_out_data.clear();
    h_Norm_data.clear();
    #endif

}


void
c_TransportSolver:: Perform_SelfConsistencyAlgorithm() 
{
   switch(c_TransportSolver::map_AlgorithmType.at(Algorithm_Type))
   {
       case s_Algorithm_Type::broyden_first:
       {
           amrex::Abort("At present, `Broyden's first' algorithm exists with only\
                         serial implementation (BROYDEN_PARALLEL=FALSE).");
           break;
       }
       case s_Algorithm_Type::broyden_second:
       {
           #ifdef BROYDEN_SKIP_GPU_OPTIMIZATION
           Execute_Broyden_Modified_Second_Algorithm_Parallel_SkipGPU();
           #else
           Execute_Broyden_Modified_Second_Algorithm_Parallel(); 
           #endif
           break;
       }
       case s_Algorithm_Type::simple_mixing:
       {
           amrex::Abort("At present, the `simple mixing' algorithm exists with only\
                         serial implementation (BROYDEN_PARALLEL=FALSE).");
           break;
       }
       default:
       {
           amrex::Abort("selfconsistency_algorithm, " + Algorithm_Type + ", is not yet defined.");
       }
   }
}


void
c_TransportSolver:: Obtain_maximum_time(amrex::Real const* total_time_counter_diff)
{

    const int num_var = 7;
    amrex::Real total_max_time_for_current_step[num_var] = {0., 0., 0., 0., 0., 0., 0.};

    MPI_Reduce(total_time_counter_diff,
               total_max_time_for_current_step,
               num_var,
               MPI_DOUBLE,
               MPI_MAX,
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

    if(ParallelDescriptor::IOProcessor()) 
    {
        amrex::Real avg_curr[num_var] = {total_max_time_for_current_step[0]/max_iter, 
                                         total_max_time_for_current_step[1]/max_iter,
                                         total_max_time_for_current_step[2]/max_iter,
                                         total_max_time_for_current_step[3]/max_iter,
                                         total_max_time_for_current_step[4]/max_iter,
                                         total_max_time_for_current_step[5]/max_iter,
                                         total_max_time_for_current_step[6]/max_iter
                                        };
        amrex::Print() << "\nIterations in this step: " << max_iter << "\n";

        amrex::Print() << "Avg. max times for: \n ";
        amrex::Print() << " Electrostatics:   " 
                       << avg_curr[0] << "\n";    

        amrex::Print() << " Gather field:     " 
                       << avg_curr[1] << "\n";    

        amrex::Print() << " NEGF:             " 
                       << avg_curr[2] << "\n";    

        amrex::Print() << " NEGF (per intg pt):    "
                       << avg_curr[6] << "\n";    

        amrex::Print() << " Self-consistency: "             
                       << avg_curr[3] << "\n";    

        amrex::Print() << " Deposit:          " 
                       << avg_curr[4] << "\n";    

        amrex::Print() << " Write at iter:    "
                       << avg_curr[5] << "\n";    

        amrex::Print() << " Total time (write excluded): " << avg_curr[0] + 
                                                              avg_curr[1] + 
                                                              avg_curr[2] +
                                                              avg_curr[3] +
                                                              avg_curr[4] << "\n";

        for (int i=0; i<num_var; ++i) {
           total_max_time_across_all_steps[i] += total_max_time_for_current_step[i];
        }

        total_iter += max_iter;

        amrex::Real avg_all[num_var] = {total_max_time_across_all_steps[0]/total_iter, 
                                        total_max_time_across_all_steps[1]/total_iter,
                                        total_max_time_across_all_steps[2]/total_iter,
                                        total_max_time_across_all_steps[3]/total_iter,
                                        total_max_time_across_all_steps[4]/total_iter,
                                        total_max_time_across_all_steps[5]/total_iter,
                                        total_max_time_across_all_steps[6]/total_iter
                                       };        
      
        amrex::Real avg_total = avg_all[0] + avg_all[1] + avg_all[2] +
                                avg_all[3] + avg_all[4];

        amrex::Print() << "\nTotal iterations so far: " << total_iter << "\n";
        amrex::Print() << "Avg. time over all steps for: \n";
        amrex::Print() << " Electrostatics:   "
                       << avg_all[0] << std::setw(15) << (avg_all[0]/avg_total)*100 << " %" << "\n";    
        amrex::Print() << " Gather field:     " 
                       << avg_all[1] << std::setw(15) << (avg_all[1]/avg_total)*100 << " %" << "\n";    
        amrex::Print() << " NEGF:             " 
                       << avg_all[2] << std::setw(15) << (avg_all[2]/avg_total)*100 << " %" << "\n";    
        amrex::Print() << " NEGF (per intg pt):" << avg_all[6] << "\n";    
        amrex::Print() << " Self-Consistency: " 
                       << avg_all[3] << std::setw(15) << (avg_all[3]/avg_total)*100 << " %" << "\n";    
        amrex::Print() << " Deposit:          " 
                       << avg_all[4] << std::setw(15) << (avg_all[4]/avg_total)*100 << " %" << "\n";    
        
        amrex::Print() << " Write at iter:    " << avg_all[5]  << "\n";    

        amrex::Print() << " Total time (write excluded): " << avg_total << "\n";

    }
}

void
c_TransportSolver:: Reset_ForNextBiasStep() 
{
    Reset_Broyden_Parallel();

    //rMprop.ReInitializeMacroparam(NS_deposit_field_str);
    MPI_Barrier(ParallelDescriptor::Communicator());
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
