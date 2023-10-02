#include "Transport.H"
#include "Transport_Table_ReadWrite.H"

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
#include <AMReX_GpuContainers.H>

using namespace amrex;

c_TransportSolver::c_TransportSolver()
{
    ReadData();
}


c_TransportSolver::~c_TransportSolver()
{
    vp_CNT.clear();
    vp_Graphene.clear();
    //vp_Silicon.clear();
    
    #ifdef BROYDEN_PARALLEL
    Deallocate_Broyden_Parallel();
    #else
    Deallocate_Broyden_Serial();
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

        Broyden_Threshold_MaxStep = 200;
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

    Broyden_Original_Fraction = NS_Broyden_frac;
    Broyden_Norm_Type         = NS_Broyden_norm_type;
    amrex::Print() << "#####* Broyden_Original_Fraction: "  << Broyden_Original_Fraction  << "\n";
    amrex::Print() << "#####* Broyden_Norm_Type: "          << Broyden_Norm_Type          << "\n";

    #ifdef BROYDEN_PARALLEL
    amrex::Print() << "Setting Broyden PARALLEL\n";
    Set_Broyden_Parallel();
    #else
    amrex::Print() << "Setting Broyden SERIAL\n";
    Set_Broyden();
    #endif

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
           if (Broyden_Step > Broyden_Threshold_MaxStep)
           {
               amrex::Abort("Broyden_Step has exceeded the Broyden_Threshold_MaxStep!");
           }

           amrex::Print() << "\n\nSelf-consistent iteration: " << max_iter << "\n";

           rMprop.ReInitializeMacroparam(NS_gather_field_str);
           rMLMG.UpdateBoundaryConditions(update_surface_soln_flag);

           auto mlmg_solve_time = rMLMG.Solve_PoissonEqn();
           total_mlmg_solve_time += mlmg_solve_time;
           amrex::Print() << "\nmlmg_solve_time: " << mlmg_solve_time << "\n";

           rPostPro.Compute();

           for (int c=0; c < vp_CNT.size(); ++c)
           {
	           vp_CNT[c]->Set_IterationFilenameString(max_iter);

                #ifdef BROYDEN_PARALLEL
	            vp_CNT[c]->Write_InputInducedCharge(vp_CNT[c]->iter_filename_str, n_curr_in_glo_data); 
                #else
	            vp_CNT[c]->Write_InputInducedCharge(vp_CNT[c]->iter_filename_str, h_n_curr_in_data); 
                #endif  

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

               #ifdef BROYDEN_PARALLEL
                   /* (?) Need a strategy to gather data for multiple CNTs*/
                   /* (?) Future: if condition to check if proc was used for charge computation of this nanostructure*/
                   /* (?) Future: for multiple nanotubes, then gather the data into a global array and then partition */
                   #ifdef BROYDEN_SKIP_GPU_OPTIMIZATION
                   h_n_curr_out_data.copy(vp_CNT[c]->h_RhoInduced_loc_data); 
                   #else
                   d_n_curr_out_data.copy(vp_CNT[c]->h_RhoInduced_loc_data); 
                   #endif
               #else
	           vp_CNT[c]->Gatherv_NEGFComputed_LocalCharge(h_n_curr_out_data); 
               #endif
            }

            switch(map_AlgorithmType[Algorithm_Type])
            {
                case s_Algorithm::Type::broyden_first:
                {
		            #ifdef BROYDEN_PARALLEL	
                    amrex::Abort("At present, `Broyden's first' algorithm exists with only\
                                  serial implementation (BROYDEN_PARALLEL=FALSE).");
                    #else
	                Execute_Broyden_First_Algorithm(); 
                    #endif
                    break;
                }
                case s_Algorithm::Type::broyden_second:
                {
		            #ifdef BROYDEN_PARALLEL	
                        #ifdef BROYDEN_SKIP_GPU_OPTIMIZATION
                        Execute_Broyden_Modified_Second_Algorithm_Parallel_SkipGPU();
                        #else
	                    Execute_Broyden_Modified_Second_Algorithm_Parallel(); 
                        #endif
                    #else
    	            Execute_Broyden_Modified_Second_Algorithm(); 
                    #endif
                    break;
                }
                case s_Algorithm::Type::simple_mixing:
                {
		            #ifdef BROYDEN_PARALLEL	
                    amrex::Abort("At present, the `simple mixing' algorithm exists with only\
                                  serial implementation (BROYDEN_PARALLEL=FALSE).");
                    #else
	                Execute_Simple_Mixing_Algorithm(); 
                    #endif
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
               /*(?) Future: Need generalization for multiple CNTs*/
	       #ifdef BROYDEN_PARALLEL	
               vp_CNT[c]->Scatterv_BroydenComputed_GlobalCharge(n_curr_in_glo_data); 
               #else
               vp_CNT[c]->Scatterv_BroydenComputed_GlobalCharge(h_n_curr_in_data);
               #endif

	           if( vp_CNT[c]->write_at_iter )
	           {
	               Create_Global_Output_Data(); 
                   #ifdef BROYDEN_PARALLEL
                   vp_CNT[c]->Write_Data(vp_CNT[c]->iter_filename_str, n_curr_out_glo_data, Norm_glo_data);
                   #else
                   vp_CNT[c]->Write_Data(vp_CNT[c]->iter_filename_str, h_n_curr_out_data, h_Norm_data);
                   #endif
    	       }
               vp_CNT[c]->Deposit_AtomAttributeToMesh();
	        }
            update_surface_soln_flag = false;
            max_iter += 1;

        } while(Broyden_Norm > Broyden_max_norm);    


        #ifdef BROYDEN_PARALLEL
        Create_Global_Output_Data(); /*May need to be before & outside the forloop for multiple NS*/
        #endif
        for (int c=0; c < vp_CNT.size(); ++c)
        {
            #ifdef BROYDEN_PARALLEL
            vp_CNT[c]->Write_Data(vp_CNT[c]->step_filename_str, n_curr_out_glo_data, Norm_glo_data);
            #else
            vp_CNT[c]->Write_Data(vp_CNT[c]->step_filename_str, h_n_curr_out_data, h_Norm_data);

            if(map_AlgorithmType[Algorithm_Type] == s_Algorithm::Type::broyden_first) 
	        {
                Write_Table2D(h_Jinv_curr_data, common_step_folder_str + "/Jinv.dat", "Jinv");
	        }
            #endif

            vp_CNT[c]->Compute_Current();

            vp_CNT[c]->Write_Current(step, Vds, Vgs, Broyden_Step, max_iter, Broyden_fraction, Broyden_Scalar);
             
            if(rCode.use_electrostatic)
            {
                #ifdef BROYDEN_PARALLEL
                Reset_Broyden_Parallel();
                #else
                Reset_Broyden();
                #endif
            }
           //rMprop.ReInitializeMacroparam(NS_deposit_field_str);
        }

        #ifdef BROYDEN_PARALLEL
        Clear_Global_Output_Data();
        #endif

	MPI_Barrier(ParallelDescriptor::Communicator());

        amrex::Print() << "\nAverage mlmg time for self-consistency (s): " << total_mlmg_solve_time / max_iter << "\n";

   } //use electrostatics
   else 
   {
       for (int c=0; c < vp_CNT.size(); ++c)
       {
           vp_CNT[c]->Solve_NEGF();

           vp_CNT[c]->Write_Data(vp_CNT[c]->step_filename_str, h_n_curr_out_data, h_Norm_data); 

           vp_CNT[c]->Compute_Current();

           vp_CNT[c]->Write_Current(step, Vds, Vgs, Broyden_Step, max_iter, Broyden_fraction, Broyden_Scalar);
       }
   }

}

#ifdef BROYDEN_PARALLEL
void
c_TransportSolver:: Create_Global_Output_Data() 
{
    #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
    h_n_curr_out_data.resize({0}, {site_size_loc}, The_Pinned_Arena());
    h_Norm_data.resize({0}, {site_size_loc}, The_Pinned_Arena());

    h_n_curr_out_data.copy(d_n_curr_out_data); 
    h_Norm_data.copy(d_Norm_data); 
    amrex::Gpu::streamSynchronize();
    #endif
 
    auto const& n_curr_out = h_n_curr_out_data.table();
    auto const& Norm       = h_Norm_data.table();

    if (ParallelDescriptor::IOProcessor())
    {
        if(n_curr_out_glo_data.size() != 0) 
        {
           n_curr_out_glo_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        }
        if(Norm_glo_data.size() != 0) 
        {
           Norm_glo_data.resize({0},{num_field_sites_all_NS}, The_Pinned_Arena());
        }
    }
 
    auto const& n_curr_out_glo = n_curr_out_glo_data.table();
    auto const& Norm_glo       = Norm_glo_data.table();
 
    MPI_Gatherv(&n_curr_out(0),
                site_size_loc,
                MPI_DOUBLE,
               &n_curr_out_glo(0),
                MPI_recv_count.data(),
                MPI_recv_disp.data(),
                MPI_DOUBLE,
                ParallelDescriptor::IOProcessorNumber(),
                ParallelDescriptor::Communicator());

    MPI_Gatherv(&Norm(0),
                site_size_loc,
                MPI_DOUBLE,
               &Norm_glo(0),
                MPI_recv_count.data(),
                MPI_recv_disp.data(),
                MPI_DOUBLE,
                ParallelDescriptor::IOProcessorNumber(),
                ParallelDescriptor::Communicator());

    #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
    h_n_curr_out_data.clear();
    h_Norm_data.clear();
    #endif
}


void
c_TransportSolver:: Clear_Global_Output_Data() 
{
     if (ParallelDescriptor::IOProcessor())
     {
         n_curr_out_glo_data.clear();
         Norm_glo_data.clear();
     }
}
#endif

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
