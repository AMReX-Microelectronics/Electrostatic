#include "Transport.H"

#include <limits>

using namespace amrex;


#ifdef BROYDEN_PARALLEL
  #ifdef BROYDEN_SKIP_GPU_OPTIMIZATION
    void 
    c_TransportSolver::Execute_Broyden_Modified_Second_Algorithm_Parallel_SkipGPU()
    {
            //amrex::Print() << "Execute_Broyden_Modified_Second_Algorithm_Parallel_SkipGPU\n";    
            //amrex::Print() << "\nBroydenStep: " << Broyden_Step  
    		//               << ",  fraction: "   << Broyden_fraction 
    		//               << ",  scalar: " << Broyden_Scalar<< "\n";
    
            /*vectors*/
            auto const& n_curr_in      = h_n_curr_in_data.table();
            auto const& n_curr_out     = h_n_curr_out_data.table();
            auto const& n_prev_in      = h_n_prev_in_data.table();
            
            auto const& F_curr         = h_F_curr_data.table();
            auto const& delta_F_curr   = h_delta_F_curr_data.table();
            auto const& Norm           = h_Norm_data.table();
    
            /*matrices*/
            auto const& VmatTran       = h_VmatTran_data.table();
            auto const& Wmat           = h_Wmat_data.table();
    
            /*temp vectors*/
            auto const& sum_vector      = h_sum_vector_data.table();
            auto const& intermed_vector = h_intermed_vector_data.table();
    
            /*Initialize*/
            SetVal_RealTable1D(h_Norm_data, 0.);
            SetVal_RealTable1D(h_sum_vector_data, 0.);
            SetVal_RealTable1D(h_intermed_vector_data, 0.);
    
            /*Evaluate local (absolute or relative) and L2 norms*/
            Broyden_NormSum_Curr = 0.; 
            switch(map_NormType[Broyden_Norm_Type])
            {
                case s_Norm::Type::Absolute:
                {
                    for(int site=0; site < site_size_loc_all_NS; ++site)
                    {
                        amrex::Real Fcurr = n_curr_in(site) - n_curr_out(site);
                        Norm(site) = fabs(Fcurr);
                        Broyden_NormSum_Curr += pow(Fcurr,2);
                    }
                    break;
                }
                case s_Norm::Type::Relative:
                {
                    for(int site=0; site < site_size_loc_all_NS; ++site)
                    {
                        amrex::Real Fcurr = n_curr_in(site) - n_curr_out(site);
                        Norm(site) = fabs(Fcurr/(n_curr_in(site) + n_curr_out(site)));
                        Broyden_NormSum_Curr += pow(Norm(site),2);
                    }
                    break;
                }
                default:
                {
                    amrex::Abort("Norm Type " + Broyden_Norm_Type + " is not yet defined.");
                }
            }
    
            /*find maximum local norm*/
            Broyden_Norm = Norm(0);
            for(int site=1; site < site_size_loc_all_NS; ++site)
            {
                if(Broyden_Norm < Norm(site))
                {
                    Broyden_Norm = Norm(site);
                }
            }
            MPI_Allreduce(MPI_IN_PLACE,
    	        		 &Broyden_Norm,
    	        		  1,
                          MPI_DOUBLE,
    		              MPI_MAX,
                          ParallelDescriptor::Communicator());
    
            MPI_Allreduce(MPI_IN_PLACE,
    	        		 &Broyden_NormSum_Curr,
    	        		  1,
                          MPI_DOUBLE,
    		              MPI_SUM,
                          ParallelDescriptor::Communicator());
            Broyden_NormSum_Curr = sqrt(Broyden_NormSum_Curr);

            /*Evaluate denom = delta_F_curr^T * delta_F_curr */
            amrex::Real Broyden_Denom = 0.;
            for(int site=0; site < site_size_loc_all_NS; ++site)
            {
                amrex::Real Fcurr = n_curr_in(site) - n_curr_out(site);
                delta_F_curr(site) = Fcurr - F_curr(site);
                F_curr(site) = Fcurr;
                Broyden_Denom += pow(delta_F_curr(site),2.);
            }
    
            MPI_Allreduce(MPI_IN_PLACE,
    	        		 &Broyden_Denom,
    	        		  1,
                          MPI_DOUBLE,
    		              MPI_SUM,
                          ParallelDescriptor::Communicator());
    
            //amrex::Print() << "n_curr_in, n_prev_in: " << n_curr_in(0) << " " << n_prev_in(0) << "\n";
            amrex::Print() << "\n Broyden_NormSum_Curr: " << std::setw(20) << Broyden_NormSum_Curr << "\n";
            amrex::Print() <<   " Broyden_NormSum_Prev: " << std::setw(20) << Broyden_NormSum_Prev
                           << ",   Difference: " << (Broyden_NormSum_Curr - Broyden_NormSum_Prev) << "\n";
            amrex::Print() << " Broyden max norm: " << Broyden_Norm << "\n\n";
            //amrex::Print() << "Broyden_denom: " << Broyden_Denom << "\n";

            Broyden_NormSum_Prev = Broyden_NormSum_Curr; 
    
            int m = Broyden_Step-1;
            if(m > 0)
            {
                /*First, evaluate W*(V^T*deltaF), i.e. Wmat*(VmatTran*delta_F_curr)*/
                /*Use intermed_vector to temporarily store vector (VmatTran*delta_F_curr)*/
                for(int iter=1; iter <= m-1; ++iter)
                {
                    amrex::Real sum = 0.;   
                    for(int site=0; site < site_size_loc_all_NS; ++site)
                    {
    	                sum += VmatTran(site,iter) * delta_F_curr(site);  		
    	            }
    		        intermed_vector(iter) = sum;
                }
                /*Allreduce intermed_vector for complete matrix-vector multiplication*/
                MPI_Allreduce(MPI_IN_PLACE,
    	            		 &intermed_vector(0),
    	            		  Broyden_Threshold_MaxStep,
                              MPI_Vector_Type,
    			              Vector_Add,
                              ParallelDescriptor::Communicator());
    
    	        /*Use sum_vector to temporarily store Wmat*intermed_vector */
    
                for(int site=0; site < site_size_loc_all_NS; ++site)
                {
            		amrex::Real sum = 0.;   
                    for(int iter=1; iter <= m-1; ++iter)
                    {
    	                sum += Wmat(iter,site) * intermed_vector(iter);  		
    	            }
                    sum_vector(site) = sum;
    	        }
    
                /*Evaluate Wmat and VmatTran at iteration m*/
                for(int site=0; site < site_size_loc_all_NS; ++site)
                {
                    amrex::Real delta_n = n_curr_in(site) - n_prev_in(site);
    
                    VmatTran(site,m) = delta_F_curr(site)/Broyden_Denom;
    
            		/*Access to (m,site) will be slower*/
                    Wmat(m, site)  = - Broyden_fraction*delta_F_curr(site) 
    		                         + delta_n 
    		                         - sum_vector(site); 
                }
    
                /*Next, evaluate W*(V^T*F_curr), i.e. Wmat*(VmatTran*F_curr)*/
                /*Reuse intermed_vector to temporarily store vector (VmatTran*F_curr)*/
    
                SetVal_RealTable1D(h_intermed_vector_data, 0.);
                for(int iter=1; iter <= m; ++iter)
                {
            		amrex::Real sum = 0.;   
                    for(int site=0; site < site_size_loc_all_NS; ++site)
                    {
    	                sum += VmatTran(site, iter) * F_curr(site);  		
    	            }
    		        intermed_vector(iter) = sum;
                }
    
                /*Allreduce intermed_vector for complete matrix-vector multiplication*/
                MPI_Allreduce(MPI_IN_PLACE,
    			             &intermed_vector(0),
                			  Broyden_Threshold_MaxStep,
                              MPI_Vector_Type,
    			              Vector_Add,
                              ParallelDescriptor::Communicator());
    
    
                //if (ParallelDescriptor::IOProcessor())
                //{
                //    amrex::Print() << "Printing intermed_vector (location 2): \n";
                //    for(int iter=0; iter <= m; ++iter)
                //    {
                //        amrex::Print() << iter << " "<< intermed_vector(iter) << "\n";
                //    }
                //}
    
        	    /*Reuse sum_vector to temporarily store Wmat*intermed_vector */
                SetVal_RealTable1D(h_sum_vector_data, 0.);
    
    
                for(int site=0; site < site_size_loc_all_NS; ++site)
                {
    	        	amrex::Real sum = 0.;   
                    for(int iter=1; iter <= m; ++iter)
                    {
    	                sum += Wmat(iter, site) * intermed_vector(iter);  		
    	            }
                    sum_vector(site) = sum;
        	    }
            } /*end of if(m > 0) */
    
    
    
            /*Store current n in previous n, predict next n and store it in current n*/
            for(int site=0; site < site_size_loc_all_NS; ++site)
            {
                n_prev_in(site) =  n_curr_in(site);
                n_curr_in(site) =  n_prev_in(site) 
    		                     - Broyden_Scalar * Broyden_fraction * F_curr(site) 
                    			 - Broyden_Scalar * sum_vector(site);
            }
            amrex::Print() << "n_new_in: " << n_curr_in(0) << "\n";
    
    
            /*Increment Broyden_Step*/
            Broyden_Step += 1;

    }
  #else
    void 
    c_TransportSolver::Execute_Broyden_Modified_Second_Algorithm_Parallel()
    {

            //amrex::Print() << "Execute_Broyden_Modified_Second_Algorithm_Parallel\n";    
    
            //amrex::Print() << "\nBroydenStep: " << Broyden_Step  
    	    //               << ",  fraction: "   << Broyden_fraction 
            //		       << ",  scalar: " << Broyden_Scalar<< "\n";
    
            auto const& h_n_curr_in       = h_n_curr_in_data.table();
            auto const& h_intermed_vector = h_intermed_vector_data.table();
            auto* h_Intermed_values       = h_Intermed_values_vec.dataPtr();
    
            /*vectors*/
            auto const& n_curr_in      = d_n_curr_in_data.table();
            auto const& n_curr_out     = d_n_curr_out_data.table();
            auto const& n_prev_in      = d_n_prev_in_data.table();
            auto const& F_curr         = d_F_curr_data.table();
            auto const& delta_F_curr   = d_delta_F_curr_data.table();
            auto const& Norm           = d_Norm_data.table();
    
            /*matrices*/
            auto const& VmatTran       = d_VmatTran_data.table();
            auto const& Wmat           = d_Wmat_data.table();
    
            /*temp vectors*/
            auto const& sum_vector     = d_sum_vector_data.table();
            auto const& intermed_vector= d_intermed_vector_data.table();

            auto* Intermed_values      = d_Intermed_values_vec.dataPtr();

            for (auto& v: h_Intermed_values_vec) v = 0.;
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_Intermed_values_vec.begin(), 
                                                         h_Intermed_values_vec.end(), 
                                                       d_Intermed_values_vec.begin() );
            amrex::Gpu::streamSynchronize();

            switch(map_NormType[Broyden_Norm_Type])
            {
                case s_Norm::Type::Absolute:
                {
                    const int BTM = Broyden_Threshold_MaxStep;
                    const int SSL = site_size_loc_all_NS;
                    amrex::ParallelFor(site_size_loc_all_NS, [=] AMREX_GPU_DEVICE (int site) noexcept
                    {
                        Norm(site) = 0.;
                        sum_vector(site) = 0.;
                        for(int iter=site; iter < BTM; iter += SSL) 
                        {
                            intermed_vector(iter) = 0.;
                        }
    
                        amrex::Real Fcurr = n_curr_in(site) - n_curr_out(site);
                        Norm(site) = fabs(Fcurr);
                        amrex::Real norm_sq = pow(Fcurr,2);
    
                        delta_F_curr(site) = Fcurr - F_curr(site);
                        F_curr(site) = Fcurr;
                        amrex::Real deltaF_sq = pow(delta_F_curr(site),2.);
    
                        amrex::Gpu::Atomic::Max(&(Intermed_values[0]), Norm(site));
                        amrex::HostDevice::Atomic::Add(&(Intermed_values[1]), norm_sq);
                        amrex::HostDevice::Atomic::Add(&(Intermed_values[2]), deltaF_sq);
                        /*Evaluate denom = delta_F_curr^T * delta_F_curr */
                    });
                    break;
                }
                case s_Norm::Type::Relative:
                {
                    const int BTM = Broyden_Threshold_MaxStep;
                    const int SSL = site_size_loc_all_NS;
                    amrex::ParallelFor(site_size_loc_all_NS, [=] AMREX_GPU_DEVICE (int site) noexcept
                    {
                        Norm(site) = 0.;
                        sum_vector(site) = 0.;
                        for(int iter=site; iter < BTM; iter += SSL) 
                        {
                            intermed_vector(iter) = 0.;
                        }
    
                        amrex::Real Fcurr = n_curr_in(site) - n_curr_out(site);
                        Norm(site) = fabs(Fcurr/(n_curr_in(site) + n_curr_out(site)));
                        amrex::Real norm_sq = pow(Norm(site),2);
    
                        delta_F_curr(site) = Fcurr - F_curr(site);
                        F_curr(site) = Fcurr;
                        amrex::Real deltaF_sq = pow(delta_F_curr(site),2.);
    
                        amrex::Gpu::Atomic::Max(&(Intermed_values[0]), Norm(site));
                        amrex::HostDevice::Atomic::Add(&(Intermed_values[1]), norm_sq);
                        amrex::HostDevice::Atomic::Add(&(Intermed_values[2]), deltaF_sq);
                        /*Evaluate denom = delta_F_curr^T * delta_F_curr */
                    });
                    break;
                }
                default:
                {
                    amrex::Abort("Norm Type " + Broyden_Norm_Type + " is not yet defined.");
                }
            }
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_Intermed_values_vec.begin(), 
                                                         d_Intermed_values_vec.end(), 
                                                       h_Intermed_values_vec.begin() );
            amrex::Gpu::streamSynchronize();
    
            Broyden_Norm               = h_Intermed_values[0];
            Broyden_NormSum_Curr       = h_Intermed_values[1]; 
            amrex::Real Broyden_Denom  = h_Intermed_values[2];

            MPI_Allreduce(MPI_IN_PLACE,
    	        		 &Broyden_Norm,
    	        		  1,
                          MPI_DOUBLE,
    		              MPI_MAX,
                          ParallelDescriptor::Communicator());
    
            MPI_Allreduce(MPI_IN_PLACE,
    	        		 &Broyden_NormSum_Curr,
    	        		  1,
                          MPI_DOUBLE,
    		              MPI_SUM,
                          ParallelDescriptor::Communicator());
            Broyden_NormSum_Curr = sqrt(Broyden_NormSum_Curr);
    
            MPI_Allreduce(MPI_IN_PLACE,
    	        		 &Broyden_Denom,
    	        		  1,
                          MPI_DOUBLE,
    		              MPI_SUM,
                          ParallelDescriptor::Communicator());
    
            amrex::Print() << "\n Broyden_NormSum_Curr: " << std::setw(20) << Broyden_NormSum_Curr << "\n";
            amrex::Print() <<   " Broyden_NormSum_Prev: " << std::setw(20) << Broyden_NormSum_Prev
                           << ",   Difference: " << (Broyden_NormSum_Curr - Broyden_NormSum_Prev) << "\n";
            amrex::Print() << " Broyden max norm: " << Broyden_Norm << "\n\n";
            //amrex::Print() << "n_curr_in, n_prev_in: " << n_curr_in(0) << " " << n_prev_in(0) << "\n";
            //amrex::Print() << "Broyden_denom: " << Broyden_Denom << "\n";
    
            /*Swap L2 norms*/
            Broyden_NormSum_Prev = Broyden_NormSum_Curr; 
    
            int m = Broyden_Step-1;
            if(m > 0)
            {
                /*First, evaluate W*(V^T*deltaF), i.e. Wmat*(VmatTran*delta_F_curr)*/
                /*Use intermed_vector to temporarily store vector (VmatTran*delta_F_curr)*/
    
                amrex::ParallelFor(site_size_loc_all_NS, [=] AMREX_GPU_DEVICE (int site) noexcept
                {
                    for(int iter=1; iter <= m-1; ++iter)
                    {
                        amrex::Real val = VmatTran(site,iter) * delta_F_curr(site);  		
                        amrex::HostDevice::Atomic::Add(&(intermed_vector(iter)), val);
                    }
                });
                h_intermed_vector_data.copy(d_intermed_vector_data); /*from device to host*/
                amrex::Gpu::streamSynchronize();
    
                /*Allreduce intermed_vector for complete matrix-vector multiplication*/
                MPI_Allreduce(MPI_IN_PLACE,
    	            		 &h_intermed_vector(0),
    	            		  Broyden_Threshold_MaxStep,
                              MPI_Vector_Type,
    			              Vector_Add,
                              ParallelDescriptor::Communicator());
    
                d_intermed_vector_data.copy(h_intermed_vector_data); /*from host to device*/
                amrex::Gpu::streamSynchronize();
    
                const amrex::Real BF = Broyden_fraction;
                const amrex::Real Denom = Broyden_Denom; 
                amrex::ParallelFor(site_size_loc_all_NS, [=] AMREX_GPU_DEVICE (int site) noexcept
                {
    	            /*Use sum_vector to temporarily store Wmat*intermed_vector */
            		amrex::Real sum = 0.;   
                    for(int iter=1; iter <= m-1; ++iter)
                    {
    	                sum += Wmat(iter,site) * intermed_vector(iter);  		
    	            }
                    sum_vector(site) = sum;
    
                    /*Evaluate Wmat and VmatTran at iteration m*/
                    amrex::Real delta_n = n_curr_in(site) - n_prev_in(site);
    
                    VmatTran(site,m) = delta_F_curr(site)/Denom;
    
            		/*Access to (m,site) will be slower*/
                    Wmat(m, site)  = - BF*delta_F_curr(site) 
    		                         + delta_n 
    		                         - sum_vector(site); 
    
                });
                SetVal_RealTable1D(h_intermed_vector_data, 0.);
                d_intermed_vector_data.copy(h_intermed_vector_data); /*from host to device*/
                amrex::Gpu::streamSynchronize();
    
                /*Next, evaluate W*(V^T*F_curr), i.e. Wmat*(VmatTran*F_curr)*/
                /*Reuse intermed_vector to temporarily store vector (VmatTran*F_curr)*/
    
                amrex::ParallelFor(site_size_loc_all_NS, [=] AMREX_GPU_DEVICE (int site) noexcept
                {
                    for(int iter=1; iter <= m; ++iter)
                    {
                        amrex::Real val = VmatTran(site, iter) * F_curr(site);  		
                        amrex::HostDevice::Atomic::Add(&(intermed_vector(iter)), val);
                    }
                });
                h_intermed_vector_data.copy(d_intermed_vector_data); /*from device to host*/
                amrex::Gpu::streamSynchronize();
    
                //amrex::Print() << "Printing itermed_vector (all proc): \n";
                //for(int iter=0; iter <= m; ++iter)
                //{
                //    std::cout << "rank/iter/value: " 
                //              << amrex::ParallelDescriptor::MyProc() << "  " 
                //              << iter << "  " 
                //              << h_intermed_vector(iter) << "\n";
                //}
    
                /*Allreduce intermed_vector for complete matrix-vector multiplication*/
                MPI_Allreduce(MPI_IN_PLACE,
    			             &h_intermed_vector(0),
                			  Broyden_Threshold_MaxStep,
                              MPI_Vector_Type,
    			              Vector_Add,
                              ParallelDescriptor::Communicator());
    
                d_intermed_vector_data.copy(h_intermed_vector_data); /*from host to device*/
                amrex::Gpu::streamSynchronize();
    
                //if (ParallelDescriptor::IOProcessor())
                //{
                //    amrex::Print() << "Printing itermed_vector (location 2): \n";
                //    for(int iter=0; iter <= m; ++iter)
                //    {
                //        amrex::Print() << iter << " "<< h_intermed_vector(iter) << "\n";
                //    }
                //}
    
        	    /*Reuse sum_vector to temporarily store Wmat*intermed_vector */
                amrex::ParallelFor(site_size_loc_all_NS, [=] AMREX_GPU_DEVICE (int site) noexcept
                {
    	        	amrex::Real sum = 0.;   
                    for(int iter=1; iter <= m; ++iter)
                    {
    	                sum += Wmat(iter, site) * intermed_vector(iter);  		
    	            }
                    sum_vector(site) = sum;
                });
                amrex::Gpu::streamSynchronize();
            } /*end of if(m > 0) */
    
            /*Store current n in previous n, predict next n and store it in current n*/
            const amrex::Real BS = Broyden_Scalar;
            const amrex::Real BF = Broyden_fraction;
            const int my_rank = amrex::ParallelDescriptor::MyProc();
            amrex::ParallelFor(site_size_loc_all_NS, [=] AMREX_GPU_DEVICE (int site) noexcept
            {
                n_prev_in(site) =  n_curr_in(site);
                n_curr_in(site) =  n_prev_in(site) 
    		                     - BS * BF * F_curr(site) 
                    			 - BS * sum_vector(site);
            });
            //h_n_curr_in_data.copy(d_n_curr_in_data); /*from device to host*/
            //amrex::Gpu::streamSynchronize();
    
            Broyden_Step += 1;
    }
  #endif
#endif
