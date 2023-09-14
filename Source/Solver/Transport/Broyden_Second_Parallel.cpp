#include "Transport.H"

#include <limits>

using namespace amrex;

#ifdef BROYDEN_PARALLEL
void 
c_TransportSolver::Execute_Broyden_Modified_Second_Algorithm_Parallel()
{

        amrex::Print() << "Execute_Broyden_Modified_Second_Algorithm_Parallel\n";    

        amrex::Print() << "\nBroydenStep: " << Broyden_Step  
		               << ",  fraction: "   << Broyden_fraction 
        		       << ",  scalar: " << Broyden_Scalar<< "\n";

        auto const& h_n_curr_in       = h_n_curr_in_data.table();
        auto const& n_curr_in_glo     = n_curr_in_glo_data.table();
        auto const& h_intermed_vector = h_intermed_vector_data.table();
        auto* h_Intermed_values       = h_Intermed_values_vec.dataPtr();

        #ifdef BROYDEN_SKIP_GPU_OPTIMIZATION
        auto const& n_curr_in      = h_n_curr_in_data.table();
        auto const& n_curr_out     = h_n_curr_out_data.table();
        auto const& n_prev_in      = h_n_prev_in_data.table();
        auto const& F_curr         = h_F_curr_data.table();
        auto const& delta_F_curr   = h_delta_F_curr_data.table();
        auto const& Norm           = h_Norm_data.table();

        auto const& VmatTran       = h_VmatTran_data.table();
        auto const& Wmat           = h_Wmat_data.table();

        auto const& sum_vector      = h_sum_vector_data.table();
        auto const& intermed_vector = h_intermed_vector_data.table();
        auto* Intermed_values       = h_Intermed_values_vec.dataPtr();
        #else
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
        #endif

        for (auto& v: h_Intermed_values_vec) v = 0.;
        #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_Intermed_values_vec.begin(), 
                                                     h_Intermed_values_vec.end(), 
                                                   d_Intermed_values_vec.begin() );
        amrex::Gpu::streamSynchronize();
        #endif

        const int BTM = Broyden_Threshold_MaxStep;
        const int SSL = site_size_loc;
        switch(map_NormType[Broyden_Norm_Type])
        {
            case s_Norm::Type::Absolute:
            {
                amrex::ParallelFor(site_size_loc, [=] AMREX_GPU_DEVICE (int site) noexcept
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
                amrex::ParallelFor(site_size_loc, [=] AMREX_GPU_DEVICE (int site) noexcept
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
        #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_Intermed_values_vec.begin(), 
                                                     d_Intermed_values_vec.end(), 
                                                   h_Intermed_values_vec.begin() );
        amrex::Gpu::streamSynchronize();
        #endif

        Broyden_Norm         = h_Intermed_values[0];
        Broyden_NormSum_Curr = h_Intermed_values[1]; 
        Broyden_Denom        = h_Intermed_values[2];

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

        amrex::Print() << "\nBroyden_NormSum_Curr: " << std::setw(20) << Broyden_NormSum_Curr << "\n";
        amrex::Print() <<   "Broyden_NormSum_Prev: " << std::setw(20) << Broyden_NormSum_Prev
                       << ",   Difference: " << (Broyden_NormSum_Curr - Broyden_NormSum_Prev) << "\n";
        amrex::Print() << "Broyden max norm: " << Broyden_Norm << "\n\n";
        amrex::Print() << "n_curr_in, n_prev_in: " << n_curr_in(0) << " " << n_prev_in(0) << "\n";
        amrex::Print() << "Broyden_denom: " << Broyden_Denom << "\n";

        /*Swap L2 norms*/
        Broyden_NormSum_Prev = Broyden_NormSum_Curr; 

        int m = Broyden_Step-1;
        if(m > 0)
        {
            /*First, evaluate W*(V^T*deltaF), i.e. Wmat*(VmatTran*delta_F_curr)*/
            /*Use intermed_vector to temporarily store vector (VmatTran*delta_F_curr)*/

            amrex::ParallelFor(site_size_loc, [=] AMREX_GPU_DEVICE (int site) noexcept
            {
                for(int iter=1; iter <= m-1; ++iter)
                {
                    amrex::Real val = VmatTran(site,iter) * delta_F_curr(site);  		
                    amrex::HostDevice::Atomic::Add(&(intermed_vector(iter)), val);
                }
            });
            #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
            h_intermed_vector_data.copy(d_intermed_vector_data); /*from device to host*/
            amrex::Gpu::streamSynchronize();
            #endif

            /*Allreduce intermed_vector for complete matrix-vector multiplication*/
            MPI_Allreduce(MPI_IN_PLACE,
	            		 &h_intermed_vector(0),
	            		  Broyden_Threshold_MaxStep,
                          MPI_Vector_Type,
			              Vector_Add,
                          ParallelDescriptor::Communicator());

            #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
            d_intermed_vector_data.copy(h_intermed_vector_data); /*from host to device*/
            amrex::Gpu::streamSynchronize();
            #endif

            const amrex::Real BF = Broyden_fraction;
            const amrex::Real Denom = Broyden_Denom; 
            amrex::ParallelFor(site_size_loc, [=] AMREX_GPU_DEVICE (int site) noexcept
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
            #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
            d_intermed_vector_data.copy(h_intermed_vector_data); /*from host to device*/
            amrex::Gpu::streamSynchronize();
            #endif

            /*Next, evaluate W*(V^T*F_curr), i.e. Wmat*(VmatTran*F_curr)*/
            /*Reuse intermed_vector to temporarily store vector (VmatTran*F_curr)*/

            amrex::ParallelFor(site_size_loc, [=] AMREX_GPU_DEVICE (int site) noexcept
            {
                for(int iter=1; iter <= m; ++iter)
                {
                    amrex::Real val = VmatTran(site, iter) * F_curr(site);  		
                    amrex::HostDevice::Atomic::Add(&(intermed_vector(iter)), val);
                }
            });
            #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
            h_intermed_vector_data.copy(d_intermed_vector_data); /*from device to host*/
            amrex::Gpu::streamSynchronize();
            #endif

            amrex::Print() << "Printing itermed_vector (all proc): \n";
            for(int iter=0; iter <= m; ++iter)
            {
                std::cout << "rank/iter/value: " 
                          << amrex::ParallelDescriptor::MyProc() << "  " 
                          << iter << "  " 
                          << h_intermed_vector(iter) << "\n";
            }

            /*Allreduce intermed_vector for complete matrix-vector multiplication*/
            MPI_Allreduce(MPI_IN_PLACE,
			             &h_intermed_vector(0),
            			  Broyden_Threshold_MaxStep,
                          MPI_Vector_Type,
			              Vector_Add,
                          ParallelDescriptor::Communicator());

            #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
            d_intermed_vector_data.copy(h_intermed_vector_data); /*from host to device*/
            amrex::Gpu::streamSynchronize();
            #endif

            if (ParallelDescriptor::IOProcessor())
            {
                amrex::Print() << "Printing itermed_vector (location 2): \n";
                for(int iter=0; iter <= m; ++iter)
                {
                    amrex::Print() << iter << " "<< h_intermed_vector(iter) << "\n";
                }
            }

    	    /*Reuse sum_vector to temporarily store Wmat*intermed_vector */
            amrex::ParallelFor(site_size_loc, [=] AMREX_GPU_DEVICE (int site) noexcept
            {
	        	amrex::Real sum = 0.;   
                for(int iter=1; iter <= m; ++iter)
                {
	                sum += Wmat(iter, site) * intermed_vector(iter);  		
	            }
                sum_vector(site) = sum;
            });
            #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
            amrex::Gpu::streamSynchronize();
            #endif
        } /*end of if(m > 0) */

        /*Store current n in previous n, predict next n and store it in current n*/
        const amrex::Real BS = Broyden_Scalar;
        const amrex::Real BF = Broyden_fraction;
        amrex::ParallelFor(site_size_loc, [=] AMREX_GPU_DEVICE (int site) noexcept
        {
            n_prev_in(site) =  n_curr_in(site);
            n_curr_in(site) =  n_prev_in(site) 
		                     - BS * BF * F_curr(site) 
                			 - BS * sum_vector(site);
        });
        #ifndef BROYDEN_SKIP_GPU_OPTIMIZATION
        h_n_curr_in_data.copy(d_n_curr_in_data); /*from device to host*/
        amrex::Gpu::streamSynchronize();
        #endif
        amrex::Print() << "n_new_in: " << h_n_curr_in(0) << "\n";

        Broyden_Step += 1;

        MPI_Allgatherv(&h_n_curr_in(0),
                        site_size_loc,
                        MPI_DOUBLE,
                       &n_curr_in_glo(0),
                        MPI_recv_count.data(),
                        MPI_disp.data(),
                        MPI_DOUBLE,
                        ParallelDescriptor::Communicator());
}
#endif
