#include "Nanostructure.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/WarpXConst.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"
#include "../../Code.H"
#include "../../Input/GeometryProperties/GeometryProperties.H"
#include "../../Input/MacroscopicProperties/MacroscopicProperties.H"

#include "../../Code_Definitions.H"

#include <AMReX.H>
#include <AMReX_GpuContainers.H>
//
#include <iostream>

template class c_Nanostructure<c_CNT>; 
template class c_Nanostructure<c_Graphene>; 
//template class c_Nanostructure<c_Silicon>; 

template<typename NSType>
c_Nanostructure<NSType>::c_Nanostructure (const amrex::Geometry            & geom,
                                  const amrex::DistributionMapping & dm,
                                  const amrex::BoxArray            & ba,
                                  const std::string NS_name_str,
                                  const std::string NS_gather_str,
                                  const std::string NS_deposit_str,
                                  const amrex::Real NS_initial_deposit_value,
                                  const amrex::Real NS_Broyden_frac,
                                  const std::string NS_Broyden_norm_type,
                                  const int use_negf,
				  const std::string negf_foldername_str
				  )
                 : amrex::ParticleContainer<realPD::NUM, intPD::NUM, 
                                            realPA::NUM, intPA::NUM> (geom, dm, ba)
{
    NSType::name = NS_name_str;
    amrex::Print() << "Nanostructure: " << NS_name_str << "\n";

    NSType::step_foldername_str = negf_foldername_str + "/" + NSType::name; 
    /*eg. output/negf/cnt for nanostructure named cnt */

    if (ParallelDescriptor::IOProcessor())
    {
        CreateDirectory(NSType::step_foldername_str);

        NSType::current_filename_str = NSType::step_foldername_str + "/I.dat";  /*current here means charge current, I */

        NSType::outfile_I.open(NSType::current_filename_str.c_str()); 
    	NSType::outfile_I << "'step', 'Vds' , 'Vgs', ";
        for (int k=0; k <NUM_CONTACTS; ++k)
        {
	    	NSType::outfile_I << ", 'I at contact_" << k+1 << "',";
        }
	    NSType::outfile_I << "'Broyden_Step', 'Max_Iter', 'Broyden_Fraction', 'Broyden_Scalar'" << "\n";
    	NSType::outfile_I.close();
    }

    NSType::step_filename_prefix_str = NSType::step_foldername_str + "/step"; 
    /*eg. output/negf/cnt/step */
    amrex::Print() << "step_filename_prefix_str: " << NSType::step_filename_prefix_str << "\n";
     

    auto& rCode = c_Code::GetInstance();

    _use_electrostatic            = rCode.use_electrostatic;
    _use_negf                     = use_negf;

    NSType::num_proc       = amrex::ParallelDescriptor::NProcs();
    NSType::my_rank        = amrex::ParallelDescriptor::MyProc();
    NSType::initial_charge = NS_initial_deposit_value;

    if(_use_negf) 
    {
        ReadNanostructureProperties();

        NSType::DefineMatrixPartition();

        pos_vec.resize(NSType::num_atoms);
        Read_AtomLocations();
    }

    if(_use_electrostatic) 
    {
        auto& rGprop = rCode.get_GeometryProperties();

        _n_cell = &rGprop.n_cell;
        _geom = &geom;
        auto dx =  _geom->CellSizeArray();
	    NSType::cell_volume = 1.;

	    for (int i=0; i<AMREX_SPACEDIM; ++i) NSType::cell_volume *= dx[i];
	    amrex::Print() << "cell_volume: " << NSType::cell_volume << "\n";

        auto& rMprop = rCode.get_MacroscopicProperties();

        p_mf_gather = rMprop.get_p_mf(NS_gather_str);  
        p_mf_deposit = rMprop.get_p_mf(NS_deposit_str);  

        amrex::Print() << "Fill_AtomLocations()\n";
        Fill_AtomLocations();

        amrex::Print() << "Evaluate_LocalFieldSites()\n";
        Evaluate_LocalFieldSites();

        amrex::Print() << "Define_MPISendCountAndDisp()\n";
        NSType::Define_MPISendCountAndDisp();

        amrex::Print() << "Mark_CellsWithAtoms()\n";
        Mark_CellsWithAtoms();

        amrex::Print() << "Initialize_ChargeAtFieldSites()\n";
        NSType::Initialize_ChargeAtFieldSites();

        amrex::Print() << "Deposit_AtomAttributeToMesh()\n";
        Deposit_AtomAttributeToMesh();
    }

    if(_use_negf) 
    {
        InitializeNEGF();
        pos_vec.clear(); 
    }
    
}


template<typename NSType>
void
c_Nanostructure<NSType>:: ReadNanostructureProperties ()
{
    amrex::ParmParse pp_ns(NSType::name);
    read_filename = "NONE";
    pp_ns.query("read_filename", read_filename);
    amrex::Print() << "##### read_filename: " << read_filename << "\n";

    NSType::ReadNanostructureProperties();
}


template<typename NSType>
void 
c_Nanostructure<NSType>::Fill_AtomLocations() 
{

    if (ParallelDescriptor::IOProcessor()) 
    {
        auto get_1D_site_id = NSType::get_1D_site_id();

        for(int i=0; i < NSType::num_atoms; ++i) 
        {
            ParticleType p;
            p.id() = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

            int site_id = get_1D_site_id(p.id());

            for(int j=0; j < AMREX_SPACEDIM; ++j) 
            {
                p.pos(j)  = pos_vec[i].dir[j];
            } 

            std::array<int,intPA::NUM> int_attribs;
            int_attribs[intPA::cid]  = 0;

            std::array<ParticleReal,realPA::NUM> real_attribs;
            real_attribs[realPA::gather]  = 0.0;
            real_attribs[realPA::deposit]  = 0.0;

            std::pair<int,int> key {0,0}; //{grid_index, tile index}
            int lev=0;
            auto& particle_tile = GetParticles(lev)[key];
          
            particle_tile.push_back(p);
            particle_tile.push_back_int(int_attribs);
            particle_tile.push_back_real(real_attribs);
        }
    }

    Redistribute(); //This function is in amrex::ParticleContainer
    
}


template<typename NSType>
void 
c_Nanostructure<NSType>::Evaluate_LocalFieldSites() 
{
    int lev=0;
    //int num_field_sites = NSType::num_field_sites;
    //#ifdef AMREX_USE_GPU
    //amrex::Gpu::HostVector<amrex::Real> h_intermed_values_vec = {std::numeric_limits<int>::max(), 0};

    //amrex::Gpu::DeviceVector<amrex::Real> d_intermed_values_vec(2);
    //amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_intermed_values_vec.begin(),
    //                                           h_intermed_values_vec.end(),
    //                                           d_intermed_values_vec.begin() );
    //amrex::Gpu::streamSynchronize();

    //auto* intermed_values = d_intermed_values_vec.dataPtr();

    //RealTable1D d_site_id_vec({0}, {num_field_sites}, The_Device_Arena());
    //amrex::ParallelFor(num_field_sites, [=] AMREX_GPU_DEVICE (int s) noexcept 
    //{
    //    d_site_id_vec(s) = 0;
    //});
    //amrex::Gpu::streamSynchronize();
    //#else
    std::unordered_map<int, int> map;
    int counter=0;
    int min_site_id = std::numeric_limits<int>::max();
    //#endif

    int hasValidBox = false;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        hasValidBox = true;
        auto np = pti.numParticles();

        const auto& particles = pti.GetArrayOfStructs();
        const auto p_par = particles().data();
        auto get_1D_site_id = NSType::get_1D_site_id();

        //#ifdef AMREX_USE_GPU
        //amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
        //{
        //    int global_id = p_par[p].id();
        //    int site_id   = get_1D_site_id(global_id);
        //    amrex::Gpu::Atomic::Min(&(intermed_values[0]), site_id);
    
        //    if(d_site_id_vec(site_id) == 0) 
        //    {
        //        //amrex::HostDevice::Atomic::Add(&(d_site_id_vec(site_id)), 1);
        //        amrex::Gpu::Atomic::Add(&(d_site_id_vec(site_id)), site_id);
        //    }
        //});
        //amrex::Gpu::streamSynchronize();

        //#else
        //hashmap solution
        for(int p=0; p < np; ++p) 
        {   
            int global_id = p_par[p].id();
            int site_id   = get_1D_site_id(global_id);

            if(map.find(site_id) == map.end()) 
            {
                map[site_id] = counter++;
                min_site_id = std::min(min_site_id, site_id);
            }
        }
        //#endif
    }
    //#ifdef AMREX_USE_GPU
    //amrex::ParallelFor(num_field_sites, [=] AMREX_GPU_DEVICE (int s) noexcept 
    //{
    //    if(d_site_id_vec(s) == 1) 
    //    {
    //        amrex::HostDevice::Atomic::Add(&(intermed_values[1]), 1);
    //    }
    //});
    //amrex::Gpu::streamSynchronize();

    //amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_intermed_values_vec.begin(),
    //                                           d_intermed_values_vec.end(),
    //                                           h_intermed_values_vec.begin() );
    //amrex::Gpu::streamSynchronize();

    //if(hasValidBox) {
    //    NSType::site_id_offset = h_intermed_values_vec[0];
    //} else {
    //    NSType::site_id_offset = 0;   
    //}
    //NSType::num_local_field_sites = h_intermed_values_vec[1];

    //#else
    //obtained from hashmap solution
    if(hasValidBox) {
        NSType::site_id_offset = min_site_id;
    } else {
        NSType::site_id_offset = 0;   
    }
    NSType::num_local_field_sites = counter;
    //#endif
    

    //std::cout << " process/site_id_offset/num_local_field_sites: " 
    //          << NSType::my_rank << " " 
    //          << NSType::site_id_offset << " "
    //          << NSType::num_local_field_sites << "\n";

    /*define local charge density 1D table, h_n_curr_in_loc_data */
    NSType::h_n_curr_in_loc_data.resize({0},{NSType::num_local_field_sites}, The_Pinned_Arena()); 
    #if AMREX_USE_GPU
    NSType::d_n_curr_in_loc_data.resize({0},{NSType::num_local_field_sites}, The_Arena()); 
    #endif
}

template<typename NSType>
void 
c_Nanostructure<NSType>::Read_AtomLocations() 
{

    if (ParallelDescriptor::IOProcessor()) 
    {

        std::string str_null = "NONE";
        if(strcmp(read_filename.c_str(),str_null.c_str()) == 0)
        {
            NSType::Generate_AtomLocations(pos_vec);
        }
        else 
        {
            std::ifstream infile;
            infile.open(read_filename.c_str());
             
            if(infile.fail())
            {
                amrex::Abort("Failed to read file " + read_filename);
            }
            else
            {
                int filesize=0;
                std::string line; 
                while(infile.peek()!=EOF)
                {
                    std::getline(infile, line);
                    filesize++;
                }
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(filesize == NSType::num_atoms,
                "Number of atoms, " + std::to_string(NSType::num_atoms) 
                 + ", are not equal to the filesize, " + std::to_string(filesize) + " !");

                infile.seekg(0, std::ios_base::beg);
               
                std::string id[2];

                for(int i=0; i < NSType::num_atoms; ++i) 
                {
                    infile >> id[0] >> id[1];

                    for(int j=0; j < AMREX_SPACEDIM; ++j) 
                    {
                        infile >> pos_vec[i].dir[j];
                        pos_vec[i].dir[j] += NSType::offset[j];
                    }
                }
                infile.close();
            }
        } 
    }

}


template<typename NSType>
void 
c_Nanostructure<NSType>::Mark_CellsWithAtoms() 
{
    auto& rCode = c_Code::GetInstance();
    auto& rPost = rCode.get_PostProcessor();
    auto& rMprop = rCode.get_MacroscopicProperties();
    
    const auto& plo = _geom->ProbLoArray();
    const auto dx =_geom->CellSizeArray();

    auto& mf = rMprop.get_mf("atom_locations"); //define in macroscopic.fields_to_define
    
    amrex::Print() << "Marking atom locations\n";
    int lev = 0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) 
    { 
        auto np = pti.numParticles();

        const auto& particles = pti.GetArrayOfStructs();
        const auto p_par = particles().data();

        auto mf_arr = mf.array(pti);
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
        {
            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> 
                   l = { AMREX_D_DECL((p_par[p].pos(0) - plo[0])/dx[0], 
                                      (p_par[p].pos(1) - plo[1])/dx[1], 
                                      (p_par[p].pos(2) - plo[2])/dx[2]) };

            amrex::GpuArray<int,AMREX_SPACEDIM> 
               index = { AMREX_D_DECL(static_cast<int>(amrex::Math::floor(l[0])), 
                                      static_cast<int>(amrex::Math::floor(l[1])), 
                                      static_cast<int>(amrex::Math::floor(l[2]))) };

            mf_arr(index[0],index[1],index[2]) = 1;                            
        });
    }
    mf.FillBoundary(_geom->periodicity());
   
}


template<typename NSType>
void 
c_Nanostructure<NSType>::Gather_MeshAttributeAtAtoms() 
{
    
    const auto& plo = _geom->ProbLoArray();
    const auto dx =_geom->CellSizeArray();

    int lev = 0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) 
    { 
        auto np = pti.numParticles();

        const auto& particles = pti.GetArrayOfStructs();
        const auto p_par = particles().data();

        auto& par_gather  = pti.get_realPA_comp(realPA::gather);
        auto p_par_gather = par_gather.data();

        auto phi = p_mf_gather->array(pti);
        auto get_1D_site_id = NSType::get_1D_site_id();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
        {
            int global_id = p_par[p].id();

	        amrex::Real lx = (p_par[p].pos(0) - plo[0] - dx[0]*0.5)/dx[0];
	        amrex::Real ly = (p_par[p].pos(1) - plo[1] - dx[1]*0.5)/dx[1];
	        amrex::Real lz = (p_par[p].pos(2) - plo[2] - dx[2]*0.5)/dx[2];

            int i = static_cast<int>(amrex::Math::floor(lx));
            int j = static_cast<int>(amrex::Math::floor(ly));
            int k = static_cast<int>(amrex::Math::floor(lz));

            amrex::Real wx_hi = lx - i;
            amrex::Real wy_hi = ly - j;
            amrex::Real wz_hi = lz - k;

            amrex::Real wx_lo = amrex::Real(1.0) - wx_hi;
            amrex::Real wy_lo = amrex::Real(1.0) - wy_hi;
            amrex::Real wz_lo = amrex::Real(1.0) - wz_hi;

            p_par_gather[p] = wx_lo*wy_lo*wz_lo*phi(i  , j  , k  , 0)

        		            + wx_hi*wy_lo*wz_lo*phi(i+1, j  , k  , 0)
		            	    + wx_lo*wy_hi*wz_lo*phi(i  , j+1, k  , 0)
            			    + wx_lo*wy_lo*wz_hi*phi(i  , j  , k+1, 0)

		                    + wx_hi*wy_hi*wz_lo*phi(i+1, j+1, k  , 0)
            			    + wx_lo*wy_hi*wz_hi*phi(i  , j+1, k+1, 0)
            			    + wx_hi*wy_lo*wz_hi*phi(i+1, j  , k+1, 0)

			                + wx_hi*wy_hi*wz_hi*phi(i+1, j+1, k+1, 0);
        });

    }

    Obtain_PotentialAtSites();
}


template<typename NSType>
void 
c_Nanostructure<NSType>::Deposit_AtomAttributeToMesh() 
{
    const auto& plo = _geom->ProbLoArray();
    const auto dx =_geom->CellSizeArray();
    int lev = 0;

    #ifdef AMREX_USE_GPU
    NSType::d_n_curr_in_loc_data.copy(NSType::h_n_curr_in_loc_data);
    auto const& n_curr_in_loc = NSType::d_n_curr_in_loc_data.table();

    amrex::Gpu::streamSynchronize();
    #else
    auto const& n_curr_in_loc = NSType::h_n_curr_in_loc_data.table();
    #endif


    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) 
    { 
        auto np = pti.numParticles();

        const auto& particles = pti.GetArrayOfStructs();
        const auto p_par = particles().data();

        const auto& par_deposit  = pti.get_realPA_comp(realPA::deposit);
        const auto p_par_deposit = par_deposit.data();

        auto rho = p_mf_deposit->array(pti);
        auto get_1D_site_id = NSType::get_1D_site_id();

    	amrex::Real unit_charge = PhysConst::q_e;
        int atoms_per_field_site =  NSType::num_atoms_per_field_site;
        int SIO = NSType::site_id_offset;

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
        {
    	    amrex::Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

            int global_id = p_par[p].id();
            int site_id   = get_1D_site_id(global_id); 
           
            amrex::Real qp = n_curr_in_loc(site_id - SIO) * unit_charge/vol/atoms_per_field_site; 

            amrex::Real lx = (p_par[p].pos(0) - plo[0] - dx[0]*0.5)/dx[0];
    	    amrex::Real ly = (p_par[p].pos(1) - plo[1] - dx[1]*0.5)/dx[1];
	        amrex::Real lz = (p_par[p].pos(2) - plo[2] - dx[2]*0.5)/dx[2];

	        int i = static_cast<int>(amrex::Math::floor(lx)); 
    	    int j = static_cast<int>(amrex::Math::floor(ly)); 
	        int k = static_cast<int>(amrex::Math::floor(lz));

            amrex::Real wx_hi = lx - i;
            amrex::Real wy_hi = ly - j;
            amrex::Real wz_hi = lz - k;

            amrex::Real wx_lo = amrex::Real(1.0) - wx_hi;
            amrex::Real wy_lo = amrex::Real(1.0) - wy_hi;
            amrex::Real wz_lo = amrex::Real(1.0) - wz_hi;

    	    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j  , k  , 0), wx_lo*wy_lo*wz_lo*qp);

	        amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j  , k  , 0), wx_hi*wy_lo*wz_lo*qp);
	        amrex::Gpu::Atomic::AddNoRet(&rho(i  , j+1, k  , 0), wx_lo*wy_hi*wz_lo*qp);
    	    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j  , k+1, 0), wx_lo*wy_lo*wz_hi*qp);

	        amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j+1, k  , 0), wx_hi*wy_hi*wz_lo*qp);
    	    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j+1, k+1, 0), wx_lo*wy_hi*wz_hi*qp);
	        amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j  , k+1, 0), wx_hi*wy_lo*wz_hi*qp);

            amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j+1, k+1, 0), wx_hi*wy_hi*wz_hi*qp);

        });
    }
    p_mf_deposit->SumBoundary(_geom->periodicity());


}


template<typename NSType>
void 
c_Nanostructure<NSType>::Obtain_PotentialAtSites() 
{
    const int num_field_sites          = NSType::num_field_sites;
    const int num_atoms_per_field_site = NSType::num_atoms_per_field_site;
    const int num_atoms_to_avg_over    = NSType::num_atoms_to_avg_over;
    const int average_field_flag       = NSType::average_field_flag;
    const int blkCol_size_loc          = NSType::blkCol_size_loc;

    amrex::Gpu::DeviceVector<amrex::Real> d_vec_V(num_field_sites);
    std::fill(d_vec_V.begin(), d_vec_V.end(), 0);

    amrex::Real* p_V   = d_vec_V.dataPtr();  
    
    int lev = 0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto np = pti.numParticles();

        const auto& particles = pti.GetArrayOfStructs();
        const auto p_par = particles().data();

        auto& par_gather    = pti.get_realPA_comp(realPA::gather);
        auto p_par_gather   = par_gather.data();
        auto get_1D_site_id = NSType::get_1D_site_id();
          
        if(average_field_flag) 
        {
            if(NSType::avg_type == s_AVG_TYPE::ALL) 
            {
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
                {
                    int global_id = p_par[p].id();
                    int site_id = get_1D_site_id(global_id); 

                    amrex::HostDevice::Atomic::Add(&(p_V[site_id]), p_par_gather[p]);
                });
            }
            else if(NSType::avg_type == s_AVG_TYPE::SPECIFIC) 
            {
                auto get_atom_id_at_site = NSType::get_atom_id_at_site();
                #ifdef AMREX_USE_GPU
    	        auto avg_indices_ptr = NSType::gpuvec_avg_indices.dataPtr();
                #else
    	        auto avg_indices_ptr = NSType::vec_avg_indices.dataPtr();
                #endif

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
                {
                    int global_id = p_par[p].id();
                    int site_id = get_1D_site_id(global_id); 

                    int atom_id_at_site = get_atom_id_at_site(global_id); 
                    int remainder = atom_id_at_site%num_atoms_per_field_site;

	                for(int m=0; m < num_atoms_to_avg_over; ++m)
	                {
                        if(remainder == avg_indices_ptr[m]) 
	                    {
                            amrex::HostDevice::Atomic::Add(&(p_V[site_id]), p_par_gather[p]);
                        }
                    } 
                });
            }
        }
        else 
        { 
            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
            {
                int global_id      = p_par[p].id();
                int site_id        = get_1D_site_id(global_id); 
                p_V[site_id]       = p_par_gather[p];
            });
        } 
    }

    for (int l=0; l<num_field_sites; ++l)
    {
        ParallelDescriptor::ReduceRealSum(p_V[l]);
    }

    //MPI_Allreduce( MPI_IN_PLACE,
    //               &(p_V[0]),
    //               num_field_sites,
    //               MPI_DOUBLE,
    //               MPI_SUM,
    //               ParallelDescriptor::Communicator());

    //for(int i=0; i < num_field_sites; ++i)
    //{
    //    amrex::Print() << "proc/i/h_vec_V: " << NSType::my_rank << " " << i << "  " << p_V[i] << "\n";
    //}
    auto const& h_U_loc = NSType::h_U_loc_data.table();
    for (int l=0; l < blkCol_size_loc; ++l) 
    {
        int gid = NSType::vec_blkCol_gids[l];

        h_U_loc(l) = -p_V[gid] / num_atoms_to_avg_over;
    }
    for(int c=0; c < NUM_CONTACTS; ++c)
    {
        NSType::U_contact[c] = -p_V[NSType::global_contact_index[c]] / num_atoms_to_avg_over;
    }

    //for (int l=0; l < blkCol_size_loc; ++l) 
    //{
    //    int gid = vec_blkCol_gids(l);
    //    if(h_U_loc(i) - (-1*p_V[gid]/num_atoms_to_avg_over) > 1e-8 )
    //    {
    //        amrex::Abort("h_U_loc != p_V: "
    //                + std::to_string(i) + " " + std::to_string(h_U_loc(i)) + " "
    //                + std::to_string(p_V[gid]));
    //    }
    //}
}


//template<typename NSType>
//void 
//c_Nanostructure<NSType>::Obtain_PotentialAtSites() 
//{
//    const int num_field_sites          = NSType::num_field_sites;
//    const int num_atoms_per_field_site = NSType::num_atoms_per_field_site;
//    const int num_atoms_to_avg_over    = NSType::num_atoms_to_avg_over;
//    const int average_field_flag       = NSType::average_field_flag;
//
//    amrex::Gpu::DeviceVector<amrex::Real> vec_V(num_field_sites);
//
//    for (int l=0; l < num_field_sites; ++l) 
//    {
//        vec_V[l] = 0.;
//    }
//
//    amrex::Real* p_V   = vec_V.dataPtr();  
//    
//    int lev = 0;
//    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
//    {
//        auto np = pti.numParticles();
//
//        const auto& particles = pti.GetArrayOfStructs();
//        const auto p_par = particles().data();
//
//        auto& par_gather  = pti.get_realPA_comp(realPA::gather);
//        auto p_par_gather = par_gather.data();
//        auto get_1D_site_id = NSType::get_1D_site_id();
//          
//        if(average_field_flag) 
//        {
//            if(NSType::avg_type == s_AVG_TYPE::ALL) 
//            {
//                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
//                {
//                    int global_id = p_par[p].id();
//                    int site_id = get_1D_site_id(global_id); 
//
//                    amrex::HostDevice::Atomic::Add(&(p_V[site_id]), p_par_gather[p]);
//                });
//            }
//            else if(NSType::avg_type == s_AVG_TYPE::SPECIFIC) 
//            {
//                auto get_atom_id_at_site = NSType::get_atom_id_at_site();
//                #ifdef AMREX_USE_GPU
//    	        auto avg_indices_ptr = NSType::gpuvec_avg_indices.dataPtr();
//                #else
//    	        auto avg_indices_ptr = NSType::vec_avg_indices.dataPtr();
//                #endif
//
//                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
//                {
//                    int global_id = p_par[p].id();
//                    int site_id = get_1D_site_id(global_id); 
//
//                    int atom_id_at_site = get_atom_id_at_site(global_id); 
//                    int remainder = atom_id_at_site%num_atoms_per_field_site;
//
//	                for(int m=0; m < num_atoms_to_avg_over; ++m)
//	                {
//                        if(remainder == avg_indices_ptr[m]) 
//	                    {
//                            amrex::HostDevice::Atomic::Add(&(p_V[site_id]), p_par_gather[p]);
//                        }
//                    } 
//                });
//            }
//        }
//        else 
//        { 
//            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
//            {
//                int global_id  = p_par[p].id();
//                int site_id    = get_1D_site_id(global_id); 
//                p_V[site_id]   = p_par_gather[p];
//            });
//        } 
//    }
//
//    for (int l=0; l<num_field_sites; ++l) 
//    {
//        ParallelDescriptor::ReduceRealSum(p_V[l]);
//    }
//
//    for (int l=0; l<num_field_sites; ++l) 
//    {
//        NSType::Potential[l]   = -p_V[l] / num_atoms_to_avg_over;
//	    //amrex::Print() << "site_id, avg_potential: " << l << "  " << NSType::Potential[l] << "\n";
//        /*minus because Potential is experienced by electrons*/
//    }
//
//}


template<typename NSType>
void
c_Nanostructure<NSType>:: InitializeNEGF ()
{

    NSType::AllocateArrays();

    NSType::ConstructHamiltonian();

    NSType::Define_ContactInfo();

    NSType::Define_EnergyLimits();

    NSType::Define_IntegrationPaths();

    BL_PROFILE_VAR("Compute_DOS", compute_dos);

    NSType::Compute_DensityOfStates();

    BL_PROFILE_VAR_STOP(compute_dos);


    BL_PROFILE_VAR("Compute_Rho0", compute_rho0);

    NSType::Compute_Rho0();

    BL_PROFILE_VAR_STOP(compute_rho0);

    if(!_use_electrostatic)
    {
	    NSType::Define_PotentialProfile();
    }

}


template<typename NSType>
void
c_Nanostructure<NSType>:: Solve_NEGF ()
{

    BL_PROFILE_VAR("Other", compute_other);

    NSType::AddPotentialToHamiltonian();
    NSType::Update_ContactElectrochemicalPotential(); 
    NSType::Define_EnergyLimits();
    NSType::Update_IntegrationPaths();

    BL_PROFILE_VAR_STOP(compute_other);

    BL_PROFILE_VAR("Compute_RhoInduced", compute_rho_ind);

    NSType::Compute_InducedCharge();

    BL_PROFILE_VAR_STOP(compute_rho_ind);

}


template<typename NSType>
void
c_Nanostructure<NSType>:: Write_Data (const std::string filename_prefix, 
                                      RealTable1D& n_curr_out_data, 
                                      RealTable1D& Norm_data)
{
    BL_PROFILE_VAR("Write_Data", compute_write_data);

    NSType::Write_PotentialAtSites(filename_prefix);
    if (ParallelDescriptor::IOProcessor())
    {
        NSType::Write_InducedCharge(filename_prefix, n_curr_out_data);
        NSType::Write_ChargeNorm(filename_prefix, Norm_data);
    }
    BL_PROFILE_VAR_STOP(compute_write_data);
}
