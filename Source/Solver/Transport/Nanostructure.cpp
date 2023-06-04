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
                                  const int use_selfconsistent_potential,
                                  const int use_negf)
                 : amrex::ParticleContainer<realPD::NUM, intPD::NUM, 
                                            realPA::NUM, intPA::NUM> (geom, dm, ba)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Nanostructure() Constructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    NSType::name = NS_name_str;
    amrex::Print() << "Nanostructure: " << NS_name_str << "\n";
     
    auto& rCode = c_Code::GetInstance();

    _use_electrostatic            = rCode.use_electrostatic;
    _use_selfconsistent_potential = use_selfconsistent_potential;
    _use_negf                     = use_negf;

    if(_use_negf) 
    {
        NSType::Broyden_fraction = NS_Broyden_frac;

        ReadNanostructureProperties();

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

        Fill_AtomLocations();
        Redistribute(); //This function is in amrex::ParticleContainer

        MarkCellsWithAtoms();


        NSType::Initialize_ChargeAtFieldSites(NS_initial_deposit_value);
        Deposit_AtomAttributeToMesh();
    }

    if(_use_negf) 
    {
        InitializeNEGF();
        pos_vec.clear(); 
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_Nanostructure() Constructor************************\n";
#endif
}

template<typename NSType>
c_Nanostructure<NSType>::~c_Nanostructure ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Nanostructure() Destructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    NSType::Deallocate();

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_Nanostructure() Destructor************************\n";
#endif
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
        for(int i=0; i < NSType::num_atoms; ++i) 
        {
            ParticleType p;
            p.id() = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

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
c_Nanostructure<NSType>::MarkCellsWithAtoms() 
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
            //amrex::Print() << p << "  " << p_par[p].pos(0) << " "<< p_par[p].pos(1) << " " << p_par[p].pos(2) << " " << index[0] << " " << index[1] << " " << index[2] << "\n";
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
            int site_id = get_1D_site_id(global_id); 

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

	    //if(p == 0) {
            //    amrex::Print() << "\np: " << p << "\n";
            //    amrex::Print() << "dx: "   << dx[0] << " " << dx[1] << " " << dx[2] << "\n";
            //    amrex::Print() << "pos: " << p_par[p].pos(0) << " " << p_par[p].pos(1) << " " << p_par[p].pos(2) << "\n";
            //    amrex::Print() << "plo: " << plo[0] << " " << plo[1] << " " << plo[2] << "\n";
            //    amrex::Print() << "l: "   << lx << " " << ly << " " << lz << "\n";
            //    amrex::Print() << "ijk: " << i << " " << j << " " << k << "\n";
            //    amrex::Print() << "w_hi: " << wx_hi << " " << wy_hi << " " << wz_hi << "\n";
            //    amrex::Print() << "w_lo: " << wx_lo << " " << wy_lo << " " << wz_lo << "\n";
	    //}

            p_par_gather[p] = wx_lo*wy_lo*wz_lo*phi(i  , j  , k  , 0)

		            + wx_hi*wy_lo*wz_lo*phi(i+1, j  , k  , 0)
			    + wx_lo*wy_hi*wz_lo*phi(i  , j+1, k  , 0)
			    + wx_lo*wy_lo*wz_hi*phi(i  , j  , k+1, 0)

		            + wx_hi*wy_hi*wz_lo*phi(i+1, j+1, k  , 0)
			    + wx_lo*wy_hi*wz_hi*phi(i  , j+1, k+1, 0)
			    + wx_hi*wy_lo*wz_hi*phi(i+1, j  , k+1, 0)

			    + wx_hi*wy_hi*wz_hi*phi(i+1, j+1, k+1, 0);
	    //if(p == 0) {
	    //    amrex::Print() << "phi(i  , j  , k  , 0): " << phi(i  , j  , k  , 0) << "\n";
	    //    amrex::Print() << "phi(i+1, j  , k  , 0): " << phi(i+1, j  , k  , 0) << "\n";
	    //    amrex::Print() << "phi(i  , j+1, k  , 0): " << phi(i  , j+1, k  , 0)<< "\n";
	    //    amrex::Print() << "phi(i  , j  , k+1, 0): " << phi(i  , j  , k+1, 0) << "\n";
	    //    amrex::Print() << "phi(i+1, j+1, k  , 0) "  << phi(i+1, j+1, k  , 0)<< "\n";
	    //    amrex::Print() << "phi(i  , j+1, k+1, 0): " << phi(i  , j+1, k+1, 0)<< "\n";
	    //    amrex::Print() << "phi(i+1, j  , k+1, 0): " << phi(i+1, j  , k+1, 0) << "\n";
	    //    amrex::Print() << "phi(i+1, j+1, k+1, 0): " << phi(i+1, j+1, k+1, 0) << "\n";
            //    amrex::Print() << "\npar_phi: " <<  p_par_gather[p] << "\n";
	    //}	

        });

    }
    p_mf_gather->setVal(0.);
    p_mf_deposit->FillBoundary(_geom->periodicity());

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
    auto const& n_curr_in = NSType::d_n_curr_in_data.table();

    NSType::d_n_curr_in_data.copy(NSType::h_n_curr_in_data);
    amrex::Gpu::streamSynchronize();
    #else
    auto const& n_curr_in = NSType::h_n_curr_in_data.table();
    #endif

    p_mf_deposit->setVal(0.);
    p_mf_deposit->FillBoundary(_geom->periodicity());

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

        //amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
        //{
        //    amrex::Real lx = (p_par[p].pos(0) - plo[0] - dx[0]*0.5)/dx[0];
	//    amrex::Real ly = (p_par[p].pos(1) - plo[1] - dx[1]*0.5)/dx[1];
	//    amrex::Real lz = (p_par[p].pos(2) - plo[2] - dx[2]*0.5)/dx[2];

	//    int i = static_cast<int>(amrex::Math::floor(lx)); 
	//    int j = static_cast<int>(amrex::Math::floor(ly)); 
	//    int k = static_cast<int>(amrex::Math::floor(lz));

	//    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j  , k  , 0), 0.);
	//    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j  , k  , 0), 0.);
	//    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j+1, k  , 0), 0.);
	//    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j  , k+1, 0), 0.);
	//    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j+1, k  , 0), 0.);
	//    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j+1, k+1, 0), 0.);
	//    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j  , k+1, 0), 0.);
        //    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j+1, k+1, 0), 0.);

        //});

        //#ifdef AMREX_USE_GPU
	//amrex::Gpu::streamSynchronize();
        //#endif

	amrex::Print() << "np: " << np << "\n";
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
        {
	    amrex::Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

            int global_id = p_par[p].id();
            int site_id = get_1D_site_id(global_id); 

            amrex::Real qp = n_curr_in(site_id)*unit_charge/vol/atoms_per_field_site;


            amrex::Real lx = (p_par[p].pos(0) - plo[0] - dx[0]*0.5)/dx[0];
	    amrex::Real ly = (p_par[p].pos(1) - plo[1] - dx[1]*0.5)/dx[1];
	    amrex::Real lz = (p_par[p].pos(2) - plo[2] - dx[2]*0.5)/dx[2];

	    int i = static_cast<int>(amrex::Math::floor(lx)); 
	    int j = static_cast<int>(amrex::Math::floor(ly)); 
	    int k = static_cast<int>(amrex::Math::floor(lz));

            //rho(i,j,k) = qp;

            amrex::Real wx_hi = lx - i;
            amrex::Real wy_hi = ly - j;
            amrex::Real wz_hi = lz - k;

            amrex::Real wx_lo = amrex::Real(1.0) - wx_hi;
            amrex::Real wy_lo = amrex::Real(1.0) - wy_hi;
            amrex::Real wz_lo = amrex::Real(1.0) - wz_hi;

	    //if(p == 0) {
            //    amrex::Print() << "\np: " << p << "\n";
            //    amrex::Print() << "\nn_curr_in: " << n_curr_in(site_id) << "\n";
            //    amrex::Print() << "atoms_per_field_site: " << atoms_per_field_site << "\n";
            //    amrex::Print() << "dx: "   << dx[0] << " " << dx[1] << " " << dx[2] << "\n";
            //    amrex::Print() << "vol: " << vol << "\n";
            //    amrex::Print() << "qp: " << qp << "\n";
            //    amrex::Print() << "pos: " << p_par[p].pos(0) << " " << p_par[p].pos(1) << " " << p_par[p].pos(2) << "\n";
            //    amrex::Print() << "plo: " << plo[0] << " " << plo[1] << " " << plo[2] << "\n";
            //    amrex::Print() << "l: "   << lx << " " << ly << " " << lz << "\n";
            //    amrex::Print() << "ijk: " << i << " " << j << " " << k << "\n";
            //    amrex::Print() << "w_hi: " << wx_hi << " " << wy_hi << " " << wz_hi << "\n";
            //    amrex::Print() << "w_lo: " << wx_lo << " " << wy_lo << " " << wz_lo << "\n";
	    //}
 
	    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j  , k  , 0), wx_lo*wy_lo*wz_lo*qp);

	    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j  , k  , 0), wx_hi*wy_lo*wz_lo*qp);
	    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j+1, k  , 0), wx_lo*wy_hi*wz_lo*qp);
	    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j  , k+1, 0), wx_lo*wy_lo*wz_hi*qp);

	    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j+1, k  , 0), wx_hi*wy_hi*wz_lo*qp);
	    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j+1, k+1, 0), wx_lo*wy_hi*wz_hi*qp);
	    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j  , k+1, 0), wx_hi*wy_lo*wz_hi*qp);

            amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j+1, k+1, 0), wx_hi*wy_hi*wz_hi*qp);
	    //if(p == 0) {
	    //    amrex::Print() << "rho(i  , j  , k  , 0): " << rho(i  , j  , k  , 0) << "\n";
	    //    amrex::Print() << "rho(i+1, j  , k  , 0): " << rho(i+1, j  , k  , 0) << "\n";
	    //    amrex::Print() << "rho(i  , j+1, k  , 0): " << rho(i  , j+1, k  , 0)<< "\n";
	    //    amrex::Print() << "rho(i  , j  , k+1, 0): " << rho(i  , j  , k+1, 0) << "\n";
	    //    amrex::Print() << "rho(i+1, j+1, k  , 0) "  << rho(i+1, j+1, k  , 0)<< "\n";
	    //    amrex::Print() << "rho(i  , j+1, k+1, 0): " << rho(i  , j+1, k+1, 0)<< "\n";
	    //    amrex::Print() << "rho(i+1, j  , k+1, 0): " << rho(i+1, j  , k+1, 0) << "\n";
	    //    amrex::Print() << "rho(i+1, j+1, k+1, 0): " << rho(i+1, j+1, k+1, 0) << "\n";
            //    amrex::Print() << "\nrho_sum: " <<  rho(i  , j  , k  , 0) 
	    //    	                          + rho(i+1, j  , k  , 0) + rho(i  , j+1, k  , 0) + rho(i  , j  , k+1, 0)
	    //    				  + rho(i+1, j+1, k  , 0) + rho(i  , j+1, k+1, 0) + rho(i+1, j  , k+1, 0)
	    //    				  + rho(i+1, j+1, k+1, 0) << "\n";
	    //}	

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

    amrex::Gpu::DeviceVector<amrex::Real> vec_U(num_field_sites);

    for (int l=0; l < num_field_sites; ++l) 
    {
        vec_U[l] = 0.;
    }

    amrex::Real* p_U   = vec_U.dataPtr();  
    
    int lev = 0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto np = pti.numParticles();

        const auto& particles = pti.GetArrayOfStructs();
        const auto p_par = particles().data();

        auto& par_gather  = pti.get_realPA_comp(realPA::gather);
        auto p_par_gather = par_gather.data();
        auto get_1D_site_id = NSType::get_1D_site_id();
          
        if(average_field_flag) 
        {
            if(NSType::avg_type == s_AVG_TYPE::ALL) 
            {
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
                {
                    int global_id = p_par[p].id();
                    int site_id = get_1D_site_id(global_id); 
                    amrex::HostDevice::Atomic::Add(&(p_U[site_id]), p_par_gather[p]);
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
                            amrex::HostDevice::Atomic::Add(&(p_U[site_id]), p_par_gather[p]);
                        }
                    } 
                });
            }
        }
        else 
        { 
            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
            {
                int global_id  = p_par[p].id();
                int site_id    = get_1D_site_id(global_id); 
                p_U[site_id]   = p_par_gather[p];
            });
        } 
    }

    for (int l=0; l<num_field_sites; ++l) 
    {
        ParallelDescriptor::ReduceRealSum(p_U[l]);
    }

    for (int l=0; l<num_field_sites; ++l) 
    {
        NSType::Potential[l]   = -p_U[l] / num_atoms_to_avg_over;
	//amrex::Print() << "site_id, avg_potential: " << l << "  " << NSType::Potential[l] << "\n";
        /*minus because Potential is experienced by electrons*/
    }

}


template<typename NSType>
void 
c_Nanostructure<NSType>::Write_PotentialAtSites() 
{
    if (ParallelDescriptor::IOProcessor()) 
    {
        std::ofstream outfile;
        outfile.open("Potential.dat");
        
        for (int l=0; l<NSType::num_field_sites; ++l)
        {
            outfile << l << std::setw(15) << NSType::PTD[l] << std::setw(15) << NSType::Potential[l] << "\n";
        }  

        outfile.close();
    }
}


template<typename NSType>
void
c_Nanostructure<NSType>:: InitializeNEGF ()
{

    NSType::num_proc = amrex::ParallelDescriptor::NProcs();
    NSType::my_rank = amrex::ParallelDescriptor::MyProc();

    NSType::DefineMatrixPartition();
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

}


template<typename NSType>
void
c_Nanostructure<NSType>:: Solve_NEGF ()
{
    //if(!_use_electrostatic) 
    //{
    //    amrex::Array<amrex::Real,2> QD_loc = {0., 1}; //1nm away in z
    //    for (int l=0; l<NSType::num_field_sites; ++l) 
    //    {
    //        amrex::Real r = sqrt(pow((NSType::PTD[l] - QD_loc[0]),2.) + pow(QD_loc[1],2))*1e-9;
    //        NSType::Potential[l]   = -1*(1./(4.*MathConst::pi*PhysConst::ep0*1.)*(PhysConst::q_e/r));
    //        /*-1 is multiplied because potential is experienced by electrons*/
    //    }
    //}

    NSType::AddPotentialToHamiltonian();
    NSType::Update_ContactPotential(); 
    NSType::Define_EnergyLimits();
    NSType::Update_IntegrationPaths();

    BL_PROFILE_VAR("Compute_RhoInduced", compute_rho_ind);

    NSType::Compute_InducedCharge();

    BL_PROFILE_VAR_STOP(compute_rho_ind);

}


template<typename NSType>
void
c_Nanostructure<NSType>:: Reset ()
{
    if(_use_electrostatic) 
    {
        NSType::Reset_Broyden();
    }
}
