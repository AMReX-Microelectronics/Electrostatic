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
        ReadNanostructureProperties();
        InitializeNEGF();
    } 

    if(_use_selfconsistent_potential) 
    {
        auto& rGprop = rCode.get_GeometryProperties();

        _n_cell = &rGprop.n_cell;
        _geom = &geom;

        auto& rMprop = rCode.get_MacroscopicProperties();

        p_mf_gather = rMprop.get_p_mf(NS_gather_str);  
        p_mf_deposit = rMprop.get_p_mf(NS_deposit_str);  

        ReadAtomLocations();
        Redistribute(); //This function is in amrex::ParticleContainer

        MarkCellsWithAtoms();
        Initialize_AttributeToDeposit(NS_initial_deposit_value);
        Deposit_AtomAttributeToMesh();
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

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_Nanostructure() Destructor************************\n";
#endif
}


template<typename NSType>
void
c_Nanostructure<NSType>:: ReadNanostructureProperties ()
{
    amrex::ParmParse pp_ns(NSType::name);
    read_filename = "__NONE__";
    pp_ns.query("read_filename", read_filename);
    amrex::Print() << "##### read_filename: " << read_filename << "\n";

    NSType::ReadNanostructureProperties();
}


template<typename NSType>
void 
c_Nanostructure<NSType>::ReadAtomLocations() 
{

    if (ParallelDescriptor::IOProcessor()) 
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
                ParticleType p;
                p.id() = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                
                infile >> id[0] >> id[1];

        		//auto get_1D_site_id = NSType::get_1D_site_id();
                //int site_id         = get_1D_site_id(p.id());

	        	//auto get_atom_id_at_site = NSType::get_atom_id_at_site();
                //int atom_id_at_site      = get_atom_id_at_site(p.id());

                for(int j=0; j < AMREX_SPACEDIM; ++j) {
                    infile >> p.pos(j);
                    p.pos(j) += NSType::offset[j];
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
            infile.close();
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

            mf_arr(index[0],index[1],index[2]) = -1;                            
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

        auto mesh_mf_arr = p_mf_gather->array(pti);

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

            p_par_gather[p] = mesh_mf_arr(index[0],index[1],index[2]);
        });
    }
}


template<typename NSType>
void 
c_Nanostructure<NSType>::Deposit_AtomAttributeToMesh() 
{
    const auto& plo = _geom->ProbLoArray();
    const auto dx =_geom->CellSizeArray();

    int lev = 0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) 
    { 
        auto np = pti.numParticles();

        const auto& particles = pti.GetArrayOfStructs();
        const auto p_par = particles().data();

        const auto& par_deposit  = pti.get_realPA_comp(realPA::deposit);
        const auto p_par_deposit = par_deposit.data();

        auto mesh_mf_arr = p_mf_deposit->array(pti);

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

            mesh_mf_arr(index[0],index[1],index[2]) = p_par_deposit[p];
        });
    }
}


template<typename NSType>
void 
c_Nanostructure<NSType>::Initialize_AttributeToDeposit(const amrex::Real value) 
{
    const auto& plo = _geom->ProbLoArray();
    const auto dx =_geom->CellSizeArray();

    int lev = 0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) 
    { 
        auto np = pti.numParticles();

        auto& par_deposit  = pti.get_realPA_comp(realPA::deposit);
        auto p_par_deposit = par_deposit.data();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int p) noexcept 
        {
            p_par_deposit[p] = value;
        });
    }
}


template<typename NSType>
void 
c_Nanostructure<NSType>::Obtain_PotentialAtSites() 
{
    const int num_field_sites          = NSType::num_field_sites;
    const int num_atoms_per_field_site = NSType::num_atoms_per_field_site;
    const int num_atoms_to_avg_over    = NSType::num_atoms_to_avg_over;
    const int average_field_flag       = NSType::average_field_flag;
    const int PTD_id                   = NSType::primary_transport_dir;

    amrex::Gpu::DeviceVector<amrex::Real> vec_U(num_field_sites);
    amrex::Gpu::DeviceVector<amrex::Real> vec_PTD(num_field_sites);

    for (int l=0; l < num_field_sites; ++l) 
    {
        vec_U[l]   = 0.;
        vec_PTD[l] = 0.;
    }

    amrex::Real* p_U   = vec_U.dataPtr();  
    amrex::Real* p_PTD = vec_PTD.dataPtr();  
    
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
                    amrex::HostDevice::Atomic::Add(&(p_PTD[site_id]), p_par[p].pos(PTD_id));
                });
            }
            else if(NSType::avg_type == s_AVG_TYPE::SPECIFIC) 
            {
                auto get_atom_id_at_site = NSType::get_atom_id_at_site();
    	        auto avg_indices_ptr = gpuvec_avg_indices.dataPtr();

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
                            amrex::HostDevice::Atomic::Add(&(p_PTD[site_id]), p_par[p].pos(PTD_id));
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
                p_PTD[site_id] = p_par[p].pos(PTD_id);
            });
        } 
    }

    for (int l=0; l<num_field_sites; ++l) 
    {
        ParallelDescriptor::ReduceRealSum(p_U[l]);
        ParallelDescriptor::ReduceRealSum(p_PTD[l]);
    }

    for (int l=0; l<num_field_sites; ++l) 
    {
        NSType::Potential[l]   = p_U[l]   / num_atoms_to_avg_over;
        NSType::PTD[l] = p_PTD[l] / num_atoms_to_avg_over;
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

    if(_use_electrostatic) 
    {
        NSType::AddPotentialToHamiltonian();
        NSType::Update_ContactPotential(); 
    }

    NSType::Define_ContactInfo();
    NSType::Define_EnergyLimits();
    NSType::Define_IntegrationPaths();
    NSType::Compute_DensityOfStates();
    NSType::Compute_Rho0();

}
