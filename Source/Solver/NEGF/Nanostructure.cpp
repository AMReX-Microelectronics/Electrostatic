#include "Nanostructure.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"
#include "../../Code.H"
#include "../../Input/GeometryProperties/GeometryProperties.H"
#include "../../Input/MacroscopicProperties/MacroscopicProperties.H"

#include "../../Code_Definitions.H"

#include <iostream>

template class c_Nanostructure<c_CNT>; 
template class c_Nanostructure<c_Graphene>; 
template class c_Nanostructure<c_Silicon>; 

template<typename NSType>
c_Nanostructure<NSType>::c_Nanostructure (const amrex::Geometry            & geom,
                                  const amrex::DistributionMapping & dm,
                                  const amrex::BoxArray            & ba,
                                  const std::string ns_name)
                 : amrex::ParticleContainer<realPD::NUM, intPD::NUM, 
                                            realPA::NUM, intPA::NUM> (geom, dm, ba)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_Nanostructure() Constructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    NSType::name = ns_name;

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    _n_cell = &rGprop.n_cell;
    _geom = &geom;

    ReadNanostructureProperties();
    ReadAtomLocations();
  
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
            const auto dx =_geom->CellSizeArray();
            amrex::Print() << "dx: " << dx[0] << "  " 
                                     << dx[1] << "  " 
                                     << dx[2] << "\n";

            amrex::Print() << "_n_cell: " << (*_n_cell)[0] << "  " 
                                          << (*_n_cell)[1] << "  " 
                                          << (*_n_cell)[2] << "\n";

            for(int i=0; i < NSType::num_atoms; ++i) 
            {
                ParticleType p;
                p.id() = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                
                infile >> id[0] >> id[1];
                //if(p.id() < 20) {  
                //    amrex::Print() << "IDs: " << p.id() << "  " << id[0] << "  " << id[1] << "\n";
                //    NSType::Print_IDs(p.id());
                //} 
                  
                for(int j=0; j < AMREX_SPACEDIM; ++j) {
                    infile >> p.pos(j);
                    p.pos(j) += NSType::offset[j];
                }
                
                amrex::GpuArray<int,AMREX_SPACEDIM> 
                   index_arr = { AMREX_D_DECL(int(p.pos(0)/dx[0]), 
                                              int(p.pos(1)/dx[1]), 
                                              int(p.pos(2)/dx[2])) };
//
//                int cell_id = Compute_Cell_ID<int>(index_arr, *_n_cell); 
//                if(p.id() < 20) {  
//                    amrex::Print() << "index_arr: " << index_arr[0] << "  " 
//                                                    << index_arr[1] << "  " 
//                                                    << index_arr[2] << "  "
//                                   << "cell_id: "   << cell_id << "\n";
//                } 

                std::array<int,intPA::NUM> int_attribs;
                int_attribs[intPA::cid]  = 0;

                std::array<ParticleReal,realPA::NUM> real_attribs;
                real_attribs[realPA::phi]  = 0.0;
                real_attribs[realPA::charge]  = 0.0;

                std::pair<int,int> key {0,0}; //{grid_index, tile index}
                int lev=0;
                auto& particle_tile = GetParticles(lev)[key];

                particle_tile.push_back(p);
                particle_tile.push_back_int(int_attribs);
                particle_tile.push_back_real(real_attribs);
            }
            infile.close();
            //for(int i=0; i < 3; ++i) 
            //{
            //     amrex::Print() << p_atom_data[i].id[0] << "  "
            //                    << p_atom_data[i].id[1] << "  ";

            //     for(int j=0; j < AMREX_SPACEDIM; ++j) {
            //         amrex::Print() << p_atom_data[i].pos[j] << "  ";
            //     }
            //     amrex::Print() << "\n";
            //}
            
        } 
    }
    Redistribute();
    MarkCellsWithAtoms();
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

    auto& mf = rMprop.get_mf("atom_locations");  
    
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
c_Nanostructure<NSType>::ObtainFieldAtAtomLocations() 
{
    auto& rCode = c_Code::GetInstance();
    auto& rPost = rCode.get_PostProcessor();
    auto& rMprop = rCode.get_MacroscopicProperties();
    
    const auto& plo = _geom->ProbLoArray();
    const auto dx =_geom->CellSizeArray();

    auto& phi = rMprop.get_mf("phi");  
    
    int lev = 0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) 
    { 
        auto np = pti.numParticles();

        const auto& particles = pti.GetArrayOfStructs();
        const auto p_par = particles().data();

        const auto& par_phi = pti.get_realPA_comp(realPA::phi);
        const auto p_par_phi = par_phi.data();

        auto phi_arr = phi.array(pti);

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

            par_phi[p] = phi_arr(index[0],index[1],index[2]);
        });
    }
}


//template<typename NSType>
//void 
//c_Nanostructure<NSType>::FieldGather() 
//{
////    auto& rCode = c_Code::GetInstance();
////    auto& rPost = rCode.get_PostProcessor();
////    auto& rMprop = rCode.get_MacroscopicProperties();
////    auto& atom_locations = rMprop.get_mf("atom_locations");  
////
////    int lev = 0;
////    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) 
////    { 
////        auto np = pti.numParticles();
////        const auto& par_phi = pti.get_realPA_comp(realPA::phi);
////        const auto& particles = pti.GetArrayOfStructs();
////        
////        amrex::Real 
////        for (int i = 0; i < pti.numParticles; ++i) {
////            // do stuff with your SoA data...
////        }
////    } 
//}
