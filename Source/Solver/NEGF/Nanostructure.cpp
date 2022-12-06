#include "Nanostructure.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"

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
                }

                std::array<int,intPA::NUM> int_attribs;
                int_attribs[intPA::cid]  = 0.0;
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

}
