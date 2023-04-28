#include "CNT.H"

#include "../../Utils/CodeUtils/CodeUtil.H"

#include <AMReX_Particles.H>

#include <cmath>
#include <math.h>
#include<stdlib.h>

amrex::Array<int,2> c_CNT::type_id; //very important

/*Next, member class functions are defined*/

void
c_CNT:: ReadNanostructureProperties ()
{
    amrex::Print() << "\n##### NANOSTRUCTURE PROPERTIES #####\n\n";


    amrex::ParmParse pp_ns(name);

    amrex::Print() << "##### Properties Specific to CNT: \n";

    getWithParser(pp_ns,"acc", acc);
    amrex::Print() << "##### acc: " << acc << "\n";

    amrex::Vector<int> vec_type_id;
    getArrWithParser(pp_ns, "type_id", vec_type_id, 0, 2);
    type_id = vecToArr_Templated<int, 2>(vec_type_id);
    amrex::Print() << "##### type_id: ";
    for (int i=0; i<2; ++i) amrex::Print() << type_id[i] << "  ";
    amrex::Print() << "\n";

    //getWithParser(pp_ns,"num_unitcells", num_unitcells);
    //amrex::Print() << "##### num_unitcells: " << num_unitcells << "\n";

    getWithParser(pp_ns,"gamma", gamma); 
    amrex::Print() << "##### gamma: " << gamma << "\n";

    //amrex::Vector<amrex::Real> vec_offset;
    //getArrWithParser(pp_ns, "offset", vec_offset, 0, AMREX_SPACEDIM);
    //offset = vecToArr(vec_offset);

    //amrex::Print() << "##### offset: ";
    //for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << offset[i] << "  ";
    //amrex::Print() << "\n";

    if(type_id[1] == 0) {
        rings_per_unitcell = 4;
        atoms_per_ring = type_id[0];
    } else {
        //need to verify this 
        rings_per_unitcell = 2;
        atoms_per_ring = type_id[0]-1;
    }
    Define_SortedModeVector();
    amrex::Print() << "#####* rings_per_unitcell: " << rings_per_unitcell << "\n";
    amrex::Print() << "#####* atoms_per_ring: " << atoms_per_ring << "\n";

    c_NEGF_Common<BlkType>::ReadNanostructureProperties();

}


void 
c_CNT::set_material_specific_parameters()
{
    num_atoms                = num_unitcells*rings_per_unitcell*atoms_per_ring;
    num_atoms_per_unitcell   = atoms_per_ring*rings_per_unitcell;

    num_field_sites          = num_unitcells*rings_per_unitcell;
    average_field_flag       = 1;
    if(average_field_flag) 
    {
        num_atoms_per_field_site = atoms_per_ring;
        if(avg_type == s_AVG_TYPE::ALL) 
        {
            num_atoms_to_avg_over    = num_atoms_per_field_site;
        }
        else if(avg_type == s_AVG_TYPE::SPECIFIC)
        {
            num_atoms_to_avg_over    = vec_avg_indices.size();
        }
    }
    primary_transport_dir    = 1; /*Y*/

    block_degen_vec.resize(NUM_MODES);
    for(int m=0; m<NUM_MODES; ++m) 
    {
         block_degen_vec[m] = mode_degen_vec[m];
    }

    #if AMREX_USE_GPU
    block_degen_gpuvec.resize(NUM_MODES);

    amrex::Gpu::copy(amrex::Gpu::hostToDevice, 
                     block_degen_vec.begin(), block_degen_vec.end(), 
                     block_degen_gpuvec.begin());
    #endif
}


void 
c_CNT::Define_SortedModeVector()
{
    int num_double_degenerate_modes = int(atoms_per_ring/2);
    int num_singly_degenerate_modes = atoms_per_ring%2;
    int total_modes = num_double_degenerate_modes + num_singly_degenerate_modes;
    mode_vec.resize(total_modes);
    mode_degen_vec.resize(total_modes);

    for(int m=0; m<total_modes; ++m) 
    {
        if(m < num_double_degenerate_modes) 
        {
            mode_degen_vec[m] = 2;
        }
        else 
        {
            mode_degen_vec[m] = 1;
        }
    }

    /*Hard-coding modes for 17,0 nanotube*/
    mode_vec[0] = 6;
    mode_vec[1] = 5;
    mode_vec[2] = 7;
    mode_vec[3] = 4;
    mode_vec[4] = 8;
    mode_vec[5] = 3;
    mode_vec[6] = 2;
    mode_vec[7] = 1;
    mode_vec[8] = 17;
}

ComplexType 
c_CNT::get_beta(int J)
{
   ComplexType arg(0., -MathConst::pi*J/type_id[0]);
   return 2. * gamma * cos(-1*arg.imag());// * exp(arg); 
}

void
c_CNT::Define_MPI_BlkType ()           
{

    MPI_Type_vector(1, NUM_MODES, NUM_MODES, MPI_DOUBLE_COMPLEX, &MPI_BlkType);
    MPI_Type_commit(&MPI_BlkType);

}


void 
c_CNT::AllocateArrays () 
{
    c_NEGF_Common<BlkType>:: AllocateArrays (); 

}
//
//
//void 
//c_CNT::DeallocateArrays () 
//{
//
//}
//
//


void 
c_CNT::ConstructHamiltonian() 
{
    /*Here we define -H0 where H0 is Hamiltonian of flat bands*/

    auto const& h_Ha = h_Ha_loc_data.table();
    auto const& h_Hb = h_Hb_loc_data.table();
    auto const& h_Hc = h_Hc_loc_data.table();

    for (std::size_t i = 0; i < blkCol_size_loc; ++i)
    {
        h_Ha(i) = 0.;
    }

    for (int j=0; j<NUM_MODES; ++j) 
    {
        int J = mode_vec[j];
        beta.block[j] = get_beta(J);
    }
    amrex::Print() << "\n Printing beta: "<< "\n";
    amrex::Print() << beta << "\n";

    for (std::size_t i = 0; i < offDiag_repeatBlkSize; ++i)
    {
       if(i%offDiag_repeatBlkSize == 0) {
          h_Hb(i) = -1*beta; /*negative sign because (E[I] - [H]) will have negative B and C*/
          h_Hc(i) = -1*beta;
       }
       else {
          h_Hb(i) = -gamma;
          h_Hc(i) = -gamma;  
       }
    }

}

void
c_CNT:: Define_ContactInfo ()
{

    /*define arrays depending on Hsize_glo*/
    global_contact_index[0] = 0;
    global_contact_index[1] = Hsize_glo-1;
    contact_transmission_index[0] = Hsize_glo-1;
    contact_transmission_index[1] = 0;

    /*define tau*/
    auto const& h_tau = h_tau_glo_data.table();
    for (std::size_t c = 0; c < NUM_CONTACTS; ++c)
    {
        h_tau(c) = gamma;
    }
}

AMREX_GPU_HOST_DEVICE
void
c_CNT:: Compute_SurfaceGreensFunction(MatrixBlock<BlkType>& gr, const ComplexType EmU)
{

   auto EmU_sq = pow(EmU,2.);
   auto gamma_sq = pow(gamma,2.);

   for(int i=0; i < NUM_MODES; ++i) 
   {
       auto Factor = EmU_sq + gamma_sq - pow(beta.block[i],2);

       auto Sqrt = sqrt(pow(Factor,2) - 4. * EmU_sq * gamma_sq);

       auto Denom = 2. * gamma_sq * EmU;

       auto val1 = (Factor + Sqrt) / Denom;
       auto val2 = (Factor - Sqrt) / Denom;

       if(val1.imag() < 0.) gr.block[i] = val1;
       else if(val2.imag() < 0.) gr.block[i] = val2; 

       //amrex::Print() << "EmU: " << EmU << "\n";
       //amrex::Print() << "Factor: " << Factor << "\n";
       //amrex::Print() << "Sqrt: "  << Sqrt << "\n";
       //amrex::Print() << "Denom: " << Denom << "\n";
       //amrex::Print() << "Numerator: " << Factor+Sqrt << "\n";
       //amrex::Print() << "Value: " << (Factor+Sqrt)/Denom << "\n";
   }
}


void 
c_CNT::Define_EnergyLimits ()
{

    /*set in the input*/
    E_f = -1.2; 
    E_valence_min = -10.; 
    E_pole_max = 3; 
    c_NEGF_Common<BlkType>:: Define_EnergyLimits ();

}


void 
c_CNT::Define_IntegrationPaths ()
{
    c_NEGF_Common<BlkType>:: Define_IntegrationPaths ();
}

//void 
//c_CNT::ComputeChargeDensity () 
//{
//
//}
