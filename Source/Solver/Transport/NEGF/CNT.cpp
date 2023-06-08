#include "CNT.H"

#include "../../Utils/CodeUtils/CodeUtil.H"

#include <AMReX_Particles.H>

#include <cmath>
#include <math.h>
#include<stdlib.h>

amrex::Array<int,2> c_CNT::type_id; //very important

/*Next, member class functions are defined*/

void
c_CNT::Deallocate ()
{
    c_NEGF_Common<BlkType>::Deallocate();
}

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

    R_cnt = acc*sqrt(3.*(pow(type_id[0],2) 
                   + pow(type_id[1],2) 
                   + type_id[0]*type_id[1])) / (2.*MathConst::pi);
    amrex::Print() << "#####* R_cnt / (nm): " << R_cnt/1.e-9 << "\n";
    amrex::Print() << "#####* D_cnt / (nm): " << 2*R_cnt/1.e-9 << "\n";
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


    E_f = 0.;
    queryWithParser(pp_ns,"E_f", E_f); 
    amrex::Print() << "#####* Fermi level, E_f (eV): " << E_f << "\n";

    E_valence_min = -10.; 
    queryWithParser(pp_ns,"E_valence_min", E_valence_min); 
    amrex::Print() << "#####* valence_band_lower_limit, E_valence_min (eV): " << E_valence_min << "\n";

    E_pole_max = 3; 
    queryWithParser(pp_ns,"E_pole_max", E_pole_max); 
    amrex::Print() << "#####* pole_energy_upper_limit, E_pole_max (eV): " << E_pole_max << "\n";


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


s_Position3D
c_CNT:: get_AtomPosition_ZigZag_CNT(int ring_id, int atom_id)
{
     amrex::Real theta = 0;  
     amrex::Real m = static_cast<amrex::Real>(num_atoms_per_field_site);

     s_Position3D center_offset;

     center_offset.dir[0] = 0.;
     center_offset.dir[1] = -(num_unitcells/2.)*3.*acc - acc/2.;
     center_offset.dir[2] = 0.;
  
     s_Position3D pos;
  
     int unit_cell_id = ring_id/rings_per_unitcell;
     amrex::Real unit_cell_offset = 3*acc*unit_cell_id;

     if(ring_id % rings_per_unitcell == 0) pos.dir[1] = acc + unit_cell_offset;
     if(ring_id % rings_per_unitcell == 1) pos.dir[1] = 1.5*acc + unit_cell_offset;
     if(ring_id % rings_per_unitcell == 2) pos.dir[1] = 2.5*acc + unit_cell_offset;
     if(ring_id % rings_per_unitcell == 3) pos.dir[1] = 3*acc + unit_cell_offset;
   
     pos.dir[1] += center_offset.dir[1];

     if(ring_id %4 == 0 or ring_id %4 == 3) 
     {
         theta = (atom_id/m)*2.*MathConst::pi;
     }
     else 
     {
         theta = (atom_id/m)*2.*MathConst::pi + (sqrt(3.)/2.)*(acc/R_cnt);
     }
     pos.dir[2] = R_cnt*cos(theta) + center_offset.dir[2];
     pos.dir[0] = R_cnt*sin(theta) + center_offset.dir[0];
  
     return pos; 
}


void
c_CNT:: Generate_AtomLocations (amrex::Vector<s_Position3D>& pos)
{
    int m = num_atoms_per_field_site;
    int counter = 0;
    for(int i = 0; i<num_field_sites; ++i) 
    {
       for(int j = 0; j<num_atoms_per_field_site; ++j) 
       {
           pos[i*m + j] = get_AtomPosition_ZigZag_CNT(i, j);
           counter +=1;
           //amrex::Print() << counter << " " << pos[i*m+j].dir[0]/1.e-9 << "  " 
           //                                 << pos[i*m+j].dir[1]/1.e-9 << "  "
           //                                 << pos[i*m+j].dir[2]/1.e-9 << "\n";
       }
    }
    c_NEGF_Common<BlkType>::Generate_AtomLocations(pos);
}


amrex::Real
c_CNT::Get_Bandgap_Of_Mode(int p)
{
    int m = type_id[0];
    int n = type_id[1];
    amrex::Array<amrex::Real,6> Eg_case;
    Eg_case[0] = (acc*gamma/R_cnt)*fabs( 3*p - (2*m + n  ));
    Eg_case[1] = (acc*gamma/R_cnt)*fabs( 3*p - (  m + 2*n));
    Eg_case[2] = (acc*gamma/R_cnt)*fabs( 3*p - (  m - n  ));
    Eg_case[3] = (acc*gamma/R_cnt)*fabs(-3*p + (2*m + n  ));
    Eg_case[4] = (acc*gamma/R_cnt)*fabs(-3*p + (  m + 2*n));
    Eg_case[5] = (acc*gamma/R_cnt)*fabs( 3*p - (  m - n  ));

    amrex::Real Eg_min = Eg_case[0];
    int min_index = 0;
    for(int i=1; i<6; ++i) 
    {
        if(Eg_min > Eg_case[i]) 
	{
            Eg_min = Eg_case[i];
	    min_index = i;
	}
    }
    return Eg_min;
}

void swap(amrex::Vector<amrex::Real>& vec, amrex::Vector<int>& mode_vec, int pos1, int pos2){
	amrex::Real temp;
	temp = vec[pos1];
	vec[pos1] = vec[pos2];
	vec[pos2] = temp;
	
	int temp_id;
	temp_id = mode_vec[pos1];
	mode_vec[pos1] = mode_vec[pos2];
	mode_vec[pos2] = temp_id;
}

int partition(amrex::Vector<amrex::Real>& vec, amrex::Vector<int>& mode_vec, int low, int high, amrex::Real pivot){
	int i = low;
	int j = low;
	while( i <= high){
		if(vec[i] > pivot){
			i++;
		}
		else{
			swap(vec, mode_vec, i, j);
			i++;
			j++;
		}
	}
	return j-1;
}

void quickSort(amrex::Vector<amrex::Real>& vec, amrex::Vector<int>& mode_vec, int low, int high)
{
    //Credit: https://favtutor.com/blogs/quick-sort-cpp
    if(low < high)
    {
        amrex::Real pivot = vec[high];
        int pos = partition(vec, mode_vec, low, high, pivot);
        
        quickSort(vec, mode_vec, low, pos-1);
        quickSort(vec, mode_vec, pos+1, high);
    }
}


void 
c_CNT::Define_SortedModeVector()
{
    amrex::Vector<amrex::Real> Eg_vec;
    Eg_vec.resize(atoms_per_ring);

    amrex::Vector<int> mode_index_vec;
    mode_index_vec.resize(atoms_per_ring);

    amrex::Print() << "All modes: \n";
    for(int p=1; p<atoms_per_ring+1; ++p) 
    {
	amrex::Real Eg = Get_Bandgap_Of_Mode(p);
	
	Eg_vec[p-1] = Eg;
	mode_index_vec[p-1] = p;

	amrex::Print() << "mode: " << p << " bandgap (nm): " << Eg << "\n";
    }

    quickSort(Eg_vec, mode_index_vec, 0, atoms_per_ring - 1) ;


    mode_vec.resize(0);
    bandgap_vec.resize(0);
    mode_degen_vec.resize(0);

    mode_vec.push_back(mode_index_vec[0]);
    bandgap_vec.push_back(Eg_vec[0]);
    mode_degen_vec.push_back(1);

    int counter = 1;
    for(int p=1; p < atoms_per_ring; ++p) 
    {
	 if((Eg_vec[p] - bandgap_vec[counter-1]) > 1e-6) 
	 {
             mode_vec.push_back(mode_index_vec[p]);		 
             bandgap_vec.push_back(Eg_vec[p]);
             mode_degen_vec.push_back(1);
	     counter++;
	 }
	 else 
	 { 
             mode_degen_vec[counter-1] += 1;
	 }
    }

    int total_modes = mode_vec.size();

    amrex::Print() << "\nSorted, distinct modes: " << total_modes << "\n";
    for(int i=0; i<total_modes; ++i) 
    {
        amrex::Print() << "mode: " << std::setw(5) << mode_vec[i] << std::setw(5) 
		       << ", bandgap / (nm): " << std::setw(15)  << bandgap_vec[i] << std::setw(15) 
		       << ", degeneracy: " << std::setw(5) << mode_degen_vec[i]  << std::setw(5) << "\n";

    }

    Eg_vec.clear();
    mode_index_vec.clear();
  
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

    auto const& h_minusHa = h_minusHa_loc_data.table();
    auto const& h_Hb = h_Hb_loc_data.table();
    auto const& h_Hc = h_Hc_loc_data.table();

    for (std::size_t i = 0; i < blkCol_size_loc; ++i)
    {
        h_minusHa(i) = 0.;
    }

    for (int j=0; j<NUM_MODES; ++j) 
    {
        int J = mode_vec[j];
        beta.block[j] = get_beta(J);
    }
    //amrex::Print() << "\n Printing beta: "<< "\n";
    //amrex::Print() << beta << "\n";

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
