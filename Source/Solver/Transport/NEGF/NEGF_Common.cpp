#include "NEGF_Common.H"
#include "Matrix_Block_Util.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/WarpXConst.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"


/*Explicit specializations*/
template class c_NEGF_Common<ComplexType[NUM_MODES]>;            //of c_CNT
template class c_NEGF_Common<ComplexType[NUM_MODES][NUM_MODES]>; //c_Graphene
								 

template<typename T>
void
c_NEGF_Common<T>:: Deallocate ()
{
    Potential.clear();
    PTD.clear();

    /*Broyden*/
    h_n_curr_in_data.clear();
    d_n_curr_in_data.clear();
    n_curr_out_data.clear();
    n_prev_in_data.clear();
    F_curr_data.clear();

}

template<typename T>
void
c_NEGF_Common<T>:: Initialize_ChargeAtFieldSites(const amrex::Real value)
{

    SetVal_Table1D(h_n_curr_in_data, value);

}

template<typename T>
void
c_NEGF_Common<T>:: Reset_Broyden ()
{
    int size = W_Broyden.size();	
    for(int j=0; j<size; ++j) 
    {
        W_Broyden[j]->clear();    
        V_Broyden[j]->clear();    
    }
    W_Broyden.clear();
    V_Broyden.clear();

    /* At present, we set n_curr_in as the converged charge density for previous gate voltage*/
    SetVal_Table1D(n_curr_out_data,0.);
    SetVal_Table1D(n_prev_in_data,0.);
    SetVal_Table1D(F_curr_data,0.);

    SetVal_Table2D(Jinv_curr_data,0.);

    Broyden_Step = 1;
    Broyden_Norm = 1;
}

template<typename T>
void
c_NEGF_Common<T>:: ReadNanostructureProperties ()
{

    amrex::ParmParse pp_ns(name);

    getWithParser(pp_ns,"num_unitcells", num_unitcells);
    amrex::Print() << "##### num_unitcells: " << num_unitcells << "\n";

    amrex::Vector<amrex::Real> vec_offset;
    getArrWithParser(pp_ns, "offset", vec_offset, 0, AMREX_SPACEDIM);
    offset = vecToArr(vec_offset);

    amrex::Print() << "##### offset: ";
    for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << offset[i] << "  ";
    amrex::Print() << "\n";

    std::string avg_type_str = "all";
    pp_ns.query("field_averaging_type", avg_type_str);
    if(avg_type_str == "all" or avg_type_str == "All" or avg_type_str == "ALL")
    {
        avg_type = s_AVG_TYPE::ALL;
    }
    if(avg_type_str == "specific" or avg_type_str == "Specific" or avg_type_str == "SPECIFIC")
    {
        avg_type = s_AVG_TYPE::SPECIFIC;
    }
    amrex::Print() << "##### field_averaging_type: " << avg_type_str << " enum id: " << avg_type << "\n";
   
    if(avg_type == s_AVG_TYPE::SPECIFIC) {
       pp_ns.getarr("atom_indices_for_averaging", vec_avg_indices);
       amrex::Print() << "##### atom_indices_for_averaging: \n";
       for (int i=0; i<vec_avg_indices.size(); ++i) 
       {
           amrex::Print() << vec_avg_indices[i] << "  ";
       }

       #if AMREX_USE_GPU
       gpuvec_avg_indices.resize(vec_avg_indices.size());
       amrex::Gpu::copy(amrex::Gpu::hostToDevice, 
                        vec_avg_indices.begin(), 
                        vec_avg_indices.end(), 
                        gpuvec_avg_indices.begin());
       #endif
    }

    amrex::Vector<amrex::Real> vec_Temperature;
    auto is_specified = queryArrWithParser(pp_ns, "contact_temperature", vec_Temperature, 0, NUM_CONTACTS);
    if(is_specified) {
        for(int c=0; c<NUM_CONTACTS; ++c) {
            Contact_Temperature[c] = vec_Temperature[c];
        }
    }

    amrex::Vector<amrex::Real> vec_contact_potential;
    is_specified = queryArrWithParser(pp_ns, "contact_potential", vec_contact_potential, 0, NUM_CONTACTS);
    if(is_specified) {
        for(int c=0; c<NUM_CONTACTS; ++c) {
            Contact_Potential[c] = vec_contact_potential[c];
        }
    }

    pp_ns.query("impose_potential", flag_impose_potential);
    if(flag_impose_potential) 
    {
        pp_ns.query("potential_profile_type", potential_profile_type_str);
    }
    
    set_material_specific_parameters();

    amrex::Print() << "\n#####* num_atoms: "                << num_atoms                  << "\n";
    amrex::Print() << "#####* num_atoms_per_unitcell: "   << num_atoms_per_unitcell     << "\n";
    amrex::Print() << "#####* num_unitcells: "            << num_unitcells              << "\n";
    amrex::Print() << "#####* num_field_sites: "          << num_field_sites            << "\n";
    amrex::Print() << "#####* average_field_flag: "       << average_field_flag         << "\n";
    amrex::Print() << "#####* num_atoms_per_field_site: " << num_atoms_per_field_site   << "\n";
    amrex::Print() << "#####* num_atoms_to_avg_over: "    << num_atoms_to_avg_over      << "\n";
    amrex::Print() << "#####* primary_transport_dir: "    << primary_transport_dir      << "\n";
    amrex::Print() << "#####* Broyden_fraction: "         << Broyden_fraction           << "\n";
    amrex::Print() << "#####* Contact_Temperature: \n";
    for(int c=0; c<NUM_CONTACTS; ++c) {
        amrex::Print() << "#####*   contact, temperature (K): " << c << "  " << Contact_Temperature[c] <<"\n";
    }
    amrex::Print() << "#####* Contact_Potential: \n";
    for(int c=0; c<NUM_CONTACTS; ++c) {
        amrex::Print() << "#####*   contact, potential (V): " << c << "  " << Contact_Potential[c] <<"\n";
    }
    amrex::Print() << "##### flag_impose_potential: " << flag_impose_potential << "\n";
    amrex::Print() << "##### potential_profile_type: " << potential_profile_type_str << "\n";

    Potential.resize(num_field_sites);
    PTD.resize(num_field_sites);

    /*Broyden*/
    Broyden_Step = 1;

    h_n_curr_in_data.resize({0},{num_field_sites}, The_Pinned_Arena());
    
    #if AMREX_USE_GPU
    d_n_curr_in_data.resize({0},{num_field_sites}, The_Arena());
    #endif 

    n_curr_out_data.resize({0},{num_field_sites}, The_Pinned_Arena());
    n_prev_in_data.resize({0},{num_field_sites}, The_Pinned_Arena());
    F_curr_data.resize({0},{num_field_sites}, The_Pinned_Arena());

    SetVal_Table1D(h_n_curr_in_data,0.);
    SetVal_Table1D(n_curr_out_data,0.);
    SetVal_Table1D(n_prev_in_data,0.);
    SetVal_Table1D(F_curr_data,0.);


    Jinv_curr_data.resize({0,0},{num_field_sites, num_field_sites}, The_Pinned_Arena());
    SetVal_Table2D(Jinv_curr_data,0.);

    auto const& Jinv_curr    = Jinv_curr_data.table();
    for(int a=0; a < num_field_sites; ++a) 
    {
        Jinv_curr(a,a) = Broyden_fraction;
    }

    for(int l=0; l < num_field_sites; ++l) 
    {
        Potential[l] = 0.;
        PTD[l]       = 0.;
    }

}


template<typename T>
void 
c_NEGF_Common<T>::set_material_specific_parameters() 
{
    /*set the following in the specialization*/
    /*num_atoms: number of total atoms*/
    /*num_atoms_per_unitcell: number of atoms per unitcell*/
    /*num_unitcells*/
    /*num_field_sites*/
    /*average_field_flag*/
    /*num_atoms_per_field_site*/
    /*num_atoms_to_avg_over*/
    /*primary_transport_dir*/
    /*block_degen_vec: this is vector of degeneracy factors of size equal to that of a block element*/ 
    /*block_degen_gpuvec: device vector copy of block_degen_vec*/
    /*set the size ofprimary_transport_dir*/
}


template<typename T>
void 
c_NEGF_Common<T>:: Generate_AtomLocations (amrex::Vector<s_Position3D>& pos)
{
    /*First code atom locations for particular material in the specialization*/
    /*Then call this function and apply global offset specified by the user*/
    /*Also, specify PrimaryTransportDirection PTD array*/

    for(int i = 0; i<num_atoms; ++i)
    {
       for(int j = 0; j < AMREX_SPACEDIM; ++j)
       {
           pos[i].dir[j] += offset[j];
       }  
    }

    //amrex::Print() << "PTD is specified in nm! \n";
    for(int l=0; l< num_field_sites; ++l) 
    {
       PTD[l] = pos[l*num_atoms_per_field_site].dir[primary_transport_dir] / 1.e-9;
       //amrex::Print() << l << "  " << PTD[l] << "\n";
    }
}


template<typename T>
void
c_NEGF_Common<T>:: Define_PotentialProfile()
{

    amrex::Real U1 = -Contact_Potential[0];
    amrex::Real U2 = -Contact_Potential[1];

    switch(map_PotentialProfile[potential_profile_type_str])
    {
        case s_Potential_Profile::Type::CONSTANT:
        {
            for (int l=0; l<num_field_sites; ++l)
            {
                Potential[l]   = U1;
            }
            break;
        }
        case s_Potential_Profile::Type::LINEAR:
        {
		
            for (int l=0; l<num_field_sites; ++l)
            {
                Potential[l]   = U1 + (static_cast<amrex::Real>(l)/(num_field_sites-1.))*(U2-U1);
            }
            break;
        }
        case s_Potential_Profile::Type::POINT_CHARGE:
        {
            //amrex::Array<amrex::Real,2> QD_loc = {0., 1}; //1nm away in z
            //for (int l=0; l<NSType::num_field_sites; ++l)
            //{
            //    amrex::Real r = sqrt(pow((PTD[l] - QD_loc[0]),2.) + pow(QD_loc[1],2))*1e-9;
            //    NSType::Potential[l]   = -1*(1./(4.*MathConst::pi*PhysConst::ep0*1.)*(PhysConst::q_e/r));
            //}
            break;
        }
    }
}


template<typename T>
void 
c_NEGF_Common<T>::DefineMatrixPartition() 
{

    Hsize_glo = get_Hsize(); 
    amrex::Print() << "\nHsize_glo: " << Hsize_glo << "\n";

    bool flag_fixed_blk_size = false;
    const int THRESHOLD_BLKCOL_SIZE = 40000; /*matrix size*/

    if(flag_fixed_blk_size) amrex::Print() << "max_blkCol_perProc is fixed by user\n";
    else
    {
         amrex::Print() << "max_blkCol_perProc is computed at run-time\n";

         max_blkCol_perProc = ceil(static_cast<amrex::Real>(Hsize_glo)/num_proc);
         if(max_blkCol_perProc > THRESHOLD_BLKCOL_SIZE)
         {
            max_blkCol_perProc = THRESHOLD_BLKCOL_SIZE;
            /*assert that use larger number of procs*/
         }
    }
    amrex::Print() << "max block columns per proc: " << max_blkCol_perProc << "\n";


    num_proc_with_blkCol = ceil(static_cast<amrex::Real>(Hsize_glo)/max_blkCol_perProc);
    amrex::Print() << "number of procs with block columns: " << num_proc_with_blkCol<< "\n";
    /*if num_proc_with_blk >= num_proc, assert.*/


    vec_cumu_blkCol_size.resize(num_proc_with_blkCol + 1);
    vec_cumu_blkCol_size[0] = 0;
    for(int p=1; p < num_proc_with_blkCol; ++p)
    {
        vec_cumu_blkCol_size[p] = vec_cumu_blkCol_size[p-1] + max_blkCol_perProc;
        /*All proc except the last one contains max_blkCol_perProc number of column blks.*/
    }
    vec_cumu_blkCol_size[num_proc_with_blkCol] = Hsize_glo;


    blkCol_size_loc = 0;
    if(my_rank < num_proc_with_blkCol)
    {
        int blk_gid = my_rank;
        blkCol_size_loc = vec_cumu_blkCol_size[blk_gid+1] - vec_cumu_blkCol_size[blk_gid];

        /*check later:setting vec_blkCol_gids may not be necessary*/
        for (int c=0; c < blkCol_size_loc; ++c)
        {
            int col_gid = vec_cumu_blkCol_size[blk_gid] + c;
            vec_blkCol_gids.push_back(col_gid);
        }
    }

    MPI_recv_count.resize(num_proc);
    MPI_disp.resize(num_proc);

    for(int p=0; p < num_proc; ++p) {

       if(p < num_proc_with_blkCol) {
          MPI_recv_count[p] = vec_cumu_blkCol_size[p+1] - vec_cumu_blkCol_size[p];
          MPI_disp[p] = vec_cumu_blkCol_size[p];
       }
       else {
          MPI_recv_count[p] = 0;
          MPI_disp[p] = 0;
       }
       //amrex::Print() << "p,recv, disp: " << p << "  " 
       //               << MPI_recv_count[p] << "  " 
       //               << MPI_disp[p] << "\n";
    }

    Define_MPI_BlkType();

}


template<typename T>
void 
c_NEGF_Common<T>::Define_MPI_BlkType () 
{
    /*define MPI_BlkType in overriden class*/
}


template<typename T>
void 
c_NEGF_Common<T>::Define_ContactInfo () 
{
    /*define the following in overriden functions:
     *global_contact_index
     *contact_transmission_index
     *h_tau
     */
}

template<typename T>
void 
c_NEGF_Common<T>::AllocateArrays () 
{
    ComplexType zero(0.,0.);   
    h_minusHa_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    SetVal_Table1D(h_minusHa_loc_data,zero);

    offDiag_repeatBlkSize = get_offDiag_repeatBlkSize();
    h_Hb_loc_data.resize({0},{offDiag_repeatBlkSize}, The_Pinned_Arena());
    SetVal_Table1D(h_Hb_loc_data,zero);

    h_Hc_loc_data.resize({0},{offDiag_repeatBlkSize}, The_Pinned_Arena());
    SetVal_Table1D(h_Hc_loc_data,zero);

    h_tau_glo_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    SetVal_Table1D(h_tau_glo_data,zero);

    #if AMREX_USE_GPU
    d_GR_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    d_A_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    d_Rho0_loc_data.resize({0},{blkCol_size_loc}, The_Arena());
    d_RhoEq_loc_data.resize({0},{blkCol_size_loc}, The_Arena());
    d_RhoNonEq_loc_data.resize({0},{blkCol_size_loc}, The_Arena());
    d_GR_atPoles_loc_data.resize({0},{blkCol_size_loc}, The_Arena());
    #else
    h_GR_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    SetVal_Table2D(h_GR_loc_data, zero);

    h_A_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    SetVal_Table2D(h_A_loc_data, zero);
    #endif

    h_Rho0_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    SetVal_Table1D(h_Rho0_loc_data,zero);

    h_RhoEq_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    SetVal_Table1D(h_RhoEq_loc_data,zero);

    h_RhoNonEq_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    SetVal_Table1D(h_RhoNonEq_loc_data,zero);

    h_GR_atPoles_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    SetVal_Table1D(h_GR_atPoles_loc_data,zero);

    h_E_RealPath_data.resize({0},{NUM_ENERGY_PTS_REAL},The_Pinned_Arena());

}


template<typename T>
void 
c_NEGF_Common<T>:: ConstructHamiltonian () {}


template<typename T>
void 
c_NEGF_Common<T>:: AddPotentialToHamiltonian () 
{
    auto const& h_minusHa = h_minusHa_loc_data.table();
    int c=0;
    for(auto& col_gid: vec_blkCol_gids)
    {
        h_minusHa(c) = -1*(Potential[col_gid]); /*Note: Ha = H0a (=0) + U, so -Ha = -U */
        ++c;
    }
}


template<typename T>
void 
c_NEGF_Common<T>:: Update_ContactPotential () 
{
    amrex::Print() <<  "Updated contact potential: \n";
    for(int c=0; c < NUM_CONTACTS; ++c)
    {
        U_contact[c] = Potential[global_contact_index[c]];
        amrex::Print() << "  contact, U: " <<  c << " " << U_contact[c] << "\n";
    }
}


template<typename T>
void 
c_NEGF_Common<T>:: Define_EnergyLimits ()
{
    for (int c=0; c<NUM_CONTACTS; ++c)
    {
        //mu_contact[c] = E_f + U_contact[c];
        mu_contact[c] = E_f ;
        kT_contact[c] = PhysConst::kb_eVperK*Contact_Temperature[c]; /*set Temp in the input*/
    }

    mu_min = mu_contact[0];
    mu_max = mu_contact[0];
    kT_min = kT_contact[0];
    kT_max = kT_contact[0];

    flag_noneq_exists = false;

    for (int c=1; c<NUM_CONTACTS; ++c)
    {
       if(mu_min > mu_contact[c])
       {
           mu_min = mu_contact[c];
       }
       if(mu_max < mu_contact[c])
       {
           mu_max = mu_contact[c];
       }
       if(kT_min > kT_contact[c])
       {
           kT_min = kT_contact[c];
       }
       if(kT_max < kT_contact[c])
       {
           kT_max = kT_contact[c];
       }
    }
    if(fabs(mu_min - mu_max) > 1e-8) flag_noneq_exists = true;
    if(fabs(kT_min - kT_max) > 0.01) flag_noneq_exists = true;

    ComplexType val(0.,1e-8);
    E_zPlus = val;
    E_contour_left  = E_valence_min + E_zPlus; /*set in the input*/
    E_rightmost = mu_max + Fermi_tail_factor*kT_max + E_zPlus;
    //E_rightmost = E_f + Fermi_tail_factor*kT_max + E_zPlus;

    if(flag_noneq_exists)
    {
        amrex::Print() << "\nnonequilibrium exists between terminals\n";
        E_contour_right = mu_min - Fermi_tail_factor*kT_max + E_zPlus;
        num_enclosed_poles = 0;
    }
    else
    {
        E_contour_right = E_rightmost;
        num_enclosed_poles = int((E_pole_max-MathConst::pi*kT_min)/(2.*MathConst::pi*kT_min) + 1);
        E_poles_vec.resize(num_enclosed_poles);

        for(int p=0; p<num_enclosed_poles; ++p) 
        {
            ComplexType pole(mu_min, MathConst::pi*kT_min*(2*p+1));
            E_poles_vec[p] = pole;
        }
    }
    ComplexType val2(E_contour_right.real(), 2*num_enclosed_poles*MathConst::pi*kT_max);
    E_zeta = val2;
    //E_eta =  E_zeta - 2*Fermi_tail_factor*kT_max;
    ComplexType val3(mu_max - Fermi_tail_factor*kT_max, E_zeta.imag());
    //ComplexType val3(E_f - Fermi_tail_factor*kT_max, E_zeta.imag());
    E_eta =  val3;

    //amrex::Print() << "U_contact: ";
    //for (int c=0; c<NUM_CONTACTS; ++c)
    //{
    //    amrex::Print() <<  U_contact[c] << " ";
    //}
    //amrex::Print() << "\n";
    //amrex::Print() << "\nE_f: " << E_f << "\n";
    //amrex::Print() << "mu_min/max: " << mu_min << " " << mu_max << "\n";
    //amrex::Print() << "kT_min/max: " << kT_min << " " << kT_max << "\n";
    //amrex::Print() << "E_zPlus: "  << E_zPlus << "\n";
    //amrex::Print() << "E_contour_left/E_contour_right/E_rightmost: " << E_contour_left      <<  "  "
    //                                                                 << E_contour_right << "  "
    //                                                                 << E_rightmost     << "\n";
    //amrex::Print() << "E_pole_max: " << E_pole_max << ", number of poles: " << num_enclosed_poles << "\n";
    //amrex::Print() << "E_zeta: " << E_zeta << "\n";
    //amrex::Print() << "E_eta: "  << E_eta << "\n";

}


template<typename T>
ComplexType
c_NEGF_Common<T>::FermiFunction(ComplexType E_minus_Mu, const amrex::Real kT)
{
    ComplexType one(1., 0.);
    return one / (exp(E_minus_Mu / kT) + one);
}


template<typename T>
void 
c_NEGF_Common<T>:: Define_IntegrationPaths ()
{
    /* Define_ContourPath_Rho0 */
    ContourPath_Rho0.resize(3); 
    ContourPath_Rho0[0].Define_GaussLegendrePoints(E_zPlus, E_zeta, 100, 0); 
    ContourPath_Rho0[1].Define_GaussLegendrePoints(E_zeta,  E_eta, 100, 0); 
    ContourPath_Rho0[2].Define_GaussLegendrePoints(E_eta, E_contour_left, 100, 1); 

    /* Define_ContourPath_DOS */
    ComplexType minus_one(-0.2,E_zPlus.imag());
    ComplexType one(0.2,E_zPlus.imag());
    ContourPath_DOS.Define_GaussLegendrePoints(minus_one, one, 10,0);
    //ContourPath_DOS.Define_GaussLegendrePoints(E_contour_left, -E_contour_left, 400,0);

    /* Write Fermi function */
    auto* Fermi_path = &ContourPath_DOS;
    if(flag_noneq_exists) Fermi_path = &ContourPath_RhoNonEq[0];

    //amrex::Print() << "Writing Fermi Function over points" << Fermi_path->num_pts; 
    //std::string filename = "FermiFunction";
    //std::ofstream outfile;
    //outfile.open(filename.c_str());

    //for(int i=0; i<Fermi_path->num_pts; ++i) 
    //{
    //    auto E = Fermi_path->E_vec[i];
    //    outfile << std::setw(20) << E.real()
    //            << std::setw(20) << FermiFunction(E-mu_contact[0], kT_contact[0]).real() 
    //            << std::setw(20) << FermiFunction(E-mu_contact[1], kT_contact[1]).real()
    //            << "\n";
    //}

    //outfile.close();
}


template<typename T>
void 
c_NEGF_Common<T>:: Update_IntegrationPaths ()
{
    ContourPath_RhoEq.clear();
    ContourPath_RhoNonEq.clear();

    /* Define_ContourPath_RhoEq */
    ContourPath_RhoEq.resize(3); 
    ContourPath_RhoEq[0].Define_GaussLegendrePoints(E_contour_right, E_zeta, 100, 0); 
    ContourPath_RhoEq[1].Define_GaussLegendrePoints(E_zeta, E_eta, 100, 0); 
    ContourPath_RhoEq[2].Define_GaussLegendrePoints(E_eta, E_contour_left, 100, 1); 

    /* Define_ContourPath_RhoNonEq */
    if(flag_noneq_exists) 
    {
        ContourPath_RhoNonEq.resize(1); 
        ContourPath_RhoNonEq[0].Define_GaussLegendrePoints(E_contour_right, E_rightmost, 400, 0);
    }

}

template<typename T>
void 
c_NEGF_Common<T>:: Allocate_TemporaryArraysForGFComputation ()
{
    ComplexType zero(0.,0.);   
    h_Alpha_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    SetVal_Table1D(h_Alpha_loc_data,zero);

    h_Alpha_glo_data.resize({0},{Hsize_glo}, The_Pinned_Arena());
    SetVal_Table1D(h_Alpha_glo_data,zero);

    h_Xtil_glo_data.resize({0},{Hsize_glo}, The_Pinned_Arena());
    SetVal_Table1D(h_Xtil_glo_data,zero);

    h_Ytil_glo_data.resize({0},{Hsize_glo}, The_Pinned_Arena());
    SetVal_Table1D(h_Ytil_glo_data,zero);

    h_X_glo_data.resize({0},{Hsize_glo}, The_Pinned_Arena());
    SetVal_Table1D(h_X_glo_data,zero);

    h_Y_glo_data.resize({0},{Hsize_glo}, The_Pinned_Arena());
    SetVal_Table1D(h_Y_glo_data,zero);

    h_X_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    SetVal_Table1D(h_X_loc_data,zero);

    h_Y_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    SetVal_Table1D(h_Y_loc_data,zero);

    h_Sigma_contact_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    SetVal_Table1D(h_Sigma_contact_data,zero);

    h_Fermi_contact_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    SetVal_Table1D(h_Fermi_contact_data,zero);

    h_Alpha_contact_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    SetVal_Table1D(h_Alpha_contact_data,zero);

    h_X_contact_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    SetVal_Table1D(h_X_contact_data,zero);

    h_Y_contact_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    SetVal_Table1D(h_Y_contact_data,zero);

    h_Trace_r.resize(num_traces);
    h_Trace_i.resize(num_traces);

    #ifdef AMREX_USE_GPU
    d_Alpha_loc_data.resize({0},{blkCol_size_loc}, The_Arena());
    d_X_loc_data.resize({0},{blkCol_size_loc}, The_Arena());
    d_Y_loc_data.resize({0},{blkCol_size_loc}, The_Arena());
    d_Xtil_glo_data.resize({0},{Hsize_glo}, The_Arena());
    d_Ytil_glo_data.resize({0},{Hsize_glo}, The_Arena());
    d_Sigma_contact_data.resize({0},{NUM_CONTACTS}, The_Arena());
    d_Fermi_contact_data.resize({0},{NUM_CONTACTS}, The_Arena());

    d_Alpha_contact_data.resize({0},{NUM_CONTACTS}, The_Arena());
    d_X_contact_data.resize({0},{NUM_CONTACTS}, The_Arena());
    d_Y_contact_data.resize({0},{NUM_CONTACTS}, The_Arena());

    d_Trace_r.resize(num_traces);
    d_Trace_i.resize(num_traces);
    #endif 
}

template<typename T>
void 
c_NEGF_Common<T>:: Deallocate_TemporaryArraysForGFComputation ()
{
    h_Alpha_glo_data.clear();
    h_X_glo_data.clear();
    h_Y_glo_data.clear();

    h_Alpha_loc_data.clear();
    h_Xtil_glo_data.clear();
    h_Ytil_glo_data.clear();
    h_X_loc_data.clear();
    h_Y_loc_data.clear();
    h_Sigma_contact_data.clear();
    h_Fermi_contact_data.clear();

    h_Alpha_contact_data.clear();
    h_X_contact_data.clear();
    h_Y_contact_data.clear();

    h_Trace_r.clear();
    h_Trace_i.clear();

    #ifdef AMREX_USE_GPU
    d_Alpha_loc_data.clear();
    d_X_loc_data.clear();
    d_Y_loc_data.clear();
    d_Xtil_glo_data.clear();
    d_Ytil_glo_data.clear();
    d_Sigma_contact_data.clear();
    d_Fermi_contact_data.clear();

    d_Alpha_contact_data.clear();
    d_X_contact_data.clear();
    d_Y_contact_data.clear();

    d_Trace_r.clear();
    d_Trace_i.clear();
    #endif 
}


template<typename T>
void 
c_NEGF_Common<T>:: Compute_DensityOfStates ()
{
    int E_pts = ContourPath_DOS.num_pts;
    h_LDOS_loc_data.resize({0}, {E_pts}, The_Pinned_Arena());
    SetVal_Table1D(h_LDOS_loc_data,0.);

    h_Transmission_loc_data.resize({0}, {E_pts}, The_Pinned_Arena());
    SetVal_Table1D(h_Transmission_loc_data,0.);


    auto const& h_minusHa_loc  = h_minusHa_loc_data.table();
    auto const& h_Hb_loc  = h_Hb_loc_data.table();
    auto const& h_Hc_loc  = h_Hc_loc_data.table();
    auto const& h_tau     = h_tau_glo_data.table();
    auto const& h_LDOS_loc = h_LDOS_loc_data.table();
    auto const& h_Transmission_loc = h_Transmission_loc_data.table();

    Allocate_TemporaryArraysForGFComputation();

    auto const& h_Alpha_loc = h_Alpha_loc_data.table();
    auto const& h_Alpha_glo = h_Alpha_glo_data.table();
    auto const& h_Xtil_glo = h_Xtil_glo_data.table();
    auto const& h_Ytil_glo = h_Ytil_glo_data.table();
    auto const& h_X_glo = h_X_glo_data.table();
    auto const& h_Y_glo = h_Y_glo_data.table();
    auto const& h_X_loc = h_X_loc_data.table();
    auto const& h_Y_loc = h_Y_loc_data.table();
    auto const& h_Alpha_contact = h_Alpha_contact_data.table();
    auto const& h_Y_contact = h_Y_contact_data.table();
    auto const& h_X_contact = h_X_contact_data.table();
    auto const& h_Sigma_contact = h_Sigma_contact_data.table();

    #ifdef AMREX_USE_GPU
    auto const& GR_loc          = d_GR_loc_data.table();
    auto const& A_loc           = d_A_loc_data.table();
    /*constant references*/
    auto const& Alpha           = d_Alpha_loc_data.const_table();
    auto const& Xtil_glo        = d_Xtil_glo_data.const_table();
    auto const& Ytil_glo        = d_Ytil_glo_data.const_table();
    auto const& X               = d_X_loc_data.const_table();
    auto const& Y               = d_Y_loc_data.const_table();

    auto const& Alpha_contact   = d_Alpha_contact_data.const_table();
    auto const& X_contact       = d_X_contact_data.const_table();
    auto const& Y_contact       = d_Y_contact_data.const_table();
    auto const& Sigma_contact   = d_Sigma_contact_data.const_table();

    auto* trace_r               = d_Trace_r.dataPtr();
    auto* trace_i               = d_Trace_i.dataPtr();
    auto& degen_vec             = block_degen_gpuvec;
    #else
    auto const& GR_loc          = h_GR_loc_data.table();
    auto const& A_loc           = h_A_loc_data.table();
    /*constant references*/
    auto const& Alpha           = h_Alpha_loc_data.const_table();
    auto const& Xtil_glo        = h_Xtil_glo_data.const_table();
    auto const& Ytil_glo        = h_Ytil_glo_data.const_table();
    auto const& X               = h_X_loc_data.const_table();
    auto const& Y               = h_Y_loc_data.const_table();

    auto const& Alpha_contact   = h_Alpha_contact_data.const_table();
    auto const& X_contact       = h_X_contact_data.const_table();
    auto const& Y_contact       = h_Y_contact_data.const_table();
    auto const& Sigma_contact   = h_Sigma_contact_data.const_table();

    auto* trace_r               = h_Trace_r.dataPtr();
    auto* trace_i               = h_Trace_i.dataPtr();
    auto& degen_vec             = block_degen_vec;
    #endif
  
    for(int e=0; e < ContourPath_DOS.num_pts; ++e) 
    {

        ComplexType E = ContourPath_DOS.E_vec[e];
        ComplexType weight = ContourPath_DOS.weight_vec[e];
        ComplexType mul_factor = ContourPath_DOS.mul_factor_vec[e];


        for(int n=0; n<blkCol_size_loc; ++n)
        {
            h_Alpha_loc(n) = E + h_minusHa_loc(n); 
            /*+ because h_minusHa is defined previously as -(H0+U)*/
        }

        get_Sigma_at_contacts(h_Sigma_contact_data, E);

        for (int c=0; c<NUM_CONTACTS; ++c)
        {
            int n_glo = global_contact_index[c];

            if(n_glo >= vec_cumu_blkCol_size[my_rank] && n_glo < vec_cumu_blkCol_size[my_rank+1])
            {
                int n = n_glo - vec_cumu_blkCol_size[my_rank];
                h_Alpha_loc(n) = h_Alpha_loc(n) - h_Sigma_contact(c);
            }
        }

        /*MPI_Allgather*/
        MPI_Allgatherv(&h_Alpha_loc(0),
                        blkCol_size_loc,
                        MPI_BlkType,
                       &h_Alpha_glo(0),
                        MPI_recv_count.data(),
                        MPI_disp.data(),
                        MPI_BlkType,
                        ParallelDescriptor::Communicator());

        for (int c=0; c<NUM_CONTACTS; ++c)
        {
            int n_glo = global_contact_index[c];
            h_Alpha_contact(c) = h_Alpha_glo(n_glo);
        }

        h_Y_glo(0) = 0;
        for (int n = 1; n < Hsize_glo; ++n)
        {
        int p = (n-1)%offDiag_repeatBlkSize;
            h_Ytil_glo(n) = h_Hc_loc(p) / ( h_Alpha_glo(n-1) - h_Y_glo(n-1) );
            h_Y_glo(n) = h_Hb_loc(p) * h_Ytil_glo(n);
        }

        h_X_glo(Hsize_glo-1) = 0;
        for (int n = Hsize_glo-2; n > -1; n--)
        {
        int p = n%offDiag_repeatBlkSize;
            h_Xtil_glo(n) = h_Hb_loc(p)/(h_Alpha_glo(n+1) - h_X_glo(n+1));
            h_X_glo(n) = h_Hc_loc(p)*h_Xtil_glo(n);
        }

        for (int c = 0; c < blkCol_size_loc; ++c)
        {
            int n = c + vec_cumu_blkCol_size[my_rank];
            h_Y_loc(c) = h_Y_glo(n);
            h_X_loc(c) = h_X_glo(n);
        }
 
        for (int c = 0; c < NUM_CONTACTS; ++c)
        {
            int n = global_contact_index[c];
            h_Y_contact(c) = h_Y_glo(n);
            h_X_contact(c) = h_X_glo(n);
        }

        for (int t=0; t< num_traces; ++t)
        {
            h_Trace_r[t] = 0.;
            h_Trace_i[t] = 0.;
        }

        #ifdef AMREX_USE_GPU
        d_Alpha_loc_data.copy(h_Alpha_loc_data);
        d_Xtil_glo_data.copy(h_Xtil_glo_data);
        d_Ytil_glo_data.copy(h_Ytil_glo_data);
        d_X_loc_data.copy(h_X_loc_data);
        d_Y_loc_data.copy(h_Y_loc_data);

        d_Alpha_contact_data.copy(h_Alpha_contact_data);
        d_X_contact_data.copy(h_X_contact_data);
        d_Y_contact_data.copy(h_Y_contact_data);
        d_Sigma_contact_data.copy(h_Sigma_contact_data);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_Trace_r.begin(), h_Trace_r.end(), d_Trace_r.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_Trace_i.begin(), h_Trace_i.end(), d_Trace_i.begin());
        amrex::Gpu::streamSynchronize();
        #endif        

    	/*following is for lambda capture*/
        int cumulative_columns = vec_cumu_blkCol_size[my_rank];
        int Hsize = Hsize_glo;  
     	auto& GC_ID = global_contact_index;
    	auto& CT_ID = contact_transmission_index;
        auto* degen_vec_ptr = degen_vec.dataPtr();

        amrex::ParallelFor(blkCol_size_loc, [=] AMREX_GPU_DEVICE (int n) noexcept
        {
            int n_glo = n + cumulative_columns; /*global column number*/
            ComplexType one(1., 0.);
            ComplexType minus_one(-1., 0.);
	        ComplexType imag(0., 1.);

            GR_loc(n_glo,n) =  one/(Alpha(n) - X(n) - Y(n));

            for (int m = n_glo; m > 0; m--)
            {
                GR_loc(m-1,n) =  -1*Ytil_glo(m)*GR_loc(m,n);
            }
            for (int m = n_glo; m < Hsize-1; ++m)
            {
                GR_loc(m+1,n) = -1*Xtil_glo(m)*GR_loc(m,n);
            }

            MatrixBlock<T> A_tk[NUM_CONTACTS];
            MatrixBlock<T> Gamma[NUM_CONTACTS];
            for (int m=0; m < Hsize; ++m)
            {
                A_loc(m, n) = 0.;
            }
            for (int k=0; k < NUM_CONTACTS; ++k)
            {
                int k_glo = GC_ID[k];
                MatrixBlock<T> G_contact_kk =  
                               one/(Alpha_contact(k) - X_contact(k) - Y_contact(k));

                MatrixBlock<T> temp = G_contact_kk;
                for (int m = k_glo; m < n_glo; ++m)
                {
                    temp = -1*Xtil_glo(m)*temp;
                }
                for (int m = k_glo; m > n_glo; m--)
                {
                    temp = -1*Ytil_glo(m)*temp;
                }
                MatrixBlock<T> G_contact_nk = temp;

                Gamma[k] = imag*(Sigma_contact(k) - Sigma_contact(k).Dagger());

                MatrixBlock<T> A_kn  = G_contact_kk * Gamma[k] *  G_contact_nk.Dagger();

                A_loc(k_glo, n) = A_loc(k_glo, n) + A_kn;
                for (int m = k_glo+1; m < Hsize; ++m)
                {
                    A_kn = -1*Xtil_glo(m-1)*A_kn;
                    A_loc(m,n) = A_loc(m,n) + A_kn;
                }
                for (int m = k_glo-1; m >= 0; m--)
                {
                    A_kn = -1*Ytil_glo(m+1)*A_kn;
                    A_loc(m,n) = A_loc(m,n) + A_kn;
                }
                A_tk[k] = 0.;
                if(n_glo == CT_ID[k])
                {
                    A_tk[k] = A_kn;
                }
            }

            /*LDOS*/ 
	     
            ComplexType val = A_loc(n_glo,n).DiagDotSum(degen_vec_ptr)/(2.*MathConst::pi);

            amrex::HostDevice::Atomic::Add(&(trace_r[0]), val.real());
            amrex::HostDevice::Atomic::Add(&(trace_i[0]), val.imag());

            /*Transmission*/ 
            auto T12 = Gamma[0]* A_tk[1];

            ComplexType T12_blksum = T12.DiagDotSum(degen_vec_ptr); 

            amrex::HostDevice::Atomic::Add(&(trace_r[1]), T12_blksum.real());
            amrex::HostDevice::Atomic::Add(&(trace_i[1]), T12_blksum.imag());
        }); 

        #ifdef AMREX_USE_GPU 
        amrex::Gpu::streamSynchronize();
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_Trace_r.begin(), d_Trace_r.end(), h_Trace_r.begin());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_Trace_i.begin(), d_Trace_i.end(), h_Trace_i.begin());
        #endif

        for (int t=0; t< num_traces; ++t)
        {
            amrex::ParallelDescriptor::ReduceRealSum(h_Trace_r[t]);
            amrex::ParallelDescriptor::ReduceRealSum(h_Trace_i[t]);
        }

        h_LDOS_loc(e) = spin_degen*h_Trace_r[0]/num_atoms_per_unitcell;
        h_Transmission_loc(e) = h_Trace_r[1];

        //if(e==0) 
    	//{
        //    h_GR_loc_data.copy(d_GR_loc_data); //copy from cpu to gpu 
        //    h_A_loc_data.copy(d_A_loc_data); //copy from cpu to gpu

        //    amrex::Gpu::streamSynchronize();

        //    amrex::Print() << "Printing GR_loc: \n";
        //    Print_Table2D_loc(h_GR_loc_data);

        //    amrex::Print() << "Printing A_loc: \n";
        //    Print_Table2D_loc(h_A_loc_data);

        //}
    }

    //amrex::Print() << "Printing LDOS: \n";
    //Write_Table1D(ContourPath_DOS.E_vec, 
    //              h_LDOS_loc_data, 
    //              "LDOS.dat",  "E_r LDOS_r");

    //amrex::Print() << "Printing Transmission: \n";
    //Write_Table1D(ContourPath_DOS.E_vec, 
    //              h_Transmission_loc_data, 
    //              "Transmission.dat",  "E_r T_r");

    Deallocate_TemporaryArraysForGFComputation();

}


template<typename T>
void 
c_NEGF_Common<T>:: Compute_InducedCharge ()
{

    Compute_RhoEq();

    if(flag_noneq_exists) Compute_RhoNonEq();
   

    auto const& h_Rho0_loc     = h_Rho0_loc_data.const_table();
    auto const& h_RhoEq_loc    = h_RhoEq_loc_data.const_table();
    auto const& h_RhoNonEq_loc = h_RhoNonEq_loc_data.const_table();

    RealTable1D h_RhoInduced_loc_data;
    h_RhoInduced_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    auto const& h_RhoInduced_loc = h_RhoInduced_loc_data.table();

    //amrex::Print() << "Differential Charge per Atom: \n";
    for (int n=0; n <blkCol_size_loc; ++n) 
    {
        h_RhoInduced_loc(n) = ( h_RhoEq_loc(n).DiagSum().imag() 
                              + h_RhoNonEq_loc(n).DiagSum().real() 
                             - h_Rho0_loc(n).DiagSum().imag() );

	//if(n < 5) {
        //amrex::Print() << "n/rho_Eq/Rho_NonEq/Rho0/RhoInduced: " << n << " " <<  h_RhoEq_loc(n).DiagSum().imag() 
        //        	                                                  << " " <<  h_RhoNonEq_loc(n).DiagSum().real() 
        //        							  << " " <<  h_Rho0_loc(n).DiagSum().imag() 
        //        							  << " " <<  h_RhoInduced_loc(n) << "\n";
	//}
    }
    amrex::Print() << "induced charge at site 0: " << h_RhoInduced_loc(0) << "\n";
    auto const& n_curr_out = n_curr_out_data.table();

    MPI_Allgatherv(&h_RhoInduced_loc(0),
                    blkCol_size_loc,
                    MPI_DOUBLE,
                   &n_curr_out(0),
                    MPI_recv_count.data(),
                    MPI_disp.data(),
                    MPI_DOUBLE,
                    ParallelDescriptor::Communicator());

    MPI_Barrier(ParallelDescriptor::Communicator());

    h_RhoInduced_loc_data.clear();

    //std::string filename = "Qm_out_" + std::to_string(Broyden_Step) + ".dat";
    //Write_Table1D(PTD, n_curr_out_data, filename.c_str(), 
    //              "'axial location / (nm)', 'Induced charge per site / (e)'");

}


template<typename T>
void 
c_NEGF_Common<T>:: GuessNewCharge_Broyden_FirstAlg ()
{

    amrex::Print() << "BroydenStep: " << Broyden_Step << "\n";
    auto const& n_curr_in  = h_n_curr_in_data.table();
    auto const& n_curr_out = n_curr_out_data.table();
    auto const& n_prev_in  = n_prev_in_data.table();
    auto const& F_curr     = F_curr_data.table();

    RealTable1D sum_Fcurr_data({0},{num_field_sites}, The_Pinned_Arena());
    RealTable1D sum_deltaFcurr_data({0},{num_field_sites}, The_Pinned_Arena());
    RealTable1D delta_F_curr_data({0},{num_field_sites}, The_Pinned_Arena());
    RealTable1D delta_n_Jinv_data({0},{num_field_sites}, The_Pinned_Arena());

    RealTable1D Norm_data({0},{num_field_sites}, The_Pinned_Arena());

    auto const& sum_Fcurr      = sum_Fcurr_data.table();
    auto const& sum_deltaFcurr = sum_deltaFcurr_data.table();
    auto const& delta_F_curr   = delta_F_curr_data.table();
    auto const& delta_n_Jinv   = delta_n_Jinv_data.table();
    auto const& Jinv_curr    = Jinv_curr_data.table();
    auto const& Norm     = Norm_data.table();

    amrex::Real denom = 0.;
    amrex::Real total_diff = 0.;
    int m = Broyden_Step-1;

    for(int l=0; l < num_field_sites; ++l) 
    {
        amrex::Real Fcurr = n_curr_in(l) - n_curr_out(l);

	delta_F_curr(l) = Fcurr - F_curr(l);    

        F_curr(l) = Fcurr;

	Norm(l) = fabs(Fcurr);

	total_diff += pow(Fcurr,2);

        sum_deltaFcurr(l) = 0;		 
        sum_Fcurr(l) = 0;		 
        delta_n_Jinv(l) = 0.;
    }
    total_diff = sqrt(total_diff);

    if(m > 0 ) {
        for(int a=0; a < num_field_sites; ++a) 
        {
            amrex::Real sum = 0.;		
            for(int b=0; b < num_field_sites; ++b) 
            {
                sum += Jinv_curr(a,b)*delta_F_curr(b);
            }	
            sum_deltaFcurr(a) += sum;
        }
        
        for(int l=0; l < num_field_sites; ++l) 
        {
            denom += (n_curr_in(l) - n_prev_in(l)) * sum_deltaFcurr(l);
        }

        for(int b=0; b < num_field_sites; ++b) 
        {
            amrex::Real sum = 0.;		
            for(int a=0; a < num_field_sites; ++a) 
            {
                sum += (n_curr_in(a) - n_prev_in(a))*Jinv_curr(a,b);
            }	
            delta_n_Jinv(b) += sum;
        }

        for(int a=0; a < num_field_sites; ++a) 
        {
            for(int b=0; b < num_field_sites; ++b) 
            {
                Jinv_curr(a,b) += ( (n_curr_in(a) - n_prev_in(a)) - sum_deltaFcurr(a) )  * delta_n_Jinv(b) / denom;
            }	
        }
    }
    for(int a=0; a < num_field_sites; ++a) 
    {
        amrex::Real sum = 0.;		
        for(int b=0; b < num_field_sites; ++b) 
        {
            sum += Jinv_curr(a,b)*F_curr(b);
        }	
        sum_Fcurr(a) += sum;
    }

    for(int l=0; l < num_field_sites; ++l) 
    {
	n_prev_in(l) = n_curr_in(l); 
        n_curr_in(l) = n_prev_in(l) - sum_Fcurr(l);
    }

    sum_Fcurr_data.clear();
    sum_deltaFcurr_data.clear();
    delta_F_curr_data.clear();
    delta_n_Jinv_data.clear();
    
    for(int l=0; l < 4; ++l) 
    {
        amrex::Print() << "Qm_in, Qm_out, Qm+1_in, sumFcrr, F_curr, deltaF: " << l 
		                                                          << "  "  << n_prev_in(l) 
		                                                          << "  " << n_curr_out(l)  
									  << "  " << n_curr_in(l) 
									  << "  " << sum_Fcurr(l)
									  << "  " << F_curr(l)
									  << "  " << delta_F_curr(l) 
									  << "\n";
    }


    Broyden_Norm = Norm(0);
    int norm_index = 0;
    for(int l=1; l < num_field_sites; ++l)
    {
        if(Broyden_Norm < Norm(l))
        {
            Broyden_Norm = Norm(l);
	    norm_index = l;
        }
    }

    amrex::Print() << "\nL2 norm: " << total_diff << " Max norm: " << Broyden_Norm << " location: " <<  norm_index << "\n";

    //std::string filename = "norm_" + std::to_string(Broyden_Step) + ".dat";
    //Write_Table1D(PTD, Norm_data, filename.c_str(), 
    //              "'axial location / (nm)', 'norm");
    //Norm_data.clear();

    Broyden_Step += 1;
}

template<typename T>
void 
c_NEGF_Common<T>:: GuessNewCharge_ModifiedBroydenSecondAlg ()
{
    /*update h_RhoInduced_glo*/

    amrex::Print() << "BroydenStep: " << Broyden_Step << "\n";

    auto const& n_curr_in  = h_n_curr_in_data.table();
    auto const& n_curr_out = n_curr_out_data.table();
    auto const& n_prev_in  = n_prev_in_data.table();
    auto const& F_curr     = F_curr_data.table();

    RealTable1D sum_Fcurr_data({0},{num_field_sites}, The_Pinned_Arena());
    RealTable1D sum_deltaFcurr_data({0},{num_field_sites}, The_Pinned_Arena());
    RealTable1D delta_F_curr_data({0},{num_field_sites}, The_Pinned_Arena());
    RealTable1D W_curr_data({0},{num_field_sites}, The_Pinned_Arena());
    RealTable1D V_curr_data({0},{num_field_sites}, The_Pinned_Arena());
    RealTable1D Norm_data({0},{num_field_sites}, The_Pinned_Arena());

    auto const& sum_Fcurr      = sum_Fcurr_data.table();
    auto const& sum_deltaFcurr = sum_deltaFcurr_data.table();
    auto const& delta_F_curr   = delta_F_curr_data.table();
    auto const& W_curr   = W_curr_data.table();
    auto const& V_curr   = V_curr_data.table();
    auto const& Norm     = Norm_data.table();

    amrex::Real denom = 0.;
    amrex::Real total_diff = 0.;
    int m = Broyden_Step-1;
    

    /*maintain arrays: n_curr_in, n_curr_out, n_prev_in, F_curr, W_Broyden, V_Broyden*/
    /*initialize arrays before Broyden_Step=1: n_curr_in=small number, n_prev_in=0, F_curr=0, W_Broyden=0, V_Broyden=0*/
    //std::string filename = "n_in_" + std::to_string(Broyden_Step) + ".dat";
    //Write_Table1D(PTD, h_n_curr_in_data, filename.c_str(), 
    //              "'axial location / (nm)', 'Induced charge density / (m^3)");

    //if(m<4) {
    //    for(int l=0; l < 4; ++l) 
    //    {
    //        amrex::Print() << "l, n_curr_in(l), n_curr_out(l), F_prev(l): " << l << " "  << n_curr_in(l) << "  " << n_curr_out(l) << "  "<< F_curr(l) << "\n";
    //    }    
    //}

    for(int l=0; l < num_field_sites; ++l) 
    {
        amrex::Real Fcurr = n_curr_in(l) - n_curr_out(l);

	delta_F_curr(l) = Fcurr - F_curr(l);    

        F_curr(l) = Fcurr;

	denom += pow(delta_F_curr(l),2.);

	Norm(l) = fabs(Fcurr);

	total_diff += pow(Fcurr,2);

        sum_deltaFcurr(l) = 0;		 
        sum_Fcurr(l) = 0;		 
	W_curr(l) = 0;
	V_curr(l) = 0;

    }
    total_diff = sqrt(total_diff);

    //if(m<4) {
    //    for(int l=0; l < 4; ++l) 
    //    {
    //        amrex::Print() << "l, F_curr(l), delta_F_curr(l): " << l << " "  << F_curr(l) << "  " << n_curr_out(l) << "\n";
    //    }    
    //}
    //amrex::Print() << "denom: " << denom<< "\n";

    W_Broyden.push_back(new RealTable1D({0},{num_field_sites}, The_Pinned_Arena()));
    V_Broyden.push_back(new RealTable1D({0},{num_field_sites}, The_Pinned_Arena()));

    if(m > 0) 
    {
        for(int j=1; j <= m-1; ++j) 
        {
            auto const& W_j = W_Broyden[j]->table();
            auto const& V_j = V_Broyden[j]->table();

            for(int a=0; a < num_field_sites; ++a) 
            {
                amrex::Real sum = 0.;		
                for(int b=0; b < num_field_sites; ++b) 
                {
            	sum += W_j(a)*V_j(b)*delta_F_curr(b);
                }	
                sum_deltaFcurr(a) += sum;
            }
        }

        for(int l=0; l < num_field_sites; ++l) 
        {
            amrex::Real delta_n = n_curr_in(l) - n_prev_in(l);

              V_curr(l) = delta_F_curr(l)/denom;
              W_curr(l) = -Broyden_fraction*delta_F_curr(l) + delta_n - sum_deltaFcurr(l);
        }
    
	W_Broyden[m]->copy(W_curr_data);
	V_Broyden[m]->copy(V_curr_data);
        
        auto const& W_m = W_Broyden[m]->table();
        auto const& V_m = V_Broyden[m]->table();
	//amrex::Print() << "W_curr/W_Broyden_m: " << W_curr(0) << " " << W_m(0) << "\n";
	//amrex::Print() << "V_curr/V_Broyden_m: " << V_curr(10) << " " << V_m(10) << "\n";
	//if(m-1 > 0) {
        //auto const& W_mMinus1 = W_Broyden[m-1]->table();
        //auto const& V_mMinus1 = V_Broyden[m-1]->table();
	//amrex::Print() << "W_Broyden_m-1: " << W_mMinus1(0) << "\n";
	//amrex::Print() << "V_Broyden_m-1: " << V_mMinus1(10) << "\n";
	//}

        for(int j=1; j <= m; ++j) 
        {
            auto const& W_j = W_Broyden[j]->table();
            auto const& V_j = V_Broyden[j]->table();

            for(int a=0; a < num_field_sites; ++a) 
            {
                amrex::Real sum = 0.;		
                for(int b=0; b < num_field_sites; ++b) 
                {
                    sum += W_j(a)*V_j(b)*F_curr(b);
                }	
                sum_Fcurr(a) += sum;
            }
        }
    }

    for(int l=0; l < num_field_sites; ++l) 
    {
	n_prev_in(l) = n_curr_in(l); 
        n_curr_in(l) = n_prev_in(l) - Broyden_fraction*F_curr(l) - sum_Fcurr(l);
    }

    //filename = "n_next_" + std::to_string(Broyden_Step) + ".dat";
    //Write_Table1D(PTD, h_n_curr_in_data, filename.c_str(), 
    //              "'axial location / (nm)', 'Induced charge density / (m^3)");
    
    sum_Fcurr_data.clear();
    sum_deltaFcurr_data.clear();
    delta_F_curr_data.clear();
    W_curr_data.clear();
    V_curr_data.clear();
    
    for(int l=0; l < 1; ++l) 
    {
        amrex::Print() << "Qm_in, Qm_out, Qm+1_in, F_curr: " << l << " "  << n_prev_in(l) 
		                                                  << "  " << n_curr_out(l)  
								  << "  " << n_curr_in(l) 
								  << "  " << F_curr(l)
								  << "\n";
    }
    for(int l=num_field_sites-1; l < num_field_sites; ++l) 
    {
        amrex::Print() << "Qm_in, Qm_out, Qm+1_in, F_curr: " << l << " "  << n_prev_in(l) 
		                                                  << "  " << n_curr_out(l)  
								  << "  " << n_curr_in(l) 
								  << "  " << F_curr(l)
								  << "\n";
    }


    Broyden_Norm = Norm(0);
    int norm_index = 0;
    for(int l=1; l < num_field_sites; ++l)
    {
        if(Broyden_Norm < Norm(l))
        {
            Broyden_Norm = Norm(l);
	    norm_index = l;
        }
    }

    amrex::Print() << "\nL2 norm: " << total_diff << " Max norm: " << Broyden_Norm << " location: " <<  norm_index << "\n";

    std::string filename = "norm_" + std::to_string(Broyden_Step) + ".dat";
    Write_Table1D(PTD, Norm_data, filename.c_str(), 
                  "'axial location / (nm)', 'norm");
    Norm_data.clear();

    Broyden_Step += 1;

}

template<typename T>
void 
c_NEGF_Common<T>:: GuessNewCharge_SimpleMixingAlg ()
{
    /*update h_RhoInduced_glo*/

    amrex::Print() << "BroydenStep: " << Broyden_Step << "\n";

    auto const& n_curr_in  = h_n_curr_in_data.table();
    auto const& n_curr_out = n_curr_out_data.table();
    auto const& n_prev_in  = n_prev_in_data.table();
    auto const& F_curr     = F_curr_data.table();

    RealTable1D sum_Fcurr_data({0},{num_field_sites}, The_Pinned_Arena());
    RealTable1D sum_deltaFcurr_data({0},{num_field_sites}, The_Pinned_Arena());
    RealTable1D delta_F_curr_data({0},{num_field_sites}, The_Pinned_Arena());
    RealTable1D Norm_data({0},{num_field_sites}, The_Pinned_Arena());

    auto const& sum_Fcurr      = sum_Fcurr_data.table();
    auto const& sum_deltaFcurr = sum_deltaFcurr_data.table();
    auto const& delta_F_curr   = delta_F_curr_data.table();
    auto const& Norm     = Norm_data.table();

    amrex::Real denom = 0.;
    amrex::Real total_diff = 0.;
    int m = Broyden_Step-1;
    

    for(int l=0; l < num_field_sites; ++l) 
    {
        amrex::Real Fcurr = n_curr_in(l) - n_curr_out(l);

	delta_F_curr(l) = Fcurr - F_curr(l);    

        F_curr(l) = Fcurr;

	denom += pow(delta_F_curr(l),2.);

	Norm(l) = fabs(Fcurr);

	total_diff += pow(Fcurr,2);

        sum_deltaFcurr(l) = 0;		 
        sum_Fcurr(l) = 0;		 

    }
    total_diff = sqrt(total_diff);

    for(int l=0; l < num_field_sites; ++l) 
    {
	n_prev_in(l) = n_curr_in(l); 
        n_curr_in(l) = n_prev_in(l) - Broyden_fraction*F_curr(l);
    }

    sum_Fcurr_data.clear();
    sum_deltaFcurr_data.clear();
    delta_F_curr_data.clear();
    
    for(int l=0; l < num_field_sites; ++l) 
    {
        amrex::Print() << "l, Qm_in, Qm_out, Qm+1_in, F_curr, Norm: " << l << " "  << n_prev_in(l) 
		                                                          << "  " << n_curr_out(l)  
									  << "  " << n_curr_in(l) 
									  << "  " << F_curr(l)
									  << "  " << Norm(l)
									  << "\n";
    }
    //int half_sites = int(num_field_sites/2);
    //for(int l=half_sites; l < half_sites+1; ++l) 
    //{
    //    amrex::Print() << "l, Qm_in, Qm_out, Qm+1_in, F_curr, Norm: " << l << " "  << n_prev_in(l) 
    //    	                                                          << "  " << n_curr_out(l)  
    //    								  << "  " << n_curr_in(l) 
    //    								  << "  " << F_curr(l)
    //    								  << "  " << Norm(l)
    //    								  << "\n";
    //}
    //for(int l=num_field_sites-1; l < num_field_sites; ++l) 
    //{
    //    amrex::Print() << "l, Qm_in, Qm_out, Qm+1_in, F_curr, Norm: " << l << " "  << n_prev_in(l) 
    //    	                                                          << "  " << n_curr_out(l)  
    //    								  << "  " << n_curr_in(l) 
    //    								  << "  " << F_curr(l)
    //    								  << "  " << Norm(l)
    //    								  << "\n";
    //}


    Broyden_Norm = Norm(0);
    int norm_index = 0;
    for(int l=1; l < num_field_sites; ++l)
    {
        if(Broyden_Norm < Norm(l))
        {
            Broyden_Norm = Norm(l);
	    norm_index = l;
        }
    }

    amrex::Print() << "\nL2 norm: " << total_diff << " Max norm: " << Broyden_Norm << " location: " <<  norm_index << "\n";

    std::string filename = "norm_" + std::to_string(Broyden_Step) + ".dat";
    Write_Table1D(PTD, Norm_data, filename.c_str(), 
                  "'axial location / (nm)', 'norm");
    Norm_data.clear();

    Broyden_Step += 1;

}


template<typename T>
void 
c_NEGF_Common<T>:: Compute_RhoNonEq ()
{

    auto const& h_minusHa_loc  = h_minusHa_loc_data.table();
    auto const& h_Hb_loc  = h_Hb_loc_data.table();
    auto const& h_Hc_loc  = h_Hc_loc_data.table();
    auto const& h_tau     = h_tau_glo_data.table();

    Allocate_TemporaryArraysForGFComputation();

    auto const& h_Alpha_loc = h_Alpha_loc_data.table();
    auto const& h_Alpha_glo = h_Alpha_glo_data.table();
    auto const& h_Xtil_glo = h_Xtil_glo_data.table();
    auto const& h_Ytil_glo = h_Ytil_glo_data.table();
    auto const& h_X_glo = h_X_glo_data.table();
    auto const& h_Y_glo = h_Y_glo_data.table();
    auto const& h_X_loc = h_X_loc_data.table();
    auto const& h_Y_loc = h_Y_loc_data.table();
    auto const& h_Alpha_contact = h_Alpha_contact_data.table();
    auto const& h_Y_contact = h_Y_contact_data.table();
    auto const& h_X_contact = h_X_contact_data.table();
    auto const& h_Sigma_contact = h_Sigma_contact_data.table();
    auto const& h_Fermi_contact = h_Fermi_contact_data.table();

    auto const& h_RhoNonEq_loc  = h_RhoNonEq_loc_data.table();
    ComplexType zero(0.,0.); 
    SetVal_Table1D(h_RhoNonEq_loc_data,zero);

    #ifdef AMREX_USE_GPU
    auto const& RhoNonEq_loc    = d_RhoNonEq_loc_data.table();
    auto const& GR_loc          = d_GR_loc_data.table();
    auto const& A_loc           = d_A_loc_data.table();
    /*constant references*/
    auto const& Alpha           = d_Alpha_loc_data.const_table();
    auto const& Xtil_glo        = d_Xtil_glo_data.const_table();
    auto const& Ytil_glo        = d_Ytil_glo_data.const_table();
    auto const& X               = d_X_loc_data.const_table();
    auto const& Y               = d_Y_loc_data.const_table();

    auto const& Alpha_contact   = d_Alpha_contact_data.const_table();
    auto const& X_contact       = d_X_contact_data.const_table();
    auto const& Y_contact       = d_Y_contact_data.const_table();
    auto const& Sigma_contact   = d_Sigma_contact_data.const_table();
    auto const& Fermi_contact   = d_Fermi_contact_data.const_table();

    auto& degen_vec             = block_degen_gpuvec;

    d_RhoNonEq_loc_data.copy(h_RhoNonEq_loc_data);

    #else
    auto const& RhoNonEq_loc    = h_RhoNonEq_loc_data.table();
    auto const& GR_loc          = h_GR_loc_data.table();
    auto const& A_loc           = h_A_loc_data.table();
    /*constant references*/
    auto const& Alpha           = h_Alpha_loc_data.const_table();
    auto const& Xtil_glo        = h_Xtil_glo_data.const_table();
    auto const& Ytil_glo        = h_Ytil_glo_data.const_table();
    auto const& X               = h_X_loc_data.const_table();
    auto const& Y               = h_Y_loc_data.const_table();

    auto const& Alpha_contact   = h_Alpha_contact_data.const_table();
    auto const& X_contact       = h_X_contact_data.const_table();
    auto const& Y_contact       = h_Y_contact_data.const_table();
    auto const& Sigma_contact   = h_Sigma_contact_data.const_table();
    auto const& Fermi_contact   = h_Fermi_contact_data.const_table();

    auto& degen_vec             = block_degen_vec;
    #endif
  
    for(int p=0; p < ContourPath_RhoNonEq.size(); ++p) 
    {
        for(int e=0; e < ContourPath_RhoNonEq[p].num_pts; ++e) 
        {
            
            ComplexType E = ContourPath_RhoNonEq[p].E_vec[e];
            ComplexType weight = ContourPath_RhoNonEq[p].weight_vec[e];
            ComplexType mul_factor = ContourPath_RhoNonEq[p].mul_factor_vec[e];


            for(int n=0; n<blkCol_size_loc; ++n)
            {
                h_Alpha_loc(n) = E + h_minusHa_loc(n); 
                /*+ because h_minusHa is defined previously as -(H0+U)*/
            }

            get_Sigma_at_contacts(h_Sigma_contact_data, E);

            for (int c=0; c<NUM_CONTACTS; ++c)
            {
                int n_glo = global_contact_index[c];
                int n = n_glo - vec_cumu_blkCol_size[my_rank];

                if(n_glo >= vec_cumu_blkCol_size[my_rank] && n_glo < vec_cumu_blkCol_size[my_rank+1])
                {
                    h_Alpha_loc(n) = h_Alpha_loc(n) - h_Sigma_contact(c);
                }
                h_Fermi_contact(c) = FermiFunction(E-mu_contact[c], kT_contact[c]);
            }

            /*MPI_Allgather*/
            MPI_Allgatherv(&h_Alpha_loc(0),
                            blkCol_size_loc,
                            MPI_BlkType,
                           &h_Alpha_glo(0),
                            MPI_recv_count.data(),
                            MPI_disp.data(),
                            MPI_BlkType,
                            ParallelDescriptor::Communicator());

            for (int c=0; c<NUM_CONTACTS; ++c)
            {
                int n_glo = global_contact_index[c];
                h_Alpha_contact(c) = h_Alpha_glo(n_glo);
            }

            h_Y_glo(0) = 0;
            for (int n = 1; n < Hsize_glo; ++n)
            {
            int p = (n-1)%offDiag_repeatBlkSize;
                h_Ytil_glo(n) = h_Hc_loc(p) / ( h_Alpha_glo(n-1) - h_Y_glo(n-1) );
                h_Y_glo(n) = h_Hb_loc(p) * h_Ytil_glo(n);
            }

            h_X_glo(Hsize_glo-1) = 0;
            for (int n = Hsize_glo-2; n > -1; n--)
            {
            int p = n%offDiag_repeatBlkSize;
                h_Xtil_glo(n) = h_Hb_loc(p)/(h_Alpha_glo(n+1) - h_X_glo(n+1));
                h_X_glo(n) = h_Hc_loc(p)*h_Xtil_glo(n);
            }

            for (int c = 0; c < blkCol_size_loc; ++c)
            {
                int n = c + vec_cumu_blkCol_size[my_rank];
                h_Y_loc(c) = h_Y_glo(n);
                h_X_loc(c) = h_X_glo(n);
            }
 
            for (int c = 0; c < NUM_CONTACTS; ++c)
            {
                int n = global_contact_index[c];
                h_Y_contact(c) = h_Y_glo(n);
                h_X_contact(c) = h_X_glo(n);
            }

            #ifdef AMREX_USE_GPU
            d_Alpha_loc_data.copy(h_Alpha_loc_data);
            d_Xtil_glo_data.copy(h_Xtil_glo_data);
            d_Ytil_glo_data.copy(h_Ytil_glo_data);
            d_X_loc_data.copy(h_X_loc_data);
            d_Y_loc_data.copy(h_Y_loc_data);

            d_Alpha_contact_data.copy(h_Alpha_contact_data);
            d_X_contact_data.copy(h_X_contact_data);
            d_Y_contact_data.copy(h_Y_contact_data);
            d_Sigma_contact_data.copy(h_Sigma_contact_data);
            d_Fermi_contact_data.copy(h_Fermi_contact_data);

            amrex::Gpu::streamSynchronize();
            #endif        

    	    /*following is for lambda capture*/
            int cumulative_columns = vec_cumu_blkCol_size[my_rank];
            int Hsize = Hsize_glo;  
     	    auto& GC_ID = global_contact_index;
    	    auto& CT_ID = contact_transmission_index;
            auto* degen_vec_ptr = degen_vec.dataPtr();

	    amrex::Real const_multiplier = -1*spin_degen/(2*MathConst::pi);

            amrex::ParallelFor(blkCol_size_loc, [=] AMREX_GPU_DEVICE (int n) noexcept
            {
                int n_glo = n + cumulative_columns; /*global column number*/
                ComplexType one(1., 0.);
                ComplexType minus_one(-1., 0.);
	            ComplexType imag(0., 1.);

                GR_loc(n_glo,n) =  one/(Alpha(n) - X(n) - Y(n));

                for (int m = n_glo; m > 0; m--)
                {
                    GR_loc(m-1,n) =  -1*Ytil_glo(m)*GR_loc(m,n);
                }
                for (int m = n_glo; m < Hsize-1; ++m)
                {
                    GR_loc(m+1,n) = -1*Xtil_glo(m)*GR_loc(m,n);
                }

                MatrixBlock<T> A_tk[NUM_CONTACTS];
                MatrixBlock<T> Gamma[NUM_CONTACTS];
                MatrixBlock<T> AnF_sum;
                AnF_sum = 0.;
                for (int m=0; m < Hsize; ++m)
                {
                    A_loc(m, n) = 0.;
                }
                for (int k=0; k < NUM_CONTACTS; ++k)
                {
                    int k_glo = GC_ID[k];
                    MatrixBlock<T> G_contact_kk =  
                                   one/(Alpha_contact(k) - X_contact(k) - Y_contact(k));

                    MatrixBlock<T> temp = G_contact_kk;
                    for (int m = k_glo; m < n_glo; ++m)
                    {
                        temp = -1*Xtil_glo(m)*temp;
                    }
                    for (int m = k_glo; m > n_glo; m--)
                    {
                        temp = -1*Ytil_glo(m)*temp;
                    }
                    MatrixBlock<T> G_contact_nk = temp;

                    Gamma[k] = imag*(Sigma_contact(k) - Sigma_contact(k).Dagger());

                    MatrixBlock<T> A_nn  = G_contact_nk * Gamma[k] *  G_contact_nk.Dagger();
                    MatrixBlock<T> A_kn  = G_contact_kk * Gamma[k] *  G_contact_nk.Dagger();

                    A_loc(k_glo, n) = A_loc(k_glo, n) + A_kn;
                    for (int m = k_glo+1; m < Hsize; ++m)
                    {
                        A_kn = -1*Xtil_glo(m-1)*A_kn;
                        A_loc(m,n) = A_loc(m,n) + A_kn;
                    }
                    for (int m = k_glo-1; m >= 0; m--)
                    {
                        A_kn = -1*Ytil_glo(m+1)*A_kn;
                        A_loc(m,n) = A_loc(m,n) + A_kn;
                    }
                    A_tk[k] = 0.;
                    if(n_glo == CT_ID[k])
                    {
                        A_tk[k] = A_kn;
                    }

                    AnF_sum = AnF_sum + A_nn*Fermi_contact(k);
                }

                /*RhoNonEq*/
                MatrixBlock<T> RhoNonEq_n = const_multiplier*AnF_sum*weight*mul_factor;
                RhoNonEq_loc(n) = RhoNonEq_loc(n) + RhoNonEq_n.DiagMult(degen_vec_ptr);

            }); 
            #ifdef AMREX_USE_GPU
            amrex::Gpu::streamSynchronize();
            #endif
        } 
    }

    #ifdef AMREX_USE_GPU
    h_RhoNonEq_loc_data.copy(d_RhoNonEq_loc_data); 
    amrex::Gpu::streamSynchronize();
    #endif

    //amrex::Print() << "RhoNonEq_loc: \n";
    //for (int n=0; n <blkCol_size_loc; ++n) 
    //{
    //    amrex::Print() << n << "  " <<std::setprecision(3)<< h_RhoNonEq_loc(n)  << "\n";
    //}

    Deallocate_TemporaryArraysForGFComputation();

}

template<typename T>
void 
c_NEGF_Common<T>:: Compute_RhoEq ()
{
    auto const& h_minusHa_loc  = h_minusHa_loc_data.table();
    auto const& h_Hb_loc  = h_Hb_loc_data.table();
    auto const& h_Hc_loc  = h_Hc_loc_data.table();
    auto const& h_tau     = h_tau_glo_data.table();

    Allocate_TemporaryArraysForGFComputation();

    auto const& h_Alpha_loc = h_Alpha_loc_data.table();
    auto const& h_Alpha_glo = h_Alpha_glo_data.table();
    auto const& h_Xtil_glo = h_Xtil_glo_data.table();
    auto const& h_Ytil_glo = h_Ytil_glo_data.table();
    auto const& h_X_glo = h_X_glo_data.table();
    auto const& h_Y_glo = h_Y_glo_data.table();
    auto const& h_X_loc = h_X_loc_data.table();
    auto const& h_Y_loc = h_Y_loc_data.table();
    auto const& h_Sigma_contact = h_Sigma_contact_data.table();

    auto const& h_RhoEq_loc      = h_RhoEq_loc_data.table();
    ComplexType zero(0.,0.); 
    SetVal_Table1D(h_RhoEq_loc_data,zero);

    #ifdef AMREX_USE_GPU
    auto const& RhoEq_loc        = d_RhoEq_loc_data.table();
    /*constant references*/
    auto const& Alpha           = d_Alpha_loc_data.const_table();
    auto const& Xtil_glo        = d_Xtil_glo_data.const_table();
    auto const& Ytil_glo        = d_Ytil_glo_data.const_table();
    auto const& X               = d_X_loc_data.const_table();
    auto const& Y               = d_Y_loc_data.const_table();
    auto const& Sigma_contact   = d_Sigma_contact_data.const_table();
    auto& degen_vec             = block_degen_gpuvec;

    d_RhoEq_loc_data.copy(h_RhoEq_loc_data);

    #else
    auto const& RhoEq_loc        = h_RhoEq_loc_data.table();
    /*constant references*/
    auto const& Alpha           = h_Alpha_loc_data.const_table();
    auto const& Xtil_glo        = h_Xtil_glo_data.const_table();
    auto const& Ytil_glo        = h_Ytil_glo_data.const_table();
    auto const& X               = h_X_loc_data.const_table();
    auto const& Y               = h_Y_loc_data.const_table();
    auto const& Sigma_contact   = h_Sigma_contact_data.const_table();
    auto& degen_vec             = block_degen_vec;
    #endif
  
    for(int p=0; p < ContourPath_RhoEq.size(); ++p) 
    {
        for(int e=0; e < ContourPath_RhoEq[p].num_pts; ++e) 
        {

            ComplexType E = ContourPath_RhoEq[p].E_vec[e];
            ComplexType weight = ContourPath_RhoEq[p].weight_vec[e];
            ComplexType mul_factor = ContourPath_RhoEq[p].mul_factor_vec[e];

            for(int n=0; n<blkCol_size_loc; ++n)
            {
                h_Alpha_loc(n) = E + h_minusHa_loc(n); 
                /*+ because -Ha = -H0a (=0) -U = -U as defined previously*/
            }

            get_Sigma_at_contacts(h_Sigma_contact_data, E);

            for (int c=0; c<NUM_CONTACTS; ++c)
            {
                int n_glo = global_contact_index[c];

                if(n_glo >= vec_cumu_blkCol_size[my_rank] && n_glo < vec_cumu_blkCol_size[my_rank+1])
                {
                    int n = n_glo - vec_cumu_blkCol_size[my_rank];
                    h_Alpha_loc(n) = h_Alpha_loc(n) - h_Sigma_contact(c);
                }
		//if(e==0) {
                //     amrex::Print() << "contact, Sigma, E, U, EmU: " << c << "  " << h_Sigma_contact(c) << "  " << E << "  " << U_contact[c] 
		//	                                             << "  " << E - U_contact[c]<< "\n";
		//}
            }

            /*MPI_Allgather*/
            MPI_Allgatherv(&h_Alpha_loc(0),
                            blkCol_size_loc,
                            MPI_BlkType,
                           &h_Alpha_glo(0),
                            MPI_recv_count.data(),
                            MPI_disp.data(),
                            MPI_BlkType,
                            ParallelDescriptor::Communicator());

            h_Y_glo(0) = 0;
            for (int n = 1; n < Hsize_glo; ++n)
            {
            int p = (n-1)%offDiag_repeatBlkSize;
                h_Ytil_glo(n) = h_Hc_loc(p) / ( h_Alpha_glo(n-1) - h_Y_glo(n-1) );
                h_Y_glo(n) = h_Hb_loc(p) * h_Ytil_glo(n);
            }

            h_X_glo(Hsize_glo-1) = 0;
            for (int n = Hsize_glo-2; n > -1; n--)
            {
            int p = n%offDiag_repeatBlkSize;
                h_Xtil_glo(n) = h_Hb_loc(p)/(h_Alpha_glo(n+1) - h_X_glo(n+1));
                h_X_glo(n) = h_Hc_loc(p)*h_Xtil_glo(n);
            }

            for (int c = 0; c < blkCol_size_loc; ++c)
            {
                int n = c + vec_cumu_blkCol_size[my_rank];
                h_Y_loc(c) = h_Y_glo(n);
                h_X_loc(c) = h_X_glo(n);
            }
 
            #ifdef AMREX_USE_GPU
            d_Alpha_loc_data.copy(h_Alpha_loc_data);
            d_X_loc_data.copy(h_X_loc_data);
            d_Y_loc_data.copy(h_Y_loc_data);

            amrex::Gpu::streamSynchronize();
            #endif        

	    /*following is for lambda capture*/
            int cumulative_columns = vec_cumu_blkCol_size[my_rank];
            auto* degen_vec_ptr = degen_vec.dataPtr();
	    //amrex::Print() << "mu_min: " << mu_min << "\n";
            ComplexType nF_eq = FermiFunction(E-mu_min, kT_min);
	    amrex::Real const_multiplier = -1.*spin_degen/MathConst::pi;

            amrex::ParallelFor(blkCol_size_loc, [=] AMREX_GPU_DEVICE (int n) noexcept
            {
                ComplexType one(1., 0.);
                MatrixBlock<T> G_nn =  one/(Alpha(n) - X(n) - Y(n));
                MatrixBlock<T> RhoEq_n = const_multiplier*G_nn*weight*mul_factor*nF_eq;
                RhoEq_loc(n) = RhoEq_loc(n) + RhoEq_n.DiagMult(degen_vec_ptr);

            }); 
            #ifdef AMREX_USE_GPU
            amrex::Gpu::streamSynchronize();
            #endif

        } /*Energy loop*/
    } /*Path loop*/
 
    #ifdef AMREX_USE_GPU
    h_RhoEq_loc_data.copy(d_RhoEq_loc_data); 
    amrex::Gpu::streamSynchronize();
    #endif

    Deallocate_TemporaryArraysForGFComputation();

    //amrex::Print() << "RhoEq_loc: \n";
    if(!flag_noneq_exists) 
    {
        Compute_GR_atPoles();
        auto const& h_GR_atPoles_loc = h_GR_atPoles_loc_data.const_table();

        for (int n=0; n <blkCol_size_loc; ++n) 
        {
            //if(n < 5) {		
            //   amrex::Print() << "RhoEq and GR: " << n << "  " <<std::setprecision(3)<< h_RhoEq_loc(n) << "  " << h_GR_atPoles_loc(n)  << "\n";
	    //}
            h_RhoEq_loc(n) = h_RhoEq_loc(n) + h_GR_atPoles_loc(n);
        }
    }

}



template<typename T>
void 
c_NEGF_Common<T>:: Compute_GR_atPoles ()
{
    auto const& h_minusHa_loc  = h_minusHa_loc_data.table();
    auto const& h_Hb_loc  = h_Hb_loc_data.table();
    auto const& h_Hc_loc  = h_Hc_loc_data.table();
    auto const& h_tau     = h_tau_glo_data.table();

    Allocate_TemporaryArraysForGFComputation();

    auto const& h_Alpha_loc = h_Alpha_loc_data.table();
    auto const& h_Alpha_glo = h_Alpha_glo_data.table();
    auto const& h_Xtil_glo = h_Xtil_glo_data.table();
    auto const& h_Ytil_glo = h_Ytil_glo_data.table();
    auto const& h_X_glo = h_X_glo_data.table();
    auto const& h_Y_glo = h_Y_glo_data.table();
    auto const& h_X_loc = h_X_loc_data.table();
    auto const& h_Y_loc = h_Y_loc_data.table();
    auto const& h_Sigma_contact = h_Sigma_contact_data.table();

    auto const& h_GR_atPoles_loc      = h_GR_atPoles_loc_data.table();

    ComplexType zero(0.,0.); 
    SetVal_Table1D(h_GR_atPoles_loc_data,zero);

    #ifdef AMREX_USE_GPU
    auto const& GR_atPoles_loc        = d_GR_atPoles_loc_data.table();
    /*constant references*/
    auto const& Alpha           = d_Alpha_loc_data.const_table();
    auto const& Xtil_glo        = d_Xtil_glo_data.const_table();
    auto const& Ytil_glo        = d_Ytil_glo_data.const_table();
    auto const& X               = d_X_loc_data.const_table();
    auto const& Y               = d_Y_loc_data.const_table();
    auto const& Sigma_contact   = d_Sigma_contact_data.const_table();
    auto& degen_vec             = block_degen_gpuvec;

    d_GR_atPoles_loc_data.copy(h_GR_atPoles_loc_data);
    #else
    auto const& GR_atPoles_loc        = h_GR_atPoles_loc_data.table();
    /*constant references*/
    auto const& Alpha           = h_Alpha_loc_data.const_table();
    auto const& Xtil_glo        = h_Xtil_glo_data.const_table();
    auto const& Ytil_glo        = h_Ytil_glo_data.const_table();
    auto const& X               = h_X_loc_data.const_table();
    auto const& Y               = h_Y_loc_data.const_table();
    auto const& Sigma_contact   = h_Sigma_contact_data.const_table();
    auto& degen_vec             = block_degen_vec;
    #endif
  
    for(int e=0; e < E_poles_vec.size(); ++e) 
    {

        ComplexType E = E_poles_vec[e];
	//amrex::Print() << "e, E: " << e << " " << E << "\n";

        for(int n=0; n<blkCol_size_loc; ++n)
        {
            h_Alpha_loc(n) = E + h_minusHa_loc(n); 
            /*+ because h_minusHa is defined previously as -(H0+U)*/
        }

        get_Sigma_at_contacts(h_Sigma_contact_data, E);

        for (int c=0; c<NUM_CONTACTS; ++c)
        {
            int n_glo = global_contact_index[c];

            if(n_glo >= vec_cumu_blkCol_size[my_rank] && n_glo < vec_cumu_blkCol_size[my_rank+1])
            {
                int n = n_glo - vec_cumu_blkCol_size[my_rank];
                h_Alpha_loc(n) = h_Alpha_loc(n) - h_Sigma_contact(c);
            }
        }

        /*MPI_Allgather*/
        MPI_Allgatherv(&h_Alpha_loc(0),
                        blkCol_size_loc,
                        MPI_BlkType,
                       &h_Alpha_glo(0),
                        MPI_recv_count.data(),
                        MPI_disp.data(),
                        MPI_BlkType,
                        ParallelDescriptor::Communicator());

        h_Y_glo(0) = 0;
        for (int n = 1; n < Hsize_glo; ++n)
        {
        int p = (n-1)%offDiag_repeatBlkSize;
            h_Ytil_glo(n) = h_Hc_loc(p) / ( h_Alpha_glo(n-1) - h_Y_glo(n-1) );
            h_Y_glo(n) = h_Hb_loc(p) * h_Ytil_glo(n);
        }

        h_X_glo(Hsize_glo-1) = 0;
        for (int n = Hsize_glo-2; n > -1; n--)
        {
        int p = n%offDiag_repeatBlkSize;
            h_Xtil_glo(n) = h_Hb_loc(p)/(h_Alpha_glo(n+1) - h_X_glo(n+1));
            h_X_glo(n) = h_Hc_loc(p)*h_Xtil_glo(n);
        }

        for (int c = 0; c < blkCol_size_loc; ++c)
        {
            int n = c + vec_cumu_blkCol_size[my_rank];
            h_Y_loc(c) = h_Y_glo(n);
            h_X_loc(c) = h_X_glo(n);
        }
 
        #ifdef AMREX_USE_GPU
        d_Alpha_loc_data.copy(h_Alpha_loc_data);
        d_X_loc_data.copy(h_X_loc_data);
        d_Y_loc_data.copy(h_Y_loc_data);

        amrex::Gpu::streamSynchronize();
        #endif        

	/*following is for lambda capture*/
        int cumulative_columns = vec_cumu_blkCol_size[my_rank];

        ComplexType pole_const(0., -2*kT_min*spin_degen);
        auto* degen_vec_ptr = degen_vec.dataPtr();

        amrex::ParallelFor(blkCol_size_loc, [=] AMREX_GPU_DEVICE (int n) noexcept
        {
            ComplexType one(1., 0.);
            MatrixBlock<T> G_nn =  one/(Alpha(n) - X(n) - Y(n));
            MatrixBlock<T> GR_atPoles_n  = pole_const*G_nn;

	    GR_atPoles_loc(n) = GR_atPoles_loc(n) + GR_atPoles_n.DiagMult(degen_vec_ptr);

        }); 
        #ifdef AMREX_USE_GPU
        amrex::Gpu::streamSynchronize();
        #endif

    } /*Energy loop*/
 
    #ifdef AMREX_USE_GPU
    h_GR_atPoles_loc_data.copy(d_GR_atPoles_loc_data); 
    amrex::Gpu::streamSynchronize();
    #endif

    Deallocate_TemporaryArraysForGFComputation();

}

template<typename T>
void 
c_NEGF_Common<T>:: Compute_Rho0 ()
{
    auto const& h_minusHa_loc  = h_minusHa_loc_data.table();
    auto const& h_Hb_loc  = h_Hb_loc_data.table();
    auto const& h_Hc_loc  = h_Hc_loc_data.table();
    auto const& h_tau     = h_tau_glo_data.table();

    Allocate_TemporaryArraysForGFComputation();

    auto const& h_Alpha_loc = h_Alpha_loc_data.table();
    auto const& h_Alpha_glo = h_Alpha_glo_data.table();
    auto const& h_Xtil_glo = h_Xtil_glo_data.table();
    auto const& h_Ytil_glo = h_Ytil_glo_data.table();
    auto const& h_X_glo = h_X_glo_data.table();
    auto const& h_Y_glo = h_Y_glo_data.table();
    auto const& h_X_loc = h_X_loc_data.table();
    auto const& h_Y_loc = h_Y_loc_data.table();
    auto const& h_Sigma_contact = h_Sigma_contact_data.table();

    auto const& h_Rho0_loc      = h_Rho0_loc_data.table();

    ComplexType zero(0.,0.); 
    SetVal_Table1D(h_Rho0_loc_data,zero);

    #ifdef AMREX_USE_GPU
    auto const& Rho0_loc        = d_Rho0_loc_data.table();
    /*constant references*/
    auto const& Alpha           = d_Alpha_loc_data.const_table();
    auto const& Xtil_glo        = d_Xtil_glo_data.const_table();
    auto const& Ytil_glo        = d_Ytil_glo_data.const_table();
    auto const& X               = d_X_loc_data.const_table();
    auto const& Y               = d_Y_loc_data.const_table();
    auto const& Sigma_contact   = d_Sigma_contact_data.const_table();
    auto& degen_vec             = block_degen_gpuvec;

    d_Rho0_loc_data.copy(h_Rho0_loc_data);

    #else
    auto const& Rho0_loc        = h_Rho0_loc_data.table();
    /*constant references*/
    auto const& Alpha           = h_Alpha_loc_data.const_table();
    auto const& Xtil_glo        = h_Xtil_glo_data.const_table();
    auto const& Ytil_glo        = h_Ytil_glo_data.const_table();
    auto const& X               = h_X_loc_data.const_table();
    auto const& Y               = h_Y_loc_data.const_table();
    auto const& Sigma_contact   = h_Sigma_contact_data.const_table();
    auto& degen_vec             = block_degen_vec;
    #endif
  
    for(int p=0; p < ContourPath_Rho0.size(); ++p) 
    {
        for(int e=0; e < ContourPath_Rho0[p].num_pts; ++e) 
        {

            ComplexType E = ContourPath_Rho0[p].E_vec[e];
            ComplexType weight = ContourPath_Rho0[p].weight_vec[e];
            ComplexType mul_factor = ContourPath_Rho0[p].mul_factor_vec[e];

            for(int n=0; n<blkCol_size_loc; ++n)
            {
                h_Alpha_loc(n) = E + h_minusHa_loc(n); 
                /*+ because h_minusHa is defined previously as -(H0+U)*/
            }

            get_Sigma_at_contacts(h_Sigma_contact_data, E);

            for (int c=0; c<NUM_CONTACTS; ++c)
            {
                int n_glo = global_contact_index[c];

                if(n_glo >= vec_cumu_blkCol_size[my_rank] && n_glo < vec_cumu_blkCol_size[my_rank+1])
                {
                    int n = n_glo - vec_cumu_blkCol_size[my_rank];
                    h_Alpha_loc(n) = h_Alpha_loc(n) - h_Sigma_contact(c);
                }
            }

            /*MPI_Allgather*/
            MPI_Allgatherv(&h_Alpha_loc(0),
                            blkCol_size_loc,
                            MPI_BlkType,
                           &h_Alpha_glo(0),
                            MPI_recv_count.data(),
                            MPI_disp.data(),
                            MPI_BlkType,
                            ParallelDescriptor::Communicator());

            h_Y_glo(0) = 0;
            for (int n = 1; n < Hsize_glo; ++n)
            {
            int p = (n-1)%offDiag_repeatBlkSize;
                h_Ytil_glo(n) = h_Hc_loc(p) / ( h_Alpha_glo(n-1) - h_Y_glo(n-1) );
                h_Y_glo(n) = h_Hb_loc(p) * h_Ytil_glo(n);
            }

            h_X_glo(Hsize_glo-1) = 0;
            for (int n = Hsize_glo-2; n > -1; n--)
            {
            int p = n%offDiag_repeatBlkSize;
                h_Xtil_glo(n) = h_Hb_loc(p)/(h_Alpha_glo(n+1) - h_X_glo(n+1));
                h_X_glo(n) = h_Hc_loc(p)*h_Xtil_glo(n);
            }

            for (int c = 0; c < blkCol_size_loc; ++c)
            {
                int n = c + vec_cumu_blkCol_size[my_rank];
                h_Y_loc(c) = h_Y_glo(n);
                h_X_loc(c) = h_X_glo(n);
            }
 
            #ifdef AMREX_USE_GPU
            d_Alpha_loc_data.copy(h_Alpha_loc_data);
            d_X_loc_data.copy(h_X_loc_data);
            d_Y_loc_data.copy(h_Y_loc_data);
            amrex::Gpu::streamSynchronize();
            #endif        

	    /*following is for lambda capture*/
            int cumulative_columns = vec_cumu_blkCol_size[my_rank];
            auto* degen_vec_ptr = degen_vec.dataPtr();
	    amrex::Real const_multiplier = -1*spin_degen/(MathConst::pi);

            amrex::ParallelFor(blkCol_size_loc, [=] AMREX_GPU_DEVICE (int n) noexcept
            {
                ComplexType one(1., 0.);
                MatrixBlock<T> G_nn =  one/(Alpha(n) - X(n) - Y(n));
                MatrixBlock<T> Rho0_n = const_multiplier*G_nn*weight*mul_factor;
                Rho0_loc(n) = Rho0_loc(n) + Rho0_n.DiagMult(degen_vec_ptr);

            }); 
            #ifdef AMREX_USE_GPU
            amrex::Gpu::streamSynchronize();
	    #endif 
        } /*Energy loop*/
    } /*Path loop*/
 
    #ifdef AMREX_USE_GPU
    h_Rho0_loc_data.copy(d_Rho0_loc_data); 
    amrex::Gpu::streamSynchronize();
    #endif

    //amrex::Print() << "Rho0_loc: \n";
    //for (int n=0; n <blkCol_size_loc; ++n) 
    //{
    //   amrex::Print() << n << "  " <<std::setprecision(4) << h_Rho0_loc(n)  << "\n";
    //}

    Deallocate_TemporaryArraysForGFComputation();

}

template<typename T>
AMREX_GPU_HOST_DEVICE
void 
c_NEGF_Common<T>:: Compute_SurfaceGreensFunction (MatrixBlock<T>& gr, const ComplexType EmU) {}


template<typename T>
void 
c_NEGF_Common<T>:: get_Sigma_at_contacts 
                              (BlkTable1D& h_Sigma_contact_data, 
                               ComplexType E)
{

    auto const& h_tau   = h_tau_glo_data.table();
    auto const& h_Sigma = h_Sigma_contact_data.table();
    
    for (std::size_t c = 0; c < NUM_CONTACTS; ++c)
    {
        MatrixBlock<T> gr;
        Compute_SurfaceGreensFunction(gr, E - U_contact[c]);
        //amrex::Print() << "c, E, gr: " << c << " " << E << " " << gr << "\n";
        h_Sigma(c) = h_tau(c)*gr*h_tau(c).Dagger();
    }

}


//template<typename T>
//AMREX_GPU_HOST_DEVICE
//MatrixBlock<T>
//c_NEGF_Common<T>::get_Gamma(const MatrixBlock<T>& Sigma)
//{
//    //ComplexType imag(0., 1.);
//    //return imag*(Sigma - Sigma.Dagger());
//    return Sigma;
//}


template<typename T>
void
c_NEGF_Common<T>::Write_BlkTable1D_asaf_E(const amrex::Vector<ComplexType>& E_vec,
                                          const BlkTable1D& Arr_data, 
                                                std::string filename, 
                                                std::string header)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "Root Writing " << filename << "\n";
        std::ofstream outfile;
        outfile.open(filename.c_str());

        auto const& Arr = Arr_data.const_table();
        auto thi = Arr_data.hi();

        outfile << header  << "\n";
        if(E_vec.size() == thi[0]) {   
            for (int e=0; e< thi[0]; ++e)
            {
                outfile << std::setw(10) << E_vec[e]
                        << std::setw(15) << Arr(e) << "\n";
            }
        }
        else {
            outfile << "Mismatch in the size of E_vec and Arr_data!"  << "\n";
        }
        outfile.close();
    }
}


template<typename T>
template<typename VectorType, typename TableType>
void
c_NEGF_Common<T>::Write_Table1D(const amrex::Vector<VectorType>& Vec,
                                const TableType& Arr_data, 
                                std::string filename, 
                                std::string header)
{ 
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\nRoot Writing " << filename << "\n";
        std::ofstream outfile;
        outfile.open(filename.c_str());

        auto const& Arr = Arr_data.const_table();
        auto thi = Arr_data.hi();

        outfile << header  << "\n";
        if(Vec.size() == thi[0]) {   
            for (int e=0; e< thi[0]; ++e)
            {
                outfile << std::setw(10) << Vec[e] 
                        << std::setw(15) << Arr(e) << "\n";
            }
        }
        else {
            outfile << "Mismatch in the size of Vec and Table1D_data!"  << "\n";
        }
        outfile.close();
    }
}


template<typename T>
void
c_NEGF_Common<T>::Write_BlkTable2D_asaf_E(const BlkTable2D& Arr_data, 
                                                std::string filename, 
                                                std::string header)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\n Root Writing " << filename << "\n";
        std::ofstream outfile;
        outfile.open(filename.c_str());

        auto const& E = h_E_RealPath_data.const_table();
        auto thi = h_E_RealPath_data.hi();
        auto const& Arr = Arr_data.const_table();

        outfile << header << "\n";
        for (int e=0; e< thi[0]; ++e)
        {
            outfile << std::setw(10) << E(e)
                    << std::setw(15) << Arr(0,e) 
                    << std::setw(15) << Arr(1,e) << "\n";
        }

        outfile.close();
    }
}


template<typename T>
template<typename U>
void
c_NEGF_Common<T>::Print_Table1D_loc(const U& Tab1D_data)
{
    auto tlo = Tab1D_data.lo();
    auto thi = Tab1D_data.hi();

    auto const& Tab1D = Tab1D_data.table();

    for (int i = tlo[0]; i < thi[0]; ++i)
    {   
        std::cout << std::setw(15) << std::setprecision(2) << Tab1D(i);
        std::cout << "\n";
    }
}




template<typename T>
template<typename U>
void
c_NEGF_Common<T>::Print_Table2D_loc(const U& Tab2D_data)
{
    auto tlo = Tab2D_data.lo();
    auto thi = Tab2D_data.hi();

    auto const& Tab2D = Tab2D_data.table();

    std::cout << "[\n";
    for (int i = tlo[0]; i < thi[0]; ++i)
    {   
        for (int j = tlo[1]; j < thi[1]; ++j) //slow moving index. printing slow
        {       
            std::cout << std::setw(5) << std::setprecision(2) << Tab2D(i,j);
        }
        std::cout << "\n";
    }
    std::cout << "]\n";
}

//template<typename T>
//void 
//c_NEGF_Common<T>::DeallocateArrays () 
//{
//
//}
//
//
//template<typename T>
//void 
//c_NEGF_Common<T>::ComputeChargeDensity () 
//{
//
//}


//
//template<typename T>
//AMREX_GPU_HOST_DEVICE
//ComplexType 
//c_NEGF_Common<T>::conjugate(ComplexType a)
//{
//   ComplexType a_conj(a.real(), -1.*a.imag());
//   return a_conj;
//}

