#include "NEGF_Common.H"
#include "Matrix_Block_Util.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"


/*Explicit specializations*/
template class c_NEGF_Common<ComplexType[NUM_MODES]>;            //of c_CNT
template class c_NEGF_Common<ComplexType[NUM_MODES][NUM_MODES]>; //c_Graphene

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
    }
    
    set_material_specific_parameters();

    amrex::Print() << "#####* num_atoms: " << num_atoms << "\n";
    amrex::Print() << "#####* num_atoms_per_unitcell: " << num_atoms_per_unitcell << "\n";

}


template<typename T>
void 
c_NEGF_Common<T>::set_material_specific_parameters() 
{
    /*set the following in the specialization*/
    /*num_atoms: number of total atoms*/
    /*num_atoms_per_unitcell: number of atoms per unitcell*/
    /*block_degen_vec: this is vector of degeneracy factors of size equal to that of a block element*/ 
    /*block_degen_gpuvec: device vector copy of block_degen_vec*/
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
       amrex::Print() << "p,recv, disp: " << p << "  " 
                      << MPI_recv_count[p] << "  " 
                      << MPI_disp[p] << "\n";
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
    h_Ha_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    SetVal_Table1D(h_Ha_loc_data,zero);

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
    #else
    h_GR_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    SetVal_Table2D(h_GR_loc_data, zero);

    h_A_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    SetVal_Table2D(h_A_loc_data, zero);
    #endif

    h_Rho0_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    SetVal_Table1D(h_Rho0_loc_data,zero);

    h_E_RealPath_data.resize({0},{NUM_ENERGY_PTS_REAL},The_Pinned_Arena());

}


template<typename T>
void 
c_NEGF_Common<T>:: ConstructHamiltonian () {}


template<typename T>
void 
c_NEGF_Common<T>:: AddPotentialToHamiltonian () 
{
    auto const& h_Ha = h_Ha_loc_data.table();
    int c=0;
    for(auto& col_gid: vec_blkCol_gids)
    {
        h_Ha(c) = h_Ha(c) - avg_gatherField[col_gid];
        ++c;
    }
}


template<typename T>
void 
c_NEGF_Common<T>:: Update_ContactPotential () 
{
    for(int c=0; c<NUM_CONTACTS; ++c)
    {
        U_contact[c] = avg_gatherField[global_contact_index[c]];
        ++c;
    }
}


template<typename T>
void 
c_NEGF_Common<T>:: Define_EnergyLimits ()
{
    for (int c=0; c<NUM_CONTACTS; ++c)
    {
        mu_contact[c] = E_f + U_contact[c];
        kT_contact[c] = PhysConst::kb_eVperK*298.; /*set in the input*/
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
           flag_noneq_exists=true;
       }
       if(mu_max < mu_contact[c])
       {
           mu_max = mu_contact[c];
           flag_noneq_exists=true;
       }
       if(kT_min > kT_contact[c])
       {
           kT_min = kT_contact[c];
           flag_noneq_exists=true;
       }
       if(kT_max < kT_contact[c])
       {
           kT_max = kT_contact[c];
           flag_noneq_exists=true;
       }
    }

    ComplexType val(0.,1e-8);
    E_zPlus = val;
    E_contour_left  = E_valence_min + E_zPlus; /*set in the input*/
    E_rightmost = mu_max + Fermi_tail_factor*kT_max + E_zPlus;

    if(flag_noneq_exists)
    {
        amrex::Print() << "\nnonequilibrium exists between terminals\n";
        E_contour_right = mu_min - Fermi_tail_factor*kT_max + E_zPlus;
    }
    else E_contour_right = E_rightmost;

    num_enclosed_poles = int((E_pole_max-MathConst::pi*kT_max)/(2.*MathConst::pi*kT_max) + 1);
    ComplexType val2(0.,2*num_enclosed_poles*MathConst::pi*kT_max);
    E_zeta = val2;
    E_eta =  E_contour_right.real() - Fermi_tail_factor*kT_max + E_zeta;

    amrex::Print() << "\nE_f: " << E_f << "\n";
    amrex::Print() << "U_contact: ";
    for (int c=0; c<NUM_CONTACTS; ++c)
    {
        amrex::Print() <<  U_contact[c] << " ";
    }
    amrex::Print() << "\n";
    amrex::Print() << "mu_min/max: " << mu_min << " " << mu_max << "\n";
    amrex::Print() << "kT_min/max: " << kT_min << " " << kT_max << "\n";
    amrex::Print() << "E_zPlus: "  << E_zPlus << "\n";
    amrex::Print() << "E_contour_left/E_contour_right/E_rightmost: " << E_contour_left      <<  "  "
                                                                 << E_contour_right << "  "
                                                                 << E_rightmost     << "\n";
    amrex::Print() << "E_pole_max: " << E_pole_max << ", number of poles: " << num_enclosed_poles << "\n";
    amrex::Print() << "E_zeta: " << E_zeta << "\n";
    amrex::Print() << "E_eta: "  << E_eta << "\n";

}



template<typename T>
void 
c_NEGF_Common<T>:: Define_IntegrationPaths ()
{
    /* Define_ContourPath_Rho0 */
    //ContourPath_Rho0[0].Define_GaussLegendrePoints(E_zPlus, E_zeta, 50, "line"); 
    //ContourPath_Rho0[1].Define_GaussLegendrePoints(E_zeta,  E_eta, 50, "line"); 
    //ContourPath_Rho0[2].Define_GaussLegendrePoints(E_eta, E_contour_left, 50, "circle"); 

    //ContourPath_RhoEq_Direct.Define_RegularPoints(E_contour_left, E_contour_right, 2000);
    
    //ContourPath_RhoEq_Direct.Define_RegularPoints(E_contour_left, E_zPlus, 20000);
    ComplexType minus_one(-1.,0.);
    ComplexType one(1.,0.);
    ContourPath_RhoEq_Direct.Define_RegularPoints(minus_one + E_zPlus, one + E_zPlus, 200);

    //for(int e=0; e < ContourPath_RhoEq_Direct.num_pts; ++e) 
    //{
    //   amrex::Print() << e << std::setprecision(4)
    //                       << "  " << ContourPath_RhoEq_Direct.E_vec[e]
    //                       << "  " << ContourPath_RhoEq_Direct.weight_vec[e] << "\n" ;
    //}
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

    h_Alpha_contact_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    SetVal_Table1D(h_Alpha_contact_data,zero);

    h_X_contact_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    SetVal_Table1D(h_X_contact_data,zero);

    h_Y_contact_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    SetVal_Table1D(h_Y_contact_data,zero);

    h_Sigma_contact_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    SetVal_Table1D(h_Sigma_contact_data,zero);

    //h_Gamma_contact_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    //SetVal_Table1D(h_Gamma_contact_data,zero);

    h_Trace_r.resize(num_traces);
    h_Trace_i.resize(num_traces);

    #ifdef AMREX_USE_GPU
    d_Alpha_loc_data.resize({0},{blkCol_size_loc}, The_Arena());
    d_X_loc_data.resize({0},{blkCol_size_loc}, The_Arena());
    d_Y_loc_data.resize({0},{blkCol_size_loc}, The_Arena());
    d_Xtil_glo_data.resize({0},{Hsize_glo}, The_Arena());
    d_Ytil_glo_data.resize({0},{Hsize_glo}, The_Arena());

    d_Alpha_contact_data.resize({0},{NUM_CONTACTS}, The_Arena());
    d_X_contact_data.resize({0},{NUM_CONTACTS}, The_Arena());
    d_Y_contact_data.resize({0},{NUM_CONTACTS}, The_Arena());
    d_Sigma_contact_data.resize({0},{NUM_CONTACTS}, The_Arena());

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

    h_Alpha_contact_data.clear();
    h_X_contact_data.clear();
    h_Y_contact_data.clear();
    h_Sigma_contact_data.clear();

    h_Trace_r.clear();
    h_Trace_i.clear();

    #ifdef AMREX_USE_GPU
    d_Alpha_loc_data.clear();
    d_X_loc_data.clear();
    d_Y_loc_data.clear();
    d_Xtil_glo_data.clear();
    d_Ytil_glo_data.clear();

    d_Alpha_contact_data.clear();
    d_X_contact_data.clear();
    d_Y_contact_data.clear();
    d_Sigma_contact_data.clear();

    d_Trace_r.clear();
    d_Trace_i.clear();
    #endif 
}


template<typename T>
void 
c_NEGF_Common<T>:: Compute_DensityOfStates ()
{
    int E_pts = ContourPath_RhoEq_Direct.num_pts;
    h_LDOS_loc_data.resize({0}, {E_pts}, The_Pinned_Arena());
    SetVal_Table1D(h_LDOS_loc_data,0.);

    h_Transmission_loc_data.resize({0}, {E_pts}, The_Pinned_Arena());
    SetVal_Table1D(h_Transmission_loc_data,0.);

    //#ifdef AMREX_USE_GPU
    //#else
    //auto const& Rho0_loc   = h_Rho0_loc_data.table();
    //#endif 

    auto const& h_Ha_loc  = h_Ha_loc_data.table();
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
    auto const& Alpha         = d_Alpha_loc_data.const_table();
    auto const& Xtil_glo      = d_Xtil_glo_data.const_table();
    auto const& Ytil_glo      = d_Ytil_glo_data.const_table();
    auto const& X             = d_X_loc_data.const_table();
    auto const& Y             = d_Y_loc_data.const_table();

    auto const& Alpha_contact = d_Alpha_contact_data.const_table();
    auto const& X_contact     = d_X_contact_data.const_table();
    auto const& Y_contact     = d_Y_contact_data.const_table();
    auto const& Sigma_contact = d_Sigma_contact_data.const_table();

    auto* trace_r             = d_Trace_r.dataPtr();
    auto* trace_i             = d_Trace_i.dataPtr();
    auto& degen_vec           = block_degen_gpuvec;
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
  

    for(int e=0; e < ContourPath_RhoEq_Direct.num_pts; ++e) 
    {

        ComplexType E = ContourPath_RhoEq_Direct.E_vec[e];
        ComplexType weight = ContourPath_RhoEq_Direct.weight_vec[e];
        ComplexType mul_factor = ContourPath_RhoEq_Direct.mul_factor_vec[e];


        for(int n=0; n<blkCol_size_loc; ++n)
        {
            h_Alpha_loc(n) = E + h_Ha_loc(n); 
            /*+ because h_Ha is defined previously as -(H0+U)*/
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
        //for(int n=0; n<blkCol_size_loc; ++n)
        //{
        //    ComplexType small(0., 1e-8);
        //    if(n ==0 or n==3) {
        //      h_Alpha_loc(n) = -0.1432818278 + small;
        //    } 
        //    else if (n == 1 or n==2){
        //      h_Alpha_loc(n) = -0.4137707681 + small;
        //    } 
        //}
        //amrex::Print() << "Alpha_loc: " << "\n";
        //for(int n=0; n< blkCol_size_loc; ++n) 
        //{
        //    amrex::Print() << n << " "<< h_Alpha_loc(n) << "\n";
        //}

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
        auto degen_vec_ptr = degen_vec.dataPtr();

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

//            MatrixBlock<T> G_nn =  one/(Alpha(n) - X(n) - Y(n));
//            Rho0_loc(n) = Rho0_loc(n) + G_nn*weight*mul_factor;

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

    //amrex::Print() << "Rho0_loc: \n";
    //for (int n=0; n <blkCol_size_loc; ++n) 
    //{
    //    Rho0_loc(n) = Rho0_loc(n)/(2.*MathConst::pi);
    //    amrex::Print() << n << "  " <<std::setprecision(3)<< Rho0_loc(n)  << "\n";
    //}

    amrex::Print() << "Printing LDOS: \n";
    Write_Table1D_asaf_E(ContourPath_RhoEq_Direct.E_vec, 
                         h_LDOS_loc_data, 
                        "LDOS.dat",  "E_r LDOS_r");

    amrex::Print() << "Printing Transmission: \n";
    Write_Table1D_asaf_E(ContourPath_RhoEq_Direct.E_vec, 
                         h_Transmission_loc_data, 
                        "Transmission.dat",  "E_r T_r");

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
        Compute_SurfaceGreensFunction(gr, E-U_contact[c]);
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
        amrex::Print() << "\n Root Writing " << filename << "\n";
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
template<typename U>
void
c_NEGF_Common<T>::Write_Table1D_asaf_E(const amrex::Vector<ComplexType>& E_vec,
                                       const U& Arr_data, 
                                       std::string filename, 
                                       std::string header)
{ 
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\n Root Writing " << filename << "\n";
        std::ofstream outfile;
        outfile.open(filename.c_str());

        auto const& Arr = Arr_data.const_table();
        auto thi = Arr_data.hi();

        outfile << header  << "\n";
        if(E_vec.size() == thi[0]) {   
            for (int e=0; e< thi[0]; ++e)
            {
                outfile << std::setw(10) << E_vec[e].real() 
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

