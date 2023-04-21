#include "Common_Properties.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"


/*Explicit specializations*/
template class c_Common_Properties<ComplexType[NUM_MODES]>;            //of c_CNT
template class c_Common_Properties<ComplexType[NUM_MODES][NUM_MODES]>; //c_Graphene


template<typename T>
AMREX_GPU_HOST_DEVICE
ComplexType 
c_Common_Properties<T>::conjugate(ComplexType a)
{
   ComplexType a_conj(a.real(), -1.*a.imag());
   return a_conj;
}


template<typename T>
void
c_Common_Properties<T>:: ReadNanostructureProperties ()
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
    
}


template<typename T>
void 
c_Common_Properties<T>::DefineMatrixPartition() 
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

}


template<typename T>
void 
c_Common_Properties<T>::AllocateArrays () 
{
    h_Ha_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());
    SetZero_BlkTable1D(h_Ha_loc_data);

    offDiag_repeatBlkSize = get_offDiag_repeatBlkSize();
    h_Hb_loc_data.resize({0},{offDiag_repeatBlkSize}, The_Pinned_Arena());
    SetZero_BlkTable1D(h_Hb_loc_data);

    h_Hc_loc_data.resize({0},{offDiag_repeatBlkSize}, The_Pinned_Arena());
    SetZero_BlkTable1D(h_Hc_loc_data);

    h_tau_glo_data.resize({0},{NUM_CONTACTS},The_Pinned_Arena());
    SetZero_BlkTable1D(h_tau_glo_data);

    h_Sigma_glo_data.resize({0,0},{NUM_CONTACTS,NUM_ENERGY_PTS_REAL},The_Pinned_Arena());
    SetZero_BlkTable2D(h_Sigma_glo_data);

    #if AMREX_USE_GPU
    d_GR_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    d_A_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    #else
    h_GR_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    SetZero_BlkTable2D(h_GR_loc_data);

    h_A_loc_data.resize({0,0}, {Hsize_glo, blkCol_size_loc}, The_Arena());
    SetZero_BlkTable2D(h_A_loc_data);
    #endif

    h_E_RealPath_data.resize({0},{NUM_ENERGY_PTS_REAL},The_Pinned_Arena());

    Print_BlkTable2D_loc(h_A_loc_data);
}


template<typename T>
void 
c_Common_Properties<T>:: ConstructHamiltonian () {}


template<typename T>
void 
c_Common_Properties<T>:: AddPotentialToHamiltonian () 
{
    auto const& h_Ha = h_Ha_loc_data.table();
    int c=0;
    for(auto& col_gid: vec_blkCol_gids)
    {
        h_Ha(c) = h_Ha(c) + avg_gatherField[col_gid];
        ++c;
    }
}


template<typename T>
void 
c_Common_Properties<T>:: Update_ContactPotential () 
{
    for(int c=0; c<NUM_CONTACTS; ++c)
    {
        U_contact[c] = avg_gatherField[global_contact_index[c]];
        ++c;
    }
}

//template<typename T>
//void 
//c_Common_Properties<T>:: DefineIntegrationPaths ()
//{
//
//
//}
//
//

template<typename T>
void 
c_Common_Properties<T>:: Compute_SurfaceGreensFunction (MatrixBlock<T> gr, const ComplexType EmU) {}


template<typename T>
void 
c_Common_Properties<T>:: Define_tau () {}


template<typename T>
void 
c_Common_Properties<T>:: Define_SelfEnergy ()
{
    Define_tau();
    auto const& h_tau   = h_tau_glo_data.table();
    auto const& h_Sigma = h_Sigma_glo_data.table();
    auto const& h_E     = h_E_RealPath_data.table();
    
    amrex::Print() << "Printing tau:\n";
    Print_BlkTable1D_loc(h_tau_glo_data);
    for (std::size_t c = 0; c < NUM_CONTACTS; ++c)
    {
        for (std::size_t e = 0; e < NUM_ENERGY_PTS_REAL; ++e)
        {
            MatrixBlock<T> gr;
            Compute_SurfaceGreensFunction(gr, h_E(e)-U_contact[c]);
            h_Sigma(c,e) = h_tau(c)*gr*h_tau(c);
        }
    }
    Write_BlkTable2D_asaf_E(h_Sigma_glo_data, "Sigma", "Er Ei Sigma_s_r Sigma_s_i Sigma_d_r Sigma_d_i");
    //amrex::Print() << "Printing Sigma in common: \n";
    //amrex::Print() << h_Sigma(0,0).block[0] << "\n";

}


template<typename T>
AMREX_GPU_HOST_DEVICE
ComplexType 
c_Common_Properties<T>::get_Gamma(ComplexType Sigma)
{
    ComplexType val(-2.*Sigma.imag(), 0.);
    return val;
}


template<typename T>
void
c_Common_Properties<T>::Write_BlkTable1D_asaf_E(const BlkTable1D& Arr_data, 
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

        outfile << header  << "\n";
        for (int e=0; e< thi[0]; ++e)
        {
            outfile << std::setw(10) << E(e)
                    << std::setw(15) << Arr(e) << "\n";
        }

        outfile.close();
    }
}


template<typename T>
void
c_Common_Properties<T>::Write_BlkTable2D_asaf_E(const BlkTable2D& Arr_data, 
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
void
c_Common_Properties<T>::Print_BlkTable1D_loc(const BlkTable1D& Tab1D_data)
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
void
c_Common_Properties<T>::Print_BlkTable2D_loc(const BlkTable2D& Tab2D_data)
{
    auto tlo = Tab2D_data.lo();
    auto thi = Tab2D_data.hi();

    auto const& Tab2D = Tab2D_data.table();

    for (int i = tlo[0]; i < thi[0]; ++i)
    {   
        for (int j = tlo[1]; j < thi[1]; ++j) //slow moving index. printing slow
        {       
            std::cout << std::setw(15) << std::setprecision(2) << Tab2D(i,j);
        }
        std::cout << "\n";
    }
}


template<typename T>
void 
c_Common_Properties<T>:: SetZero_BlkTable1D (BlkTable1D& Tab1D_data) 
{
    auto tlo = Tab1D_data.lo();
    auto thi = Tab1D_data.hi();

    auto const& Tab1D = Tab1D_data.table();

    for (int i = tlo[0]; i < thi[0]; ++i)
    {   
        ComplexType val(0.,0.);   
        Tab1D(i) = val; 
    }
}

template<typename T>
void 
c_Common_Properties<T>:: SetZero_BlkTable2D (BlkTable2D& Tab2D_data) 
{
    auto tlo = Tab2D_data.lo();
    auto thi = Tab2D_data.hi();

    auto const& Tab2D = Tab2D_data.table();

    for (int i = tlo[0]; i < thi[0]; ++i)
    {   
        for (int j = tlo[1]; j < thi[1]; ++j) //slow moving index. printing slow
        {   
            ComplexType val(0.,0.);   
            Tab2D(i,j) = val; 
        }
    }
}
//template<typename T>
//void 
//c_Common_Properties<T>::DeallocateArrays () 
//{
//
//}
//
//
//template<typename T>
//void 
//c_Common_Properties<T>::ComputeChargeDensity () 
//{
//
//}


