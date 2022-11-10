#include <CodeUtil.H>

using namespace amrex;

void 
Multifab_Manipulation::InitializeMacroMultiFabUsingParser_3vars (amrex::MultiFab *macro_mf,
                                                                 amrex::ParserExecutor<3> const& macro_parser,
                                                                 amrex::Geometry& geom)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************Multifab_Manipulation::InitializeMacroMultiFabUsingParser************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    auto dx = geom.CellSizeArray();
    //amrex::Print() << "dx: " << dx[0] << " " << dx[1] << " " << dx[2] << "\n";

    auto& real_box = geom.ProbDomain();
    //amrex::Print() << "real_box_lo: " << real_box.lo(0) << " " << real_box.lo(1) << " " << real_box.lo(2) << "\n";

    auto iv = macro_mf->ixType().toIntVect();
    //amrex::Print() << "iv: " << iv[0] << " " << iv[1] << " " << iv[2] << "\n";

    for ( amrex::MFIter mfi(*macro_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const auto& tb = mfi.tilebox( iv, macro_mf->nGrowVect() ); /** initialize ghost cells in addition to valid cells.
                                                                       auto = amrex::Box
                                                                    */
        auto const& mf_array =  macro_mf->array(mfi); //auto = amrex::Array4<amrex::Real>

        amrex::ParallelFor (tb,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                amrex::Real fac_x = (1._rt - iv[0]) * dx[0] * 0.5_rt;
                amrex::Real x = i * dx[0] + real_box.lo(0) + fac_x;

                amrex::Real fac_y = (1._rt - iv[1]) * dx[1] * 0.5_rt;
                amrex::Real y = j * dx[1] + real_box.lo(1) + fac_y;

                amrex::Real fac_z = (1._rt - iv[2]) * dx[2] * 0.5_rt;
                amrex::Real z = k * dx[2] + real_box.lo(2) + fac_z;

                mf_array(i,j,k) = macro_parser(x,y,z);
        });
    }
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************Multifab_Manipulation::InitializeMacroMultiFabUsingParser************************\n";
#endif
}


void 
Multifab_Manipulation::AverageCellCenteredMultiFabToCellFaces(const amrex::MultiFab& cc_arr,
                                       std::array< amrex::MultiFab, 
                                       AMREX_SPACEDIM >& face_arr)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************Multifab_Manipulation::AverageCellCenteredMultiFabToCellFaces()************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif
    for (MFIter mfi(cc_arr, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Array4<Real const> & cc = cc_arr.array(mfi);
        AMREX_D_TERM(const Array4<Real> & facex = face_arr[0].array(mfi);,
                     const Array4<Real> & facey = face_arr[1].array(mfi);,
                     const Array4<Real> & facez = face_arr[2].array(mfi););

        AMREX_D_TERM(const Box & nodal_x = mfi.nodaltilebox(0);,
                     const Box & nodal_y = mfi.nodaltilebox(1);,
                     const Box & nodal_z = mfi.nodaltilebox(2););

        amrex::ParallelFor(nodal_x, nodal_y, nodal_z,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            facex(i,j,k) = 0.5*(cc(i,j,k)+cc(i-1,j,k));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            facey(i,j,k) = 0.5*(cc(i,j,k)+cc(i,j-1,k));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            facez(i,j,k) = 0.5*(cc(i,j,k)+cc(i,j,k-1));
        });
    }
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************Multifab_Manipulation::AverageCellCenteredMultiFabToCellFaces()************************\n";
#endif

}

void
Multifab_Manipulation::SpecifyValueOnlyOnCutcells(amrex::MultiFab* mf, amrex::Real value) 
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************Multifab_Manipulation::SpecifyValueOnlyOnCutcells************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    auto factory  = dynamic_cast<amrex::EBFArrayBoxFactory const*>(&(mf->Factory()));

    auto const &flags = factory->getMultiEBCellFlagFab();
    auto const &vfrac = factory->getVolFrac();

    auto iv = mf->ixType().toIntVect();

    for ( amrex::MFIter mfi(flags, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi ) 
    {
        const auto& box = mfi.tilebox( iv, mf->nGrowVect() ); 

        auto const& mf_array =  mf->array(mfi); //auto = amrex::Array4<amrex::Real>

        amrex::FabType fab_type = flags[mfi].getType(box);

        if(fab_type == amrex::FabType::regular) 
        {
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
               mf_array(i, j, k) = amrex::Real(0.);
            });
        }
        else if (fab_type == amrex::FabType::covered) 
        {
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
               mf_array(i, j, k) = amrex::Real(0.);
            });
        }
        else //box contains some cutcells
        {
            auto const &vfrac_array = vfrac.const_array(mfi);

            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
               if(vfrac_array(i,j,k) > 0 and vfrac_array(i,j,k) < 1) 
               {
                   mf_array(i, j, k) = value;
               } 
            });
        }
    }
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************Multifab_Manipulation::SpecifyValueOnlyOnCutcells************************\n";
#endif
}

//AMREX_GPU_DEVICE AMREX_FORCE_INLINE
//void
////Multifab_Manipulation::GetXYZ(const int i, const int j, const int k, 
////		              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx, 
////		              const amrex::RealBox& real_box, const amrex::IntVect& iv, 
////		              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& coord_vec)
//Multifab_Manipulation::GetXYZ(const int i, const int j, const int k, 
//		              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& coord_vec)
//{
//#ifdef PRINT_NAME
//    amrex::Print() << "\n\n\t\t\t\t\t{************************Multifab_Manipulation::GetXYZ()************************\n";
//    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
//#endif
//      coord_vec[0] = i * dx[0];
//      coord_vec[1] = j * dx[1];
//      coord_vec[2] = k * dx[2];
////    amrex::Real fac_x = (1._rt - iv[0]) * dx[0] * 0.5_rt;
////    coord_vec[0] = i * dx[0] + real_box.lo(0) + fac_x;
////
////    amrex::Real fac_y = (1._rt - iv[1]) * dx[1] * 0.5_rt;
////    coord_vec[1] = j * dx[1] + real_box.lo(1) + fac_y;
////
////    amrex::Real fac_z = (1._rt - iv[2]) * dx[2] * 0.5_rt;
////    coord_vec[2] = k * dx[2] + real_box.lo(2) + fac_z;
//
//#ifdef PRINT_NAME
//    amrex::Print() << "\t\t\t\t\t}************************Multifab_Manipulation::GetXYZ()************************\n";
//#endif
//}

void PrintRunDiagnostics(amrex::Real initial_time) 
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t{************************PrintRunDiagnostics()************************\n";
    amrex::Print() << "\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    // MultiFab memory usage
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    amrex::Long min_fab_megabytes  = amrex::TotalBytesAllocatedInFabsHWM()/1048576;
    amrex::Long max_fab_megabytes  = min_fab_megabytes;

    ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
    ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

    amrex::Print() << "\n\tHigh-water FAB megabyte spread across MPI nodes: ["
                   << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";

    min_fab_megabytes  = amrex::TotalBytesAllocatedInFabs()/1048576;
    max_fab_megabytes  = min_fab_megabytes;

    ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
    ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

    amrex::Print() << "\tCurent     FAB megabyte spread across MPI nodes: ["
                   << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";


    amrex::Real total_run_time = ParallelDescriptor::second() - initial_time;
    ParallelDescriptor::ReduceRealMax(total_run_time);

    amrex::Print() << "\tTotal run time " << total_run_time << " seconds\n\n";

#ifdef PRINT_NAME
    amrex::Print() << "\t}************************PrintRunDiagnostics()************************\n";
#endif
}
