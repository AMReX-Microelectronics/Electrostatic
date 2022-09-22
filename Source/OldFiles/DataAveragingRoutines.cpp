#include "DataAveragingRoutines.H"

void AverageCellCenteredMultiFabToCellFaces(const MultiFab& cc_arr,
                   std::array< MultiFab, AMREX_SPACEDIM >& face_arr)
{
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
}


void ComputeSumOfFluxes(amrex::Vector<Array< MultiFab, AMREX_SPACEDIM >>& minus_epsGradPhi)
{
//     auto reduction_function_parser = m_parser->compile<m_nvars>();
//     amrex::ReduceOps<ReduceOp> reduceFluxX_op;
//     amrex::ReduceData<amrex::Real> reduceFluxX_data(reduceFluxX_op);
//     using ReduceTuple = typename decltype(reduceFluxX_data)::Type;
//
//     for ( amrex::MFIter mfi(Ex, false); mfi.isValid(); ++mfi)
//     {
//         const amrex::Box& tx = mfi.tilebox(Ex_nodalType);
//         const amrex::Box& ty = mfi.tilebox(Ey_nodalType);
//         const amrex::Box& tz = mfi.tilebox(Ez_nodalType);
//         reduceEx_op.eval(tx, reduceEx_data,
//         [=] AMREX_GPU_DEVICE (int i, int j, int k) ->ReduceTuple
//         {
//           
//           return flux_arr(i,j,k);
//         });
//     }
//     amrex::Real reducedFluxX_value = amrex::get<0>(reduceFluxX_data.value());
//     //MPI Reduce
//     if (std::is_same<ReduceOp, amrex::ReduceOpSum>::value)
//     {
//          amrex::ParallelDescriptor::ReduceRealSum(reducedFluxX_value);
//
//          if (integral_type == 0) {
//              amrex::Real area = dx[1]*dx[2];
//              reducedFluxX_value *= area;
//          }
//     }
}
