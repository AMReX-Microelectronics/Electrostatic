#ifndef MACROSCOPIC_PROPERTIES_IMPL_H_
#define MACROSCOPIC_PROPERTIES_IMPL_H_

#include "MacroscopicProperties.H"
#include "Code.H"
#include "GeometryProperties.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
using namespace amrex;

//Template<typename SourceType>
//Void
//C_MacroscopicProperties::Obtain_PotentialAtChargeDensitySources(
//        std::unique_ptr<ChargeDensitySource<SourceType>>& p_ChargeDensitySource)
//{
//    if (!p_ChargeDensitySource) {
//        throw std::invalid_argument("Null pointer provided for charge source");
//    }
//    amrex::MultiFab* const p_potential_mf = get_p_mf("phi");
//    
//    auto& rCode  = c_Code::GetInstance();
//    auto& rGprop = rCode.get_GeometryProperties();
//    const auto& plo = rGprop.geom.ProbLoArray();
//    auto iv = p_potential_mf->ixType().toIntVect();
//
//    int num_sources = p_ChargeDensitySource->get_num_sources();
//
//    amrex::Gpu::DeviceVector<amrex::Real> d_vec_V(num_sources);
//    std::fill(d_vec_V.begin(), d_vec_V.end(), 0);
//    amrex::Real* p_dV   = d_vec_V.dataPtr();
//
//    for ( amrex::MFIter mfi(*p_potential_mf, amrex::TilingIfNotGPU());
//                        mfi.isValid();
//                        ++mfi )
//    {
//        const auto& tb = mfi.tilebox(iv);
//        auto const& potential_array = p_potential_mf->array(mfi);
//
//        const auto* p_vec_source = p_ChargeDensitySource->get_p_sources();
//
//        const amrex::IntVect lo = tx.smallEnd();
//        const amrex::IntVect hi = tx.bigEnd();
//
//        amrex::ParallelFor (num_sources, [=] AMREX_GPU_DEVICE (int s) 
//        {
//            amrex::Real lx = (p_vec_sources[s].pos[0] - plo[0] - dx[0]*0.5)/dx[0];
//            amrex::Real ly = (p_vec_sources[s].pos[1] - plo[1] - dx[1]*0.5)/dx[1];
//            amrex::Real lz = (p_vec_sources[s].pos[2] - plo[2] - dx[2]*0.5)/dx[2];
//
//            int i_par = static_cast<int>(amrex::Math::floor(lx));
//            int j_par = static_cast<int>(amrex::Math::floor(ly));
//            int k_par = static_cast<int>(amrex::Math::floor(lz));
//
//            if(i_par >= lo[0] && i_par < hi[0] &&
//               j_par >= lo[1] && j_par < hi[1] &&
//               k_par >= lo[2] && k_par < hi[2]) 
//            {
//                p_dV[s] = potential_array(i_par,j_par,k_par);
//            }
//        });
//    }
//    amrex::Gpu::streamSynchronize();
//
//    MPI_Allreduce(MPI_IN_PLACE,
//                  p_dV.data(),
//                  num_sources,
//                  MPI_DOUBLE,
//                  MPI_SUM,
//                  ParallelDescriptor::Communicator());
//
//    amrex::ParallelFor (num_sources, [=] AMREX_GPU_DEVICE (int s) 
//    {
//        p_vec_source[s].current_potential = p_dV[s];
//    });
//    amrex::Gpu::streamSynchronize();
//}


template<typename SourceType>
void
c_MacroscopicProperties::Deposit_ExternalChargeDensitySources(
        std::unique_ptr<ChargeDensitySource<SourceType>>& p_ChargeDensitySource)
{
    if (!p_ChargeDensitySource) {
        throw std::invalid_argument("Null pointer provided for charge source");
    }

    auto& rCode  = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    //auto& rMprop = rCode.get_MacroscopicProperties();
    auto& geom   = rGprop.geom;
    auto dx      = geom.CellSizeArray();
    auto& real_box = geom.ProbDomain();

    amrex::MultiFab* const p_charge_density_mf = get_p_mf("charge_density");
    auto iv = p_charge_density_mf->ixType().toIntVect();

    for ( amrex::MFIter mfi(*p_charge_density_mf, amrex::TilingIfNotGPU());
                        mfi.isValid();
                        ++mfi )
    {
        //const auto& tb = mfi.tilebox( iv, p_charge_density_mf->nGrowVect() );
        const auto& tb = mfi.tilebox(iv);
        auto const& mf_array = p_charge_density_mf->array(mfi);

        int num_sources = p_ChargeDensitySource->get_num_sources();
        const auto* p_vec_source = p_ChargeDensitySource->get_p_sources();
        auto charge_density_calculator = p_ChargeDensitySource->get_charge_density_calculator();
        amrex::ParallelFor (tb,
            [=]
            AMREX_GPU_DEVICE (int i, int j, int k)
        {
                amrex::Real fac_x = (1._rt - iv[0]) * dx[0] * 0.5_rt;
                amrex::Real x = i * dx[0] + real_box.lo(0) + fac_x;

                amrex::Real fac_y = (1._rt - iv[1]) * dx[1] * 0.5_rt;
                amrex::Real y = j * dx[1] + real_box.lo(1) + fac_y;

                amrex::Real fac_z = (1._rt - iv[2]) * dx[2] * 0.5_rt;
                amrex::Real z = k * dx[2] + real_box.lo(2) + fac_z;

                for(int s=0; s< num_sources; ++s) {
                    mf_array(i,j,k) += charge_density_calculator(p_vec_source[s],x,y,z);
                }
        });
    }

    amrex::Gpu::streamSynchronize();
    p_charge_density_mf->FillBoundary(geom.periodicity());

}
#endif
