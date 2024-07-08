#ifndef MACROSCOPIC_PROPERTIES_IMPL_H_
#define MACROSCOPIC_PROPERTIES_IMPL_H_

#include "MacroscopicProperties.H"
#include "Code.H"
#include "GeometryProperties.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
using namespace amrex;


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
