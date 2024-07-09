#ifndef CHARGE_DENSITY_SOURCE_H_
#define CHARGE_DENSITY_SOURCE_H_

#include "SourceTraits.H"
#include <AMReX_GpuContainers.H>


template<typename SourceType=PointCharge>
class ChargeDensitySource {

    amrex::Vector<SourceType> h_vec_source;
    amrex::Gpu::DeviceVector<SourceType> d_vec_source;
    
public:
    using ChargeDensityCalculator = typename SourceTraits<SourceType>::DensityCalculator;

    ChargeDensitySource(amrex::Vector<SourceType>&& vec_source, bool PrintFlag=false) 
        : h_vec_source(std::move(vec_source))
    {
        SourceTraits<SourceType>::CheckSource(h_vec_source);

        #ifdef AMREX_USE_GPU
        if (!h_vec_source.empty()) {
            d_vec_source.resize(h_vec_source.size());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_vec_source.begin(),
                                                       h_vec_source.end(),
                                                       d_vec_source.begin() );
        }
        #endif

        if(PrintFlag) Print_ChargeDensity();
    }

    int get_num_sources() const {
        return h_vec_source.size();
    }

    const SourceType* get_p_sources() const {
        #ifdef AMREX_USE_GPU
        return d_vec_source.dataPtr();
        #else
        return h_vec_source.dataPtr();
        #endif
    }

    ChargeDensityCalculator get_charge_density_calculator() {
       return ChargeDensityCalculator();
    }

    void Print_ChargeDensity() {
        for(int s=0; s<h_vec_source.size(); ++s) {
            h_vec_source[s].Print();
        }
    }
};

using PointChargeSource = ChargeDensitySource<PointCharge>;
#endif
