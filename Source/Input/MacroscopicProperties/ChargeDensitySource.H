#ifndef CHARGE_DENSITY_SOURCE_H_
#define CHARGE_DENSITY_SOURCE_H_

#include "../../Code.H"
#include "../GeometryProperties/GeometryProperties.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuQualifiers.H>

#include <algorithm>  // For std::copy
#include <type_traits>  // For std::is_trivially_copyable



class PointCharge {
public:
    amrex::Real pos[AMREX_SPACEDIM];  
    amrex::Real sigma;
    int charge_unit;

    PointCharge() : sigma(2e-10), charge_unit(1) {
        std::fill_n(pos, AMREX_SPACEDIM, 0.0);
    }
    
    PointCharge(const amrex::Real p[AMREX_SPACEDIM], 
               amrex::Real sigma_val = 2e-10, 
               int charge_unit_val = 1)
    : sigma(sigma_val), charge_unit(charge_unit_val) {
        std::copy(p, p + AMREX_SPACEDIM, pos);
    }    

    void Print() const {
        amrex::Print() << " Point Charge location: [";
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            amrex::Print() << std::scientific << std::setprecision(3) << pos[i];
            if (i < AMREX_SPACEDIM - 1) std::cout << ", ";
        }
        amrex::Print() << "], ";
        amrex::Print() << "  sigma: " << std::scientific << std::setprecision(3)  
                       << sigma 
                       << ", charge_unit: " << charge_unit << "\n";
    }
};


// Static assertion to ensure the class is trivially copyable
static_assert(std::is_trivially_copyable<PointCharge>::value, 
              "PointCharge must be trivially copyable");    

struct PointChargeDensityCalculator {
    AMREX_GPU_HOST_DEVICE
    double operator()(const PointCharge& pc, double x, double y, double z) const {

        static constexpr amrex::Real pi = static_cast<amrex::Real>(3.14159265358979323846);
        static constexpr amrex::Real q_e  = static_cast<amrex::Real>(1.602176634e-19);
        return q_e * pc.charge_unit *
                   ( 1./pow((sqrt(2*pi)*pc.sigma),3.)) *
                     exp(-1.*( (x-pc.pos[0])*(x-pc.pos[0]) +
                               (y-pc.pos[1])*(y-pc.pos[1])+
                               (z-pc.pos[2])*(z-pc.pos[2])
                            ) / (2 * pc.sigma*pc.sigma)
                        );
    }
};


//define a traits class
template<typename SourceType>
struct SourceTraits;

//specialization for PointCharge
template<>
struct SourceTraits<PointCharge>{
    using DensityCalculator = PointChargeDensityCalculator;

    static void CheckSource(const amrex::Vector<PointCharge>& vec_source) {
        auto& rCode = c_Code::GetInstance();
        auto& rGprop = rCode.get_GeometryProperties();

        std::string assert_msg = "Point charge is out of the physical domain! See pos: [";

        for(auto& source : vec_source) {

            assert_msg += std::to_string(source.pos[0]) 
                        + std::to_string(source.pos[1]) 
                        + std::to_string(source.pos[2]) + "]\n";

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    rGprop.Is_Point_Inside_Physical_Domain(source.pos), assert_msg);
        }
    }    
};


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
