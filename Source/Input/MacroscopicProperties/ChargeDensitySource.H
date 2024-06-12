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

    // default
    PointCharge() : sigma(2e-10), charge_unit(1) {
        // Initialize pos array to zero or another suitable default value
        std::fill_n(pos, AMREX_SPACEDIM, 0.0);
    }
    
    PointCharge(const amrex::Real p[AMREX_SPACEDIM], 
               amrex::Real sigma_val = 2e-10, 
               int charge_unit_val = 1)
    : sigma(sigma_val), charge_unit(charge_unit_val) {
        std::copy(p, p + AMREX_SPACEDIM, pos);
    }    

    void Print() const {
        amrex::Print() << " location: [";
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

    AMREX_GPU_DEVICE
    double operator()(const PointCharge& pc, double x, double y, double z) const {

        //static constexpr amrex::Real pi = static_cast<amrex::Real>(3.14159265358979323846);
        //static constexpr amrex::Real q_e  = static_cast<amrex::Real>(1.602176634e-19);
        amrex::Real pi = static_cast<amrex::Real>(3.14159265358979323846);
        amrex::Real q_e  = static_cast<amrex::Real>(1.602176634e-19);
        return q_e * pc.charge_unit *
                   ( 1./pow((sqrt(2*pi)*pc.sigma),3.)) *
                     exp(-1.*( (x-pc.pos[0])*(x-pc.pos[0]) +
                               (y-pc.pos[1])*(y-pc.pos[1])+
                               (z-pc.pos[2])*(z-pc.pos[2])
                            ) / (2 * pc.sigma*pc.sigma)
                        );
    }
};

template<typename SourceType=PointCharge,
         typename ChargeDensityCalculator=PointChargeDensityCalculator>
class ChargeDensitySource {

    int num_sources;  

    SourceType* p_vec_source;

    amrex::Vector<SourceType> h_vec_source;

    #ifdef AMREX_USE_GPU
    int * d_num_sources;
    amrex::Gpu::DeviceVector<SourceType> d_vec_source;
    #endif

    ChargeDensityCalculator charge_density_calculator;

public:

    ChargeDensitySource(int _num_sources,
                        const ChargeDensityCalculator& calculator = ChargeDensityCalculator())
        : num_sources(_num_sources), 
          charge_density_calculator(calculator) 
    {
        amrex::Print() << " num_sources: " << num_sources << "\n";
        h_vec_source.resize(num_sources); 

        #ifdef AMREX_USE_GPU
        d_vec_source.resize(num_sources); 
        p_vec_source = d_vec_source.dataPtr();

        amrex::Gpu::htod_memcpy((void*)d_num_sources, &num_sources, sizeof(int));
        #else
        p_vec_source = h_vec_source.dataPtr();
        #endif
    }

    void Copy_HostVectorToGPU() 
    {
        if (!h_vec_source.empty()) {
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_vec_source.begin(),
                                                       h_vec_source.end(),
                                                       d_vec_source.begin() );
        }
    }

    AMREX_GPU_DEVICE
    double get_charge_sum(double x, double y, double z) const
    {
        double sum=0;
        amrex::Real pi = static_cast<amrex::Real>(3.14159265358979323846);
        amrex::Real q_e  = static_cast<amrex::Real>(1.602176634e-19);

        //for(int s=0; s< num_sources; ++s)
        //{
        //    //sum += charge_density_calculator(p_vec_source[s], x, y, z);

        //}

        if(p_vec_source != nullptr) {
            const PointCharge& pc = p_vec_source[0];
            sum += *d_num_sources;
        }
        //sum += q_e * pc.charge_unit *
        //       ( 1./pow((sqrt(2*pi)*pc.sigma),3.)) *
        //         exp(-1.*( (x-pc.pos[0])*(x-pc.pos[0]) +
        //                   (y-pc.pos[1])*(y-pc.pos[1])+
        //                   (z-pc.pos[2])*(z-pc.pos[2])
        //                ) / (2 * pc.sigma*pc.sigma)
        //            );
        return sum;
    }

    void Define_PointCharge(size_t index, PointCharge&& pc)
    {
        std::string assert_msg = "index: " + std::to_string(index)
                   + ", h_vec_source.size(): " + std::to_string(h_vec_source.size()) + "\n";
    
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(index < h_vec_source.size(), assert_msg);
    
        assert_msg = "point charge: " + std::to_string(index)
                   + " is out of the physical domain!\n";
    
        auto& rCode = c_Code::GetInstance();
        auto& rGprop = rCode.get_GeometryProperties();
    
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rGprop.Is_Point_Inside_Physical_Domain(pc.pos),
                                         assert_msg);
    
        h_vec_source[index] = std::move(pc);
    }

    void Print_PointCharge()
    {
        amrex::Print() << "Point Charges:\n";
        for(int s=0; s<h_vec_source.size(); ++s) {
            h_vec_source[s].Print();
        }
    }
};


using PointChargeSource = 
      ChargeDensitySource<PointCharge, PointChargeDensityCalculator>;

//#include "ChargeDensitySource.tpp"  // Include template definitions

#endif
