#ifndef SOURCE_TRAITS_H_
#define SOURCE_TRAITS_H_

#include "PointCharge.H"

#include "../../Code.H"
#include "../GeometryProperties/GeometryProperties.H"

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
#endif
