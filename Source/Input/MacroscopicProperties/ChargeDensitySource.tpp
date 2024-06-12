template<typename SourceType, typename ChargeDensityCalculator>
void 
ChargeDensitySource<SourceType, ChargeDensityCalculator>::
Define_PointCharge(size_t index, PointCharge&& pc);

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
