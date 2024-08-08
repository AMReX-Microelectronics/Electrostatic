#ifndef POINT_CHARGE_CONTAINER_H_
#define POINT_CHARGE_CONTAINER_H_

#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <mutex>    

struct PointChargeAttributes {
    struct realPA { enum { charge_unit, occupation, potential, NUM }; };
    struct intPA { enum { id, NUM }; };
    struct realPD { enum {  NUM }; };
    struct intPD { enum { NUM }; };
};
using PCA = PointChargeAttributes;

class c_PointChargeContainer
    : public amrex::ParticleContainer<PCA::realPD::NUM,
                                      PCA::intPD::NUM,
                                      PCA::realPA::NUM,
                                      PCA::intPA::NUM>
{

public:
    using ParticleContainerType = 
        amrex::ParticleContainer<PCA::realPD::NUM, PCA::intPD::NUM, 
                                 PCA::realPA::NUM, PCA::intPA::NUM>;
    
    //using ParticleType = typename ParticleContainerType::ParticleType;

    c_PointChargeContainer(const amrex::Geometry& geom,
                           const amrex::DistributionMapping& dm,
                           const amrex::BoxArray& ba);

    void Define(const amrex::Vector<amrex::Vector<amrex::Real>>& particles,
                const amrex::Vector<amrex::Real>& charge_units = {},
                const amrex::Vector<amrex::Real>& occupations = {});

    static void Set_V0(amrex::Real value);

    static void Set_Et(amrex::Real value);

    struct GetV0 {
       AMREX_GPU_HOST_DEVICE int operator() () const { return V0; }
       amrex::Real V0;
    };
    static GetV0 Get_V0() { return GetV0{V0}; }    

    struct GetEt {
       AMREX_GPU_HOST_DEVICE int operator() () const { return Et; }
       amrex::Real Et;
    };
    static GetEt Get_Et() { return GetEt{Et}; }    

    void Print_Container(bool print_positions = false);

    void Compute_Occupation();

    bool Is_Occupation_Variable() {
        return flag_vary_occupation;
    }

    // Define the iterator class to iterate over the particles
    class ParIter
        : public amrex::ParIter<PCA::realPD::NUM, PCA::intPD::NUM, 
                                PCA::realPA::NUM, PCA::intPA::NUM>
    {
    public:
        using amrex::ParIter<PCA::realPD::NUM, PCA::intPD::NUM, 
                             PCA::realPA::NUM, PCA::intPA::NUM>::ParIter;

        using RealVector = amrex::ParIter<PCA::realPD::NUM,
                                          PCA::intPD::NUM,
                                          PCA::realPA::NUM,
                                          PCA::intPA::NUM>::ContainerType::RealVector;

        // Get realPA
        const std::array<RealVector, PCA::realPA::NUM>& get_realPA() const {
            return this->GetStructOfArrays().GetRealData();
        }

        std::array<RealVector, PCA::realPA::NUM>& get_realPA() {
            return this->GetStructOfArrays().GetRealData();
        }

        const RealVector& get_charge_unit() const {
            return this->GetStructOfArrays().GetRealData(PCA::realPA::charge_unit);
        }

        const RealVector& get_potential() const {
            return this->GetStructOfArrays().GetRealData(PCA::realPA::potential);
        }

        const RealVector& get_occupation() const {
            return this->GetStructOfArrays().GetRealData(PCA::realPA::occupation);
        }

        RealVector& get_charge_unit() {
            return this->GetStructOfArrays().GetRealData(PCA::realPA::charge_unit);
        }

        RealVector& get_potential() {
            return this->GetStructOfArrays().GetRealData(PCA::realPA::potential);
        }

        RealVector& get_occupation() {
            return this->GetStructOfArrays().GetRealData(PCA::realPA::occupation);
        }
    };

private:
   static amrex::Real V0;
   static amrex::Real Et;
   static std::mutex mtx;  

   bool flag_vary_occupation = false;

   void Read_PointCharges();
};

#endif
