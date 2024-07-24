#ifndef POINT_CHARGE_H_
#define POINT_CHARGE_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Print.H>

#include <algorithm>  
#include <type_traits>  

class PointCharge {
public:
    amrex::Real pos[AMREX_SPACEDIM]; //m
    amrex::Real sigma; 
    amrex::Real current_potential; //eV
    amrex::Real threshold_potential; //eV
    int charge_unit;
    int boxID;

    PointCharge() : boxID(-1), sigma(2e-10), charge_unit(1), 
                    threshold_potential(0.), current_potential(0.) {
        std::fill_n(pos, AMREX_SPACEDIM, 0.0);
    }

    PointCharge(const amrex::Real p[AMREX_SPACEDIM],
                int boxID_val,
               amrex::Real sigma_val = 2e-10,
               int charge_unit_val = 1,
               amrex::Real threshold_val = 0.)
    : boxID(boxID_val), 
      sigma(sigma_val), charge_unit(charge_unit_val), threshold_potential(threshold_val)
    {
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
                       << ", charge_unit: " 
                       << charge_unit 
                       << ", threshold_potential: "
                       << threshold_potential << "\n";
    }

    RealVect pos () const & {
        return RealVect(AMREX_D_DECL(this->pos[0], this->pos[1], this->pos[2]));
    }
};
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
#endif
