#include "ChargeDensitySource.H"

#include "../../Utils/CodeUtils/CloudInCell.H"
#include "../../Utils/SelectWarpXUtils/WarpXConst.H"
#include "Code.H"
#include "GeometryProperties.H"

//#include <AMReX_ParmParse.H>
//#include <AMReX_Parser.H>

using namespace amrex;

template class c_ChargeDensitySource<c_PointChargeContainer>;

template <typename ParticleContainerType>
void c_ChargeDensitySource<ParticleContainerType>::Deposit(
    amrex::MultiFab *const p_rho_mf)
{
    auto &rCode = c_Code::GetInstance();
    auto &rGprop = rCode.get_GeometryProperties();

    const auto &plo = rGprop.geom.ProbLoArray();
    const auto &dx = rGprop.geom.CellSizeArray();
    amrex::Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

    int lev = 0;
    for (auto pti = typename ParticleContainerType::ParIter(m_container, lev);
         pti.isValid(); ++pti)
    {
        auto const &rho = p_rho_mf->array(pti);

        auto np = pti.numParticles();
        const auto &particles = pti.GetArrayOfStructs();
        const auto *p_par = particles().data();

        const auto &charge_unit = pti.get_charge_unit();
        const auto *p_charge_unit = charge_unit.data();

        const auto &occupation = pti.get_occupation();
        const auto *p_occupation = occupation.data();

        amrex::Real unit_charge = PhysConst::q_e;

        amrex::ParallelFor(np,
                           [=] AMREX_GPU_DEVICE(int p)
                           {
                               amrex::Real qp = p_occupation[p] *
                                                p_charge_unit[p] * unit_charge /
                                                vol;

                               CloudInCell::Deposit_Trilinear(p_par[p].pos(),
                                                              plo, dx, rho, qp);
                           });
    }
}

template <typename ParticleContainerType>
void c_ChargeDensitySource<ParticleContainerType>::Gather(
    amrex::MultiFab *const p_potential_mf)
{
    auto &rCode = c_Code::GetInstance();
    auto &rGprop = rCode.get_GeometryProperties();

    const auto &plo = rGprop.geom.ProbLoArray();
    const auto &dx = rGprop.geom.CellSizeArray();

    int lev = 0;
    for (auto pti = typename ParticleContainerType::ParIter(m_container, lev);
         pti.isValid(); ++pti)
    {
        auto const &phi = p_potential_mf->array(pti);

        auto np = pti.numParticles();
        const auto &particles = pti.GetArrayOfStructs();
        const auto *p_par = particles().data();

        auto &par_potential = pti.get_potential();
        auto *p_par_potential = par_potential.data();

        amrex::ParallelFor(np,
                           [=] AMREX_GPU_DEVICE(int p)
                           {
                               p_par_potential[p] =
                                   CloudInCell::Gather_Trilinear(p_par[p].pos(),
                                                                 plo, dx, phi);
                           });
    }
}
