#ifndef MATH_FUNCTIONS_H_
#define MATH_FUNCTIONS_H_

#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>

namespace MathFunctions
{
AMREX_GPU_HOST_DEVICE
amrex::Real Sigmoid(amrex::Real x) { return 1.0 / (1.0 + exp(-x)); }
}  // namespace MathFunctions

namespace MathFunctionsHost
{
void GenerateRandomParticles(
    amrex::Vector<amrex::Vector<amrex::Real>> &particles, int num_particles,
    unsigned int seed,
    amrex::Vector<amrex::Real> vec_scaling =
        amrex::Vector<amrex::Real>({1.0, 1.0, 1.0}),
    amrex::Vector<amrex::Real> vec_offset =
        amrex::Vector<amrex::Real>({0.0, 0.0, 0.0}))
{
    amrex::Real epsilon = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Real adjusted_min = 0. + epsilon;
    amrex::Real adjusted_max = 1. - epsilon;

    // Create a random engine and set the seed
    std::mt19937 random_engine(seed);
    std::uniform_real_distribution<amrex::Real> distribution(adjusted_min,
                                                             adjusted_max);

    for (int i = 0; i < num_particles; ++i)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            particles[i][d] =
                distribution(random_engine) * vec_scaling[d] + vec_offset[d];
        }
    }
}

}  // namespace MathFunctionsHost
#endif
