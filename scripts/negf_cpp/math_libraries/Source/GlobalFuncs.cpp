
#include "GlobalFuncs.H"

#include <random>

ComplexType randomComplex(amrex::Real scale)
{
    static std::mt19937 rng(std::random_device{}());
    static std::uniform_real_distribution<double> dist(-1.0, 1.0);

    return ComplexType(scale * dist(rng), scale * dist(rng));
}
