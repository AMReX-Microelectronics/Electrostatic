#ifndef PARTICLE_STRUCTURE_H_
#define PARTICLE_STRUCTURE_H_

#include <AMReX_Particles.H>
#include <AMReX_REAL.H>

// Internal variables: x,y,z, id, cpu. I think id is the local id (not sure).

// Extra Particle data in Array-of-Structs form
// real: 0
// int: 0

// Extra Particle attributes in Struct-of-Arrays form
// real: phi, charge
// int: 2 //cell_id, atom_id (from atom_id we can obtain (ring_number,
// azimuthal_number)

struct realPA  // real Particle Attribute
{
    enum
    {
        gather,
        deposit,
        NUM
    };
};

struct intPA  // integer Particle Attribute
{
    enum
    {
        cid,  // cell id
        NUM
    };
};

struct realPD  // real Particle Data
{
    enum
    {
        NUM
    };
};

struct intPD  // integer Particle Data
{
    enum
    {
        NUM
    };
};

class MyParIter
    : public amrex::ParIter<realPD::NUM, intPD::NUM, realPA::NUM, intPA::NUM>
{
   public:
    using amrex::ParIter<realPD::NUM, intPD::NUM, realPA::NUM,
                         intPA::NUM>::ParIter;
    using RealVector = amrex::ParIter<realPD::NUM, intPD::NUM, realPA::NUM,
                                      intPA::NUM>::ContainerType::RealVector;
    using IntVector = amrex::ParIter<realPD::NUM, intPD::NUM, realPA::NUM,
                                     intPA::NUM>::ContainerType::IntVector;

    // get realPA
    const std::array<RealVector, realPA::NUM> &get_realPA() const
    {
        return GetStructOfArrays().GetRealData();
    }

    std::array<RealVector, realPA::NUM> &get_realPA()
    {
        return GetStructOfArrays().GetRealData();
    }
    // get realPA component
    const RealVector &get_realPA_comp(int comp) const
    {
        return GetStructOfArrays().GetRealData(comp);
    }

    RealVector &get_realPA_comp(int comp)
    {
        return GetStructOfArrays().GetRealData(comp);
    }

    // get intPA
    const std::array<IntVector, intPA::NUM> &get_intPA() const
    {
        return GetStructOfArrays().GetIntData();
    }

    std::array<IntVector, intPA::NUM> &get_intPA()
    {
        return GetStructOfArrays().GetIntData();
    }
    // get intPA component
    const IntVector &get_intPA_comp(int comp) const
    {
        return GetStructOfArrays().GetIntData(comp);
    }

    IntVector &get_intPA_comp(int comp)
    {
        return GetStructOfArrays().GetIntData(comp);
    }

    // get realPA
    // std::array<> (int comp) {
    //     return GetArrayOfStructs();
    // }
};
#endif
