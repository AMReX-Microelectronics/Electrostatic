#ifndef NANOSTRUCTURE_H_
#define NANOSTRUCTURE_H_

#include "Nanostructure_fwd.H"
#include "../../Code_definitions.H"

#include <AMReX_REAL.H>
//#include <AMReX_GpuContainers.H>
//#include <AMReX_Vector.H>
//#include <AMReX_IntVect.H>
//#include <AMReX_BoxArray.H>
//#include <AMReX_Geometry.H>
//#include <AMReX_MLLinOp.H>
//#include <AMReX_MLABecLaplacian.H>
//#include <AMReX_MLEBABecLap.H>
//#include <AMReX_MLMG.H>
//#include <AMReX_BaseFwd.H>

#include <string>

//Internal variables: x,y,z, id, cpu. I think id is the local id (not sure).

//Extra Particle data in Array-of-Structs form 
//real: 0
//int: 0

//Extra Particle attributes in Struct-of-Arrays form
//real: phi, charge
//int: 2 //cell_id, atom_id (from atom_id we can obtain (ring_number, azimuthal_number)

struct realPA //real Particle Attribute
{
    enum { 
          phi,
          charge,
          NUM 
    };
};
struct intPA //integer Particle Attribute
{
    enum { 
          cid, //cell id
          aid, //atom id
          NUM 
    };
};

struct realPD //real Particle Data
{
    enum {
          NUM 
    };
};
struct intPD //integer Particle Data
{
    enum {
          NUM 
    };
};


class MyParIter
    : public amrex::ParIter<realPD::NUM, intPD::NUM, realPA::NUM, intPA::NUM>
{
public:
    using amrex::ParIter<realPD::NUM, intPD::NUM, realPA::NUM, intPA::NUM>::ParIter;
    using RealVector = ParIter::ContainerType::RealVector;
    using IntVector = ParIter::ContainerType::IntVector;

    //get realPA
    const std::array<RealVector, realPA::NUM>& get_realPA () const {
        return GetStructOfArrays().GetRealData();
    }
    
    std::array<RealVector, realPA::NUM>& get_realPA () {
        return GetStructOfArrays().GetRealData();
    }

    //get intPA
    const std::array<IntVector, intPA::NUM>& get_intPA () const {
        return GetStructOfArrays().GetIntData();
    }
    
    std::array<IntVector, intPA::NUM>& get_intPA () {
        return GetStructOfArrays().GetIntData();
    }
    //get realPA component
    const RealVector& get_realPA_comp (int comp) const {
        return GetStructOfArrays().GetRealData(comp);
    }
    
    RealVector& get_realPA_comp (int comp) {
        return GetStructOfArrays().GetRealData(comp);
    }
    //get intPA component
    const IntVector& get_intPA_comp (int comp) const {
        return GetStructOfArrays().GetRealData(comp);
    }
    
    IntVector& get_intPA_comp (int comp) {
        return GetStructOfArrays().GetRealData(comp);
    }

};



class 
c_Nanostructure
               : public amrex::ParticleContainer<realPD::NUM, 
                                                 intPD::NUM, 
                                                 realPA::NUM, 
                                                 intPA::NUM>
{
public:

private:
    
    c_Nanostructure (const amrex::Geometry            & geom,
                     const amrex::DistributionMapping & dmap,
                     const amrex::BoxArray            & ba)
        : amrex::ParticleContainer<realPD::NUM, intPD::NUM, 
                                   realPA::NUM, intPA::NUM> (geom, dmap, ba);
    //~c_Nanostructure_Atom_Container();

    void ReadAtomLocations ();
    void ReadNanostructureProperties ();

    //void pushX (const amrex::Real& dt);

    //void DepositCharge (const amrex::Vector<amrex::MultiFab*>& rho);

    //void FieldGather (const amrex::Vector<std::array<const amrex::MultiFab*, AMREX_SPACEDIM> >& E,
    //                  const amrex::Vector<const amrex::iMultiFab*>& masks);

    //void Evolve (const amrex::Vector<std::array<const amrex::MultiFab*, AMREX_SPACEDIM> >& E,
    //             const amrex::Vector<const amrex::MultiFab*>& rho, const amrex::Real& dt);

    //void writeParticles (int n);

    //amrex::Real charge, mass;

    amrex::Real acc; //bond length

};

#endif
