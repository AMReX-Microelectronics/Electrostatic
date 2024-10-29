#ifndef NANOSTRUCTURE_H_
#define NANOSTRUCTURE_H_

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>

#include <string>

#include "../../Utils/CodeUtils/ParticleStructure.H"
#include "NEGF/CNT.H"
#include "NEGF/Graphene.H"
#include "Nanostructure_fwd.H"

template <typename NSType>
class c_Nanostructure
    : private amrex::ParticleContainer<realPD::NUM, intPD::NUM, realPA::NUM,
                                       intPA::NUM>,
      public NSType
{
    using RealTable1D = TableData<amrex::Real, 1>;

    int _use_electrostatic = 0;
    int _use_negf = 0;
    const amrex::Geometry *const _geom;
    amrex::MultiFab *p_mf_gather = nullptr;
    amrex::MultiFab *p_mf_deposit = nullptr;
    const amrex::GpuArray<int, AMREX_SPACEDIM> *_n_cell;
    amrex::Vector<s_Position3D> pos_vec;
    int particle_id_offset = 0;

    void Fill_AtomLocations();
    void Read_AtomLocations();

    void Set_GatherAndDepositMultiFabs();
    amrex::Real Compute_CellVolume();

   public:
    c_Nanostructure(const amrex::Geometry &geom,
                    const amrex::DistributionMapping &dm,
                    const amrex::BoxArray &ba, const std::string NS_name,
                    const int NS_id, const std::string NS_gather_str,
                    const std::string NS_deposit_str,
                    const amrex::Real NS_initial_deposit_value,
                    const int use_negf, const std::string negf_foldername_str);

    void Evaluate_LocalFieldSites();
    void Gather_MeshAttributeAtAtoms();
    void Deposit_AtomAttributeToMesh();
    void Deposit_ZeroToMesh();
    void Obtain_PotentialAtSites();
    void Mark_CellsWithAtoms();
};
#endif
