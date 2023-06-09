#ifndef CNT_H_
#define CNT_H_

#include "NEGF_Common.H"
#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Array.H>

#include <string>
#include <array>

class c_CNT
           :public c_NEGF_Common<ComplexType[NUM_MODES]>
{

public:
    using BlkType = ComplexType[NUM_MODES];

    struct Get1DSiteID 
    {
       AMREX_GPU_HOST_DEVICE int operator() (int global_par_id) const 
       { 
	       return static_cast<int>(amrex::Math::floor((global_par_id-1)/type_id[0])); 
       }
       amrex::Array<int,2> type_id;
    };

    static Get1DSiteID get_1D_site_id() 
    {
       return Get1DSiteID{type_id}; 
    }

    struct GetAtomIDAtSite 
    {
       AMREX_GPU_HOST_DEVICE int operator() (int global_par_id) const 
       { 
	       return static_cast<int>((global_par_id-1) % type_id[0]); 
       }
       amrex::Array<int,2> type_id;
    };

    static GetAtomIDAtSite get_atom_id_at_site() 
    {
        return GetAtomIDAtSite{type_id}; 
    }

    static amrex::Array<int,2> type_id; 
    //type_id[0]: m, type_id[1]: n for (m,n) tube

    amrex::Real acc;
    amrex::Real R_cnt;
    int rings_per_unitcell;
    int atoms_per_ring;
    amrex::Real gamma;

    amrex::Vector<int> mode_vec;
    amrex::Vector<int> mode_degen_vec;
    amrex::Vector<amrex::Real> bandgap_vec;
    MatrixBlock<BlkType> beta;

    void Define_SortedModeVector();
    virtual int get_Hsize() override { return rings_per_unitcell*num_unitcells; }
    virtual int get_offDiag_repeatBlkSize() override { return 2; }
    virtual void set_material_specific_parameters() override;
    virtual void Define_PotentialProfile() override;
    virtual void Generate_AtomLocations(amrex::Vector<s_Position3D>& pos) override;
    s_Position3D get_AtomPosition_ZigZag_CNT(int ring_id, int atom_id);

    virtual void ReadNanostructureProperties () override;
    virtual void Define_MPI_BlkType () override;
    virtual void AllocateArrays() override;
    virtual void ConstructHamiltonian () override;
    virtual void Define_EnergyLimits () override;
    virtual void Define_IntegrationPaths () override;

    virtual void Define_ContactInfo() override;
    virtual AMREX_GPU_HOST_DEVICE
    void Compute_SurfaceGreensFunction (MatrixBlock<BlkType>& gr, 
                                        const ComplexType EmU) override;

    amrex::Real Get_Bandgap_Of_Mode(int p);

    //virtual void DeallocateArrays() override;
    //virtual void ComputeChargeDensity() override;

    ComplexType get_beta(int J);

    virtual void Deallocate() override;

};

#endif
