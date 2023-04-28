#ifndef GRAPHENE_H_
#define GRAPHENE_H_

#include "NEGF_Common.H"

#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Array.H>

#include <string>
#include <array>


class c_Graphene
           :public c_NEGF_Common<ComplexType[NUM_MODES][NUM_MODES]>
{

public:
    using BlkType = ComplexType[NUM_MODES][NUM_MODES];

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

    //template
    virtual int get_Hsize() override { return 0; }
    virtual int get_offDiag_repeatBlkSize() override { return 2; }

    //virtual void ReadNanostructureProperties () override;
    //virtual void AllocateArrays() override;
    //virtual void ConstructHamiltonian () override;
    //virtual void AddPotentialToHamiltonian () override;
    //virtual void DefineIntegrationPaths () override;
    //virtual void DeallocateArrays() override;
    //virtual void ComputeChargeDensity() override;

};


#endif
