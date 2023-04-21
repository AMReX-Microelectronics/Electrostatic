#ifndef GRAPHENE_H_
#define GRAPHENE_H_

#include "Common_Properties.H"

#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Array.H>

#include <string>
#include <array>


class c_Graphene
           :public c_Common_Properties<ComplexType[NUM_MODES][NUM_MODES]>
{

public:
    using BlkType = ComplexType[NUM_MODES][NUM_MODES];


    std::array<int,2> Get_Ring_Atom_IDs(int cell_id)
    {    
         std::array<int,2> id;
         return id;
    };

    struct Get1DLayerID {
       AMREX_GPU_HOST_DEVICE int operator() (int global_par_id) const
       {
           return static_cast<int>(amrex::Math::floor((global_par_id-1)/type_id[0]));
       }
       amrex::Array<int,2> type_id;
    };

    static Get1DLayerID get_1Dlayer_id() {
         return Get1DLayerID{type_id};
    }

    struct GetAtomIn1DLayerID {
       AMREX_GPU_HOST_DEVICE int operator() (int global_par_id) const
       {
           return static_cast<int>((global_par_id-1) % type_id[0]);                                  
       }
       amrex::Array<int,2> type_id;
    };

    static GetAtomIn1DLayerID get_atom_in_1Dlayer_id() {
         return GetAtomIn1DLayerID{type_id};            
    }

    const int get_num_layers() 
    {
         int num_layers;
         return num_layers;
    };

    const int get_num_atoms_per_layer()
    {
         int atoms_per_layer;
         return atoms_per_layer;
    };

    void Print_IDs(int cell_id)
    {
        //
    };

    static amrex::Array<int,2> type_id;


    //template
    virtual int get_Hsize() override { return 0; }
    virtual int get_offDiag_repeatBlkSize() override { return 2; }

    //virtual void ReadNanostructureProperties () override;
    //virtual void AllocateArrays() override;
    //virtual void ConstructHamiltonian () override;
    //virtual void AddPotentialToHamiltonian () override;
    //virtual void DefineIntegrationPaths () override;
    virtual void Define_SelfEnergy() override;
    //virtual void DeallocateArrays() override;
    //virtual void ComputeChargeDensity() override;

};


#endif
