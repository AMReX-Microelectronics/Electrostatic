#ifndef SILICON_H_
#define SILICON_H_

#include "Common_Properties.H"

#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>

#include <AMReX_Array.H>

#include <string>
#include <array>


class c_Silicon
           :public c_Common_Properties
{

public:

    amrex::Real blah;
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
    { //
    };

    static amrex::Array<int,2> type_id;

//    virtual void ReadNanostructureProperties () override;
};

#endif
