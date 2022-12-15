#ifndef GRAPHENE_H_
#define GRAPHENE_H_

#include "Common_Properties.H"

#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>

#include <AMReX_Array.H>

#include <string>
#include <array>

class c_Graphene
           :public c_Common_Properties
{

public:

    amrex::Real blah;

    std::array<int,2> Get_Ring_Atom_IDs(int cell_id)
    {    
         std::array<int,2> id;
         return id;
    };

    const int get_1Dlayer_id(int global_par_id)
    {
          return 0;    
    };

    const int get_atom_in_1Dlayer_id(int global_par_id)
    {
         return 0;
    };

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

//    virtual void ReadNanostructureProperties ();
};

#endif
