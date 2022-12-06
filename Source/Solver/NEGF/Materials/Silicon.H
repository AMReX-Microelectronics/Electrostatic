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
    void Print_IDs(int cell_id)                 
    {
        //
    };
//    virtual void ReadNanostructureProperties () override;
};

#endif
