#ifndef COMMON_PROPERTIES_H_
#define COMMON_PROPERTIES_H_

#include "../../../Utils/SelectWarpXUtils/WarpXUtil.H"

#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_REAL.H>

#include <string>

class c_Common_Properties
{

public:

    std::string name;
    int num_atoms;
    int num_unitcells;
    amrex::RealArray offset;

    virtual void ReadNanostructureProperties ();

};

#endif
