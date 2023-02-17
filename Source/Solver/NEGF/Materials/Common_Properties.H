#ifndef COMMON_PROPERTIES_H_
#define COMMON_PROPERTIES_H_

#include "../../../Utils/SelectWarpXUtils/WarpXUtil.H"

#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_REAL.H>
#include <AMReX.H>
//#include <AMReX_GpuContainers.H>

using namespace amrex;
struct s_AVG_TYPE 
{
    enum { 
          ALL,
          SPECIFIC
    };
};

class c_Common_Properties
{

public:

    std::string name;
    int num_atoms;
    int num_unitcells;
    amrex::RealArray offset;
    int avg_type;
    amrex::Vector<int> vec_avg_indices;

    virtual void ReadNanostructureProperties ();

};

#endif
