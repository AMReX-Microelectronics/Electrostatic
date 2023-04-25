#include "Graphene.H"
#include "../../Utils/CodeUtils/CodeUtil.H"

#include <AMReX_Particles.H>

#include <cmath>
#include <math.h>
#include<stdlib.h>


amrex::Array<int,2> c_Graphene::type_id;

void
c_Graphene::Define_SelfEnergy ()
{

    c_NEGF_Common<BlkType>:: Define_SelfEnergy ();
//    auto const& h_Sigma = h_Sigma_glo_data.table();
//
//    amrex::Print() << "Printing Sigma: \n";
//    std::cout<< h_Sigma(0,0) << "\n";
//
}
