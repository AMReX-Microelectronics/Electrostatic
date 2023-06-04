#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "Transport_fwd.H"
#include "Nanostructure.H"

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>

#include <string>

struct s_NS {
    enum Type { CNT, Graphene, Silicon };
};

class 
c_TransportSolver
{

public:
    amrex::Vector<std::unique_ptr<c_Nanostructure<c_CNT>>> vp_CNT;
    amrex::Vector<std::unique_ptr<c_Nanostructure<c_Graphene>>> vp_Graphene;
    //amrex::Vector<std::unique_ptr<c_Nanostructure<c_Silicon>>> vp_Silicon;

    c_TransportSolver();
    ~c_TransportSolver();

    void InitData();
    void Solve();

private:
    int num_NS=0;
    amrex::Geometry  dummy_geom;
    amrex::BoxArray  dummy_ba;
    amrex::DistributionMapping  dummy_dm;

    amrex::Geometry const * _geom = &dummy_geom;
    amrex::BoxArray const * _ba = &dummy_ba;
    amrex::DistributionMapping const * _dm = &dummy_dm;
    //const amrex::GpuArray<int, AMREX_SPACEDIM> * _n_cell;
    amrex::Vector< std::string > vec_NS_names;

    void ReadData();

    std::map<std::string, s_NS::Type>
    map_NSType_enum = 
        {
           {"carbon nanotube" , s_NS::Type::CNT},
           {"Carbon Nanotube" , s_NS::Type::CNT},
           {"CNT"             , s_NS::Type::CNT},
           {"graphene"        , s_NS::Type::Graphene},
           {"Graphene"        , s_NS::Type::Graphene},
           {"silicon"         , s_NS::Type::Silicon},
           {"Silicon"         , s_NS::Type::Silicon}
        };

    std::string NS_gather_field_str;
    std::string NS_deposit_field_str;
    amrex::Real NS_initial_deposit_value;
    amrex::Real Broyden_max_norm = 1.e-5;

    int use_selfconsistent_potential;
    int use_negf;

};

#endif
