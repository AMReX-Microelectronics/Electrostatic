#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "Transport_fwd.H"
#include "Nanostructure.H"

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>

#include<AMReX_TableData.H>

#include <string>

struct s_NS {
    enum Type { CNT, Graphene, Silicon };
};

struct s_Algorithm { 
    enum Type { broyden_second, simple_mixing };
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
    void Solve(const int step, const amrex::Real time);

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

    std::string negf_foldername_str = "output/negf";




    /*For Broyden*/
    int num_field_sites_all_NS = 0;

    using RealTable1D    = TableData<amrex::Real, 1>;

    int Broyden_Step = 1;
    amrex::Real Broyden_fraction = 0.1;
    amrex::Real Broyden_Original_Fraction = 0.1;
    amrex::Real Broyden_Norm = 1.0;
    int Broyden_Correction_Step = 0;
    int Broyden_NormSumIsIncreasing_Step = 0;
    int Broyden_Reset_Step = 0;
    amrex::Real Broyden_Scalar = 1.0;
    amrex::Real Broyden_NormSum_Curr = 1.e10;
    amrex::Real Broyden_NormSum_Prev = 1.e100;


    int Broyden_Threshold_MaxStep = 200;
    int Broyden_Threshold_NormSumIncreaseStep = 8;
    int Broyden_Threshold_CorrectionStep = 1;
    amrex::Real Broyden_Threshold_MinFraction = 0.01;
    amrex::Real Broyden_Fraction_Decrease_Factor = 4.;
    amrex::Real Broyden_Scalar_Decrease_Factor   = 2.;

    std::string Algorithm_Type = "broyden_second";
    std::map<std::string, s_Algorithm::Type>
    map_AlgorithmType =
        {
           {"broyden_second"        , s_Algorithm::Type::broyden_second},
           {"Broyden_second"        , s_Algorithm::Type::broyden_second},
           {"simple_mixing"         , s_Algorithm::Type::simple_mixing},
        };
    
    std::string Broyden_Norm_Type = "relative";
    std::map<std::string, s_Norm::Type>
    map_NormType =
        {
           {"absolute"        , s_Norm::Type::Absolute},
           {"relative"        , s_Norm::Type::Relative},
           {"Absolute"        , s_Norm::Type::Absolute},
           {"Relative"        , s_Norm::Type::Relative}
        };
    
    RealTable1D n_curr_in_data;
    RealTable1D n_start_in_data;
    RealTable1D n_curr_out_data;
    RealTable1D n_prev_in_data;
    RealTable1D F_curr_data;
    RealTable1D Norm_data;

    RealTable1D delta_F_curr_data;
    RealTable1D delta_n_curr_data;

    amrex::Vector<RealTable1D*> W_Broyden;
    amrex::Vector<RealTable1D*> V_Broyden;

    void Set_Broyden();
    void Reset_Broyden();
    void Restart_Broyden();
    void Execute_Broyden_Modified_Second_Algorithm();
    void Execute_Simple_Mixing_Algorithm();
    void Deallocate_Broyden();

    void SetVal_RealTable1D (RealTable1D& Tab1D_data, amrex::Real val);    

};

#endif
