#ifndef INTEGRATION_PATH_H_
#define INTEGRATION_PATH_H_

#include "../../../Code_Definitions.H"

class c_IntegrationPath
{
    private:
    ComplexType E_min  = 0.;
    ComplexType E_max  = 0.;
    int type_id     = 0; /*0=line, 1=circle*/
    //ComplexType E_center  = 0.; 
    //amrex::Real E_radius  = 0.;

    public:
    int num_pts = 0;
    amrex::Vector<ComplexType> weight_vec;
    amrex::Vector<ComplexType> mul_factor_vec;
    amrex::Vector<ComplexType> E_vec;

    void Define_GaussLegendrePoints(const ComplexType min, const ComplexType max,
                                    const int degree,
                                    const int id);

    void Define_RegularPoints(ComplexType min, ComplexType max, int pts);
};

#endif
