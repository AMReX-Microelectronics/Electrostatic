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
    
//    c_IntegrationPath() = default;
//    c_IntegrationPath(const c_IntegrationPath& rhs);
//    c_IntegrationPath& operator=(const c_IntegrationPath& rhs);
//    void swap(c_IntegrationPath& other) noexcept;

    public:
    int num_pts = 0;
    amrex::Vector<ComplexType> weight_vec;
    amrex::Vector<ComplexType> mul_factor_vec;
    amrex::Vector<ComplexType> E_vec;

    void Define_GaussLegendrePoints(const ComplexType min, const ComplexType max,
                                    const int degree,
                                    const int id);

    void Define_RegularPoints(ComplexType min, ComplexType max, int pts);

    void Reset(); 
};
//void swap(c_IntegrationPath& first, c_IntegrationPath& second) noexcept 
//{
//    first.swap(second)
//}
#endif
