#include "IntegrationPath.H"

#include "../../../Utils/CodeUtils/CodeUtil.H"

#include <AMReX_Print.H>

#include <fstream>

void
c_IntegrationPath::Define_GaussLegendrePoints(const ComplexType min, 
                                              const ComplexType max,
                                              const int degree,
                                              const int id) 
{
    type_id = id;
    num_pts = degree;

    amrex::Vector<amrex::Real> x, w;
    x.resize(degree);
    w.resize(degree); 
    Quadrature::Gauss_Legendre(x, w, degree);
 
    if(id==0) /*line*/
    {
           
    }
    else if (type_id == 1) /*circle*/
    {
  
    }

    //amrex::Print() << "Gauss-Legendre roots for degree " << degree << "\n";
    //std::string filename = "GL_Points";
    //std::ofstream outfile;
    //outfile.open(filename.c_str());

    //for (int i=0; i<degree; ++i) 
    //{
    //    outfile   << std::setprecision(15) 
    //              << std::setw(25) << x[i] 
    //              << std::setw(25) << w[i] << "\n";
    //}
    //outfile.close();
}


void
c_IntegrationPath:: Define_RegularPoints(ComplexType min, 
                                         ComplexType max, int pts)
{
    E_min = min;
    E_max = max;
    num_pts = pts;

    ComplexType dE = (E_max-E_min)/static_cast<amrex::Real>(num_pts-1);

    weight_vec.resize(num_pts);
    mul_factor_vec.resize(num_pts);
    E_vec.resize(num_pts);

    for (int i=0; i<num_pts; ++i)
    {
        E_vec[i] = E_min + static_cast<amrex::Real>(i)*dE;
        weight_vec[i] = dE;
        mul_factor_vec[i] = 1.;
    }
}
