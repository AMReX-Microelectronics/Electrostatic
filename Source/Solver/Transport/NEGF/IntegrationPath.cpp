#include "IntegrationPath.H"

#include "../../../Utils/CodeUtils/CodeUtil.H"

#include <AMReX_Print.H>

#include <fstream>

//c_IntegrationPath::c_IntegrationPath(const c_IntegrationPath& src)
//   :c_IntegrationPath {} //delegating constructor 
//{
//    //deep copy
//}
//
//void 
//c_IntegrationPath::swap(c_IntegrationPath& other) noexcept
//{
//    std::swap(m_var, other.m_var);
//}
//
//c_IntegrationPath& 
//c_IntegrationPath::operator=(const c_IntegrationPath& rhs)
//{
//    if(this == &rhs) { return *this }; //checking for self-assignment
//
//    c_IntegrationPath temp{rhs};
//    swap(temp);
//    return *this;
//}

void 
c_IntegrationPath::Reset()
{
    E_min = 0.;
    E_max = 0.;
    type_id = 0;
    weight_vec.clear();
    mul_factor_vec.clear();
    E_vec.clear();
}

void
c_IntegrationPath::Define_GaussLegendrePoints(const ComplexType min, 
                                              const ComplexType max,
                                              const int degree,
                                              const int id) 
{
    weight_vec.clear();
    mul_factor_vec.clear();
    E_vec.clear();

    amrex::Vector<amrex::Real> x, w;
    x.resize(degree);
    w.resize(degree); 
    Quadrature::Gauss_Legendre(x, w, degree);

    type_id = id;
    num_pts = degree;
    E_vec.resize(num_pts);
    mul_factor_vec.resize(num_pts);
    weight_vec.resize(num_pts);
 
    if(id==0) /*line*/
    {
        for(int i=0; i<num_pts; ++i) 
        {
            E_vec[i]          = (max-min)*0.5*x[i] + (max+min)*0.5;
            mul_factor_vec[i] = (max-min)*0.5;
            weight_vec[i]     = w[i]; 
        }           
        //amrex::Print() << "Gauss-Legendre line points for degree " << num_pts << "\n";
        //std::string filename = "GL_Line_Points";
        //std::ofstream outfile;
        //outfile.open(filename.c_str());

        //for (int i=0; i<num_pts; ++i) 
        //{
        //    outfile   << std::setprecision(15) 
        //              << std::setw(25) << E_vec[i] 
        //              << std::setw(25) << mul_factor_vec[i] 
        //              << std::setw(25) << weight_vec[i] << "\n";
        //}
        //outfile.close();
    }
    else if (type_id == 1) /*circle*/
    {
        amrex::Real E_center_real = (max.real()*max.real() - min.real()*min.real() - min.imag()*min.imag()) /
                                    (2*(max.real()-min.real()));

        amrex::Real E_center_imag = max.imag(); /*this is assumed*/

        amrex::Real R = E_center_real - max.real();

        ComplexType E_center(E_center_real, max.imag());
        //amrex::Print() << "E_center: " << E_center << "\n";
        //amrex::Print() << "R: "        << R << "\n";

        amrex::Real theta1 = asin(fabs(min.imag())/R);
        if(E_center_real > min.real()) theta1 = MathConst::pi - theta1;

        amrex::Real theta2 = MathConst::pi;
  
        //amrex::Print() << "theta1: "<< theta1 << "\n";

        for(int i=0; i<num_pts; ++i) 
        {
            amrex::Real theta_i  = (theta2-theta1)*0.5*x[i] + (theta2+theta1)*0.5;

            amrex::Real real_part = E_center_real + R*cos(theta_i);
            amrex::Real imag_part = E_center_imag + R*sin(theta_i);

            ComplexType E_i(real_part, imag_part);

            E_vec[i] = E_i;  

            ComplexType arg(0., theta_i);
            ComplexType iota(0., 1);
            mul_factor_vec[i] = (theta2-theta1)*0.5*R*exp(arg)*iota;

            weight_vec[i]     = w[i]; 
        }           
        //amrex::Print() << "Gauss-Legendre circ points for degree " << num_pts << "\n";
        //std::string filename = "GL_Circ_Points";
        //std::ofstream outfile;
        //outfile.open(filename.c_str());

        //for (int i=0; i<num_pts; ++i) 
        //{
        //    outfile   << std::setprecision(15) 
        //              << std::setw(25) << E_vec[i] 
        //              << std::setw(25) << mul_factor_vec[i] 
        //              << std::setw(25) << weight_vec[i] << "\n";
        //}
        //outfile.close();
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
