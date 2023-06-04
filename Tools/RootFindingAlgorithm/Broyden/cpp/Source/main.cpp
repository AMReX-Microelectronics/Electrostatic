
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_MFParallelForC.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuComplex.H>
#include<AMReX_Print.H>
#include<AMReX_TableData.H>
#include <AMReX_TinyProfiler.H>
#include <AMReX_GpuUtility.H>

#include <cmath>
#include <iomanip>
#include <math.h> 
using namespace amrex;
using RealTable1D    = TableData<amrex::Real, 1>;
using RealTable2D    = TableData<amrex::Real, 2>;

#define NUM_ROOTS 2

template<typename T>
class TD;

namespace MathConst
{
    static constexpr amrex::Real pi = static_cast<amrex::Real>(3.14159265358979323846);
}

namespace PhysConst
{
    static constexpr auto q_e   = static_cast<amrex::Real>( 1.602176634e-19 );
    static constexpr auto ep0   = static_cast<amrex::Real>( 8.8541878128e-12 );
//    static constexpr auto hbar  = static_cast<amrex::Real>( 1.054571817e-34 );
}

//amrex::GpuComplex<amrex::Real> get_beta(amrex::Real gamma, int M, int J) 
//{
//   amrex::GpuComplex arg(0., -MathConst::pi*J/M);
//   return 2 * gamma * cos(-1*arg.imag()) * exp(arg); 
//}
//
//amrex::GpuComplex<amrex::Real> conjugate(amrex::GpuComplex<amrex::Real> a) 
//{
//   amrex::GpuComplex a_conj(a.real(), -1.*a.imag());
//   return a_conj;
//}
//
//void PrintTable(const TableData<MatrixDType, 2>& G)
//{
//    auto const& h_table = G.table();
//    auto tlo = G.lo();
//    auto thi = G.hi();
//    for (int i = tlo[0]; i <= thi[0]; ++i) { 
//        for (int j = tlo[1]; j <= thi[1]; ++j) { //slow moving index. printing slow
//            amrex::Print() <<std::setw(12) << std::setprecision(2) << h_table(i,j);
//        }
//        amrex::Print() << "\n";
//    }
//}
template<typename U, typename V>
void SetVal_Table1D (U& Tab1D_data, V val)
{
    auto tlo = Tab1D_data.lo();
    auto thi = Tab1D_data.hi();

    auto const& Tab1D = Tab1D_data.table();

    for (int i = tlo[0]; i < thi[0]; ++i)
    {
        Tab1D(i) = val;
    }
}

template<typename U, typename V>
void
SetVal_Table2D (U& Tab2D_data, V val)
{
    auto tlo = Tab2D_data.lo();
    auto thi = Tab2D_data.hi();

    auto const& Tab2D = Tab2D_data.table();

    for (int i = tlo[0]; i < thi[0]; ++i)
    {
        for (int j = tlo[1]; j < thi[1]; ++j) //slow moving index. printing slow
        {
            Tab2D(i,j) = val;
        }
    }
}

void Define_InitGuess(RealTable1D& h_n_curr_in_data)
{ 
    auto const& n_curr_in  = h_n_curr_in_data.table();
    for( int n=0; n < NUM_ROOTS; ++n) {
        switch(n) 
        {	   
            case 0:
                n_curr_in(0) = 1;
                break;

            case 1:
                n_curr_in(1) = 2;         
                break;
	    default:
		amrex::Print() << "equation not defined!\n";
        }
    }
}


void Compute_FofX(RealTable1D& h_n_curr_in_data,
                  RealTable1D& n_curr_out_data)
{ 
    auto const& n_curr_in  = h_n_curr_in_data.table();
    auto const& n_curr_out = n_curr_out_data.table();
    for( int n=0; n < NUM_ROOTS; ++n) {
        switch(n) 
        {	   
            case 0:
                n_curr_out(0) = pow(n_curr_in(0),2) - n_curr_in(1) - 1;         
                break;

            case 1:
                n_curr_out(1) = n_curr_in(0) - pow(n_curr_in(1),2) + 1;         
                break;
	    default:
		amrex::Print() << "equation not defined!\n";
        }
    }
}


amrex::Real Guess_BroydenFirstAlg(const int Broyden_Step, 
		                  const amrex::Real Broyden_fraction,
		                  RealTable1D& h_n_curr_in_data,
                                  RealTable1D& n_curr_out_data,
				  RealTable1D& n_prev_in_data,
				  RealTable1D& F_curr_data,
				  RealTable2D& Jinv_curr_data)
{
    amrex::Print() << "\nBroydenStep: " << Broyden_Step << "\n";

    auto const& n_curr_in  = h_n_curr_in_data.table();
    auto const& n_curr_out = n_curr_out_data.table();
    auto const& n_prev_in  = n_prev_in_data.table();
    auto const& F_curr     = F_curr_data.table();

    RealTable1D sum_Fcurr_data({0},{NUM_ROOTS}, The_Pinned_Arena());
    RealTable1D sum_deltaFcurr_data({0},{NUM_ROOTS}, The_Pinned_Arena());
    RealTable1D delta_F_curr_data({0},{NUM_ROOTS}, The_Pinned_Arena());
    RealTable1D delta_n_Jinv_data({0},{NUM_ROOTS}, The_Pinned_Arena());

    RealTable1D Norm_data({0},{NUM_ROOTS}, The_Pinned_Arena());

    auto const& sum_Fcurr      = sum_Fcurr_data.table();
    auto const& sum_deltaFcurr = sum_deltaFcurr_data.table();
    auto const& delta_F_curr   = delta_F_curr_data.table();
    auto const& delta_n_Jinv   = delta_n_Jinv_data.table();
    auto const& Jinv_curr    = Jinv_curr_data.table();
    auto const& Norm     = Norm_data.table();

    amrex::Real denom = 0.;
    int m = Broyden_Step-1;

    for(int l=0; l < NUM_ROOTS; ++l)
    {
        amrex::Real Fcurr = n_curr_out(l);

        delta_F_curr(l) = Fcurr - F_curr(l);

        F_curr(l) = Fcurr;

        Norm(l) = fabs(Fcurr);

        sum_deltaFcurr(l) = 0;
        sum_Fcurr(l) = 0;
	delta_n_Jinv(l) = 0.;
    }

    if(m > 0 ) 
    {
        for(int a=0; a < NUM_ROOTS; ++a)
        {
            amrex::Real sum = 0.;
            for(int b=0; b < NUM_ROOTS; ++b)
            {
                sum += Jinv_curr(a,b)*delta_F_curr(b);
            }
            sum_deltaFcurr(a) += sum;
        }

        for(int l=0; l < NUM_ROOTS; ++l)
        {
            denom += (n_curr_in(l) - n_prev_in(l)) * sum_deltaFcurr(l);
        }

        for(int b=0; b < NUM_ROOTS; ++b)
        {
            amrex::Real sum = 0.;
            for(int a=0; a < NUM_ROOTS; ++a)
            {
                sum += (n_curr_in(a) - n_prev_in(a))*Jinv_curr(a,b);
            }
            delta_n_Jinv(b) += sum;
        }

        for(int a=0; a < NUM_ROOTS; ++a)
        {
            for(int b=0; b < NUM_ROOTS; ++b)
            {
                Jinv_curr(a,b) += ( (n_curr_in(a) - n_prev_in(a)) - sum_deltaFcurr(a) )  * delta_n_Jinv(b) / denom;
            }
        }
    }
    amrex::Print() << "\nJ_inverse: \n";
    amrex::Print() << "  " << Jinv_curr(0,0) << "  " << Jinv_curr(0,1) << "\n";
    amrex::Print() << "  " << Jinv_curr(1,0) << "  " << Jinv_curr(1,1) << "\n";

    amrex::Real j_denom = ( Jinv_curr(0,0)*Jinv_curr(1,1) - Jinv_curr(0,1)*Jinv_curr(1,0) );
    amrex::Print() << "j_denom: " << j_denom << "\n";
    amrex::Print() << "\nJ: \n";
    amrex::Print() << "  " <<  Jinv_curr(1,1)/j_denom << "  " << -Jinv_curr(0,1)/j_denom << "\n";
    amrex::Print() << "  " << -Jinv_curr(1,0)/j_denom << "  " <<  Jinv_curr(0,0)/j_denom << "\n";

    for(int a=0; a < NUM_ROOTS; ++a)
    {
        amrex::Real sum = 0.;
        for(int b=0; b < NUM_ROOTS; ++b)
        {
            sum += Jinv_curr(a,b)*F_curr(b);
        }
        sum_Fcurr(a) += sum;
    }

    amrex::Real max_norm = fabs(n_curr_in(0) - n_curr_in(1));
    amrex::Print() << "max_norm: " << max_norm << "\n";


    for(int l=0; l < NUM_ROOTS; ++l)
    {
        n_prev_in(l) = n_curr_in(l);
        n_curr_in(l) = n_prev_in(l) - sum_Fcurr(l);
    }

    sum_Fcurr_data.clear();
    sum_deltaFcurr_data.clear();
    delta_F_curr_data.clear();
    delta_n_Jinv_data.clear();

    for(int l=0; l < NUM_ROOTS; ++l)
    {
        amrex::Print() << "Xn_in, Xn_out, Xn+1_in, norm: " << l << "  " << n_prev_in(l) 
		                                                << "  " << n_curr_out(l) 
		                                                << "  " << n_curr_in(l)
                                                                << "  " << Norm(l) << "\n";
    }

    //std::string filename = "norm_" + std::to_string(Broyden_Step) + ".dat";
    //Write_Table1D(PTD, Norm_data, filename.c_str(), 
    //              "'axial location / (nm)', 'norm");
    Norm_data.clear();

    return max_norm;
}




amrex::Real Guess_ModifiedBroydenSecondAlg(const int Broyden_Step, 
		                    const amrex::Real Broyden_fraction,
		                    RealTable1D& h_n_curr_in_data,
                                    RealTable1D& n_curr_out_data,
				    RealTable1D& n_prev_in_data,
				    RealTable1D& F_curr_data,
				    amrex::Vector<RealTable1D*>& W_Broyden,
				    amrex::Vector<RealTable1D*>& V_Broyden)
{
    /*update h_RhoInduced_perAtom_glo*/

    amrex::Print() << "\nBroydenStep: " << Broyden_Step << "\n";

    auto const& n_curr_in  = h_n_curr_in_data.table();
    auto const& n_curr_out = n_curr_out_data.table();
    auto const& n_prev_in  = n_prev_in_data.table();
    auto const& F_curr     = F_curr_data.table();

    RealTable1D sum_Fcurr_data({0},{NUM_ROOTS}, The_Pinned_Arena());
    RealTable1D sum_deltaFcurr_data({0},{NUM_ROOTS}, The_Pinned_Arena());
    RealTable1D delta_F_curr_data({0},{NUM_ROOTS}, The_Pinned_Arena());
    RealTable1D W_curr_data({0},{NUM_ROOTS}, The_Pinned_Arena());
    RealTable1D V_curr_data({0},{NUM_ROOTS}, The_Pinned_Arena());
    RealTable1D Norm_data({0},{NUM_ROOTS}, The_Pinned_Arena());

    auto const& sum_Fcurr      = sum_Fcurr_data.table();
    auto const& sum_deltaFcurr = sum_deltaFcurr_data.table();
    auto const& delta_F_curr   = delta_F_curr_data.table();
    auto const& W_curr   = W_curr_data.table();
    auto const& V_curr   = V_curr_data.table();
    auto const& Norm     = Norm_data.table();

    amrex::Real denom = 0.;
    int m = Broyden_Step-1;

    //std::string filename = "n_in_" + std::to_string(Broyden_Step) + ".dat";
    //Write_Table1D(PTD, h_n_curr_in_data, filename.c_str(), 
    //              "'axial location / (nm)', 'Induced charge density / (m^3)");


    for(int l=0; l < NUM_ROOTS; ++l)
    {
        amrex::Real Fcurr = n_curr_out(l);

        delta_F_curr(l) = Fcurr - F_curr(l);

        F_curr(l) = Fcurr;

        denom += pow(delta_F_curr(l),2.);

        Norm(l) = fabs(n_curr_in(l)-n_prev_in(l));

        sum_deltaFcurr(l) = 0;
        sum_Fcurr(l) = 0;
        W_curr(l) = 0;
        V_curr(l) = 0;

    }
    //amrex::Print() << "denom: " << denom<< "\n";

    W_Broyden.push_back(new RealTable1D({0},{NUM_ROOTS}, The_Pinned_Arena()));
    V_Broyden.push_back(new RealTable1D({0},{NUM_ROOTS}, The_Pinned_Arena()));

    if(m > 0)
    {
        for(int j=1; j <= m-1; ++j)
        {
            auto const& W_j = W_Broyden[j]->table();
            auto const& V_j = V_Broyden[j]->table();

            for(int a=0; a < NUM_ROOTS; ++a)
            {
                amrex::Real sum = 0.;
                for(int b=0; b < NUM_ROOTS; ++b)
                {
                sum += W_j(a)*V_j(b)*delta_F_curr(b);
                }
                sum_deltaFcurr(a) += sum;
            }
        }

        for(int l=0; l < NUM_ROOTS; ++l)
        {
            amrex::Real delta_n = n_curr_in(l) - n_prev_in(l);

              V_curr(l) = delta_F_curr(l)/denom;
              W_curr(l) = -Broyden_fraction*delta_F_curr(l) + delta_n - sum_deltaFcurr(l);
        }

        W_Broyden[m]->copy(W_curr_data);
        V_Broyden[m]->copy(V_curr_data);

        auto const& W_m = W_Broyden[m]->table();
        auto const& V_m = V_Broyden[m]->table();
        //amrex::Print() << "W_curr/W_Broyden_m: " << W_curr(0) << " " << W_m(0) << "\n";
        //amrex::Print() << "V_curr/V_Broyden_m: " << V_curr(10) << " " << V_m(10) << "\n";
        //if(m-1 > 0) {
        //auto const& W_mMinus1 = W_Broyden[m-1]->table();
        //auto const& V_mMinus1 = V_Broyden[m-1]->table();
        //amrex::Print() << "W_Broyden_m-1: " << W_mMinus1(0) << "\n";
        //amrex::Print() << "V_Broyden_m-1: " << V_mMinus1(10) << "\n";
        //}

        for(int j=1; j <= m; ++j)
        {
            auto const& W_j = W_Broyden[j]->table();
            auto const& V_j = V_Broyden[j]->table();

            for(int a=0; a < NUM_ROOTS; ++a)
            {
                amrex::Real sum = 0.;
                for(int b=0; b < NUM_ROOTS; ++b)
                {
                    sum += W_j(a)*V_j(b)*F_curr(b);
                }
                sum_Fcurr(a) += sum;
            }
        }
    }
    amrex::Print() << "size of W and V: " << W_Broyden.size() << "\n";


    amrex::Real max_norm = fabs(n_curr_in(0) - n_curr_in(1));
    amrex::Print() << "max_norm: " << max_norm << "\n";
    //for(int l=1; l < NUM_ROOTS; ++l)
    //{
    //    //if(fabs(max_norm - Norm(l)) > 1e-10) 
    //    if(max_norm < Norm(l))
    //    {
    //        max_norm = Norm(l);
    //    }
    //}
    for(int l=0; l < NUM_ROOTS; ++l)
    {
        n_prev_in(l) = n_curr_in(l);
        n_curr_in(l) = n_prev_in(l) - Broyden_fraction*F_curr(l) - sum_Fcurr(l);
    }

    //filename = "n_next_" + std::to_string(Broyden_Step) + ".dat";
    //Write_Table1D(PTD, h_n_curr_in_data, filename.c_str(), 
    //              "'axial location / (nm)', 'Induced charge density / (m^3)");

    sum_Fcurr_data.clear();
    sum_deltaFcurr_data.clear();
    delta_F_curr_data.clear();
    W_curr_data.clear();
    V_curr_data.clear();

    for(int l=0; l < NUM_ROOTS; ++l)
    {
        amrex::Print() << "Xn_in, Xn_out, Xn+1_in, norm: " << l << " "<< n_prev_in(l) << "  " << n_curr_out(l) << "  " << n_curr_in(l)
                                                      << "  " << Norm(l) << "\n";
    }


    //std::string filename = "norm_" + std::to_string(Broyden_Step) + ".dat";
    //Write_Table1D(PTD, Norm_data, filename.c_str(),
    //              "'axial location / (nm)', 'norm");
    Norm_data.clear();

    return max_norm;
}



int main (int argc, char* argv[])
{

    amrex::Initialize(argc,argv);

    int Broyden_Step = 1;
    amrex::Real Broyden_fraction = 0.1;
    //amrex::Real init_guess = 0.5;

    RealTable1D h_n_curr_in_data({0},{NUM_ROOTS}, The_Pinned_Arena());;
    //#if AMREX_USE_GPU
    //d_n_curr_in_data.resize({0},{NUM_ROOTS}, The_Arena());
    //#endif
    RealTable1D n_curr_out_data({0},{NUM_ROOTS}, The_Pinned_Arena());;
    RealTable1D n_prev_in_data({0},{NUM_ROOTS}, The_Pinned_Arena());;
    RealTable1D F_curr_data({0},{NUM_ROOTS}, The_Pinned_Arena());;

    amrex::Vector<RealTable1D*> W_Broyden;
    amrex::Vector<RealTable1D*> V_Broyden;

    //SetVal_Table1D(h_n_curr_in_data,init_guess);
    Define_InitGuess(h_n_curr_in_data);
    SetVal_Table1D(n_curr_out_data,0.);
    SetVal_Table1D(n_prev_in_data,0.);
    SetVal_Table1D(F_curr_data,0.);
    
    RealTable2D Jinv_curr_data({0,0},{NUM_ROOTS, NUM_ROOTS}, The_Pinned_Arena());
    SetVal_Table2D(Jinv_curr_data,0.);

    auto const& Jinv_curr    = Jinv_curr_data.table();
    /*In this example, we know J is: ([2x, -1]; [1, -2y]). We can use the inverse of this matrix as initial guess for Jinv. 
     *This is only required for the first algorithm of Broyden*/
    Jinv_curr(0,0) =  4./7.; 
    Jinv_curr(0,1) = -1./7.;
    Jinv_curr(1,0) =  1./7.;
    Jinv_curr(1,1) = -2./7.;
    /*For the first algorithm any random diagonal matrix as a guess for Jinv doesn't work*/
    //for(int a=0; a < NUM_ROOTS; ++a)
    //{
    //    Jinv_curr(a,a) = Broyden_fraction;
    //}

    int step=0;
    amrex::Real max_norm = 1;
    do
    {
       Compute_FofX(h_n_curr_in_data, n_curr_out_data);	    

       max_norm = Guess_ModifiedBroydenSecondAlg(Broyden_Step, Broyden_fraction, 
		                                 h_n_curr_in_data, n_curr_out_data, n_prev_in_data, F_curr_data,
		                                 W_Broyden, V_Broyden);

//       max_norm = Guess_BroydenFirstAlg(Broyden_Step, Broyden_fraction, 
//		                        h_n_curr_in_data, n_curr_out_data, n_prev_in_data, F_curr_data,
//		                        Jinv_curr_data);
//
       Broyden_Step += 1;
    } while(max_norm > 1e-6);

    /*Broyden*/
    h_n_curr_in_data.clear();
    //d_n_curr_in_data.clear();
    n_curr_out_data.clear();
    n_prev_in_data.clear();
    F_curr_data.clear();

    int size = W_Broyden.size();
    for(int j=0; j<size; ++j)
    {
        W_Broyden[j]->clear();
        V_Broyden[j]->clear();
    }
    W_Broyden.clear();
    V_Broyden.clear();
 
    amrex::Finalize();

}
