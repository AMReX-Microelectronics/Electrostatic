#include "CNT.H"
#include "../../Utils/CodeUtils/CodeUtil.H"

#include <AMReX_Particles.H>

#include <cmath>
#include <math.h>
#include<stdlib.h>

amrex::Array<int,2> c_CNT::type_id; //very important

/*First operator loaded functions are defined*/
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator=(const amrex::Real c)
{
   ComplexType c_complex(c, 0.);
   *this = c_complex;
   return *this;
}


template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator=(const ComplexType c_comp)
{
   for(int i=0; i < NUM_MODES; ++i)
   {
       this->block[i] = c_comp;
   }
   return *this;
}


template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator*(const amrex::Real c)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i]*c;
   }
   return result;
   /*if S = A*1; In this case &A is this, 1 is m, and S is the returned value*/
}

template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator*(const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i]*rhs.block[i];
   }
   return result;
}

template<typename T>
MatrixBlock<T> operator*(const amrex::Real c, const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = c*rhs.block[i];
   }
   return result;
}


template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator+(const ComplexType c)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] + c;
   }
   return result;
   /*if S = A*1; In this case &A is this, 1 is m, and S is the returned value*/
}


template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator+(const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] + rhs.block[i];
   }
   return result;
}


template<typename T>
MatrixBlock<T> operator+(const ComplexType c, const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = c + rhs.block[i];
   }
   return result;
}


template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator-(const amrex::Real c)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] - c;
   }
   return result;
   /*if S = A*1; In this case &A is this, 1 is m, and S is the returned value*/
}


template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator-(const ComplexType c)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] - c;
   }
   return result;
   /*if S = A*1; In this case &A is this, 1 is m, and S is the returned value*/
}


template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator-(const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] - rhs.block[i];
   }
   return result;
}


template<typename T>
MatrixBlock<T> operator-(const ComplexType c, const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = c - rhs.block[i];
   }
   return result;
}


template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator/(ComplexType c)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] / c;
   }
   return result;
   /*if S = A*1; In this case &A is this, 1 is m, and S is the returned value*/
}


template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator/(const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] / rhs.block[i];
   }
   return result;
}


template<typename T>
std::ostream&
operator<<(std::ostream& stream, const MatrixBlock<T>& rhs)
{
    stream << "[";
    for (int i = 0; i < NUM_MODES; ++i)
    {
        if (i != NUM_MODES - 1) {
            stream << rhs.block[i] << ", ";
        } else {
            stream << rhs.block[i];
        }
    }
    stream << "]";
    return stream;
}


/*Next, member class functions are defined*/

void
c_CNT:: ReadNanostructureProperties ()
{
    amrex::Print() << "\n##### NANOSTRUCTURE PROPERTIES #####\n\n";

    c_Common_Properties<BlkType>::ReadNanostructureProperties();

    amrex::ParmParse pp_ns(name);

    amrex::Print() << "##### Properties Specific to CNT: \n";

    getWithParser(pp_ns,"acc", acc);
    amrex::Print() << "##### acc: " << acc << "\n";

    amrex::Vector<int> vec_type_id;
    getArrWithParser(pp_ns, "type_id", vec_type_id, 0, 2);
    type_id = vecToArr_Templated<int, 2>(vec_type_id);
    amrex::Print() << "##### type_id: ";
    for (int i=0; i<2; ++i) amrex::Print() << type_id[i] << "  ";
    amrex::Print() << "\n";

    //getWithParser(pp_ns,"num_unitcells", num_unitcells);
    //amrex::Print() << "##### num_unitcells: " << num_unitcells << "\n";

    getWithParser(pp_ns,"gamma", gamma); 
    amrex::Print() << "##### gamma: " << gamma << "\n";

    //amrex::Vector<amrex::Real> vec_offset;
    //getArrWithParser(pp_ns, "offset", vec_offset, 0, AMREX_SPACEDIM);
    //offset = vecToArr(vec_offset);

    //amrex::Print() << "##### offset: ";
    //for (int i=0; i<AMREX_SPACEDIM; ++i) amrex::Print() << offset[i] << "  ";
    //amrex::Print() << "\n";

    if(type_id[1] == 0) {
        rings_per_unitcell = 4;
        atoms_per_ring = type_id[0];
    } else {
        //need to verify this 
        rings_per_unitcell = 2;
        atoms_per_ring = type_id[0]-1;
    }
    num_atoms = num_unitcells*rings_per_unitcell*atoms_per_ring;
    amrex::Print() << "#####* rings_per_unitcell: " << rings_per_unitcell << "\n";
    amrex::Print() << "#####* atoms_per_ring: " << atoms_per_ring << "\n";
    amrex::Print() << "#####* num_atoms: " << num_atoms << "\n";

}


void
c_CNT::DefineMatrixPartition (int num_proc) 
{
    /*Define here Hsize_glo, max_blkCol_perProc*/

    Hsize_glo = get_num_layers();

    bool flag_fixed_blk_size = false;
    const int THRESHOLD_BLKCOL_SIZE = 40000; /*matrix size*/
    
    if(flag_fixed_blk_size) amrex::Print() << "max_blkCol_perProc is fixed by user\n";
    else 
    {
         amrex::Print() << "max_blkCol_perProc is computed at run-time\n";
    
         max_blkCol_perProc = ceil(static_cast<amrex::Real>(Hsize_glo)/num_proc);
         if(max_blkCol_perProc > THRESHOLD_BLKCOL_SIZE) 
         {
            max_blkCol_perProc = THRESHOLD_BLKCOL_SIZE;
            /*assert that use larger number of procs*/
         }
    }
    
}


AMREX_GPU_HOST_DEVICE ComplexType 
c_CNT::get_beta(int J)
{
   ComplexType arg(0., -MathConst::pi*J/type_id[0]);
   return 2. * gamma * cos(-1*arg.imag());// * exp(arg); 
}


//AMREX_GPU_HOST_DEVICE
//void
//c_CNT:: Compute_SurfaceGreensFunction(MatrixBlock<BlkType>& gr, const ComplexType EmU)
//{
//
//   auto EmU_sq = pow(EmU,2.);
//   auto gamma_sq = pow(gamma,2.);
//
//   for(int i=0; i < NUM_MODES; ++i) 
//   {
//       auto Factor = EmU_sq + gamma_sq - pow(beta.block[i],2);
//
//       auto Sqrt = sqrt(pow(Factor,2) - 4. * EmU_sq * gamma_sq);
//
//       auto Denom = 2. * gamma_sq * EmU;
//
//       auto val1 = (Factor + Sqrt) / Denom;
//       auto val2 = (Factor - Sqrt) / Denom;
//
//       if(val1.imag() < 0.) gr.block[i] = val1;
//       else if(val2.imag() < 0.) gr.block[i] = val2; 
//       //amrex::Print() << "EmU: " << EmU << "\n";
//       //amrex::Print() << "Factor: " << Factor << "\n";
//       //amrex::Print() << "Sqrt: "  << Sqrt << "\n";
//       //amrex::Print() << "Denom: " << Denom << "\n";
//       //amrex::Print() << "Numerator: " << Factor+Sqrt << "\n";
//       //amrex::Print() << "Value: " << (Factor+Sqrt)/Denom << "\n";
//   }
//
//}
//
//
//
//AMREX_GPU_HOST_DEVICE 
//void
//c_CNT::get_Sigma(MatrixBlock<BlkType>& Sigma, const ComplexType EmU)
//{
//
//   MatrixBlock<BlkType> gr;
//   Compute_SurfaceGreensFunction(gr, EmU);
//
//   amrex::Print() << "Printing gr: \n";
//   amrex::Print() << gr << "\n";
//
//   Sigma = gr*pow(gamma,2.);
//}
//
//
void 
c_CNT::AllocateArrays () 
{
//    h_Ha_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());

    /*specific to mode-space approximation*/
    h_Hb_loc_data.resize({0},{2},The_Pinned_Arena());
    h_Hc_loc_data.resize({0},{2},The_Pinned_Arena());
    /*end*/
}
//
//
//void 
//c_CNT::DeallocateArrays () 
//{
//
//}
//
//
void 
c_CNT::AddPotentialToHamiltonian() 
{
    auto const& h_Ha = h_Ha_loc_data.table();
    int c=0; 
    for(auto& col_gid: vec_blkCol_gids)
    {
        ComplexType val(avg_gatherField[col_gid],0);
        h_Ha(c) = h_Ha(c) + val;
        ++c;
    }
  
}


void 
c_CNT::ConstructHamiltonian() 
{

    for (int j=0; j<NUM_MODES; ++j) 
    {
        int J = mode_arr[j];
        beta.block[j] = get_beta(J);
    }
    amrex::Print() << "\n Printing beta: "<< "\n";
    amrex::Print() << beta << "\n";

    /*specifying Hb and Hc vectors for the special case of mode-space Hamiltonian of CNT*/
    auto const& h_Hb = h_Hb_loc_data.table();
    auto const& h_Hc = h_Hc_loc_data.table();

    for (std::size_t i = 0; i < 2; ++i)
    {
       if(i%2 == 0) {
          h_Hb(i) = -1*beta; /*negative sign because (E[I] - [H]) will have negative B and C*/
          h_Hc(i) = -1*beta;
       }
       else {
          ComplexType val(-gamma,0);
          h_Hb(i) = val;
          h_Hc(i) = val;  
       }
    }
}


//void 
//c_CNT::DefineIntegrationPaths ()
//{
//
//
//}
//
//
//void 
//c_CNT::DefineSelfEnergy ()
//{
//    //auto const& h_Sigma = h_Sigma_glo_data.table();
//
//    MatrixBlock<BlkType> test_sigma;
//    ComplexType EmU(1, 0.);
//    get_Sigma(test_sigma, EmU);
//    amrex::Print() << "Printing Sigma: \n";
//    std::cout<< test_sigma << "\n";
//
//}
//
//
//void 
//c_CNT::ComputeChargeDensity () 
//{
//
//}
