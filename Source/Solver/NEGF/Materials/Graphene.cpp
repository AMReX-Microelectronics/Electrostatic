#include "Graphene.H"
#include <AMReX_Particles.H>

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"


amrex::Array<int,2> c_Graphene::type_id;

template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator*(const amrex::Real m)
{
   MatrixBlock<T> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       for(int j=0; j < NUM_MODES; ++j)
       {
           result.block[i][j] = this->block[i][j]*m;
       }
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
       for(int j=0; j < NUM_MODES; ++j)
       {
           result.block[i][j] = this->block[i][j]*rhs.block[i][j];
       }
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
        for (int j = 0; j < NUM_MODES; ++j)
        {
            stream << std::setw(10)<< rhs.block[i][j];
        } 
        if(i != NUM_MODES-1) stream << "\n ";
    }
    stream << "]";
    return stream;
}


void
c_Graphene:: ReadNanostructureProperties ()
{
}


void
c_Graphene::DefineMatrixPartition(int num_proc) {}


void
c_Graphene::AllocateArrays ()
{
//    h_Ha_loc_data.resize({0},{blkCol_size_loc}, The_Pinned_Arena());

    /*specific to mode-space approximation*/
    h_Hb_loc_data.resize({0},{2},The_Pinned_Arena());
    h_Hc_loc_data.resize({0},{2},The_Pinned_Arena());
    /*end*/
}


void
c_Graphene:: ConstructHamiltonian ()
{
}

void 
c_Graphene::AddPotentialToHamiltonian() 
{
}

//void
//c_Graphene::DefineSelfEnergy ()
//{
//
//    MatrixBlock<BlkType> test_sigma;
//    ComplexType EmU(1, 0.);
//
//    MatrixBlock<BlkType> gr;
//    for(int i=0; i< NUM_MODES; ++i) {
//        for(int j=0; j< NUM_MODES; ++j) {
//            ComplexType val(1, -1);
//            gr.block[i][j] = val;
//        }
//    }
//    amrex::Print() << "Printing gr for Graphene: \n";
//    amrex::Print() << gr << "\n";
//    test_sigma = gr*gr;
//
//    amrex::Print() << "Printing sigma for Graphene: \n";
//    amrex::Print() << test_sigma << "\n";
//
//}
