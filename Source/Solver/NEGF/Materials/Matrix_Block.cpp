#include "Matrix_Block.H"


template struct MatrixBlock<ComplexType[NUM_MODES]>;            //of c_CNT
template struct MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>; //of c_Graphene

/* Definitions */
/* Operation [R] = c_complex, i.e. a complex constant */
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator=(const ComplexType c_comp) { return *this;}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator=(const ComplexType c_comp)
{
   for(int i=0; i < NUM_MODES; ++i)
   {   
       this->block[i] = c_comp;
   }
   return *this;
}


template<>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::operator=(const ComplexType c_comp)
{
   for(int i=0; i < NUM_MODES; ++i)
   {   
       for(int j=0; j < NUM_MODES; ++j)
       {   
           this->block[i][j] = c_comp;
       }
   }
   return *this;
}



/* Operation [R] = c_real, i.e. a real constant */
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator=(const amrex::Real c) { amrex::Print() << "HEY! \n"; return *this;}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator=(const amrex::Real c)
{
   ComplexType c_complex(c, 0.);
   *this = c_complex;
   return *this;
}


template<>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::operator=(const amrex::Real c)
{
   ComplexType c_complex(c, 0.);
   *this = c_complex;
   return *this;
}

/* Operation [R] = [B]*c_real */
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator*(const amrex::Real c) 
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator*(const amrex::Real c)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {   
       result.block[i] = this->block[i]*c;
   }
   return result;
   /*if R = A*1; In this case &A is this, 1 is m, and S is the returned value*/
}


template<>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::operator*(const amrex::Real c)
{
   MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       for(int j=0; j < NUM_MODES; ++j)
       {
           result.block[i][j] = this->block[i][j]*c;
       }
   }
   return result;
   /*if S = A*1; In this case &A is this, 1 is m, and S is the returned value*/
}


/* Operation [R] = [B]*[C] */
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator*(const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator*(const MatrixBlock<ComplexType[NUM_MODES]>& rhs)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i]*rhs.block[i];
   }
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::operator*(const MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>& rhs)
{
   MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       for(int j=0; j < NUM_MODES; ++j)
       {   
           result.block[i][j] = this->block[i][j]*rhs.block[i][j];
       }
   }
   return result;
}


/* Operation [R] = c_real*[B] */
template<typename T>
MatrixBlock<T> operator*(const amrex::Real c, const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]> operator*(const amrex::Real c, const MatrixBlock<ComplexType[NUM_MODES]>& rhs)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = c*rhs.block[i];
   }
   return result;
}



/* Operation [R] = [B] + c_complex */
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator+(const ComplexType c)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator+(const ComplexType c)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] + c;
   }
   return result;
   /*if S = A*1; In this case &A is this, 1 is m, and S is the returned value*/
}



/* Operation [R] = [B] + [C] */
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator+(const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator+(const MatrixBlock<ComplexType[NUM_MODES]>& rhs)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] + rhs.block[i];
   }
   return result;
}



/* Operation [R] = c_complex + [C] */
template<typename T>
MatrixBlock<T> operator+(const ComplexType c, const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]> operator+(const ComplexType c, const MatrixBlock<ComplexType[NUM_MODES]>& rhs)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = c + rhs.block[i];
   }
   return result;
}



/* Operation [R] = [B] - c_real */
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator-(const amrex::Real c)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator-(const amrex::Real c)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] - c;
   }
   return result;
   /*if S = A*1; In this case &A is this, 1 is m, and S is the returned value*/
}



/* Operation [R] = [B] - c_complex */
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator-(const ComplexType c)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator-(const ComplexType c)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] - c;
   }
   return result; 
   /*if S = A*1; In this case &A is this, 1 is m, and S is the returned value*/
}

 
  
/* Operation [R] = [B] - [C] */
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator-(const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator-(const MatrixBlock<ComplexType[NUM_MODES]>& rhs)
{  
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] - rhs.block[i];
   }
   return result; 
}



/* Operation [R] = c_complex - [C] */
template<typename T>
MatrixBlock<T> operator-(const ComplexType c, const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]> operator-(const ComplexType c, const MatrixBlock<ComplexType[NUM_MODES]>& rhs)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = c - rhs.block[i];
   }
   return result;
}



/* Operation [R] = [B] / c_complex */
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator/(ComplexType c)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator/(ComplexType c)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] / c;
   }
   return result;
   /*if S = A*1; In this case &A is this, 1 is m, and S is the returned value*/
}



/* Operation [R] = [B] / [C] */
template<typename T>
MatrixBlock<T>
MatrixBlock<T>::operator/(const MatrixBlock<T>& rhs)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator/(const MatrixBlock<ComplexType[NUM_MODES]>& rhs)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = this->block[i] / rhs.block[i];
   }
   return result;
}



/* Operation amrex::Print() << [R] */
template<typename T>
std::ostream&
operator<<(std::ostream& stream, const MatrixBlock<T>& rhs) { return stream; }


template<>
std::ostream&
operator<<(std::ostream& stream, const MatrixBlock<amrex::Real[NUM_MODES]>& rhs)
{
    stream << "{";
    for (int i = 0; i < NUM_MODES; ++i)
    {
        if (i != NUM_MODES - 1) {
            stream << rhs.block[i] << ", ";
        } else {
            stream << rhs.block[i];
        }
    }
    stream << "}";
    return stream;
}


template<>
std::ostream&
operator<<(std::ostream& stream, const MatrixBlock<ComplexType[NUM_MODES]>& rhs)
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


template<>
std::ostream&
operator<<(std::ostream& stream, const MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>& rhs)
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


//template<>
//MatrixBlock<amrex::Real[NUM_MODES]> real(const MatrixBlock<T>& rhs)
//{ 
//    amrex::Real[NUM_MODES] block_real;
//
//    for (int i=0; i < NUM_MODES; ++i) 
//    {
//        block_real[i]  = rhs.block[i].real();
//    }
//    return block_real; 
//}
