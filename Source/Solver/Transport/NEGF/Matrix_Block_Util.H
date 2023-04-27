#include "Matrix_Block.H"


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
MatrixBlock<T>::operator*(const ComplexType c) 
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::operator*(const ComplexType c)
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
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::operator*(const ComplexType c)
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


/* Operation [R] = c_complex*[B] */
template<typename T>
MatrixBlock<T> operator*(const ComplexType c_complex, const MatrixBlock<T>& B)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]> operator*(const ComplexType c_complex, const MatrixBlock<ComplexType[NUM_MODES]>& B)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = c_complex*B.block[i];
   }
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> operator*(const ComplexType c_complex, const MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>& B)
{
   MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       for(int j=0; j < NUM_MODES; ++j)
       {
           result.block[i][j] = c_complex*B.block[i][j];
       }
   }
   return result;
}


/* Operation [R] = c_real*[B] */
template<typename T>
MatrixBlock<T> operator*(const amrex::Real c_real, const MatrixBlock<T>& B)
{
   ComplexType c_complex(c_real, 0.);
   return operator*(c_complex, B);
}


template<>
MatrixBlock<ComplexType[NUM_MODES]> operator*(const amrex::Real c_real, 
                                       const MatrixBlock<ComplexType[NUM_MODES]>& B)
{
   ComplexType c_complex(c_real, 0.);
   return operator*(c_complex,B);
}


template<>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> operator*(const amrex::Real c_real, 
                                const MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>& B)
{
   ComplexType c_complex(c_real, 0.);
   return operator*(c_complex,B);
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


template<>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::operator+(const MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>& rhs)
{
   MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       for(int j=0; j < NUM_MODES; ++j)
       {
           result.block[i][j] = this->block[i][j] + rhs.block[i][j];
       }  
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


template<>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> operator+(const ComplexType c, const MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>& rhs)
{
   MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       for(int j=0; j < NUM_MODES; ++j)
       {
           result.block[i][j] = c + rhs.block[i][j];
       }  
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


/* Operation [R] = [B] - [C] */
template<typename T>
MatrixBlock<T> operator-(const MatrixBlock<T>& B, const MatrixBlock<T>& C)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]> operator-(const MatrixBlock<ComplexType[NUM_MODES]>& B, const MatrixBlock<ComplexType[NUM_MODES]>& C)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = B.block[i] - C.block[i];
   }
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> operator-(const MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>& B, const MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>& C)
{
   MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       for(int j=0; j < NUM_MODES; ++j)
       {
           result.block[i][j] = B.block[i][j] - C.block[i][j];
       }
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

/* Operation [R] = c_complex / [C] */
template<typename T>
MatrixBlock<T> operator/(const ComplexType c_complex, const MatrixBlock<T>& C)
{
   MatrixBlock<T> result; 
   return result;
}


template<>
MatrixBlock<ComplexType[NUM_MODES]> operator/(const ComplexType c_complex, 
                                              const MatrixBlock<ComplexType[NUM_MODES]>& C)
{
   MatrixBlock<ComplexType[NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       result.block[i] = c_complex / C.block[i];
   }
   return result;
}

template<>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> operator/(const ComplexType c_complex, 
                                              const MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>& C)
{
   MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> result;
   for(int i=0; i < NUM_MODES; ++i)
   {
       for(int j=0; j < NUM_MODES; ++j)
       {
           result.block[i][j] = c_complex / C.block[i][j];
       }
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
            stream << rhs.block[i]  << ", ";
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
            stream << std::setw(10) << rhs.block[i][j];
        }
        if(i != NUM_MODES-1) stream << "\n ";
    }
    stream << "]";
    return stream;
}


template <typename T>
MatrixBlock<T> 
MatrixBlock<T>::Conj() const
{
    MatrixBlock<T> result;
    return result;
}

template <>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::Conj() const
{
    MatrixBlock<ComplexType[NUM_MODES]> result; 
    for (int i = 0; i < NUM_MODES; ++i)
    {
        amrex::GpuComplex X_conj(this->block[i].real(), -1.*this->block[i].imag());
        result.block[i] = X_conj;
    }
    return result;
}

template <>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::Conj() const
{
    MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> result; 
    for (int i = 0; i < NUM_MODES; ++i)
    {
        for (int j = 0; j < NUM_MODES; ++j)
        {
            amrex::GpuComplex X_conj(this->block[i][j].real(), -1.*this->block[i][j].imag());
            result.block[i][j] = X_conj;
        }
    }
    return result;
}


template <typename T>
MatrixBlock<T> 
MatrixBlock<T>::Tran() const
{
    return *this;
}

template <>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::Tran() const
{
    return *this;
}

template <>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::Tran() const
{
    MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> result; 
    for (int i = 0; i < NUM_MODES; ++i)
    {
        for (int j = 0; j < NUM_MODES; ++j)
        {
            if(i==j) 
            {
                result.block[i][j] = this->block[i][j];
            }
            else 
            {
                result.block[i][j] = this->block[j][i];
            } 
        }
    }
    return result;
}


template <typename T>
MatrixBlock<T> 
MatrixBlock<T>::Dagger() const
{
    return this->Conj().Tran();
}

template <>
MatrixBlock<ComplexType[NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES]>::Dagger() const
{
    return this->Conj().Tran();
}

template <>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::Dagger() const
{
    return this->Conj().Tran();
}


///* Operation AtomicAdd */
//template <typename T>
//MatrixBlock<T> 
//MatrixBlock<T>::AtomicAdd(MatrixBlock<T>& B) const
//{
//    return *this;
//}
//
//template <>
//MatrixBlock<ComplexType[NUM_MODES]>
//MatrixBlock<ComplexType[NUM_MODES]>::AtomicAdd(MatrixBlock<T>& B) const
//{
//    
//    for (int i = 0; i < NUM_MODES; ++i)
//    {
//        amrex::HostDevice::Atomic::Add(&(this->trace_r.block[i]), B.block[i].real());
//        amrex::HostDevice::Atomic::Add(&(this->trace_i.block[i]), B.block[i].imag());
//    }
//    return result;
//}
//
//template <>
//MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>
//MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::AtomicAdd(MatrixBlock<T>& B) const
//{
//    MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]> result; 
//    for (int i = 0; i < NUM_MODES; ++i)
//    {
//        for (int j = 0; j < NUM_MODES; ++j)
//        {
//            if(i==j) 
//            {
//                result.block[i][j] = this->block[i][j];
//            }
//            else 
//            {
//                result.block[i][j] = this->block[j][i];
//            } 
//        }
//    }
//    return result;
//}

//template <typename T, std::size_t N>
//MatrixBlock<amrex::Real[N]> 
//MatrixBlock<T>::Real() const
//{
//    using N = sizeof(T);
//    MatrixBlock<amrex::Real> result;
//    return result;
//}
//
//template <>
//MatrixBlock<amrex::Real[NUM_MODES]>
//MatrixBlock<ComplexType[NUM_MODES]>::Real() const
//{
//    MatrixBlock<amrex::Real[NUM_MODES]> result;
//
//    for (int i = 0; i < NUM_MODES; ++i)
//    {
//        result.block[i] = this->block[i].real();
//    }
//    return result;
//}
//
//template <>
//MatrixBlock<amrex::Real[NUM_MODES][NUM_MODES]>
//MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::Real() const
//{
//    MatrixBlock<amrex::Real[NUM_MODES][NUM_MODES]> result; 
//
//    for (int i = 0; i < NUM_MODES; ++i)
//    {
//        for (int j = 0; j < NUM_MODES; ++j)
//        {
//            result.block[i][j] = this->block[i][j].real();
//        }
//    }
//    return result;
//}


/* Operation DiagSum */
template <typename T>
ComplexType
MatrixBlock<T>::DiagSum() const
{
    ComplexType result(0.,0.);
    return *this;
}

template <>
ComplexType
MatrixBlock<ComplexType[NUM_MODES]>::DiagSum() const
{
    
    ComplexType result(0.,0.);
    for (int i = 0; i < NUM_MODES; ++i)
    {
        result += this->block[i];
    }
    return result;
}

template <>
ComplexType
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::DiagSum() const
{
    ComplexType result(0.,0.);
    for (int i = 0; i < NUM_MODES; ++i)
    {
        for (int j = 0; j < NUM_MODES; ++j)
        {
            if(i==j) 
            {
                result += this->block[i][j];
            }
        }
    }
    return result;
}


///* Operation DiagDotSum */
template <typename T>
template <typename U>
ComplexType
MatrixBlock<T>::DiagDotSum(U& vec) const
{
    ComplexType result(0.,0.);
    return *this;
}

template <>
template <typename U>
ComplexType
MatrixBlock<ComplexType[NUM_MODES]>::DiagDotSum(U& vec) const
{
    
    ComplexType result(0.,0.);
    for (int i = 0; i < NUM_MODES; ++i)
    {
        ComplexType factor(vec[i],0.);
        result += factor*this->block[i];
    }
    return result;
}

template <>
template <typename U>
ComplexType
MatrixBlock<ComplexType[NUM_MODES][NUM_MODES]>::DiagDotSum(U& vec) const
{
    ComplexType result(0.,0.);
    for (int i = 0; i < NUM_MODES; ++i)
    {
        for (int j = 0; j < NUM_MODES; ++j)
        {
            if(i==j) 
            {
                ComplexType factor(vec[i],0.);
                result += factor*this->block[i][j];
            }
        }
    }
    return result;
}
