#ifndef ROTATION_MATRIX_H_
#define ROTATION_MATRIX_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

enum class AngleType : int { Degrees, Radians };
enum class AxisType : int { X, Y, Z };

template<typename T>
class c_RotationMatrix
{
    /*X,Y,Z rotation angles*/
    amrex::Vector<T> _angles;
    amrex::Vector<amrex::Vector<T>> _rotX;
    amrex::Vector<amrex::Vector<T>> _rotY;
    amrex::Vector<amrex::Vector<T>> _rotZ;

    void Define_RotationMatrices();

    void Set_RotationAngles(amrex::Vector<T>& new_angles);

    template<typename ContainerType>
    void Rotate(ContainerType& v2,
                const ContainerType v1,
                const amrex::Vector<amrex::Vector<T>>& rotM);

    template<typename ContainerType>
    void Apply_RotationOrder(ContainerType& v,
                             const std::vector<AxisType>& order);

public:

    c_RotationMatrix(amrex::Vector<T> angles,
                     AngleType type); //default type AngleType::Degrees

    amrex::Vector<T> Get_RotationAngles() const { 
        return _angles;
    }

    void Update_RotationMatrices(amrex::Vector<T>& new_angles);

    template<typename ContainerType>
    void RotateContainer(const ContainerType& v1,
                         const std::vector<AxisType>& order); 
                      //default order = {} 

};
#endif
