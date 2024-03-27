#ifndef ROTATION_MATRIX_H_
#define ROTATION_MATRIX_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"


enum class AngleType : int { Degrees, Radians };
enum class AxisType : int { X, Y, Z };

class c_RotationMatrix
{
    /*X,Y,Z rotation angles*/
    amrex::Vector<amrex::Real> _angles;
    AngleType _angle_type;
    amrex::Vector<AxisType> _rotation_order;
    amrex::Vector<amrex::Vector<amrex::Real>> _rotX;
    amrex::Vector<amrex::Vector<amrex::Real>> _rotY;
    amrex::Vector<amrex::Vector<amrex::Real>> _rotZ;

    void Define_RotationMatrices();

    void Set_RotationAngles(amrex::Vector<amrex::Real>& new_angles);

    template<typename ContainerType>
    void Rotate(ContainerType& v2,
                const ContainerType v1,
                const amrex::Vector<amrex::Vector<amrex::Real>>& rotM);

    template<typename ContainerType>
    void Apply_RotationOrder(ContainerType& v,
                             const std::vector<AxisType>& order);

public:

    c_RotationMatrix(amrex::Vector<amrex::Real> angles,
                     amrex::Vector<AxisType> rotation_order,
                     AngleType type); //default type AngleType::Degrees

    amrex::Vector<amrex::Real> Get_RotationAngles() const { 
        return _angles;
    }

    amrex::Vector<AxisType> Get_RotationOrder() const { 
        return _rotation_order;
    }

    void Update_RotationOrder(amrex::Vector<AxisType>& new_order) {
        _rotation_order = new_order;
    }

    void Update_RotationMatrices(amrex::Vector<amrex::Real>& new_angles);


    template<typename ContainerType>
    void RotateContainer(ContainerType& v1);

};



template<typename ContainerType>
void
c_RotationMatrix::Rotate(ContainerType& v2,
                            const ContainerType v1,
                            const amrex::Vector<amrex::Vector<amrex::Real>>& rotM)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(v1.size() == AMREX_SPACEDIM, "v1 size must be AMREX_SPACEDIM");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(v1.size() == AMREX_SPACEDIM, "v2 size must be AMREX_SPACEDIM");

    for(int i=0; i<3; ++i) {
        for(int j=0; j<3; ++j) {
            v2[i] += rotM[i][j] * v1[j];
        }
    }
}


template<typename ContainerType>
void
c_RotationMatrix::Apply_RotationOrder(ContainerType& v,
                        const std::vector<AxisType>& order)
{
    for (AxisType axis : order)
    {
        switch(axis) {
            case AxisType::X:
                Rotate(v, v, _rotX);
                break;
            case AxisType::Y:
                Rotate(v, v, _rotY);
                break;
            case AxisType::Z:
                Rotate(v, v, _rotZ);
                break;
        }
    }
}


template<typename ContainerType>
void
c_RotationMatrix::RotateContainer(ContainerType& v1)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(v1.size() == AMREX_SPACEDIM, "v1 s3ize must be AMREX_SPACEDIM");

    Apply_RotationOrder(v1, _rotation_order);
}

#endif
