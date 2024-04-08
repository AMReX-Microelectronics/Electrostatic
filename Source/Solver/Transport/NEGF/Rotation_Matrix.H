#ifndef ROTATION_MATRIX_H_
#define ROTATION_MATRIX_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"


enum class AngleType : int { Degrees, Radians };
enum class AxisType : int { X, Y, Z };

struct RotationInputParams {
    amrex::Vector<amrex::Real> angles;
    amrex::Vector<AxisType> rotation_order;
    AngleType angle_type = AngleType::Degrees;
};


class c_RotationMatrix
{
    std::unique_ptr<RotationInputParams> _p;

    /*X,Y,Z rotation angles*/
    amrex::Vector<amrex::Vector<amrex::Real>> _rotX;
    amrex::Vector<amrex::Vector<amrex::Real>> _rotY;
    amrex::Vector<amrex::Vector<amrex::Real>> _rotZ;

    void Define_RotationMatrices();


    void Rotate(amrex::Vector<amrex::Real>& v2,
                const amrex::Vector<amrex::Real> v1,
                const amrex::Vector<amrex::Vector<amrex::Real>>& rotM);

    void Apply_RotationOrder(amrex::Vector<amrex::Real>& v,
                             const std::vector<AxisType>& order);

public:

    c_RotationMatrix(std::unique_ptr<RotationInputParams> rotationParams);

    void Print_RotationParams() const;

    const RotationInputParams* Get_RotationParamsPtr() const { return _p.get(); }
    amrex::Vector<amrex::Real> Get_RotationAnglesInDegrees() const;
    amrex::Vector<AxisType> Get_RotationOrder() const { return _p->rotation_order; }

    void Update_RotationAngles(amrex::Vector<amrex::Real>& new_angles)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(new_angles.size() == AMREX_SPACEDIM,
                                     "Angles vector must be of size AMREX_SPACEDIM");
        _p->angles = new_angles;
    }

    void Update_RotationOrder(amrex::Vector<AxisType>& new_order) {_p->rotation_order = new_order; }

    void Update_RotationMatrices(amrex::Vector<amrex::Real>& new_angles);

    template<typename ContainerType, typename OffsetType>
    void RotateContainer(ContainerType& v1, const OffsetType& offset=OffsetType{AMREX_D_DECL(0.,0.,0.)});

};


template<typename ContainerType, typename OffsetType = amrex::Vector<amrex::Real>>
void
c_RotationMatrix::RotateContainer(ContainerType& v1, 
                                  const OffsetType& offset)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(v1.size() == AMREX_SPACEDIM, "v1 size must be AMREX_SPACEDIM");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(offset.size() == AMREX_SPACEDIM, "v1 size must be AMREX_SPACEDIM");

    amrex::Vector<amrex::Real> v(AMREX_SPACEDIM,0.);

    for(int i=0; i<AMREX_SPACEDIM; ++i) v[i] = v1[i] - offset[i];

    Apply_RotationOrder(v, _p->rotation_order);

    for(int i=0; i<AMREX_SPACEDIM; ++i) v1[i] = v[i] + offset[i];
}

#endif
