#include "Rotation_Matrix.H"

#include "../../../Utils/SelectWarpXUtils/WarpXConst.H"

#include <cmath>


c_RotationMatrix::c_RotationMatrix(amrex::Vector<amrex::Real> angles,
                                   AngleType type = AngleType::Degrees) :
    _angles(angles.size()),
    _rotX(3, amrex::Vector<amrex::Real>(3)),
    _rotY(3, amrex::Vector<amrex::Real>(3)),
    _rotZ(3, amrex::Vector<amrex::Real>(3))
{
    
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(angles.size() == AMREX_SPACEDIM, "Angles vector must be of size AMREX_SPACEDIM");
    for (int i = 0; i < angles.size(); ++i) {
        _angles[i] = (type == AngleType::Degrees) ? angles[i] * (MathConst::pi / 180) : angles[i];
    }

    Define_RotationMatrices();
}


void 
c_RotationMatrix::Define_RotationMatrices() 
{
        _rotX[0] = {1, 0, 0};
        _rotX[1] = {0,  cos(_angles[0]), sin(_angles[0])};
        _rotX[2] = {0, -sin(_angles[0]), cos(_angles[0])};

        _rotY[0] = {cos(_angles[1]), 0, -sin(_angles[1])};
        _rotY[1] = {0, 1, 0};
        _rotY[2] = {sin(_angles[1]), 0, cos(_angles[1])};

        _rotZ[0] = { cos(_angles[2]), sin(_angles[2]), 0};
        _rotZ[1] = {-sin(_angles[2]), cos(_angles[2]), 0};
        _rotZ[2] = {0, 0, 1};
}


void 
c_RotationMatrix::Set_RotationAngles(amrex::Vector<amrex::Real>& new_angles) 
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(new_angles.size() == AMREX_SPACEDIM, 
                                     "Angles vector must be of size AMREX_SPACEDIM");
    _angles = new_angles;

}


void 
c_RotationMatrix::Update_RotationMatrices(amrex::Vector<amrex::Real>& new_angles) 
{
    Set_RotationAngles(new_angles);

    Define_RotationMatrices();
}
