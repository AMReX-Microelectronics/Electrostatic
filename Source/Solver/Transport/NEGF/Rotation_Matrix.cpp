#include "Rotation_Matrix.H"

#include "../../../Utils/SelectWarpXUtils/WarpXConst.H"

#include <cmath>

c_RotationMatrix::c_RotationMatrix(std::unique_ptr<RotationInputParams> rotationParams) :
    _p {std::move(rotationParams)},
    _rotX{3, amrex::Vector<amrex::Real>{3}},
    _rotY{3, amrex::Vector<amrex::Real>{3}},
    _rotZ{3, amrex::Vector<amrex::Real>{3}}
{
    
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(_p->angles.size() == AMREX_SPACEDIM, "Angles vector must be of size AMREX_SPACEDIM");

    if(_p->angle_type == AngleType::Degrees) 
        for (auto& angle : _p->angles) angle *= (MathConst::pi / 180);

    Define_RotationMatrices();

    Print_RotationParams();

    //test
    //amrex::Vector<amrex::Real> test_vec = {1., 1., 1.};
    //RotateContainer(test_vec);
    //amrex::Print() << "##### test_vec(1.,1.,1.) rotated to:\n";
    //for(auto val: test_vec) amrex::Print() << val<< " ";
    //amrex::Print() << "\n";
}


void 
c_RotationMatrix::Define_RotationMatrices() 
{
        _rotX[0] = {1, 0, 0};
        _rotX[1] = {0,  cos(_p->angles[0]), sin(_p->angles[0])};
        _rotX[2] = {0, -sin(_p->angles[0]), cos(_p->angles[0])};

        _rotY[0] = {cos(_p->angles[1]), 0, -sin(_p->angles[1])};
        _rotY[1] = {0, 1, 0};
        _rotY[2] = {sin(_p->angles[1]), 0, cos(_p->angles[1])};

        _rotZ[0] = { cos(_p->angles[2]), sin(_p->angles[2]), 0};
        _rotZ[1] = {-sin(_p->angles[2]), cos(_p->angles[2]), 0};
        _rotZ[2] = {0, 0, 1};
}


amrex::Vector<amrex::Real> 
c_RotationMatrix::Get_RotationAnglesInDegrees() const
{
    amrex::Vector<amrex::Real> anglesD(_p->angles.size());

    for (int i = 0; i < _p->angles.size(); ++i) {
        anglesD[i] = _p->angles[i] * (180/MathConst::pi);
    }

    return anglesD;
}


void 
c_RotationMatrix::Print_RotationParams() const
{
    amrex::Print() << "##### read angle type: "
                   << (_p->angle_type == AngleType::Degrees ? "Degrees" : "Radians")
                   << "\n";

    amrex::Vector<amrex::Real> anglesD = Get_RotationAnglesInDegrees();
    amrex::Print() << "##### rotation angles (in degrees): ";
    for (auto angle: anglesD) amrex::Print() << angle << " ";
    amrex::Print() << "\n";

    amrex::Print() << "##### rotation order vector: ";
    for (AxisType axis : _p->rotation_order)
    {
        switch(axis) {
            case AxisType::X:
                amrex::Print() << "X ";
                break;
            case AxisType::Y:
                amrex::Print() << "Y ";
                break;
            case AxisType::Z:
                amrex::Print() << "Z ";
                break;
        }
    }
    amrex::Print() << "\n";
}


void 
c_RotationMatrix::Update_RotationMatrices(amrex::Vector<amrex::Real>& new_angles) 
{
    Update_RotationAngles(new_angles);

    Define_RotationMatrices();
}


void
c_RotationMatrix::Rotate(amrex::Vector<amrex::Real>& v2,
                            const amrex::Vector<amrex::Real> v1,
                            const amrex::Vector<amrex::Vector<amrex::Real>>& rotM)
{
    for(int i=0; i<3; ++i) {
        amrex::Real sum=0.;
        for(int j=0; j<3; ++j) {
            sum += rotM[i][j] * v1[j];
        }
        v2[i] = sum;
    }
}


void
c_RotationMatrix::Apply_RotationOrder(amrex::Vector<amrex::Real>& v,
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
