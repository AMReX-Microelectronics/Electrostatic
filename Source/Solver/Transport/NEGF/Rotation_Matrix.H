#ifndef ROTATION_MATRIX_H_
#define ROTATION_MATRIX_H_

#include "matrix.hpp"
#include "nvector.hpp"
#include "ops.hpp"
#include "rotationMatrix.hpp"
#include <cmath>

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

    void Define_RotationMatrices() 
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

    void Rotate(amrex::Vector<T>& v2,
                const amrex::Vector<T> v1, 
                const amrex::Vector<amrex::Vector<T>>& rotM) 
    {
        assert(v1.size() == 3);
        assert(v2.size() == 3);

        for(int i=0; i<3; ++i) {
            for(int j=0; j<3; ++j) {
                v2[i] += rotM[i][j] * v1[j];
            }
        }
    }
 
    void ApplyRotationOrder(amrex::Vector<T>& v,
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

public:

    c_RotationMatrix(amrex::Vector<T> angles, 
                   AngleType type = AngleType::Degrees) :
    _angles(type == AngleType::Degrees ? angles * M_PI / 180 : angles),
    _rotX(3, amrex::Vector<T>(3)),
    _rotY(3, amrex::Vector<T>(3)),
    _rotZ(3, amrex::Vector<T>(3))
    {
        static_assert(angles.size() == 3, "Angles vector must be of size 3");
        Define_RotationMatrices();
    }

    void Update_RotationAngles(amrex::Vector<T> new_angles) 
    {
        static_assert(new_angles.size() == 3, "Angles vector must be of size 3");
        _angles = new_angles;

        Define_RotationMatrices();
    }

    amrex::Vector<T> RotateVector(const amrex::Vector<T>& v1, 
                                  const std::vector<AxisType>& order = {})
    {
        assert(v1.size() == 3);

        if (order.empty()) return; 

        ApplyRotationOrder(order, v1);
    }

};
