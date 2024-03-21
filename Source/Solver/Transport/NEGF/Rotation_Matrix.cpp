#include "Rotation_Matrix.H"

#include <cmath>


template<typename T>
c_RotationMatrix<T>::c_RotationMatrix(amrex::Vector<T> angles,
                                   AngleType type = AngleType::Degrees) :
    _angles(type == AngleType::Degrees ? angles * M_PI / 180 : angles),
    _rotX(3, amrex::Vector<T>(3)),
    _rotY(3, amrex::Vector<T>(3)),
    _rotZ(3, amrex::Vector<T>(3))
{
    static_assert(angles.size() == 3, "Angles vector must be of size 3");
    Define_RotationMatrices();
}


template<typename T>
void 
c_RotationMatrix<T>::Define_RotationMatrices() 
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


template<typename T>
template<typename ContainerType>
void 
c_RotationMatrix<T>::Rotate(ContainerType& v2,
                            const ContainerType v1, 
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


template<typename T>
template<typename ContainerType>
void 
c_RotationMatrix<T>::Apply_RotationOrder(ContainerType& v,
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


template<typename T>
void 
c_RotationMatrix<T>::Set_RotationAngles(amrex::Vector<T>& new_angles) 
{
    static_assert(new_angles.size() == 3, "Angles vector must be of size 3");
    _angles = new_angles;

}


template<typename T>
void 
c_RotationMatrix<T>::Update_RotationMatrices(amrex::Vector<T>& new_angles) 
{
    Set_RotationAngles(new_angles);

    Define_RotationMatrices();
}


template<typename T>
template<typename ContainerType>
void
c_RotationMatrix<T>::RotateContainer(const ContainerType& v1, 
                                     const std::vector<AxisType>& order = {})
{
    assert(v1.size() == 3);

    if (order.empty()) return; 

    ApplyRotationOrder(order, v1);
}
