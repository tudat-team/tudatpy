#ifndef TUDAT_ROTATIONREPRESENTATIONS_H
#define TUDAT_ROTATIONREPRESENTATIONS_H

#include <Eigen/Core>
#include <Eigen/Dense>

namespace tudat
{

namespace basic_mathematics
{

Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartial(
        const Eigen::Quaterniond& quaternion );

Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartialFromEulerAngles(
        const Eigen::Vector3d& eulerAngles );

Eigen::Quaterniond getQuaternionFrom313EulerAngles(
        const Eigen::Vector3d& eulerAngles );

//! Get classical 1-3-2 Euler angles set from rotation matrix
/*!
 * Get classical 1-3-2 Euler angles set from rotation matrix R. That is, the Euler angles x, y, z are returned such that
 * R = R_{1}(x)R_{3}(y)R_{2}(z)
 * \param rotationMatrix Rotation matrix for which the equivalent Euler angles are to be computed.
 * \return Euler angles x,y,z (about 1, 3 and 2 axes, respectively).
 */
Eigen::Vector3d get132EulerAnglesFromRotationMatrix(
        const Eigen::Matrix3d& rotationMatrix );

Eigen::Vector3d get313EulerAnglesFromQuaternion(
        const Eigen::Quaterniond& quaternion );

Eigen::Vector3d get313EulerAnglesFromRotationMatrix(
        const Eigen::Matrix3d& rotationMatrix );
}

}

#endif // TUDAT_ROTATIONREPRESENTATIONS_H
