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

Eigen::Vector3d get313EulerAnglesFromQuaternion(
        const Eigen::Quaterniond& quaternion );
}

}

#endif // TUDAT_ROTATIONREPRESENTATIONS_H
