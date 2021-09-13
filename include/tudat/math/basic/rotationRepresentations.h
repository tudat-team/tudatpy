/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */


#ifndef TUDAT_ROTATIONREPRESENTATIONS_H
#define TUDAT_ROTATIONREPRESENTATIONS_H

#include <Eigen/Core>
#include <Eigen/Dense>

namespace tudat
{

namespace basic_mathematics
{

Eigen::Matrix< double, 4, 3 > calculateQuaternionWrtEulerAngle313Partial(
        const Eigen::Quaterniond& quaternion );

//! Function to compute the partial derivative of 3-1-3 Euler angles w.r.t. entries of associated quaternion
/*!
 * Function to compute the partial derivative of 3-1-3 Euler angles w.r.t. entries of associated quaternion, with quaternion
 * as input
 * \param quaternion Quaternion defining rotation at which partials are to be evaluated
 * \return Partial derivative matrix of 3-1-3 Euler angles w.r.t. entries of associated quaternion
 */
Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartial(
        const Eigen::Quaterniond& quaternion );

//! Function to compute the partial derivative of 3-1-3 Euler angles w.r.t. entries of associated quaternion
/*!
 * Function to compute the partial derivative of 3-1-3 Euler angles w.r.t. entries of associated quaternion, with Euler angles
 * as input
 * \param eulerAngles Euler angles (3-1-3) defining rotation at which partials are to be evaluated
 * \return Partial derivative matrix of 3-1-3 Euler angles w.r.t. entries of associated quaternion
 */
Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartialFromEulerAngles(
        const Eigen::Vector3d& eulerAngles );

//! Get quaternion from associated 3-1-3 Euler angles
/*!
 * Get  quaternion q from 3-1-3 Euler angles set. That is, the Euler angles x, y, z are defined such that the associated
 * rotation matrix R = R_{3}(x)R_{1}(y)R_{3}(z)
 * \param eulerAngles Euler angles for which the equivalent quaternion is to be computed.
 * \return Quaternion defining same rotation as Euler angles
 */
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

//! Get classical 3-1-3 Euler angles set from quaternion
/*!
 * Get classical 3-1-3 Euler angles set from quaternion q. That is, the Euler angles x, y, z are returned such that the associated
 * rotation matrix R = R_{3}(x)R_{1}(y)R_{3}(z)
 * \param quaternion Quaternion for which the equivalent Euler angles are to be computed.
 * \return Euler angles x,y,z (about 3, 1 and 3 axes, respectively).
 */
Eigen::Vector3d get313EulerAnglesFromQuaternion(
        const Eigen::Quaterniond& quaternion );

//! Get classical 3-1-3 Euler angles set from rotation matrix
/*!
 * Get classical 3-1-3 Euler angles set from rotation matrix R. That is, the Euler angles x, y, z are returned such that
 * R = R_{3}(x)R_{1}(y)R_{3}(z)
 * \param rotationMatrix Rotation matrix for which the equivalent Euler angles are to be computed.
 * \return Euler angles x,y,z (about 3, 1 and 3 axes, respectively).
 */
Eigen::Vector3d get313EulerAnglesFromRotationMatrix(
        const Eigen::Matrix3d& rotationMatrix );
}

}

#endif // TUDAT_ROTATIONREPRESENTATIONS_H
