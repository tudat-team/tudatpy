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

#include <limits>

#include <Eigen/Geometry>

namespace tudat
{

namespace basic_mathematics
{

Eigen::Matrix< double, 4, 3 > calculateQuaternionWrtEulerAngle313Partial( const Eigen::Quaterniond& quaternion );

//! Function to compute the partial derivative of 3-1-3 Euler angles w.r.t. entries of associated quaternion
/*!
 * Function to compute the partial derivative of 3-1-3 Euler angles w.r.t. entries of associated quaternion, with quaternion
 * as input
 * \param quaternion Quaternion defining rotation at which partials are to be evaluated
 * \return Partial derivative matrix of 3-1-3 Euler angles w.r.t. entries of associated quaternion
 */
Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartial( const Eigen::Quaterniond& quaternion );

//! Function to compute the partial derivative of 3-1-3 Euler angles w.r.t. entries of associated quaternion
/*!
 * Function to compute the partial derivative of 3-1-3 Euler angles w.r.t. entries of associated quaternion, with Euler angles
 * as input
 * \param eulerAngles Euler angles (3-1-3) defining rotation at which partials are to be evaluated
 * \return Partial derivative matrix of 3-1-3 Euler angles w.r.t. entries of associated quaternion
 */
Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartialFromEulerAngles( const Eigen::Vector3d& eulerAngles );

//! Get quaternion from associated 3-1-3 Euler angles
/*!
 * Get  quaternion q from 3-1-3 Euler angles set. That is, the Euler angles x, y, z are defined such that the associated
 * rotation matrix R = R_{3}(x)R_{1}(y)R_{3}(z)
 * \param eulerAngles Euler angles for which the equivalent quaternion is to be computed.
 * \return Quaternion defining same rotation as Euler angles
 */
Eigen::Quaterniond getQuaternionFrom313EulerAngles( const Eigen::Vector3d& eulerAngles );

//! Get quaternion from associated 3-2-1 Euler angles
/*!
 * Get  quaternion q from 3-2-1 Euler angles set (Roll-Pitch-Yaw). That is, the Euler angles x, y, z are defined such that the associated
 * rotation matrix R = R_{3}(x)R_{2}(y)R_{1}(z)
 * \param eulerAngles Euler angles for which the equivalent quaternion is to be computed.
 * \return Quaternion defining same rotation as Euler angles
 */
Eigen::Quaterniond getQuaternionFrom321EulerAngles( const Eigen::Vector3d& eulerAngles );

//! Get classical 1-3-2 Euler angles set from rotation matrix
/*!
 * Get classical 1-3-2 Euler angles set from rotation matrix R. That is, the Euler angles x, y, z are returned such that
 * R = R_{1}(x)R_{3}(y)R_{2}(z)
 * \param rotationMatrix Rotation matrix for which the equivalent Euler angles are to be computed.
 * \return Euler angles x,y,z (about 1, 3 and 2 axes, respectively).
 */
Eigen::Vector3d get132EulerAnglesFromRotationMatrix( const Eigen::Matrix3d& rotationMatrix );

//! Get classical 3-1-3 Euler angles set from quaternion
/*!
 * Get classical 3-1-3 Euler angles set from quaternion q. That is, the Euler angles x, y, z are returned such that the associated
 * rotation matrix R = R_{3}(x)R_{1}(y)R_{3}(z)
 * \param quaternion Quaternion for which the equivalent Euler angles are to be computed.
 * \return Euler angles x,y,z (about 3, 1 and 3 axes, respectively).
 */
Eigen::Vector3d get313EulerAnglesFromQuaternion( const Eigen::Quaterniond& quaternion );

//! Get classical 3-1-3 Euler angles set from rotation matrix
/*!
 * Get classical 3-1-3 Euler angles set from rotation matrix R. That is, the Euler angles x, y, z are returned such that
 * R = R_{3}(x)R_{1}(y)R_{3}(z)
 * \param rotationMatrix Rotation matrix for which the equivalent Euler angles are to be computed.
 * \return Euler angles x,y,z (about 3, 1 and 3 axes, respectively).
 */
Eigen::Vector3d get313EulerAnglesFromRotationMatrix( const Eigen::Matrix3d& rotationMatrix );

//! Get quaternion from associated rotation vector (exponential map of SO(3))
/*!
 * Get quaternion q from a rotation vector, i.e. a vector whose direction defines the rotation axis and
 * whose norm defines the rotation angle (in radians). This is the exponential map Exp: so(3) -> SO(3),
 * and is the representation used for small attitude corrections that are to be estimated, since it is
 * free of the constraints and singularities of quaternion/Euler-angle parameterizations near zero.
 * For a vanishing rotation angle, the identity quaternion is returned.
 * \param rotationVector Rotation vector for which the equivalent quaternion is to be computed.
 * \return Quaternion defining same rotation as rotation vector
 */
inline Eigen::Quaterniond getQuaternionFromRotationVector( const Eigen::Vector3d& rotationVector )
{
    const double rotationAngle = rotationVector.norm( );
    if( rotationAngle <= std::numeric_limits< double >::epsilon( ) )
    {
        return Eigen::Quaterniond::Identity( );
    }
    return Eigen::Quaterniond( Eigen::AngleAxisd( rotationAngle, rotationVector / rotationAngle ) );
}

//! Get rotation vector from associated quaternion (logarithmic map of SO(3))
/*!
 * Get the rotation vector associated with quaternion q, i.e. the inverse of getQuaternionFromRotationVector.
 * The returned vector has the rotation axis as its direction and the rotation angle (in radians, in [0, pi])
 * as its norm.
 * \param quaternion Quaternion for which the equivalent rotation vector is to be computed.
 * \return Rotation vector defining same rotation as quaternion
 */
inline Eigen::Vector3d getRotationVectorFromQuaternion( const Eigen::Quaterniond& quaternion )
{
    const Eigen::AngleAxisd angleAxis( quaternion.normalized( ) );
    return angleAxis.angle( ) * angleAxis.axis( );
}
}  // namespace basic_mathematics

}  // namespace tudat

#endif  // TUDAT_ROTATIONREPRESENTATIONS_H
