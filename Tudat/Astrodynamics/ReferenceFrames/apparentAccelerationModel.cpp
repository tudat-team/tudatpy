/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      J.S. Torok, Analytical Mechanics, Wiley-Interscience, 2000.
 *
 */

#include "Tudat/Astrodynamics/ReferenceFrames/apparentAccelerationModel.h"

namespace tudat
{
namespace reference_frames
{

//! Compute apparent acceleration due to non-inertiality of reference frame.
Eigen::Vector3d computeApparentAcceleration(
        const Eigen::Vector3d& accelerationOfNonInertialReferenceFrame,
        const Eigen::Vector3d& angularVelocityOfNonInertialReferenceFrame,
        const Eigen::Vector3d& angularAccelerationOfNonInertialReferenceFrame,
        const Eigen::Vector3d& positionOfBodyInNonInertialReferenceFrame,
        const Eigen::Vector3d& velocityOfBodyInNonInertialReferenceFrame )
{
    return -accelerationOfNonInertialReferenceFrame +
            computeCentripetalAcceleration( angularVelocityOfNonInertialReferenceFrame,
                                            positionOfBodyInNonInertialReferenceFrame ) +
            computeCoriolisAcceleration( angularVelocityOfNonInertialReferenceFrame,
                                         velocityOfBodyInNonInertialReferenceFrame ) +
            computeEulerAcceleration( angularAccelerationOfNonInertialReferenceFrame,
                                      positionOfBodyInNonInertialReferenceFrame );
}

//! Compute centripetal acceleration due to non-inertiality of reference frame.
Eigen::Vector3d computeCentripetalAcceleration(
        const Eigen::Vector3d& angularVelocityOfNonInertialReferenceFrame,
        const Eigen::Vector3d& positionOfBodyInNonInertialReferenceFrame )
{
    return -angularVelocityOfNonInertialReferenceFrame.cross(
                ( angularVelocityOfNonInertialReferenceFrame.cross(
                      positionOfBodyInNonInertialReferenceFrame ) ) );
}

//! Compute Coriolis acceleration due to non-inertiality of reference frame.
Eigen::Vector3d computeCoriolisAcceleration(
        const Eigen::Vector3d& angularVelocityOfNonInertialReferenceFrame,
        const Eigen::Vector3d& velocityOfBodyInNonInertialReferenceFrame )
{
    return -2.0 * ( angularVelocityOfNonInertialReferenceFrame.cross(
                       velocityOfBodyInNonInertialReferenceFrame ) );
}

//! Compute Euler acceleration due to non-inertiality of reference frame.
Eigen::Vector3d computeEulerAcceleration(
        const Eigen::Vector3d& angularAccelerationOfNonInertialReferenceFrame,
        const Eigen::Vector3d& positionOfBodyInNonInertialReferenceFrame )
{
    return -angularAccelerationOfNonInertialReferenceFrame.cross(
                positionOfBodyInNonInertialReferenceFrame );
}

//! Update member variables used by apparent acceleration model.
void ApparentAccelerationModel::updateMembers( const double currentTime )
{
    if( !( this->currentTime_ == currentTime ) )
    {
        currentAccelerationOfNonInertialReferenceFrame_ =
                accelerationOfNonInertialReferenceFrameFunction_( );
        currentAngularVelocityOfNonInertialReferenceFrame_ =
                angularVelocityOfNonInertialReferenceFrameFunction_( );
        currentAngularAccelerationOfNonInertialReferenceFrame_ =
                angularAccelerationOfNonInertialReferenceFrameFunction_( );
        currentPositionOfBodyInNonInertialReferenceFrame_ =
                positionOfBodyInNonInertialReferenceFrameFunction_( );
        currentVelocityOfBodyInNonInertialReferenceFrame_ =
                velocityOfBodyInNonInertialReferenceFrameFunction_( );
    }
}

} // namespace reference_frames
} // namespace tudat
