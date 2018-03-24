/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include <stdexcept>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/External/SpiceInterface/spiceRotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{


//! Function to calculate the rotation quaternion from target frame to original frame.
Eigen::Quaterniond SpiceRotationalEphemeris::getRotationToBaseFrame(
        const double secondsSinceEpoch )
{
    // Get rotational quaternion from spice wrapper function
    return spice_interface::computeRotationQuaternionBetweenFrames(
                targetFrameOrientation_, baseFrameOrientation_, secondsSinceEpoch + referenceDayOffSet_ );
}

//! Function to calculate the derivative of the rotation matrix from target frame to original
//! frame.
Eigen::Matrix3d SpiceRotationalEphemeris::getDerivativeOfRotationToBaseFrame(
        const double secondsSinceEpoch )
{
    // Get rotation matrix derivative from spice wrapper function
    return spice_interface::computeRotationMatrixDerivativeBetweenFrames(
                targetFrameOrientation_, baseFrameOrientation_, secondsSinceEpoch + referenceDayOffSet_ );
}

//! Function to calculate the full rotational state at given time
void SpiceRotationalEphemeris::getFullRotationalQuantitiesToTargetFrame(
        Eigen::Quaterniond& currentRotationToLocalFrame,
        Eigen::Matrix3d& currentRotationToLocalFrameDerivative,
        Eigen::Vector3d& currentAngularVelocityVectorInGlobalFrame,
        const double secondsSinceEpoch)
{
    // Calculate rotation (and its time derivative) directly from spice.
    std::pair< Eigen::Quaterniond, Eigen::Matrix3d > fullRotation =
            spice_interface::computeRotationQuaternionAndRotationMatrixDerivativeBetweenFrames(
                baseFrameOrientation_, targetFrameOrientation_, secondsSinceEpoch + referenceDayOffSet_ );
    currentRotationToLocalFrame = fullRotation.first;
    currentRotationToLocalFrameDerivative = fullRotation.second;

    // Calculate angular velocity vector.
    currentAngularVelocityVectorInGlobalFrame = getRotationalVelocityVectorInBaseFrameFromMatrices(
                Eigen::Matrix3d( currentRotationToLocalFrame ), currentRotationToLocalFrameDerivative.transpose( ) );
}


} // namespace ephemerides

} // namespace tudat
