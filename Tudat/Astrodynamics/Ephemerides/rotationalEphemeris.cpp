/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
  
#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{
namespace ephemerides
{

//! Function to calculate the rotational velocity vector of frame B w.r.t frame A.
Eigen::Vector3d getRotationalVelocityVectorInBaseFrameFromMatrices(
        const Eigen::Matrix3d& rotationToTargetFrame,
        const Eigen::Matrix3d& rotationMatrixToGlobalFrameDerivative )
{
    Eigen::Matrix3d crossProductMatrix =
            rotationMatrixToGlobalFrameDerivative * rotationToTargetFrame;
    return ( Eigen::Vector3d( ) << crossProductMatrix( 2, 1 ),
             crossProductMatrix( 0, 2 ), crossProductMatrix( 1, 0 ) ).finished( );

}

//! Function to calculate the time derivative of rotation matrix from frame A to frame B.
Eigen::Matrix3d getDerivativeOfRotationMatrixToFrame(
        const Eigen::Matrix3d& rotationToTargetFrame,
        const Eigen::Vector3d& rotationalVelocityVectorOfTargetFrameInBaseFrame )
{
    return linear_algebra::getCrossProductMatrix(
                -rotationToTargetFrame * rotationalVelocityVectorOfTargetFrameInBaseFrame ) *
            rotationToTargetFrame;
}

//! Get rotation quaternion from target frame to base frame.
template< >
Eigen::Quaterniond RotationalEphemeris::getRotationToBaseFrameTemplated< double >(
            const double timeSinceEpoch )
{
    return getRotationToBaseFrame( timeSinceEpoch );
}

//! Get rotation quaternion from target frame to base frame.
template< >
Eigen::Quaterniond RotationalEphemeris::getRotationToBaseFrameTemplated< Time >(
            const Time timeSinceEpoch )
{
    return getRotationToBaseFrameFromExtendedTime( timeSinceEpoch );
}

//! Get rotation quaternion to target frame from base frame.
template< >
Eigen::Quaterniond RotationalEphemeris::getRotationToTargetFrameTemplated< double >(
        const double secondsSinceEpoch )
{
    return getRotationToTargetFrame( secondsSinceEpoch );
}

//! Get rotation quaternion to target frame from base frame.
template< >
Eigen::Quaterniond RotationalEphemeris::getRotationToTargetFrameTemplated< Time >(
        const Time secondsSinceEpoch )
{
    return getRotationToTargetFrameFromExtendedTime( secondsSinceEpoch );
}

//! Function to calculate the derivative of the rotation matrix from target frame to original frame.
template< >
Eigen::Matrix3d RotationalEphemeris::getDerivativeOfRotationToBaseFrameTemplated< double >(
            const double timeSinceEpoch )
{
    return getDerivativeOfRotationToBaseFrame( timeSinceEpoch );
}

//! Function to calculate the derivative of the rotation matrix from target frame to original frame.
template< >
Eigen::Matrix3d RotationalEphemeris::getDerivativeOfRotationToBaseFrameTemplated< Time >(
            const Time timeSinceEpoch )
{
    return getDerivativeOfRotationToBaseFrameFromExtendedTime( timeSinceEpoch );
}


//! Function to calculate the derivative of the rotation matrix to target frame from original frame.
template< >
Eigen::Matrix3d RotationalEphemeris::getDerivativeOfRotationToTargetFrameTemplated< double >(
        const double secondsSinceEpoch )
{
    return getDerivativeOfRotationToTargetFrame( secondsSinceEpoch );
}

//! Function to calculate the derivative of the rotation matrix to target frame from original frame.
template< >
Eigen::Matrix3d RotationalEphemeris::getDerivativeOfRotationToTargetFrameTemplated< Time >(
        const Time secondsSinceEpoch )
{
    return getDerivativeOfRotationToTargetFrameFromExtendedTime( secondsSinceEpoch );
}

//! Function to calculate the full rotational state at given time
template< >
void RotationalEphemeris::getFullRotationalQuantitiesToTargetFrameTemplated< double >(
        Eigen::Quaterniond& currentRotationToLocalFrame,
        Eigen::Matrix3d& currentRotationToLocalFrameDerivative,
        Eigen::Vector3d& currentAngularVelocityVectorInGlobalFrame,
        const double timeSinceEpoch )
{
    getFullRotationalQuantitiesToTargetFrame(
                currentRotationToLocalFrame, currentRotationToLocalFrameDerivative, currentAngularVelocityVectorInGlobalFrame,
                timeSinceEpoch );
}

//! Function to calculate the full rotational state at given time
template< >
void  RotationalEphemeris::getFullRotationalQuantitiesToTargetFrameTemplated< Time >(
        Eigen::Quaterniond& currentRotationToLocalFrame,
        Eigen::Matrix3d& currentRotationToLocalFrameDerivative,
        Eigen::Vector3d& currentAngularVelocityVectorInGlobalFrame,
        const Time timeSinceEpoch )
{
    getFullRotationalQuantitiesToTargetFrameFromExtendedTime(
                currentRotationToLocalFrame, currentRotationToLocalFrameDerivative, currentAngularVelocityVectorInGlobalFrame,
                timeSinceEpoch );
}



} // namespace tudat

} // namespace ephemerides

