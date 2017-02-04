/*    Copyright (c) 2010-2017, Delft University of Technology
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

} // namespace tudat
} // namespace ephemerides

