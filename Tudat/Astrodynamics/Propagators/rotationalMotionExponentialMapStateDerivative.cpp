/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/Astrodynamics/Propagators/rotationalMotionExponentialMapStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to obtain the time derivative of exponantial map of body-fixed to inertial frame
Eigen::Vector3d calculateExponentialMapDerivative( const Eigen::Vector3d& currentExponentialMapToBaseFrame,
                                                   const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame )
{
    // Define the tolerance of a singularity
    double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Declare eventual output vector
    Eigen::Vector3d exponentialMapDerivative;

    // Get intermediate variables
    double exponentialMapMagnitude = currentExponentialMapToBaseFrame.norm( );

    // Compute kinematic equation, i.e., derivative of exponential map (also valid for SEM)
    if ( exponentialMapMagnitude < singularityTolerance )
    {
        double exponentialMapMagnitudeSquared = std::pow( exponentialMapMagnitude, 2 );
        exponentialMapDerivative = 0.5 * ( ( ( 12.0 - exponentialMapMagnitudeSquared ) / 6.0 ) *
                                           angularVelocityVectorInBodyFixedFrame - angularVelocityVectorInBodyFixedFrame.cross(
                                               currentExponentialMapToBaseFrame ) - angularVelocityVectorInBodyFixedFrame.dot(
                                               currentExponentialMapToBaseFrame ) *
                                           ( ( 60.0 + exponentialMapMagnitudeSquared ) / 360.0 ) * currentExponentialMapToBaseFrame );
    }
    else
    {
        double cotangentHalfExponentialMapMagnitude = std::cos( 0.5 * exponentialMapMagnitude ) /
                std::sin( 0.5 * exponentialMapMagnitude );
        Eigen::Vector3d exponentialMapCrossRotationalVelocityVector = currentExponentialMapToBaseFrame.cross(
                    angularVelocityVectorInBodyFixedFrame );
        exponentialMapDerivative = angularVelocityVectorInBodyFixedFrame + 0.5 * exponentialMapCrossRotationalVelocityVector +
                ( 1 - 0.5 * exponentialMapMagnitude * cotangentHalfExponentialMapMagnitude ) /
                std::pow( exponentialMapMagnitude, 2 ) *
                currentExponentialMapToBaseFrame.cross( exponentialMapCrossRotationalVelocityVector );
    }

    // Give output
    return exponentialMapDerivative;
}

} // namespace propagators

} // namespace tudat
