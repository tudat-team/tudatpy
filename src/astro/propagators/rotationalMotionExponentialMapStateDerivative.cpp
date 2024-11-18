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

#include "tudat/astro/propagators/rotationalMotionExponentialMapStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to obtain the time derivative of exponantial map of body-fixed to inertial frame.
Eigen::Vector4d calculateExponentialMapDerivative( const Eigen::Vector4d& currentExponentialMapToBaseFrame,
                                                   const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame )
{
    // Define the tolerance of a singularity
    double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Declare eventual output vector
    Eigen::Vector4d exponentialMapDerivative = Eigen::Vector4d::Zero( );

    // Get intermediate variables
    Eigen::Vector3d exponentialMapVector = currentExponentialMapToBaseFrame.segment( 0, 3 );
    double exponentialMapMagnitude = exponentialMapVector.norm( );

    // Compute kinematic equation, i.e., derivative of exponential map (also valid for SEM)
    if ( exponentialMapMagnitude < singularityTolerance )
    {
        double exponentialMapMagnitudeSquared = std::pow( exponentialMapMagnitude, 2 );
        exponentialMapDerivative.segment( 0, 3 ) = 0.5 * (
                    ( ( 12.0 - exponentialMapMagnitudeSquared ) / 6.0 ) *
                    angularVelocityVectorInBodyFixedFrame - angularVelocityVectorInBodyFixedFrame.cross(
                        exponentialMapVector ) - angularVelocityVectorInBodyFixedFrame.dot(
                        exponentialMapVector ) *
                    ( ( 60.0 + exponentialMapMagnitudeSquared ) / 360.0 ) * exponentialMapVector );
    }
    else
    {
        double cotangentHalfExponentialMapMagnitude = std::cos( 0.5 * exponentialMapMagnitude ) /
                std::sin( 0.5 * exponentialMapMagnitude );
        Eigen::Vector3d exponentialMapCrossRotationalVelocityVector = exponentialMapVector.cross(
                    angularVelocityVectorInBodyFixedFrame );
        exponentialMapDerivative.segment( 0, 3 ) =
                angularVelocityVectorInBodyFixedFrame + 0.5 * exponentialMapCrossRotationalVelocityVector +
                ( 1 - 0.5 * exponentialMapMagnitude * cotangentHalfExponentialMapMagnitude ) /
                std::pow( exponentialMapMagnitude, 2 ) *
                exponentialMapVector.cross( exponentialMapCrossRotationalVelocityVector );
    }

    // Give output
    return exponentialMapDerivative;
}

//template class RotationalMotionExponentialMapStateDerivative< double, double >;

} // namespace propagators

} // namespace tudat
