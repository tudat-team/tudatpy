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
 *      lambertTargeterIzzo.h/.cpp source files, tudat revision 455/466.
 *      PyKEP kepler toolbox, Dario Izzo, ESA Advanced Concepts Team.
 *
 *    Notes
 *      This is a new implementation of the lambertTargeterIzzo class, for better adaptability and
 *      extension towards subclasses and future improvements/additions. Therefore, it replaces the
 *      lambertTargeterIzzo class while still providing the same functionality.
 *
 */

#include <cmath>

#include <boost/math/special_functions.hpp> // For asinh and acosh
#include <boost/exception/all.hpp> // For exceptions in sanity checks
#include <Eigen/Dense> // for cross product issues (can someone explain why, exactly?)

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/MissionSegments/zeroRevolutionLambertTargeterIzzo.h"
#include "Tudat/Mathematics/BasicMathematics/convergenceException.h"

namespace tudat
{
namespace mission_segments
{

//! Get radial velocity at departure.
double ZeroRevolutionLambertTargeterIzzo::getRadialVelocityAtDeparture( )
// Based on lambertTargeterIzzo class.
{
    // If execute has not been called yet, execute.
    if ( !solved ) execute( );

    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtDeparture = cartesianPositionAtDeparture.normalized( );

    // Compute radial velocity at departure.
    return cartesianVelocityAtDeparture.dot( radialUnitVectorAtDeparture );
}

//! Get transverse velocity at departure.
double ZeroRevolutionLambertTargeterIzzo::getTransverseVelocityAtDeparture( )
// Based on lambertTargeterIzzo class.
{
    // If execute has not been called yet, execute.
    if ( !solved ) execute( );

    // Compute angular momemtum vector.
    const Eigen::Vector3d angularMomentumVector =
            cartesianPositionAtDeparture.cross( cartesianVelocityAtDeparture );

    // Compute normalized angular momentum vector.
    const Eigen::Vector3d angularMomentumUnitVector = angularMomentumVector.normalized( );

    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtDeparture = cartesianPositionAtDeparture.normalized( );

    // Compute tangential unit vector.
    Eigen::Vector3d tangentialUnitVectorAtDeparture =
            angularMomentumUnitVector.cross( radialUnitVectorAtDeparture );

    // Compute tangential velocity at departure.
    return cartesianVelocityAtDeparture.dot( tangentialUnitVectorAtDeparture );
}

//! Get radial velocity at arrival.
double ZeroRevolutionLambertTargeterIzzo::getRadialVelocityAtArrival( )
// Based on lambertTargeterIzzo class.
{
    // If execute has not been called yet, execute.
    if ( !solved ) execute( );

    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtArrival = cartesianPositionAtArrival.normalized( );

    // Compute radial velocity at arrival.
    return cartesianVelocityAtArrival.dot( radialUnitVectorAtArrival );
}

//! Get transverse velocity at arrival.
double ZeroRevolutionLambertTargeterIzzo::getTransverseVelocityAtArrival( )
// Based on lambertTargeterIzzo class.
{
    // If execute has not been called yet, execute.
    if ( !solved ) execute( );

    // Compute angular momemtum vector.
    const Eigen::Vector3d angularMomentumVector =
            cartesianPositionAtArrival.cross( cartesianVelocityAtArrival );

    // Compute normalized angular momentum vector.
    const Eigen::Vector3d angularMomentumUnitVector = angularMomentumVector.normalized( );

    // Determine radial unit vector.
    const Eigen::Vector3d radialUnitVectorAtArrival = cartesianPositionAtArrival.normalized( );

    // Compute tangential unit vector.
    Eigen::Vector3d tangentialUnitVectorAtArrival
            = angularMomentumUnitVector.cross( radialUnitVectorAtArrival );

    // Compute tangential velocity at departure.
    return cartesianVelocityAtArrival.dot( tangentialUnitVectorAtArrival );
}

//! Get semi-major axis.
double ZeroRevolutionLambertTargeterIzzo::getSemiMajorAxis( )
// Based on lambertTargeterIzzo class.
{
    // If execute has not been called yet, execute.
    if ( !solved ) execute( );

    // Compute specific orbital energy: eps = v^2/ - mu/r.
    const double specificOrbitalEnergy = cartesianVelocityAtDeparture.squaredNorm( ) / 2.0
            - gravitationalParameter / cartesianPositionAtDeparture.norm( );

    // Compute semi-major axis: a = -mu / 2*eps.
    return -gravitationalParameter / ( 2.0 * specificOrbitalEnergy );
}

// What about get retrograde flag, tolerance and max number of iterations? Might be useful to check.

//! Execute the solving procedure.
void ZeroRevolutionLambertTargeterIzzo::execute( )
{
    // Sanity checks
    sanityCheckTimeOfFlight( );
    sanityCheckGravitationalParameter( );

    // Transform dimensional parameters to dimensionless parameters if not done already (e.g. in
    // multirevolution case
    if( !transformed ) transformDimensions( );

    // Solve root (single rev)
    double xResult = ZeroRevolutionLambertTargeterIzzo::computeRootTimeOfFlight( );

    // Reconstruct Vs
    computeVelocities( xResult );

    solved = true;
}

//! Sanity check time of flight.
void ZeroRevolutionLambertTargeterIzzo::sanityCheckTimeOfFlight( )
{
    // If time of flight is negative, throw an exception
    if ( timeOfFlight < 0 )
    {
        // Throw exception.
        throw std::runtime_error(
                    "Time-of-flight specified in Lambert problem must be strictly positive. Specified time-of-flight in days." +
                                         std::to_string( timeOfFlight ) );
    }
    // Else, do nothing and continue.
}

//! Sanity check gravitational parameter.
void ZeroRevolutionLambertTargeterIzzo::sanityCheckGravitationalParameter( )
{
    // If gravitational parameter is negative, throw an exception
    if ( gravitationalParameter < 0 )
    {
        // Throw exception.
        throw std::runtime_error(
                    "Gravitational parameter specified in Lambert problem must be strictly positive. Specified gravitational parameter: " +
                    std::to_string( gravitationalParameter ) );
    }
    // Else, do nothing and continue.
}

//! Transform input to sub results in adimensional units.
void ZeroRevolutionLambertTargeterIzzo::transformDimensions( )
// Created using theory from PyKEP toolbox and LambertTargeterIzzo class
{
    // Compute normalizing values.
    const double distanceNormalizingValue = cartesianPositionAtDeparture.norm( );
    velocityNormalizingValue = std::sqrt( gravitationalParameter /
                                          distanceNormalizingValue );
    const double timeNormalizingValue = distanceNormalizingValue / velocityNormalizingValue;

    // Compute transfer geometry parameters in adimensional units.
    // Time of Flight.
    normalizedTimeOfFlight = timeOfFlight / timeNormalizingValue;

    // Cosine of transfer angle.
    const double cosineOfTransferAngle =
            cartesianPositionAtDeparture.dot( cartesianPositionAtArrival )
            / (distanceNormalizingValue * cartesianPositionAtArrival.norm( ) );

    // Normalized Cartesian position at arrival.
    normalizedRadiusAtArrival = cartesianPositionAtArrival.norm( ) / distanceNormalizingValue;

    // Chord.
    normalizedChord = std::sqrt( 1.0 + normalizedRadiusAtArrival
                                 * ( normalizedRadiusAtArrival
                                     - 2.0 * cosineOfTransferAngle ) );

    // Semi-perimeter.
    normalizedSemiPerimeter = ( 1.0 + normalizedRadiusAtArrival + normalizedChord ) / 2.0;

    // Assuming a prograde motion, determine whether the transfer corresponds to the long- or the
    // short-way solution: longway if x1*y2 - x2*y1 < 0.
    isLongway = false;
    if ( cartesianPositionAtDeparture.x( ) * cartesianPositionAtArrival.y( )
         - cartesianPositionAtDeparture.y( ) * cartesianPositionAtArrival.x( ) < 0.0 )
    {
        isLongway = true;
    }

    // If retrograde is true, switch longway flag.
    if ( isRetrograde )
    {
        isLongway = !isLongway;
    }

    // Semi-major axis of the minimum energy ellipse.
    normalizedMinimumEnergySemiMajorAxis = normalizedSemiPerimeter / 2.0;

    // Transfer angle.
    transferAngle = std::acos( cosineOfTransferAngle );
    if ( isLongway )
    {
        transferAngle = 2.0 * mathematical_constants::PI - transferAngle;
    }

    // Lambda parameter.
    lambdaParameter = std::sqrt( normalizedRadiusAtArrival )
            * std::cos( transferAngle / 2.0 ) / normalizedSemiPerimeter;

    // Set transformed flag to true, as the dimension transformation has been performed
    transformed = true;
}

//! Compute time-of-flight using Lagrange's equation.
double ZeroRevolutionLambertTargeterIzzo::computeTimeOfFlight( const double xParameter )
// Created using theory from PyKEP toolbox and LambertTargeterIzzo class
{
    // Determine semi-major axis.
    const double semiMajorAxis = normalizedMinimumEnergySemiMajorAxis
            / ( 1.0 - xParameter * xParameter );

    // If x < 1, the solution is an ellipse.
    if ( xParameter < 1.0 )
    {
        // Alpha parameter in Lagrange's equation (no explanation available).
        const double alphaParameter = 2.0 * std::acos( xParameter );

        // Beta parameter in Lagrange's equation (no explanation available).
        double betaParameter;
        // If long transfer arc
        if ( isLongway )
        {
            betaParameter = -2.0 * std::asin(
                        std::sqrt( ( normalizedSemiPerimeter - normalizedChord )
                                   / ( 2.0 * semiMajorAxis ) ) );
        }
        // Otherwise short transfer arc
        else
        {
            betaParameter = 2.0 * std::asin(
                        std::sqrt( ( normalizedSemiPerimeter - normalizedChord )
                                   / ( 2.0 * semiMajorAxis ) ) );
        }

        // Time-of-flight according to Lagrange including multiple revolutions.
        const double timeOfFlight = semiMajorAxis * std::sqrt( semiMajorAxis ) *
                ( ( alphaParameter - std::sin( alphaParameter ) )
                  - ( betaParameter - std::sin( betaParameter ) ) );

        return timeOfFlight;
    }
    // Otherwise it is a hyperbola.
    else
    {
        // Alpha parameter in Lagrange's equation (no explanation available).
        const double alphaParameter = 2.0 * boost::math::acosh( xParameter );

        // Beta parameter in Lagrange's equation (no explanation available).
        double betaParameter;
        // If long transfer arc
        if ( isLongway )
        {
            betaParameter = -2.0 * boost::math::asinh ( std::sqrt( ( normalizedSemiPerimeter
                                                                     - normalizedChord )
                                                                   / ( -2.0 * semiMajorAxis ) ) );
        }
        // Otherwise short transfer arc
        else
        {
            betaParameter = 2.0 * boost::math::asinh ( std::sqrt( ( normalizedSemiPerimeter
                                                                    - normalizedChord )
                                                                  / ( -2.0 * semiMajorAxis ) ) );
        }

        // Time-of-flight according to Lagrange.
        const double timeOfFlightLagrange = -semiMajorAxis * std::sqrt( -semiMajorAxis ) *
                ( ( std::sinh( alphaParameter ) - alphaParameter )
                  - ( std::sinh( betaParameter ) - betaParameter ) );

        return timeOfFlightLagrange;
    }
}

//! Solve the time of flight equation for x.
double ZeroRevolutionLambertTargeterIzzo::computeRootTimeOfFlight( )
// Created using theory from PyKEP toolbox and LambertTargeterIzzo class
{
    // Find root (secant method, currently hard coded)
    // Optimize log(t_spec).
    const double logarithmOfTheSpecifiedTimeOfFlight = std::log( normalizedTimeOfFlight );

    // Define initial guesses for abcissae (x) and ordinates (y).
    double x1 = std::log( 0.5 ), x2 = std::log( 1.5 );

    double y1 = std::log( computeTimeOfFlight( -0.5 ) ) - logarithmOfTheSpecifiedTimeOfFlight;

    double y2 = std::log( computeTimeOfFlight( 0.5 ) ) - logarithmOfTheSpecifiedTimeOfFlight;

    // Declare and initialize root-finding parameters.
    double rootFindingError = 1.0, xNew = 0.0, yNew = 0.0;
    int iterator = 0;

    // Root-finding loop.
    while ( ( rootFindingError > convergenceTolerance ) && (y1 != y2)
            && ( iterator < maximumNumberOfIterations ) )
    {
        // Update iterator.
        iterator++;

        // Compute new x-value.
        xNew = ( x1 * y2 - y1 * x2 ) / ( y2 - y1 );

        // Compute corresponding y-value.
        yNew = std::log( computeTimeOfFlight( std::exp( xNew ) - 1.0 ) )
                - logarithmOfTheSpecifiedTimeOfFlight;

        // Update abcissae and ordinates.
        x1 = x2;
        y1 = y2;
        x2 = xNew;
        y2 = yNew;

        // Compute root-finding error.
        rootFindingError = std::fabs( x1 - xNew );
    }

    // Verify that root-finder has converged.
    if ( iterator == maximumNumberOfIterations )
    {
        throw std::runtime_error(
                    "Multi-Revolution Lambert targeter failed to converge to a solution. Reached the maximum number of iterations: " +
                    std::to_string( maximumNumberOfIterations ) );
    }

    // Recovering x parameter and returning it.
    double xParameter = std::exp( xNew ) - 1.0;
    return xParameter;
}

//! Compute velocities at departure and arrival.
void ZeroRevolutionLambertTargeterIzzo::computeVelocities( const double xParameter )
// Created using theory from PyKEP toolbox and LambertTargeterIzzo class
{
    // Then it is possible to retrieve a sensible value from the x-parameter computed)
    // Determine semi-major axis of the conic.
    const double semiMajorAxis = normalizedMinimumEnergySemiMajorAxis
            / ( 1.0 - xParameter * xParameter );

    // Declare variables.
    double etaParameter, etaParameterSquared, psiParameter;

    // If x < 1, the solution is an ellipse.
    if ( xParameter < 1.0 )
    {
        // Alpha parameter in Lagrange's equation (no explanation available).
        const double alphaParameter = 2.0 * std::acos( xParameter );

        // Beta parameter in Lagrange's equation (no explanation available).
        double betaParameter = 2.0 * std::asin( std::sqrt( ( normalizedSemiPerimeter
                                                             - normalizedChord )
                                                           / ( 2.0 * semiMajorAxis ) ) );

        if ( isLongway )
        {
            betaParameter = -betaParameter;
        }

        // Psi parameter in Izzo's approach (no explanation available).
        psiParameter = ( alphaParameter - betaParameter ) / 2.0;

        // Eta parameter in Izzo's approach (no explanation available).
        etaParameterSquared = 2.0 * semiMajorAxis * std::sin( psiParameter )
                * std::sin( psiParameter ) / normalizedSemiPerimeter;
        etaParameter = std::sqrt( etaParameterSquared );
    }

    // Otherwise it is a hyperbola.
    else
    {
        // Alpha parameter in Lagrange's equation (no explanation available).
        const double alphaParameter = 2.0 * boost::math::acosh( xParameter );

        // Beta parameter in Lagrange's equation (no explanation available).
        double betaParameter = 2.0 * boost::math::asinh (
                    std::sqrt( ( normalizedSemiPerimeter - normalizedChord )
                               / ( -2.0 * semiMajorAxis ) ) );

        if ( isLongway )
        {
            betaParameter = -betaParameter;
        }

        // Psi parameter in Izzo's approach (no explanation available).
        psiParameter = (alphaParameter - betaParameter ) / 2.0;

        // Eta parameter in Izzo's approach (no explanation available).
        etaParameterSquared = -2.0 * semiMajorAxis * std::sinh( psiParameter )
                * std::sinh( psiParameter ) / normalizedSemiPerimeter;
        etaParameter = std::sqrt( etaParameterSquared );
    }

    // Determine semi-latus rectum, p.
    const double semiLatusRectum = ( normalizedRadiusAtArrival
                                     / ( normalizedMinimumEnergySemiMajorAxis
                                         * etaParameterSquared ) )
            * std::sin( transferAngle / 2.0 )
            * std::sin( transferAngle / 2.0 );

    // Velocity components at departure.
    const double radialVelocityAtDeparture =
            ( 1.0 / ( etaParameter * std::sqrt( normalizedMinimumEnergySemiMajorAxis ) ) )
            * ( 2.0 * lambdaParameter * normalizedMinimumEnergySemiMajorAxis
                - ( lambdaParameter + xParameter * etaParameter ) );

    const double transverseVelocityAtDeparture = std::sqrt( semiLatusRectum );

    // Velocity components at arrival.
    const double transverseVelocityAtArrival = transverseVelocityAtDeparture
            / normalizedRadiusAtArrival;
    const double radialVelocityAtArrival = ( transverseVelocityAtDeparture
                                             - transverseVelocityAtArrival )
            / std::tan( transferAngle / 2.0 )
            - radialVelocityAtDeparture;

    // Determining inertial vectors
    // Determine radial unit vectors.
    const Eigen::Vector3d radialUnitVectorAtDeparture = cartesianPositionAtDeparture.normalized( );
    const Eigen::Vector3d radialUnitVectorAtArrival = cartesianPositionAtArrival.normalized( );

    // Determine plane of motion.
    Eigen::Vector3d angularMomentumVector;

    if ( isLongway )
    {
        angularMomentumVector = radialUnitVectorAtArrival.cross( radialUnitVectorAtDeparture );
    }
    else
    {
        angularMomentumVector = radialUnitVectorAtDeparture.cross( radialUnitVectorAtArrival );
    }

    // Compute normalized angular momentum vector.
    const Eigen::Vector3d angularMomentumUnitVector = angularMomentumVector.normalized( );

    // Compute transverse unit vectors.
    const Eigen::Vector3d transverseUnitVectorAtDeparture =
            radialUnitVectorAtDeparture.cross( angularMomentumUnitVector );
    const Eigen::Vector3d transverseUnitVectorAtArrival =
            radialUnitVectorAtArrival.cross( angularMomentumUnitVector );

    // Reconstruct non-dimensional velocity vectors.
    cartesianVelocityAtDeparture << radialVelocityAtDeparture * radialUnitVectorAtDeparture
                                    - transverseVelocityAtDeparture * transverseUnitVectorAtDeparture;

    cartesianVelocityAtArrival << radialVelocityAtArrival * radialUnitVectorAtArrival
                                  - transverseVelocityAtArrival * transverseUnitVectorAtArrival;

    // Return to dimensions of initial problem definition.
    cartesianVelocityAtDeparture *= velocityNormalizingValue;
    cartesianVelocityAtArrival *= velocityNormalizingValue;
}

} // namespace mission_segments
} // namespace tudat
