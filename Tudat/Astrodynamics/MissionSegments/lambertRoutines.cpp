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
 *      Battin, R.H. An Introduction to the Mathematics and Methods of Astrodynamics,
 *          AIAA Education Series, 1999.
 *      Izzo, D. lambert_problem.h, keptoolbox.
 *      Gooding, R.H. A procedure for the solution of Lambert's orbital boundary-value problem,
 *          Celestial Mechanics and Dynamical Astronomy, 48:145-165, 1990.
 *
 */

#include <cmath>
#include <stdexcept>

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/special_functions.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"
#include "Tudat/Mathematics/BasicMathematics/functionProxy.h"

namespace tudat
{
namespace mission_segments
{

using namespace root_finders;

//! Solve Lambert Problem using Izzo's algorithm.
void solveLambertProblemIzzo( const Eigen::Vector3d& cartesianPositionAtDeparture,
                              const Eigen::Vector3d& cartesianPositionAtArrival,
                              const double timeOfFlight,
                              const double gravitationalParameter,
                              Eigen::Vector3d& cartesianVelocityAtDeparture,
                              Eigen::Vector3d& cartesianVelocityAtArrival,
                              const bool isRetrograde,
                              const double convergenceTolerance,
                              const unsigned int maximumNumberOfIterations )
{
    // Sanity check for specified time-of-flight.
    if ( timeOfFlight <= 0.0 )
    {
        // Throw exception.
        throw std::runtime_error( "Specified time-of-flight must be strictly positive: " + std::to_string( timeOfFlight ) + " days." );
    }

    // Compute normalizing values.
    const double distanceNormalizingValue = cartesianPositionAtDeparture.norm( );
    const double velocityNormalizingValue = std::sqrt( gravitationalParameter /
                                                       distanceNormalizingValue );
    const double timeNormalizingValue = distanceNormalizingValue / velocityNormalizingValue;

    // Compute transfer geometry parameters in adimensional units.
    // Cosine of transfer angle.
    const double cosineOfTransferAngle =
            cartesianPositionAtDeparture.dot( cartesianPositionAtArrival )
            / (distanceNormalizingValue * cartesianPositionAtArrival.norm( ) );

    // Normalized Cartesian position at arrival.
    const double normalizedRadiusAtArrival = cartesianPositionAtArrival.norm( )
            / distanceNormalizingValue;

    // Chord.
    const double chord = std::sqrt( 1.0 + normalizedRadiusAtArrival
                                    * ( normalizedRadiusAtArrival - 2.0 * cosineOfTransferAngle ) );

    // Semi-perimeter.
    const double semiPerimeter = ( 1.0 + normalizedRadiusAtArrival + chord ) / 2.0;

    // Assuming a prograde motion, determine whether the transfer corresponds to the long- or the
    // short-way solution: longway if x1*y2 - x2*y1 < 0.
    bool isLongway = false;
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
    const double semiMajorAxisOfTheMinimumEnergyEllipse = semiPerimeter / 2.0;

    // Transfer angle.
    double transferAngle = std::acos( cosineOfTransferAngle );
    if ( isLongway )
    {
        transferAngle = 2.0 * mathematical_constants::PI - transferAngle;
    }

    // Lambda parameter.
    const double lambdaParameter = std::sqrt( normalizedRadiusAtArrival )
            * std::cos( transferAngle / 2.0 ) / semiPerimeter;

    // Optimize log(t_spec).
    const double normalizedSpecifiedTimeOfFlight = timeOfFlight / timeNormalizingValue;
    const double logarithmOfTheSpecifiedTimeOfFlight = std::log( normalizedSpecifiedTimeOfFlight );

    // Secant Method.
    // Define initial guesses for abcissae (x) and ordinates (y).
    double x1 = std::log( 0.5 ), x2 = std::log( 1.5 );

    double y1 = std::log( computeTimeOfFlightIzzo( -0.5, semiPerimeter, chord, isLongway,
                                                   semiMajorAxisOfTheMinimumEnergyEllipse ) )
            - logarithmOfTheSpecifiedTimeOfFlight;

    double y2 = std::log( computeTimeOfFlightIzzo( 0.5, semiPerimeter, chord, isLongway,
                                                   semiMajorAxisOfTheMinimumEnergyEllipse ) )
            - logarithmOfTheSpecifiedTimeOfFlight;

    // Declare and initialize root-finding parameters.
    double rootFindingError = 1.0, xNew = 0.0, yNew = 0.0;
    unsigned int iterator = 0;


    // Root-finding loop.
    while ( ( rootFindingError > convergenceTolerance ) && (y1 != y2)
            && ( iterator < maximumNumberOfIterations ) )
    {
        // Update iterator.
        iterator++;

        // Compute new x-value.
        xNew = ( x1 * y2 - y1 * x2 ) / ( y2 - y1 );

        // Compute corresponding y-value.
        yNew = std::log( computeTimeOfFlightIzzo( std::exp( xNew ) - 1.0, semiPerimeter, chord,
                                                  isLongway,
                                                  semiMajorAxisOfTheMinimumEnergyEllipse ) )
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
        std::string errorMessage = "Lambert Solver did not converge within the maximum number of iterations: " +
                std::to_string( maximumNumberOfIterations );
        throw std::runtime_error( errorMessage );
    }

    // Revert to x parameter.
    const double xParameter = std::exp( xNew ) - 1.0;

    // Determine semi-major axis of the conic.
    const double semiMajorAxis = semiMajorAxisOfTheMinimumEnergyEllipse
            / ( 1.0 - xParameter * xParameter );

    // Declare variables.
    double etaParameter, etaParameterSquared, psiParameter;

    // If x < 1, the solution is an ellipse.
    if ( xParameter < 1.0 )
    {
        // Alpha parameter.
        const double alphaParameter = 2.0 * std::acos( xParameter );

        // Beta parameter.
        double betaParameter = 2.0 * std::asin( std::sqrt( ( semiPerimeter - chord )
                                                           / ( 2.0 * semiMajorAxis ) ) );

        if ( isLongway )
        {
            betaParameter = -betaParameter;
        }

        // Psi parameter.
        psiParameter = ( alphaParameter - betaParameter ) / 2.0;

        // Eta parameter.
        etaParameterSquared = 2.0 * semiMajorAxis * std::sin( psiParameter )
                * std::sin( psiParameter ) / semiPerimeter;
        etaParameter = std::sqrt( etaParameterSquared );
    }

    // Otherwise it is a hyperbola.
    else
    {
        // Alpha parameter.
        const double alphaParameter = 2.0 * boost::math::acosh( xParameter );

        // Beta parameter.
        double betaParameter = 2.0 * boost::math::asinh (
                    std::sqrt( ( semiPerimeter - chord ) / ( -2.0 * semiMajorAxis ) ) );

        if ( isLongway )
        {
            betaParameter = -betaParameter;
        }

        // Psi parameter.
        psiParameter = (alphaParameter - betaParameter ) / 2.0;

        // Eta parameter.
        etaParameterSquared = -2.0 * semiMajorAxis * std::sinh( psiParameter )
                * std::sinh( psiParameter ) / semiPerimeter;
        etaParameter = std::sqrt( etaParameterSquared );
    }

    // Determine semi-latus rectum, p.
    const double semiLatusRectum = ( normalizedRadiusAtArrival
                                     / ( semiMajorAxisOfTheMinimumEnergyEllipse
                                         * etaParameterSquared ) )
            * std::sin( transferAngle / 2.0 )
            * std::sin( transferAngle / 2.0 );

    // Determine sigma.
    const double sigmaParameter =
            ( 1.0 / ( etaParameter * std::sqrt( semiMajorAxisOfTheMinimumEnergyEllipse ) ) )
            * ( 2.0 * lambdaParameter * semiMajorAxisOfTheMinimumEnergyEllipse
                - ( lambdaParameter + xParameter * etaParameter ) );

    // Velocity components at departure.
    const double radialVelocityAtDeparture = sigmaParameter;
    const double transverseVelocityAtDeparture = std::sqrt( semiLatusRectum );

    // Velocity components at arrival.
    const double transverseVelocityAtArrival = transverseVelocityAtDeparture
            / normalizedRadiusAtArrival;
    const double radialVelocityAtArrival = ( transverseVelocityAtDeparture
                                             - transverseVelocityAtArrival )
                                           / std::tan( transferAngle / 2.0 )
                                           - radialVelocityAtDeparture;

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
    cartesianVelocityAtDeparture
            << radialVelocityAtDeparture * radialUnitVectorAtDeparture.x( )
               - transverseVelocityAtDeparture * transverseUnitVectorAtDeparture.x( ),
            radialVelocityAtDeparture * radialUnitVectorAtDeparture.y( )
            - transverseVelocityAtDeparture * transverseUnitVectorAtDeparture.y( ),
            radialVelocityAtDeparture * radialUnitVectorAtDeparture.z( )
            - transverseVelocityAtDeparture * transverseUnitVectorAtDeparture.z( );

    cartesianVelocityAtArrival
            << radialVelocityAtArrival * radialUnitVectorAtArrival.x( )
               - transverseVelocityAtArrival * transverseUnitVectorAtArrival.x( ),
            radialVelocityAtArrival * radialUnitVectorAtArrival.y( )
            - transverseVelocityAtArrival * transverseUnitVectorAtArrival.y( ),
            radialVelocityAtArrival * radialUnitVectorAtArrival.z( )
            - transverseVelocityAtArrival * transverseUnitVectorAtArrival.z( );

    // Return dimensions.
    cartesianVelocityAtDeparture *= velocityNormalizingValue;
    cartesianVelocityAtArrival *= velocityNormalizingValue;

}

//! Compute time-of-flight using Lagrange's equation.
double computeTimeOfFlightIzzo( const double xParameter, const double semiPerimeter,
                                const double chord, const bool isLongway,
                                const double semiMajorAxisOfTheMinimumEnergyEllipse )
{
    // Determine semi-major axis.
    const double semiMajorAxis = semiMajorAxisOfTheMinimumEnergyEllipse
            / ( 1.0 - xParameter * xParameter );

    // If x < 1, the solution is an ellipse.
    if ( xParameter < 1.0 )
    {
        // Alpha parameter.
        const double alphaParameter = 2.0 * std::acos( xParameter );

        // Beta parameter.
        double betaParameter = 2.0 * std::asin ( std::sqrt( ( semiPerimeter - chord )
                                                            / ( 2.0 * semiMajorAxis ) ) );

        if ( isLongway )
        {
            betaParameter = -betaParameter;
        }

        // Time-of-flight according to Lagrange.
        const double timeOfFlight = semiMajorAxis * std::sqrt( semiMajorAxis ) *
                ( ( alphaParameter - std::sin( alphaParameter ) )
                  - ( betaParameter - std::sin( betaParameter ) ) );

        return timeOfFlight;
    }
    // Otherwise it is a hyperbola.
    else
    {
        // Alpha parameter.
        const double alphaParameter = 2.0 * boost::math::acosh( xParameter );

        // Beta parameter.
        double betaParameter = 2.0 * boost::math::asinh(
                    std::sqrt( ( semiPerimeter - chord ) / ( -2.0 * semiMajorAxis ) ) );

        if ( isLongway )
        {
            betaParameter = -betaParameter;
        }

        // Time-of-flight according to Lagrange.
        const double timeOfFlight = -semiMajorAxis * std::sqrt( -semiMajorAxis ) *
                ( ( std::sinh( alphaParameter ) - alphaParameter )
                  - ( std::sinh( betaParameter ) - betaParameter ) );

        return timeOfFlight;
    }


}

//! Solve Lambert Problem using Gooding's algorithm.
void solveLambertProblemGooding( const Eigen::Vector3d& cartesianPositionAtDeparture,
                                 const Eigen::Vector3d& cartesianPositionAtArrival,
                                 const double timeOfFlight,
                                 const double gravitationalParameter,
                                 Eigen::Vector3d& cartesianVelocityAtDeparture,
                                 Eigen::Vector3d& cartesianVelocityAtArrival,
                                 RootFinderPointer rootFinder )
{
    if ( !rootFinder.get( ) )
    {
        rootFinder = boost::make_shared< root_finders::NewtonRaphson >( 1.0e-12, 1000 );
    }

    // Normalize positions.
    const double radiusAtDeparture = cartesianPositionAtDeparture.norm( );
    const double radiusAtArrival = cartesianPositionAtArrival.norm( );

    // Compute angle between positions.
    Eigen::Vector3d planeNormalPosition =
            cartesianPositionAtDeparture.cross( cartesianPositionAtArrival ).normalized( );

    double reducedLambertAngle = linear_algebra::computeAngleBetweenVectors(
                cartesianPositionAtDeparture, cartesianPositionAtArrival );

    if ( planeNormalPosition.z( ) < 0.0 )
    {
        reducedLambertAngle = 2.0 * mathematical_constants::PI
                - reducedLambertAngle;
    }

    // Compute chord.
    const double chord = std::sqrt( radiusAtDeparture * radiusAtDeparture +
                                    radiusAtArrival * radiusAtArrival -
                                    2.0 * radiusAtDeparture * radiusAtArrival *
                                    std::cos( reducedLambertAngle ) );

    // Compute semi-perimeter.
    const double semiPerimeter = ( radiusAtDeparture + radiusAtArrival + chord ) / 2.0;

    // Compute normalized time of flight.
    // Formula (7) [3].
    const double normalizedTimeOfFlight = std::sqrt( 8.0 * gravitationalParameter / semiPerimeter
                                                     / semiPerimeter / semiPerimeter )
                                          * timeOfFlight;

    // Compute q-parameter.
    // Formula (5) [3].
    const double qParameter = ( std::sqrt( radiusAtDeparture * radiusAtArrival ) /
                                semiPerimeter ) * cos( reducedLambertAngle / 2.0 );

    // Compute values of parameters needed to compute the initial guess of x-parameter.
    const double zParameterInitial = std::sqrt( 1.0 - qParameter * qParameter );

    // Formula (6.11) [2].
    const double lambertEccentricAnomalyInitial = -1.0;

    // Formula (6.12) [2].
    const double yParameterInitial = std::sqrt( std::fabs( lambertEccentricAnomalyInitial ) );

    // Formula (6.14) [2].
    const double fParameterInitial = yParameterInitial * zParameterInitial;

    // Formula (6.15) [2].
    const double gParameterInitial = -1.0 * qParameter * lambertEccentricAnomalyInitial;

    // Page 47 [2].
    const double dParameterInitial = atan( fParameterInitial / gParameterInitial );

    // Value of T(x) for x=0.
    // Formula (6.9) [2].
    const double tFunctionInitial = -2.0 * ( qParameter * zParameterInitial +
            dParameterInitial / yParameterInitial ) / lambertEccentricAnomalyInitial;

    // Declare initial guess of x-parameter.
    double initialLambertGuess;

    // Determine initial Lambert guess.
    if ( tFunctionInitial > normalizedTimeOfFlight )
    {
        // Formula (11) [3].
        initialLambertGuess = tFunctionInitial *
                ( tFunctionInitial - normalizedTimeOfFlight ) /
                ( 4.0 * normalizedTimeOfFlight );
    }
    else
    {
        // Formula (13) [3].
        const double x01 = -1.0 * ( normalizedTimeOfFlight - tFunctionInitial ) /
                     ( normalizedTimeOfFlight - tFunctionInitial + 4.0 );

        // Formula (15) [3].
        const double x02 = -1.0 * std::sqrt( ( normalizedTimeOfFlight -
                     tFunctionInitial ) / ( normalizedTimeOfFlight +
                     0.5 * tFunctionInitial ) );

        // Formula (20) [3].
        const double phi = 2.0 * atan2( ( 1.0 - qParameter * qParameter ), 2.0 * qParameter );

        // Formula (16) [3].
        const double W = x01 + 1.7 * std::sqrt(
                    2.0 - phi / mathematical_constants::PI );

        // Formula (17) [3].
        double x03;
        if ( W >= 0.0 )
        {
           x03 = x01;
        }
        else
        {
           x03 = x01 + std::pow( -W , 1.0 / 16.0 ) * ( x02 - x01 );
        }

        // Formula (19) [3].
        const double lambdax = 1.0 + 0.5 * x03 * ( 1.0 + x01 ) - 0.03 * x03 * x03
                * std::sqrt( 1.0 + x01 );

        // Formula (18) [3].
        initialLambertGuess = lambdax * x03;
    }

    // Newton-Raphson method implementation.
    // Set the class that contains the functions needed for Newton-Raphson.
    LambertFunctionsGooding lambertFunctionsGooding( qParameter, normalizedTimeOfFlight );

    // Create an object containing the function of which we whish to obtain the root from.
    using basic_mathematics::UnivariateProxyPointer;
    using basic_mathematics::UnivariateProxy;
    UnivariateProxyPointer rootFunction = boost::make_shared< UnivariateProxy >(
                boost::bind( &LambertFunctionsGooding::computeLambertFunctionGooding,
                             lambertFunctionsGooding, _1 ) );

    // Add the first derivative of the root function.
    rootFunction->addBinding( -1, boost::bind( &LambertFunctionsGooding::
                                               computeFirstDerivativeLambertFunctionGooding,
                                               lambertFunctionsGooding, _1 ) );

    // Initialize the xParameter.
    double xParameter = TUDAT_NAN;

    // Set initial guess of the variable computed in Newton-Rapshon method. A patch for negative
    // initialLambertGuess is applied. This has not been stated in the paper by Gooding, but this
    // patch has been found by trial and error.
    if ( initialLambertGuess * initialLambertGuess - 1.0 < 0.0 )
    {
        // Set xParameter based on result of Newton-Raphson root-finding algorithm.
        xParameter = rootFinder->execute( rootFunction, std::fabs( initialLambertGuess ) );
    }
    else
    {
        // Set xParameter based on result of Newton-Raphson root-finding algorithm.
        xParameter = rootFinder->execute( rootFunction, initialLambertGuess );
    }

    // Compute velocities at departure and at arrival.

    // Compute gamma, rho and sigma parameters, needed to compute the velocities.
    // Formula (8) [3].
    const double zParameter = std::sqrt( 1.0 - qParameter * qParameter +
            qParameter * qParameter * xParameter * xParameter );

    // Formula (6.23) [2].
    const double lambertGamma = std::sqrt( gravitationalParameter * semiPerimeter / 2.0 );

    // Formula (6.24) [2].
    const double lambertRho = ( radiusAtDeparture - radiusAtArrival ) / chord;

    // Formula (6.25) [2].
    const double lambertSigma = 2.0 * ( std::sqrt(
                                            radiusAtDeparture * radiusAtArrival / chord / chord ) )
            * std::sin( reducedLambertAngle / 2.0 );

    // Compute radial speeds at departure and at arrival.
    // Formula (23) [3].
    const double radialSpeedAtDeparture =
            lambertGamma * ( qParameter * zParameter - xParameter -
                             lambertRho * ( qParameter * zParameter + xParameter ) ) /
            radiusAtDeparture;

    // Formula (24) [3].
    const double radialSpeedAtArrival =
            -lambertGamma * ( qParameter * zParameter - xParameter +
                              lambertRho * ( qParameter * zParameter + xParameter ) ) /
            radiusAtArrival;

    // Compute transverse speeds at departure and at arrival.
    // Compute large part of formulas (25) and (26).
    const double transverseSpeedHelper = lambertGamma * lambertSigma *
            ( zParameter + qParameter * xParameter );

    // Formula (25) [3].
    const double transverseSpeedAtDeparture = transverseSpeedHelper / radiusAtDeparture;

    // Formula (26) [3].
    const double transverseSpeedAtArrival = transverseSpeedHelper / radiusAtArrival;

    // Compute inertial velocities (those velocity are computed with respect to the central body,
    // so they can be either heliocentric or planetocentric).

    // Compute radial unit vectors at departure and at arrival.
    const Eigen::Vector3d radialUnitVectorAtDeparture = cartesianPositionAtDeparture.normalized( );
    const Eigen::Vector3d radialUnitVectorAtArrival = cartesianPositionAtArrival.normalized( );

    Eigen::Vector3d unitZVector;
    unitZVector.x( ) = 0.0;
    unitZVector.y( ) = 0.0;
    unitZVector.z( ) = 1.0;

    Eigen::Vector3d planeNormal =
            ( ( cartesianPositionAtDeparture.cross( unitZVector ) ).cross(
                    cartesianPositionAtDeparture ) ).normalized( );

    // Compute unit vector that is normal to the plane in which the trajectory takes place, and
    // points in the positive z-direction.
    if ( planeNormal.z( ) < 0.0 )
    {
        planeNormal = -planeNormal;
    }

    // Compute transverse unit vectors at departure and at arrival.
    const Eigen::Vector3d transverseUnitVectorAtDeparture = planeNormal.cross(
            cartesianPositionAtDeparture.normalized( ) );
    const Eigen::Vector3d transverseUnitVectorAtArrival = planeNormal.cross(
            cartesianPositionAtArrival.normalized( ) );

    // Compute radial inertial velocities at departure and at arrival.
    const Eigen::Vector3d radialInertialVelocityAtDeparture
            = radialSpeedAtDeparture * radialUnitVectorAtDeparture;
    const Eigen::Vector3d radialInertialVelocityAtArrival
            = radialSpeedAtArrival * radialUnitVectorAtArrival;

    // Compute transverse heliocentric velocities at departure and
    // at arrival.
    const Eigen::Vector3d transverseInertialVelocityAtDeparture
            = transverseSpeedAtDeparture * transverseUnitVectorAtDeparture;
    const Eigen::Vector3d transverseInertialVelocityAtArrival
            = transverseSpeedAtArrival * transverseUnitVectorAtArrival;

    // Compute heliocentric velocities at departure and at arrival.
    // Define output velocities.
    cartesianVelocityAtDeparture
            = radialInertialVelocityAtDeparture + transverseInertialVelocityAtDeparture;
    cartesianVelocityAtArrival
            = radialInertialVelocityAtArrival + transverseInertialVelocityAtArrival;
}

//! Define general Lambert function.
double LambertFunctionsGooding::computeLambertFunctionGooding( const double xParameter )
{
    const double lambertEccentricAnomaly = xParameter * xParameter - 1.0;

    if ( lambertEccentricAnomaly > 0.0 )
    {
        return lambertFunctionPositiveGooding( xParameter );
    }

    else
    {
        return lambertFunctionNegativeGooding( xParameter );
    }
}

//! Define first derivative of general Lambert function.
double LambertFunctionsGooding::computeFirstDerivativeLambertFunctionGooding( const double xParameter )
{
    const double lambertEccentricAnomaly = xParameter * xParameter - 1.0;

    if ( lambertEccentricAnomaly > 0.0 )
    {
        return lambertFirstDerivativeFunctionPositiveGooding( xParameter );
    }

    else
    {
        return lambertFirstDerivativeFunctionNegativeGooding( xParameter );
    }
}

//! Define Lambert function for positive lambertEccentricAnomaly.
double LambertFunctionsGooding::lambertFunctionPositiveGooding( const double xParameter )
{
    const double lambertEccentricAnomaly = xParameter * xParameter - 1.0;
    const double yParameter = std::sqrt( std::fabs( lambertEccentricAnomaly ) );
    const double zParameter = std::sqrt( 1.0 - qParameter * qParameter+
                                         qParameter * qParameter * xParameter  * xParameter );
    const double fParameter = yParameter * ( zParameter - qParameter * xParameter );
    const double gParameter = xParameter * zParameter - qParameter * lambertEccentricAnomaly;
    const double dParameter = log( fParameter + gParameter );

    return normalizedTimeOfFlight - 2.0 * ( xParameter - qParameter * zParameter - dParameter /
                                            yParameter ) / lambertEccentricAnomaly;
}

//! Define Lambert function for negative lambertEccentricAnomaly.
double LambertFunctionsGooding::lambertFunctionNegativeGooding( const double xParameter )
{
    const double lambertEccentricAnomaly = xParameter * xParameter - 1.0;
    const double yParameter = std::sqrt( std::fabs( lambertEccentricAnomaly ) );
    const double zParameter = std::sqrt( 1.0 - qParameter * qParameter +
                                         qParameter * qParameter * xParameter * xParameter );
    const double fParameter = yParameter * ( zParameter - qParameter * xParameter );
    const double gParameter = xParameter * zParameter - qParameter * lambertEccentricAnomaly;
    const double dParameter = std::atan( fParameter / gParameter );

    return normalizedTimeOfFlight - 2.0 * ( xParameter - qParameter * zParameter - dParameter
                                            / yParameter ) / lambertEccentricAnomaly;
}

//! Define first derivative of Lambert function for positive lambertEccentricAnomaly.
double LambertFunctionsGooding::lambertFirstDerivativeFunctionPositiveGooding(
        const double xParameter )
{
    const double lambertEccentricAnomaly = xParameter * xParameter - 1.0;
    const double yParameter = std::sqrt( std::fabs( lambertEccentricAnomaly ) );
    const double zParameter = std::sqrt( 1.0 - qParameter * qParameter +
                                         qParameter * qParameter * xParameter * xParameter );
    const double fParameter = yParameter * ( zParameter - qParameter * xParameter );
    const double gParameter = xParameter * zParameter - qParameter * lambertEccentricAnomaly;
    const double dParameter = std::log( fParameter + gParameter );
    const double lambertEccentricAnomalyDerivative = 2.0 * xParameter;
    const double yParameterDerivative = xParameter / std::sqrt( xParameter * xParameter - 1.0 );
    const double zParameterDerivative = qParameter * qParameter * xParameter / zParameter;
    const double fParameterDerivative = yParameterDerivative  *
            ( zParameter - xParameter * qParameter ) + yParameter *
            ( zParameterDerivative - qParameter );
    const double gParameterDerivative = zParameter + xParameter *
            zParameterDerivative - qParameter * lambertEccentricAnomalyDerivative;
    const double dParameterDerivative = ( fParameterDerivative + gParameterDerivative ) /
            ( fParameter + gParameter );

    return - 2.0 * ( lambertEccentricAnomaly
                     * ( 1.0 - qParameter * zParameterDerivative
                         - ( ( dParameterDerivative * yParameter
                               - yParameterDerivative * dParameter )
                             / ( yParameter * yParameter ) ) )
                     - ( ( xParameter - qParameter * zParameter - ( dParameter / yParameter ) )
                         * lambertEccentricAnomalyDerivative ) )
            / ( lambertEccentricAnomaly * lambertEccentricAnomaly );
}

//! Define first derivative of Lambert function for negative lambertEccentricAnomaly.
double LambertFunctionsGooding::lambertFirstDerivativeFunctionNegativeGooding(
        const double xParameter )
{
    const double lambertEccentricAnomaly = xParameter * xParameter - 1.0;
    const double yParameter = std::sqrt( std::fabs( lambertEccentricAnomaly ) );
    const double zParameter = std::sqrt( 1.0 - qParameter * qParameter +
                                         qParameter * qParameter * xParameter * xParameter );
    const double fParameter = yParameter * ( zParameter - qParameter * xParameter );
    const double gParameter = xParameter * zParameter - qParameter * lambertEccentricAnomaly;
    const double dParameter = std::atan( fParameter / gParameter );
    const double lambertEccentricAnomalyDerivative = 2.0 * xParameter;
    const double yParameterDerivative = -1.0 * xParameter
            / std::sqrt( 1.0 - xParameter * xParameter);
    const double zParameterDerivative = qParameter * qParameter * xParameter / zParameter;
    const double fParameterDerivative = yParameterDerivative *
            ( zParameter - xParameter * qParameter ) + yParameter *
            ( zParameterDerivative - qParameter );
    const double gParameterDerivative = zParameter + xParameter *
            zParameterDerivative - qParameter *  lambertEccentricAnomalyDerivative;
    const double dParameterDerivative = ( fParameterDerivative * gParameter -
                                          fParameter * gParameterDerivative ) /
            ( fParameter * fParameter + gParameter * gParameter );

    return - 2.0 * ( lambertEccentricAnomaly * ( 1.0 - qParameter * zParameterDerivative -
                                                 ( ( dParameterDerivative * yParameter -
                                                     yParameterDerivative * dParameter ) /
                                                   yParameter / yParameter ) ) -
                     ( ( xParameter - qParameter * zParameter - ( dParameter / yParameter ) ) *
                       lambertEccentricAnomalyDerivative ) ) /
            ( lambertEccentricAnomaly * lambertEccentricAnomaly );
}

} // namespace mission_segments
} // namespace tudat
