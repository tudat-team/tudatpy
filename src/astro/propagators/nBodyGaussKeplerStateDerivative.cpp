/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/propagators/nBodyGaussKeplerStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate the Gauss planetary equations for Kepler elements
Eigen::Vector6d computeGaussPlanetaryEquationsForKeplerElements(
        const Eigen::Vector6d& currentOsculatingKeplerElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double semiLatusRectum,
        const double distance,
        const double meanMotion,
        const double orbitalAngularMomentum )
{
    // Retrieve Kepler elements and precompute reused terms
    double semiMajorAxis = currentOsculatingKeplerElements( 0 );
    double eccentricity = currentOsculatingKeplerElements( 1 );
    double inclination = currentOsculatingKeplerElements( 2 );

    double meanMotionTimesASquared = meanMotion * semiMajorAxis * semiMajorAxis;
    double eccentricityTerm = std::sqrt( 1.0 - eccentricity * eccentricity );
    double sineTrueAnomaly = std::sin( currentOsculatingKeplerElements( 5 ) );
    double cosineTrueAnomaly = std::cos( currentOsculatingKeplerElements( 5 ) );
    double argumentOfLatitude = currentOsculatingKeplerElements( 5 ) + currentOsculatingKeplerElements( 3 );

    Eigen::Vector6d stateDerivative;

    // Evaluate Gauss equations.
    stateDerivative( 0 ) = 2.0 / ( meanMotion * eccentricityTerm ) * (
                eccentricity * sineTrueAnomaly * accelerationsInRswFrame( 0 ) +
                semiLatusRectum / distance * accelerationsInRswFrame( 1 ) );

    stateDerivative( 1 ) = eccentricityTerm / ( meanMotion * semiMajorAxis ) * (
                sineTrueAnomaly * accelerationsInRswFrame( 0 ) +
                ( cosineTrueAnomaly + ( eccentricity + cosineTrueAnomaly ) /
                  (  1.0 + eccentricity * cosineTrueAnomaly ) ) *  accelerationsInRswFrame( 1 ) );

    stateDerivative( 2 ) = distance * std::cos( argumentOfLatitude ) /
            ( meanMotionTimesASquared * eccentricityTerm ) *
            accelerationsInRswFrame( 2 );

    stateDerivative( 3 ) = eccentricityTerm /
            ( meanMotion * semiMajorAxis * eccentricity ) * (
                -cosineTrueAnomaly * accelerationsInRswFrame( 0 ) +
                sineTrueAnomaly * ( 1.0 + distance / semiLatusRectum ) * accelerationsInRswFrame( 1 ) ) -
            distance * std::sin( argumentOfLatitude ) /
            ( orbitalAngularMomentum * std::tan( inclination ) ) *
            accelerationsInRswFrame( 2 );

    stateDerivative( 4 ) = distance * std::sin( argumentOfLatitude ) /
            ( meanMotionTimesASquared * eccentricityTerm * std::sin( inclination ) ) *
            accelerationsInRswFrame( 2 );

    stateDerivative( 5 ) = meanMotion + 1.0 / ( meanMotionTimesASquared * eccentricity ) * (
                ( semiLatusRectum * cosineTrueAnomaly - 2.0 * eccentricity * distance ) * accelerationsInRswFrame( 0 ) -
                ( semiLatusRectum + distance ) * sineTrueAnomaly * accelerationsInRswFrame( 1 ) );

    return stateDerivative;

}

//! Function to evaluate the Gauss planetary equations for Kepler elements
Eigen::Vector6d computeGaussPlanetaryEquationsForKeplerElements(
        const Eigen::Vector6d& currentOsculatingKeplerElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter )
{
    // Compute orbit properties
    double semiLatusRectum =  orbital_element_conversions::computeSemiLatusRectum(
                currentOsculatingKeplerElements( 1 ),
                currentOsculatingKeplerElements( 0 ),
                std::numeric_limits< double >::epsilon( ) );
    double orbitalAngularMomentum = orbital_element_conversions::computeOrbitalAngularMomentumPerUnitMass(
                semiLatusRectum, centralBodyGravitationalParameter );
    double meanMotion = basic_astrodynamics::computeKeplerMeanMotion(
                currentOsculatingKeplerElements( 0 ), centralBodyGravitationalParameter );
    double distance = semiLatusRectum /
            ( 1.0 + currentOsculatingKeplerElements( 1 ) * std::cos( currentOsculatingKeplerElements( 5 ) ) );

    // Evaluate Gauss equations
    return computeGaussPlanetaryEquationsForKeplerElements(
                currentOsculatingKeplerElements, accelerationsInRswFrame, semiLatusRectum, distance, meanMotion,
                orbitalAngularMomentum );
}

//! Function to evaluate the Gauss planetary equations for Kepler elements
Eigen::Vector6d computeGaussPlanetaryEquationsForKeplerElements(
        const Eigen::Vector6d& currentOsculatingKeplerElements,
        const Eigen::Vector6d& currentCartesianState,
        const Eigen::Vector3d& accelerationsInInertialFrame,
        const double centralBodyGravitationalParameter )
{

    return  computeGaussPlanetaryEquationsForKeplerElements(
                currentOsculatingKeplerElements,
                reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                    currentCartesianState ) * accelerationsInInertialFrame, centralBodyGravitationalParameter );
}

//template class NBodyGaussKeplerStateDerivative< double, double >;

} // namespace propagators

} // namespace tudat
