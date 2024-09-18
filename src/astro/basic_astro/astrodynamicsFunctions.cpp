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

#include <cmath>

#include "tudat/astro/basic_astro/astrodynamicsFunctions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/astro/basic_astro/orbitalElementConversions.h"

namespace tudat
{

namespace basic_astrodynamics
{

using mathematical_constants::PI;

//! Compute Kepler orbital period.
double computeKeplerOrbitalPeriod( const double semiMajorAxis,
                                   const double gravitationalParameterOfCentralBody,
                                   const double massOfOrbitingBody )
{
    return 2.0 * PI * std::sqrt( std::pow( semiMajorAxis, 3.0 )
                                 / ( ( physical_constants::GRAVITATIONAL_CONSTANT
                                       * massOfOrbitingBody )
                                     + gravitationalParameterOfCentralBody ) );
}

//! Compute two-body radial distance.
double computeKeplerRadialDistance( const double semiMajorAxis,
                                    const double eccentricity,
                                    const double trueAnomaly )
{
    return ( semiMajorAxis * ( 1.0 - eccentricity * eccentricity ) ) / ( 1.0 + eccentricity * std::cos( trueAnomaly ) );
}

//! Compute two-body radial distance.
double computeKeplerRadialDistance( const Eigen::Vector6d& keplerianElements )
{
    return ( keplerianElements[ orbital_element_conversions::semiMajorAxisIndex ] *
            ( 1.0 - keplerianElements[ orbital_element_conversions::eccentricityIndex ] *
            keplerianElements[ orbital_element_conversions::eccentricityIndex ] ) ) /
            ( 1.0 + keplerianElements[ orbital_element_conversions::eccentricityIndex ] *
            std::cos( keplerianElements[ orbital_element_conversions::trueAnomalyIndex ] ) );
}

//! Compute two-body orbital velocity with vis-viva equation.
double computeKeplerOrbitalVelocity( const double semiMajorAxis,
                                     const double eccentricity,
                                     const double trueAnomaly,
                                     const double gravitationalParameterOfCentralBody,
                                     const double massOfOrbitingBody )
{
    return std::sqrt( ( ( physical_constants::GRAVITATIONAL_CONSTANT * massOfOrbitingBody ) +
                        gravitationalParameterOfCentralBody ) * (
                          2.0 / computeKeplerRadialDistance( semiMajorAxis, eccentricity, trueAnomaly ) -
                          1.0 / semiMajorAxis ) );
}

//! Compute two-body orbital velocity with vis-viva equation.
double computeKeplerOrbitalVelocity( const Eigen::Vector6d& keplerianElements,
                                     const double gravitationalParameterOfCentralBody,
                                     const double massOfOrbitingBody )
{
    return std::sqrt( ( ( physical_constants::GRAVITATIONAL_CONSTANT * massOfOrbitingBody ) +
                        gravitationalParameterOfCentralBody ) * (
                          2.0 / computeKeplerRadialDistance( keplerianElements ) -
                          1.0 / keplerianElements[ orbital_element_conversions::semiMajorAxisIndex ] ) );
}

//! Compute Kepler angular momentum.
double computeKeplerAngularMomentum( const double semiMajorAxis, const double eccentricity,
                                     const double gravitationalParameterOfCentralBody,
                                     const double massOfOrbitingBody )
{
    return massOfOrbitingBody * std::sqrt( gravitationalParameterOfCentralBody * semiMajorAxis
                                           * ( 1.0 - std::pow( eccentricity, 2.0 ) ) );
}

//! Compute Kepler mean motion.
double computeKeplerMeanMotion( const double semiMajorAxis,
                                const double gravitationalParameterOfCentralBody,
                                const double massOfOrbitingBody )
{
    return std::sqrt( ( ( physical_constants::GRAVITATIONAL_CONSTANT * massOfOrbitingBody )
                        + gravitationalParameterOfCentralBody ) / std::pow( semiMajorAxis, 3.0 ) );
}

//! Compute Kepler orbital energy.
double computeKeplerEnergy( const double semiMajorAxis,
                            const double gravitationalParameterOfCentralBody,
                            const double massOfOrbitingBody )
{
    return -massOfOrbitingBody * gravitationalParameterOfCentralBody / ( 2.0 * semiMajorAxis );
}

//! Compute synodic period.
double computeSynodicPeriod( const double orbitalPeriodBody1, const double orbitalPeriodBody2 )
{
    return 1.0 / std::fabs( 1.0 / orbitalPeriodBody1 - 1.0 / orbitalPeriodBody2 );
}


//! Compute periapsis altitude from Keplerian state for spherical central body.
double computePeriapsisAltitudeFromKeplerianState( const Eigen::Vector6d& state,
                                                   const double centralBodyRadius )
{
    return state( orbital_element_conversions::semiMajorAxisIndex ) *
            ( 1.0 - state( orbital_element_conversions::eccentricityIndex ) ) - centralBodyRadius;
}

//! Compute periapsis altitude from Cartesian state for spherical central body.
double computePeriapsisAltitudeFromCartesianState( const Eigen::Vector6d& state,
                                                   const double centralBodyGravitationalParameter,
                                                   const double centralBodyRadius )
{
    return computePeriapsisAltitudeFromKeplerianState(
                orbital_element_conversions::convertCartesianToKeplerianElements(
                    state, centralBodyGravitationalParameter ), centralBodyRadius );
}

//! Compute apoapsis altitude from Keplerian state for spherical central body.
double computeApoapsisAltitudeFromKeplerianState( const Eigen::Vector6d& state,
                                                   const double centralBodyRadius )
{
    return state( orbital_element_conversions::semiMajorAxisIndex ) *
            ( 1.0 + state( orbital_element_conversions::eccentricityIndex ) ) - centralBodyRadius;
}

//! Compute apoapsis altitude from Cartesian state for spherical central body.
double computeApoapsisAltitudeFromCartesianState( const Eigen::Vector6d& state,
                                                   const double centralBodyGravitationalParameter,
                                                   const double centralBodyRadius )
{
    return computeApoapsisAltitudeFromKeplerianState(
                orbital_element_conversions::convertCartesianToKeplerianElements(
                    state, centralBodyGravitationalParameter ), centralBodyRadius );
}

} // namespace basic_astrodynamics

} // namespace tudat
