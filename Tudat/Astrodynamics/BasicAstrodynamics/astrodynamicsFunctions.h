/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      Backwards compatibility of namespaces is implemented for Tudat Core 2 in this file. The
 *      code block marked "DEPRECATED!" at the end of the file should be removed in Tudat Core 3.
 *
 */

#ifndef TUDAT_ASTRODYNAMICS_FUNCTIONS_H
#define TUDAT_ASTRODYNAMICS_FUNCTIONS_H

#include <Eigen/Eigen>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Compute two-body orbital period.
/*!
 * Computes the two-body orbital period of an orbiting body that follows a closed conic section
 * (circle or ellipse Kepler orbit). The mass of the orbiting body is set to that of a test
 * particle by default.
 * \param semiMajorAxis Semi-major axis of Kepler orbit (circle or ellipse).
 * \param gravitationalParameterOfCentralBody Gravitational parameter of central body.
 * \param massOfOrbitingBody Mass of orbiting body.
 * \return Two-body orbital period.
 */
double computeKeplerOrbitalPeriod( const double semiMajorAxis,
                                   const double gravitationalParameterOfCentralBody,
                                   const double massOfOrbitingBody = 0.0 );

//! Compute two-body angular momentum.
/*!
 * Computes the angular momentum of an orbiting body that follows a conic section (Kepler orbit), 
 * relative to the center-of-mass of the central body. The default mass value is for the angular 
 * momentum per unit mass.
 * \param semiMajorAxis Semi-major axis of Kepler orbit.
 * \param eccentricity Eccentricity of Kepler orbit.
 * \param gravitationalParameterOfCentralBody Gravitational parameter of central body.
 * \param massOfOrbitingBody Mass of orbiting body.
 * \return Two-body angular momentum.
 */
double computeKeplerAngularMomentum( const double semiMajorAxis, const double eccentricity,
                                     const double gravitationalParameterOfCentralBody,
                                     const double massOfOrbitingBody = 1.0 );

//! Compute two-body mean motion.
/*!
 * Computes the two-body mean motion of an orbiting body that follows a conic section
 * (Kepler orbit). The mass of the orbiting body is set to that of a test particle by default.
 * \param semiMajorAxis Semi-major axis of Kepler orbit.
 * \param gravitationalParameterOfCentralBody Gravitational parameter of central body.
 * \param massOfOrbitingBody Mass of orbiting body.
 * \return Two-body mean motion.
 */
double computeKeplerMeanMotion( const double semiMajorAxis,
                                 const double gravitationalParameterOfCentralBody,
                                 const double massOfOrbitingBody = 0.0 );

//! Compute Kepler energy.
/*!
 * Computes the energy of an orbiting body that follows a conic section (Kepler orbit). The 
 * default mass value is for the two-body orbital energy per unit mass. For closed conic sections 
 * (circles, ellipses), the semi-major axis is positive, and for open sections (hyperbolas) the 
 * semi-major axis is negative.
 * \param semiMajorAxis Semi-major axis of Kepler orbit.
 * \param gravitationalParameterOfCentralBody Gravitational parameter of central body.
 * \param massOfOrbitingBody Mass of orbiting body.
 * \return Kepler orbital energy.
 */
double computeKeplerEnergy( const double semiMajorAxis,
                            const double gravitationalParameterOfCentralBody,
                            const double massOfOrbitingBody = 1.0 );

//! Compute synodic period.
/*!
 * Computes synodic period between two bodies in different Kepler orbits (closed conic sections).
 * The orbital periods must be positive values for the synodic period to be sensible.
 * \param orbitalPeriodBody1 Orbital period of Body 1.
 * \param orbitalPeriodBody2 Orbital period of Body 2.
 * \return Synodic period.
 */
double computeSynodicPeriod( const double orbitalPeriodBody1, const double orbitalPeriodBody2 );


//! Compute periapsis altitude from Keplerian state for spherical central body.
/*!
 * Compute periapsis altitude from Keplerian state for spherical central body.
 * \param state Keplerian state of the propagated body.
 * \param centralBodyRadius Radius of the central body (assumed spherical).
 * \return The distance from the propagated body to the central body's spherical surface at periapsis.
 */
double computePeriapsisAltitudeFromKeplerianState( const Eigen::Vector6d& state,
                                                   const double centralBodyRadius );

//! Compute periapsis altitude from Cartesian state for spherical central body.
/*!
 * Compute periapsis altitude from Cartesian state for spherical central body.
 * \param state Cartesian state of the propagated body.
 * \param centralBodyGravitationalParameter Gravitational parameter of the central body.
 * \param centralBodyRadius Radius of the central body (assumed spherical).
 * \return The distance from the propagated body to the central body's spherical surface at periapsis.
 */
double computePeriapsisAltitudeFromCartesianState( const Eigen::Vector6d& state,
                                                   const double centralBodyGravitationalParameter,
                                                   const double centralBodyRadius );


} // namespace basic_astrodynamics
} // namespace tudat


#endif // TUDAT_ASTRODYNAMICS_FUNCTIONS_H
