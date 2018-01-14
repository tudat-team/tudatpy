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

#ifndef TUDAT_KEPLER_PROPAGATOR_H
#define TUDAT_KEPLER_PROPAGATOR_H

#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"


namespace tudat
{

namespace orbital_element_conversions
{

//! Propagate Kepler orbit.
/*!
 * Propagates Kepler orbit. This function essentially takes a state in classical Keplerian elements
 * at an initial epoch and propagates it to a final state at a given final epoch. Currently both
 * elliptic and hyperbolic orbits are supported. Parabolic orbits are not supported and will result
 * in an error message.
 * IMPORTANT! Note that the true anomaly is returned within the -PI to PI spectrum. If the user
 * desires a different spectrum (possibly including the number of revolutions), these should be
 * added by the user a posteriori.
 * \param initialStateInKeplerianElements Initial state vector in classical Keplerian elements.
 *          Order is important!
 *          initialStateInKeplerianElements( 0 ) = semiMajorAxis,                               [m]
 *          initialStateInKeplerianElements( 1 ) = eccentricity,                                [-]
 *          initialStateInKeplerianElements( 2 ) = inclination,                               [rad]
 *          initialStateInKeplerianElements( 3 ) = argument of periapsis,                     [rad]
 *          initialStateInKeplerianElements( 4 ) = longitude of ascending node,               [rad]
 *          initialStateInKeplerianElements( 5 ) = true anomaly.                              [rad]
 * \param propagationTime Propagation time.                                                     [s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body      [m^3 s^-2]
 * \param aRootFinder Shared-pointer to the root-finder that is used to solve the conversion from
 *          mean to eccentric anomaly. Default is Newton-Raphson using 5.0e-14 absolute X-tolerance
 *          and 1000 iterations as maximum. Higher precision may invoke machine precision
 *          problems for some values.
 * \return finalStateInKeplerianElements Final state vector in classical Keplerian elements.
 *          Order is important!
 *          finalStateInKeplerianElements( 0 ) = semiMajorAxis,                                 [m]
 *          finalStateInKeplerianElements( 1 ) = eccentricity,                                  [-]
 *          finalStateInKeplerianElements( 2 ) = inclination,                                 [rad]
 *          finalStateInKeplerianElements( 3 ) = argument of periapsis,                       [rad]
 *          finalStateInKeplerianElements( 4 ) = longitude of ascending node,                 [rad]
 *          finalStateInKeplerianElements( 5 ) = true anomaly.                                [rad]
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > propagateKeplerOrbit(
        const Eigen::Matrix< ScalarType, 6, 1 >& initialStateInKeplerianElements,
        const ScalarType propagationTime,
        const ScalarType centralBodyGravitationalParameter,
        boost::shared_ptr< root_finders::RootFinderCore< ScalarType > > aRootFinder =
        boost::shared_ptr< root_finders::RootFinderCore< ScalarType > >( ) )
{
    // Create final state in Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > finalStateInKeplerianElements =
            initialStateInKeplerianElements;

    // Check if eccentricity is valid.
    if ( initialStateInKeplerianElements( eccentricityIndex ) <
         mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        throw std::runtime_error( "Eccentricity is invalid (smaller than 0)." );
    }

    // Check if orbit is elliptical.
    else if ( initialStateInKeplerianElements( eccentricityIndex ) <
              mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        // Convert initial true anomaly to eccentric anomaly.
        const ScalarType initialEccentricAnomaly =
                convertTrueAnomalyToEccentricAnomaly< ScalarType >(
                    initialStateInKeplerianElements( trueAnomalyIndex ),
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Convert initial eccentric anomaly to mean anomaly.
        const ScalarType initialMeanAnomaly =
                convertEccentricAnomalyToMeanAnomaly< ScalarType >(
                    initialEccentricAnomaly,
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Compute change of mean anomaly between start and end of propagation.
        const ScalarType meanAnomalyChange =
                convertElapsedTimeToEllipticalMeanAnomalyChange< ScalarType >(
                    propagationTime, centralBodyGravitationalParameter,
                    initialStateInKeplerianElements( semiMajorAxisIndex ) );

        // Compute eccentric anomaly for mean anomaly.
        const ScalarType finalEccentricAnomaly =
                convertMeanAnomalyToEccentricAnomaly< ScalarType >(
                    initialStateInKeplerianElements( eccentricityIndex ),
                    initialMeanAnomaly + meanAnomalyChange, true,
                    TUDAT_NAN, aRootFinder );

        // Compute true anomaly for computed eccentric anomaly.
        finalStateInKeplerianElements( trueAnomalyIndex ) =
                convertEccentricAnomalyToTrueAnomaly< ScalarType >(
                    finalEccentricAnomaly, finalStateInKeplerianElements( eccentricityIndex ) );
    }

    // Check if the orbit is hyperbolic.
    else if ( initialStateInKeplerianElements( eccentricityIndex ) >
              mathematical_constants::getFloatingInteger< ScalarType >( 1 ))
    {
        // Convert initial true anomaly to hyperbolic eccentric anomaly.
        const ScalarType initialHyperbolicEccentricAnomaly =
                convertTrueAnomalyToHyperbolicEccentricAnomaly(
                    initialStateInKeplerianElements( trueAnomalyIndex ),
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Convert initial hyperbolic eccentric anomaly to the hyperbolic mean anomaly.
        const ScalarType initialHyperbolicMeanAnomaly =
                convertHyperbolicEccentricAnomalyToMeanAnomaly(
                    initialHyperbolicEccentricAnomaly,
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Compute change of mean anomaly because of the propagation time.
        const ScalarType hyperbolicMeanAnomalyChange =
                convertElapsedTimeToHyperbolicMeanAnomalyChange(
                    propagationTime, centralBodyGravitationalParameter,
                    initialStateInKeplerianElements( semiMajorAxisIndex ) );

        // Compute hyperbolic eccentric anomaly for mean anomaly.
        const ScalarType finalHyperbolicEccentricAnomaly =
                convertMeanAnomalyToHyperbolicEccentricAnomaly(
                    initialStateInKeplerianElements( eccentricityIndex ),
                    initialHyperbolicMeanAnomaly + hyperbolicMeanAnomalyChange );

        // Compute true anomaly for computed hyperbolic eccentric anomaly.
        finalStateInKeplerianElements( trueAnomalyIndex ) =
                convertHyperbolicEccentricAnomalyToTrueAnomaly(
                               finalHyperbolicEccentricAnomaly,
                               finalStateInKeplerianElements( eccentricityIndex ) );
    }

    // In this case the eccentricity has to be 1.0, hence the orbit is parabolic.
    else
    {
        throw std::runtime_error( "Parabolic orbits are not (yet) supported." );
    }

    return finalStateInKeplerianElements;
}

} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_KEPLER_PROPAGATOR_H
