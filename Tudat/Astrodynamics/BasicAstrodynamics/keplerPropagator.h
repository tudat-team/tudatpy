/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110203    K. Kumar          File created.
 *      110207    E. Iorfida        Minor changes.
 *      110214    K. Kumar          Updated code to use orbital element conversion functions.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 *      120215    K. Kumar          Rewrote Kepler propagator as free function.
 *      120607    P. Musegaas       Changed interface (propagation time instead of two epochs).
 *      120813    P. Musegaas       Changed code to new root finding structure. Added option to
 *                                  specify which rootfinder and termination conditions to use.
 *      120823    P. Musegaas       Added functionality for hyperbolic and near-parabolic orbits.
 *                                  Changed some parameters to const.
 *      120903    P. Musegaas       Removed modulo option, due to errors with it. Kepler propagator
 *                                  now simply return true anomaly in -PI to PI spectrum. Added
 *                                  Comments.
 *      121205    P. Musegaas       Updated code to final version of rootfinders.
 *      130120    K. Kumar          Updated VectorXd to Vector6d.
 *      150417    D. Dirkx          Made modifications for templated element conversions.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_KEPLER_PROPAGATOR_H
#define TUDAT_KEPLER_PROPAGATOR_H

#include <boost/exception/all.hpp>
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
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid (smaller than 0)." ) ) );
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
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Parabolic orbits are not (yet) supported." ) ) );
    }

    return finalStateInKeplerianElements;
}

} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_KEPLER_PROPAGATOR_H
