/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110203    K. Kumar          File created.
 *      110207    E. Iorfida        Minor changes.
 *      110214    K. Kumar          Updated code to use orbital element conversion functions.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 *      120215    K. Kumar          Rewrote Kepler propagator as free function.
 *
 *    References
 *
 */

#ifndef TUDAT_KEPLER_PROPAGATOR_H
#define TUDAT_KEPLER_PROPAGATOR_H

#include <Eigen/Core>

namespace tudat
{
namespace orbital_element_conversions
{

//! Propagate Kepler orbit.
/*!
 * Propagates Kepler orbit. This function essentially takes a state in classical Keplerian elements
 * at an initial epoch and propagates it to a final state at a given final epoch. Currently,
 * this implementation only supports elliptical orbits ( 0 < e < 0.98 ). An error is thrown for all
 * other eccentricities.
 * \param initialStateInKeplerianElements Initial state vector in classical Keplerian elements.
 *          Order is important!
 *          initialStateInKeplerianElements( 0 ) = semiMajorAxis,                               [m]
 *          initialStateInKeplerianElements( 1 ) = eccentricity,                                [-]
 *          initialStateInKeplerianElements( 2 ) = inclination,                               [rad]
 *          initialStateInKeplerianElements( 3 ) = argument of periapsis,                     [rad]
 *          initialStateInKeplerianElements( 4 ) = longitude of ascending node,               [rad]
 *          initialStateInKeplerianElements( 5 ) = true anomaly.                              [rad]
 * \param epochOfInitialState Epoch of initial state.                                           [s]
 * \param epochOfFinalState Epoch of final state.                                               [s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body      [m^3 s^-2]
 * \param newtonRaphsonConvergenceTolerance Convergence tolerance for Newton-Raphson
 *          root-finder. This quantity represents the absolute difference in solution for the
 *          eccentric anomaly in the Newton-Raphson root-finding solution when converting
 *          eccentric to mean anomaly                                                           [-]
 * \param useModuloOption Option to propagate remainder time computed from
 *          mod( propagationTime, orbitalPeriod ). This has computational advantages when the
 *          angles in Kepler's equation become very large. The default is set to true.
 * \return finalStateInKeplerianElements Final state vector in classical Keplerian elements.
 *          Order is important!
 *          finalStateInKeplerianElements( 0 ) = semiMajorAxis,                                 [m]
 *          finalStateInKeplerianElements( 1 ) = eccentricity,                                  [-]
 *          finalStateInKeplerianElements( 2 ) = inclination,                                 [rad]
 *          finalStateInKeplerianElements( 3 ) = argument of periapsis,                       [rad]
 *          finalStateInKeplerianElements( 4 ) = longitude of ascending node,                 [rad]
 *          finalStateInKeplerianElements( 5 ) = true anomaly.                                [rad]
 */
Eigen::VectorXd propagateKeplerOrbit( const Eigen::VectorXd& initialStateInKeplerianElements,
                                      const double epochOfInitialState,
                                      const double epochOfFinalState,
                                      const double centralBodyGravitationalParameter,
                                      const double newtonRaphsonConvergenceTolerance,
                                      bool useModuloOption = true );

} // namespace orbital_element_conversions
} // namespace tudat

#endif // TUDAT_KEPLER_PROPAGATOR_H
