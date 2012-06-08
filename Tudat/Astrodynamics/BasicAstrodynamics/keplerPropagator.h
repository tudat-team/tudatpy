/*    Copyright (c) 2010-2012, Delft University of Technology
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
 * \param propagationTime Propagation time.                                                     [s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body      [m^3 s^-2]
 * \param newtonRaphsonConvergenceTolerance Convergence tolerance for Newton-Raphson
 *          root-finder. This quantity represents the absolute difference in solution for the
 *          eccentric anomaly in the Newton-Raphson root-finding solution when converting
 *          eccentric to mean anomaly.                                                          [-]
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
                                      const double propagationTime,
                                      const double centralBodyGravitationalParameter,
                                      const double newtonRaphsonConvergenceTolerance,
                                      bool useModuloOption = true );

} // namespace orbital_element_conversions
} // namespace tudat

#endif // TUDAT_KEPLER_PROPAGATOR_H
