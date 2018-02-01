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
 *      Vallado, D.A. Fundamentals of Astrodynamics and Applications. Microcosm Press, 2001.
 *
 */

#ifndef TUDAT_CLOHESSY_WILTSHIRE_PROPAGATOR_H
#define TUDAT_CLOHESSY_WILTSHIRE_PROPAGATOR_H

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! Propagate Clohessy-Wilshire equations (linearized relative motion).
/*!
 * This function propagates linearized relative motion, based on the Clohessy-Wiltshire equations.
 * It calculates the motion of a point mass A with respect to a local-vertical-local-horizontal
 * reference frame centered on a point mass B.
 *
 * The following assumptions apply:
 *     (1) Masses A and B are negligable compared to the central body mass (e.g. mass A is a
 *         daughter spacecraft, mass B is a mother spacecraft, and the central body is planet
 *         Earth).
 *     (2) Orbit perturbations are negligable (i.e. massses A and B follow Kepler orbits about
 *         the central body).
 *     (3) The orbit of mass B is circular.
 *     (4) The separation between mass A and mass B is very small compared to the circumference of
 *         the orbit of mass B.
 *
 * The Clohessy-Wilthsire equations are given by Vallado [2001] as:
 *
 * \f{eqnarray*}{
 *     x( t ) = \frac{ \dot{ x }_0 }{ n } \sin( n t )
 *         - \left( 3 x_0 + \frac{ 2 \dot{ y }_0 } { n } \right) \cos( n t )
 *         + \left( 4 x_0 + \frac{ 2 \dot{ y }_0 }{ n } \right) \\
 *     y( t ) = \left( 6 x_0 + \frac{ 4 \dot{ y_0 } }{ n } \right) \sin( n t )
 *              + \frac{ 2 \dot{ x }_0 }{ n } \cos( n t )
 *              - ( 6 x_0 n + 3 \dot{ y_0 } ) t
 *              + \left( y_0 - \frac{ 2 \dot{ x }_0 }{ n } \right) \\
 *     z( t ) = \frac{ \dot{ z }_0 }{ n } \sin( n t )
 *              + z_0 \cos ( n t ) \\
 *     \dot{ x }( t ) = ( 3 x_0 n + 2 \dot{ y }_0 ) \sin( n t )
 *                      + \dot{ x }_0 \cos( n t ) \\
 *     \dot{ y }( t ) = - 2 \dot{ x }_0 \sin( n t )
 *                      + ( 6 x_0 n + 4 \dot{ y }_0 ) \cos( n t)
 *                      - ( 6 n x_0 + 3 \dot{ y }_0 ) \\
 *     \dot{ z }( t ) = - z_0 n \sin( n t )
 *                      + \dot{ z }_0 \cos( n t )
 * \f}
 *
 * in which \f$ x( t ) \f$, \f$ y( t ) \f$ and \f$ z( t ) \f$ are the radial, along-track and
 * cross-track position respectively, of mass A. Parameters \f$ x_0 \f$, \f$ y_0 \f$ and
 * \f$ z_0 \f$ are the initial radial, along-track and cross-track position respectively.
 * Parameter \f$ n \f$ is the mean angular motion of the circular orbit of mass B. Finally,
 * parameter \f$ t \f$ is time.
 *
 * \param initialState Initial state vector in Cartesian elements.
 *          The order is important!
 *          initialState( 0 ) = radial position                                                [m],
 *          initialState( 1 ) = along-track position                                           [m],
 *          initialState( 2 ) = cross-track position                                           [m],
 *          initialState( 3 ) = radial velocity                                              [m/s],
 *          initialstate( 4 ) = along-track velocity                                         [m/s],
 *          initialState( 5 ) = cross-track velocity                                         [m/s],
 * \param propagationDuration Duration of propagation                                          [s].
 * \param centralBodyGravitationalParameter Gravitational parameter of central body     [m^3 s^-2].
 * \param referenceOrbitRadius Radius of circular orbit of mass B                              [m].
 * \return Final state vector in Cartesian elements.
 *          The order is important!
 *          finalState( 0 ) = radial position                                                  [m],
 *          finalState( 1 ) = along-track position                                             [m],
 *          finalState( 2 ) = cross-track position                                             [m],
 *          finalState( 3 ) = radial velocity                                                [m/s],
 *          finalState( 4 ) = along-track velocity                                           [m/s],
 *          finalState( 5 ) = cross-track velocity                                           [m/s].
 */

Eigen::Vector6d propagateClohessyWiltshire(
        const Eigen::Vector6d& initialState,
        const double propagationDuration,
        const double centralBodyGravitationalParameter,
        const double referenceOrbitRadius );

} // namespace basic_astrodynamics
} // namespace tudat

#endif // TUDAT_CLOHESSY_WILTSHIRE_PROPAGATOR_H
