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
 *      120530    M.I. Ganeff       Code created.
 *      121004    M.I. Ganeff       Input parameter types and variable-naming updated.
 *      121018    M.I. Ganeff       Added computeSphereOfInfluence().
 *      121123    D. Dirkx          Added computeSphereOfInfluence() function taking mass ratios.
 *      130225    D. Dirkx          Added isOrbitRetrograde(...) functions taking inclination and
 *                                  Kepler vector.
 *      130301    R.C.A. Boon       Minor textual changes.
 *      130305    R.C.A. Boon       Replaced Eigen::VectorXd by Eigen::Vector6d.
 *      131212    S. Billemont      Fixed pass-by-reference error in isOrbitRetrograde()-function.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Corrected Third Printing, Springer, 2005.
 *      Bate R. Fundamentals of Astrodynamics, Courier Dover Publications, 1971.
 *
 *    Notes
 *      In future, it might make sense to create a new MissionGeometry/ sub-directory if a lot more
 *      funtionality is added.
 *
 */

#ifndef TUDAT_MISSION_GEOMETRY_H
#define TUDAT_MISSION_GEOMETRY_H

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace mission_geometry
{

//! Compute whether an orbit is retrograde based on inclination.
/*!
 * Determines whether an orbit is retrograde, based on its inclination. The orbit is defined to
 * be retrograde if its inclination is >90 and <180 degrees, prograde if >0 and <=90 degrees.
 * Values outside [0, 180] degrees are not allowed for this function. The longitude of the ascending
 * node does not influence whether the orbit is retrograde.
 *
 * N.B: The inclination must be given in radians!
 *
 * \param inclination Inclination of the orbit [rad].
 * \return true if orbit is retrograde, false if prograde.
 */
bool isOrbitRetrograde( const double inclination );

//! Compute whether an orbit is retrograde based on Keplerian state.
/*!
 * Determines whether an orbit is retrograde, based on the kepler elements. The orbit is defined to
 * be retrograde if its inclination is >90 and <180 degrees, prograde if >0 and <=90 degrees.
 * Values for the inclination outside [0, 180] degrees are not allowed for this function. The right
 * ascension of the ascending node does not influence whether the orbit is retrograde.
 *
 * N.B: The inclination must be given in radians!
 *
 * \param keplerElements Vector of keplerian elements.
 * \return true if orbit is retrograde, false if prograde.
 */
bool isOrbitRetrograde( const Eigen::Vector6d& keplerElements );

//! Compute the shadow function.
/*!
 * Returns the value of of the shadow function. Returns 0 if the satellite is in umbra, 1 if the
 * satellite is fully exposed and a value between 0 and 1 if the satellite is in penumbra.
 *
 * The point of view is from the satellite. The occulting body (for example the Earth) is the body
 * that blocks the light from the occulting body (for example the Sun).
 *
 * Reference: Section 3.4 from ( Montebruck O, Gill E., 2005).
 *
 * \param occultedBodyPosition Vector containing Cartesian coordinates of the occulted body.
 * \param occultedBodyRadius Mean radius of occulted body.
 * \param occultingBodyPosition Vector containing Cartesian coordinates of the occulting body.
 * \param occultingBodyRadius Mean radius of occulting body.
 * \param satellitePosition Vector containing Cartesian coordinates of the satellite.
 * \return Shadow function value.
 */
double computeShadowFunction( const Eigen::Vector3d& occultedBodyPosition,
                              const double occultedBodyRadius,
                              const Eigen::Vector3d& occultingBodyPosition,
                              const double occultingBodyRadius,
                              const Eigen::Vector3d& satellitePosition );

//! Compute the radius of the sphere of influence.
/*!
 * Returns the radius of the the Sphere of Influence (SOI) for a body orbiting a central body.
 *
 * Reference: Section 7.4 from (Fundamentals of Astrodynamics, R. Bate, 1971).
 *
 * \param distanceToCentralBody Distance from the orbiting body to the central body [m].
 * \param ratioOfOrbitingToCentralBodyMass ratio of mass of orbiting body (i.e. Earth)
 *          to that of central body (i.e. Sun).
 * \return Radius of sphere of influence [m].
 */
double computeSphereOfInfluence( const double distanceToCentralBody,
                                 const double ratioOfOrbitingToCentralBodyMass );

//! Compute the radius of the sphere of influence.
/*!
 * Returns the radius of the the Sphere of Influence (SOI) for a body orbiting a central body.
 *
 * Reference: Section 7.4 from (Fundamentals of Astrodynamics, R. Bate, 1971).
 *
 * \param distanceToCentralBody Distance from the orbiting body to the central body.
 * \param massOrbitingBody Mass of orbiting body (i.e. Earth) [kg].
 * \param massCentralBody Mass of central body (i.e. Sun) [kg].
 * \return Radius of sphere of influence [m].
 */
double computeSphereOfInfluence( const double distanceToCentralBody,
                                 const double massOrbitingBody,
                                 const double massCentralBody );

} // namespace mission_geometry

} // namespace tudat

#endif // TUDAT_MISSION_GEOMETRY_H
