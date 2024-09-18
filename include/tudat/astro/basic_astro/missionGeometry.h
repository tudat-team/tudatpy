/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Corrected Third Printing, Springer, 2005.
 *      Bate R. Fundamentals of astro, Courier Dover Publications, 1971.
 *      Rui Zhang, et al. "Study of satellite shadow function model considering the overlapping parts of Earth shadow
 *          and Moon shadow and its application to GPS satellite orbit determination".
 *          Advances in Space Research 63. 9(2019): 2912–2929.
 *
 *    Notes
 *      In future, it might make sense to create a new MissionGeometry/ sub-directory if a lot more
 *      funtionality is added.
 *
 */

#ifndef TUDAT_MISSION_GEOMETRY_H
#define TUDAT_MISSION_GEOMETRY_H

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"

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
 * satellite is fully exposed and a value between 0 and 1 if the satellite is in penumbra or antumbra.
 *
 * The point of view is from the satellite. The occulting body (for example the Earth) is the body
 * that blocks the light from the occulting body (for example the Sun).
 *
 * Reference: Section 3.4 from ( Montebruck O, Gill E., 2005) and Fig. 5 from (Zhang et al., 2019).
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
 * Reference: Section 7.4 from (Fundamentals of astro, R. Bate, 1971).
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
 * Reference: Section 7.4 from (Fundamentals of astro, R. Bate, 1971).
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
