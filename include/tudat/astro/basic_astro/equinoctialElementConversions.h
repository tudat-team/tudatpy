/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and binary
 *    forms, with or without modification, are permitted exclusively under the
 *    terms of the Modified BSD license.
 */

#ifndef TUDAT_EQUINOCTIAL_ELEMENT_CONVERSIONS_H
#define TUDAT_EQUINOCTIAL_ELEMENT_CONVERSIONS_H

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace orbital_element_conversions
{

//! Convert orbital elements with mean motion and mean anomaly to nonsingular equinoctial parameters.
/*!
 * The input order is mean motion, eccentricity, inclination, right ascension of the ascending node, argument of
 * periapsis and mean anomaly. This is explicitly not Tudat's conventional Keplerian vector, whose first and final
 * entries are semi-major axis and true anomaly. The output is log mean motion, the two eccentricity components, the
 * two inclination components and mean longitude.
 *
 * The prograde definition follows Vallado and Crawford, SGP4 Orbit Determination, AIAA 2008-6770, Section II.B.
 * When useRetrogradeElements is true, cot(i/2) and the retrograde longitude of perigee are used to avoid the
 * singularity at an inclination of pi.
 */
Eigen::Vector6d convertOrbitalElementsWithMeanMotionAndMeanAnomalyToEquinoctialElements(
        const Eigen::Vector6d& orbitalElementsWithMeanMotionAndMeanAnomaly,
        const bool useRetrogradeElements );

//! Convert nonsingular equinoctial parameters to orbital elements with mean motion and mean anomaly.
/*!
 * This function is the inverse of convertOrbitalElementsWithMeanMotionAndMeanAnomalyToEquinoctialElements. Returned
 * angular elements are normalized to the interval [0, 2 pi).
 */
Eigen::Vector6d convertEquinoctialToOrbitalElementsWithMeanMotionAndMeanAnomaly( const Eigen::Vector6d& equinoctialElements,
                                                                                 const bool useRetrogradeElements );

}  // namespace orbital_element_conversions
}  // namespace tudat

#endif  // TUDAT_EQUINOCTIAL_ELEMENT_CONVERSIONS_H
