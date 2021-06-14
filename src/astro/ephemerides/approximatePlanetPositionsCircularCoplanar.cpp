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
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 */

#include <vector>
#include <cmath>

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/ephemerides/approximatePlanetPositionsCircularCoplanar.h"

namespace tudat
{
namespace ephemerides
{

//! Get state from ephemeris; circular, coplanar case
Eigen::Vector6d ApproximateJplCircularCoplanarEphemeris::
getCartesianState( const double secondsSinceEpoch )
{
    // Set Julian date.
    julianDate_ = basic_astrodynamics::convertSecondsSinceEpochToJulianDay(
                secondsSinceEpoch, referenceJulianDate_ );

    // Compute number of centuries past J2000.
    numberOfCenturiesPastJ2000_ = ( julianDate_ - 2451545.0 ) / 36525.0;

    // Compute mean longitude of planet at given Julian date.
    meanLongitudeAtGivenJulianDate_ = approximatePlanetPositionsDataContainer_.meanLongitude_
            + ( approximatePlanetPositionsDataContainer_.rateOfChangeOfMeanLongitude_
                * numberOfCenturiesPastJ2000_ );

    // Convert mean longitude at given Julian date from degrees to radians.
    meanLongitudeAtGivenJulianDate_ =
            unit_conversions::convertDegreesToRadians(
                meanLongitudeAtGivenJulianDate_);

    // Get semi-major axis at J2000 and assume constant radius of circular orbit.
    constantOrbitalRadius_ =
            unit_conversions::convertAstronomicalUnitsToMeters(
                approximatePlanetPositionsDataContainer_.semiMajorAxis_ );

    // Convert to Cartesian position.
    Eigen::VectorXd planetCartesianStateAtGivenJulianDate( 6 );
    planetCartesianStateAtGivenJulianDate.segment( 0, 3 )
            = coordinate_conversions::convertSphericalToCartesian(
                Eigen::Vector3d( constantOrbitalRadius_,
                                 0.5 * mathematical_constants::PI,
                                 meanLongitudeAtGivenJulianDate_ ) );

    // Compute orbital velocity.
    double circularOrbitalVelocity = std::sqrt( ( sunGravitationalParameter_ + planetGravitationalParameter_ ) /
                                                constantOrbitalRadius_ );

    // Convert to Cartesian velocity.
    planetCartesianStateAtGivenJulianDate( 3 ) = -sin( meanLongitudeAtGivenJulianDate_ ) *
            circularOrbitalVelocity;
    planetCartesianStateAtGivenJulianDate( 4 ) = cos( meanLongitudeAtGivenJulianDate_ ) *
            circularOrbitalVelocity;
    planetCartesianStateAtGivenJulianDate( 5 ) = 0.0;

    // Return Cartesian state of planet at given Julian date.
    return planetCartesianStateAtGivenJulianDate;
}

} // namespace ephemerides
} // namespace tudat
