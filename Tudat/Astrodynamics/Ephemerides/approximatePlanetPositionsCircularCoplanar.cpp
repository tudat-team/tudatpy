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
 *      110629    L. van der Ham    File created.
 *      110803    L. van der Ham    Seperated this code from approximatePlanetPositions.
 *      120322    D. Dirkx          Modified to new Ephemeris interfaces.
 *      130120    D. Dirkx          Updated with new Julian day + seconds since Julian day input.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 *    Notes
 *
 */

#include <cmath>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositionsCircularCoplanar.h"

namespace tudat
{
namespace ephemerides
{

//! Get state from ephemeris; circular, coplanar case
basic_mathematics::Vector6d ApproximatePlanetPositionsCircularCoplanar::
getCartesianStateFromEphemeris( const double secondsSinceEpoch, const double julianDayAtEpoch )
{
    // Set Julian date.
    julianDate_ = basic_astrodynamics::convertSecondsSinceEpochToJulianDay(
                secondsSinceEpoch, julianDayAtEpoch );

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
    double circularOrbitalVelocity = std::sqrt( sunGravitationalParameter /
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
