/*! \file approximatePlanetPositionsCircularCoplanar.cpp
 *    This source file contains the definition of an ephemeris class that makes use of the JPL
 *    "Approximate Positions of Major Planets" ( http://ssd.jpl.nasa.gov/?planet_pos ) to retrieve
 *    initial ephemeris data for a specific planet. The ephemeris is valid for circular, coplanar
 *    orbits.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : L. van der Ham
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.vanderHam@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 24 February, 2011
 *    Last modified     : 3 August, 2011
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the
 *          Major Planets, http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf,
 *          last accessed: 24 February, 2011.
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110629    L. van der Ham    First creation of code.
 *      110803    L. van der Ham    Seperated this code from approximatePlanetPositions.
 */

// Include statements.
#include <cmath>
#include "Astrodynamics/Bodies/planet.h"
#include "Astrodynamics/States/approximatePlanetPositionsCircularCoplanar.h"

//! Tudat library namespace.
namespace tudat
{

// Using declarations.
using std::cerr;
using std::endl;
using std::sin;
using std::cos;

//! Get state from ephemeris; circular, coplanar case
CartesianElements* ApproximatePlanetPositionsCircularCoplanar::getStateFromEphemeris(
        double julianDate )
{
    // Set Julian date.
    julianDate_ = julianDate;

    // Compute number of centuries past J2000.
    numberOfCenturiesPastJ2000_ = ( julianDate_ - 2451545.0 ) / 36525.0;

    // Compute mean longitude of planet at given Julian date.
    meanLongitudeAtGivenJulianDate_ = approximatePlanetPositionsDataContainer_.meanLongitude_
            + ( approximatePlanetPositionsDataContainer_.rateOfChangeOfMeanLongitude_
                * numberOfCenturiesPastJ2000_ );

    // Convert mean longitude at given Julian date from degrees to radians.
    meanLongitudeAtGivenJulianDate_ = unit_conversions::convertDegreesToRadians(
                meanLongitudeAtGivenJulianDate_);

    // Get semi-major axis at J2000 and assume constant radius of circular orbit.
    constantOrbitalRadius_ = unit_conversions::convertAstronomicalUnitsToMeters(
                approximatePlanetPositionsDataContainer_.semiMajorAxis_ );

    // Convert to Cartesian position.
    Eigen::VectorXd planetCartesianPositionAtGivenJulianDateX_( 3 );
    mathematics::convertSphericalToCartesian( constantOrbitalRadius_,
                                              meanLongitudeAtGivenJulianDate_, 0.5 * M_PI,
                                              planetCartesianPositionAtGivenJulianDateX_ );
    Eigen::Vector3d planetCartesianPositionAtGivenJulianDate_ =
            planetCartesianPositionAtGivenJulianDateX_;

    // Create predefined Sun.
    Planet predefinedSun_;
    predefinedSun_.setPredefinedPlanetSettings( Planet::sun );

    // Compute orbital velocity.
    double circularOrbitalVelocity = std::sqrt( predefinedSun_.getGravitationalParameter( ) /
                                                constantOrbitalRadius_ );

    // Convert to Cartesian velocity.
    Eigen::Vector3d planetCartesianVelocityAtGivenJulianDate_;
    planetCartesianVelocityAtGivenJulianDate_( 0 ) = -sin( meanLongitudeAtGivenJulianDate_ ) *
            circularOrbitalVelocity;
    planetCartesianVelocityAtGivenJulianDate_( 1 ) = cos( meanLongitudeAtGivenJulianDate_ ) *
            circularOrbitalVelocity;
    planetCartesianVelocityAtGivenJulianDate_( 2 ) = 0.0;

    // Set Cartesian state elements.
    planetCartesianElementsAtGivenJulianDate_.setPosition(
                planetCartesianPositionAtGivenJulianDate_ );
    planetCartesianElementsAtGivenJulianDate_.setVelocity(
                planetCartesianVelocityAtGivenJulianDate_ );

    // Return Cartesian state of planet at given Julian date.
    return &planetCartesianElementsAtGivenJulianDate_;
}

}

// End of file.

