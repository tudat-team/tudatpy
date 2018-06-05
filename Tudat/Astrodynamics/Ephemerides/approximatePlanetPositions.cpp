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
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *      http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 */

#include <cmath>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"

namespace tudat
{
namespace ephemerides
{

//! Get state from ephemeris.
Eigen::Vector6d ApproximatePlanetPositions::getCartesianState(
        const double secondsSinceEpoch)
{
    // Convert planet elements in Keplerian elements to Cartesian elements.
    return orbital_element_conversions::
            convertKeplerianToCartesianElements(
                getKeplerianStateFromEphemeris( secondsSinceEpoch ),
                sunGravitationalParameter );
}

//! Get keplerian state from ephemeris.
Eigen::Vector6d ApproximatePlanetPositions::getKeplerianStateFromEphemeris(
        const double secondsSinceEpoch )
{
    using std::pow;
    using std::sin;
    using std::cos;
    using namespace orbital_element_conversions;

    // Set Julian date.
    julianDate_ = basic_astrodynamics::convertSecondsSinceEpochToJulianDay(
                secondsSinceEpoch, referenceJulianDate_ );

    // Compute number of centuries past J2000.
    numberOfCenturiesPastJ2000_ = ( julianDate_ - 2451545.0 ) / 36525.0;

    // Compute and set semi-major axis of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_( semiMajorAxisIndex )
            = approximatePlanetPositionsDataContainer_.semiMajorAxis_
            + ( approximatePlanetPositionsDataContainer_
                .rateOfChangeOfSemiMajorAxis_ * numberOfCenturiesPastJ2000_ );

    // Compute and set eccentricity of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_( eccentricityIndex )
            = approximatePlanetPositionsDataContainer_.eccentricity_
            + ( approximatePlanetPositionsDataContainer_
                .rateOfChangeOfEccentricity_ * numberOfCenturiesPastJ2000_ );

    // Compute and set inclination of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_( inclinationIndex )
            = approximatePlanetPositionsDataContainer_.inclination_
            + ( approximatePlanetPositionsDataContainer_
                .rateOfChangeOfInclination_ * numberOfCenturiesPastJ2000_ );

    // Compute and set longitude of ascending node of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_( longitudeOfAscendingNodeIndex )
            = approximatePlanetPositionsDataContainer_.longitudeOfAscendingNode_
            + ( approximatePlanetPositionsDataContainer_
                .rateOfChangeOfLongitudeOfAscendingNode_ * numberOfCenturiesPastJ2000_ );

    // Compute longitude of perihelion of planet at given Julian date.
    longitudeOfPerihelionAtGivenJulianDate_
            = approximatePlanetPositionsDataContainer_.longitudeOfPerihelion_
            + ( approximatePlanetPositionsDataContainer_.rateOfChangeOfLongitudeOfPerihelion_
                * numberOfCenturiesPastJ2000_ );

    // Compute and set argument of periapsis of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_( argumentOfPeriapsisIndex )
            = longitudeOfPerihelionAtGivenJulianDate_
            - planetKeplerianElementsAtGivenJulianDate_( longitudeOfAscendingNodeIndex );

    // Compute mean longitude of planet at given Julian date.
    meanLongitudeAtGivenJulianDate_ = approximatePlanetPositionsDataContainer_.meanLongitude_
            + ( approximatePlanetPositionsDataContainer_.rateOfChangeOfMeanLongitude_
                * numberOfCenturiesPastJ2000_ );

    // Compute mean anomaly of planet at given Julian date.
    meanAnomalyAtGivenJulianDate_ = meanLongitudeAtGivenJulianDate_
            - longitudeOfPerihelionAtGivenJulianDate_
            + ( approximatePlanetPositionsDataContainer_.additionalTermB_
                * pow( numberOfCenturiesPastJ2000_, 2.0 ) )
            + ( approximatePlanetPositionsDataContainer_.additionalTermC_
                * cos( unit_conversions::convertDegreesToRadians(approximatePlanetPositionsDataContainer_.additionalTermF_ *
                       numberOfCenturiesPastJ2000_ )) )
            + ( approximatePlanetPositionsDataContainer_.additionalTermS_
                * sin( unit_conversions::convertDegreesToRadians(approximatePlanetPositionsDataContainer_.additionalTermF_ *
                       numberOfCenturiesPastJ2000_) ) );

    // Compute modulo of mean anomaly for interval :
    // 0 <= meanAnomalyAtGivenJulianDate_ < 360.
    meanAnomalyAtGivenJulianDate_ = basic_mathematics::computeModulo(
                meanAnomalyAtGivenJulianDate_, 360.0 );

    // Translate mean anomaly to:
    // -180 < meanAnomalyAtGivenJulianDate_ <= 180 bounds.
    if ( meanAnomalyAtGivenJulianDate_ > 180.0 )
    {
        meanAnomalyAtGivenJulianDate_ -= 360.0;
    }

    // Convert mean anomaly to eccentric anomaly.
    eccentricAnomalyAtGivenJulianDate_ = convertMeanAnomalyToEccentricAnomaly(
                planetKeplerianElementsAtGivenJulianDate_( eccentricityIndex ),
                unit_conversions::convertDegreesToRadians(
                    meanAnomalyAtGivenJulianDate_ ) );

    // Convert eccentric anomaly to true anomaly and set in planet elements.
    trueAnomalyAtGivenJulianData_
            = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                eccentricAnomalyAtGivenJulianDate_,
                planetKeplerianElementsAtGivenJulianDate_( eccentricityIndex ) );

    planetKeplerianElementsAtGivenJulianDate_( trueAnomalyIndex )
            = trueAnomalyAtGivenJulianData_;

    // Convert Keplerian elements to standard units.
    // Convert semi-major axis from AU to meters.
    planetKeplerianElementsAtGivenJulianDate_( semiMajorAxisIndex )
            = unit_conversions::convertAstronomicalUnitsToMeters(
                planetKeplerianElementsAtGivenJulianDate_( semiMajorAxisIndex ) );

    // Convert inclination from degrees to radians.
    planetKeplerianElementsAtGivenJulianDate_( inclinationIndex )
            = unit_conversions::convertDegreesToRadians(
                planetKeplerianElementsAtGivenJulianDate_( inclinationIndex ) );

    // Convert longitude of ascending node from degrees to radians.
    planetKeplerianElementsAtGivenJulianDate_( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians(
                planetKeplerianElementsAtGivenJulianDate_( longitudeOfAscendingNodeIndex ) );

    // Convert argument of periapsis from degrees to radians.
    planetKeplerianElementsAtGivenJulianDate_( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians(
                planetKeplerianElementsAtGivenJulianDate_( argumentOfPeriapsisIndex ) );

    return planetKeplerianElementsAtGivenJulianDate_;
}

} // namespace ephemerides
} // namespace tudat
