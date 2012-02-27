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
 *      110221    K. Kumar          First creation of code.
 *      110224    K. Kumar          Renamed class and file.
 *      120217    K. Kumar          Updated computeModuloForSignedValues() to computeModulo()
 *                                  from Tudat Core.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the
 *          Major Planets, http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf,
 *          last accessed: 24 February, 2011.
 *
 */

#include <cmath>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include "Tudat/Astrodynamics/Bodies/Ephemeris/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Bodies/planet.h"

// Using declarations.
using std::cerr;
using std::endl;
using std::pow;
using std::sin;
using std::cos;

namespace tudat
{

//! Get state from ephemeris.
CartesianElements* ApproximatePlanetPositions::
        getStateFromEphemeris( double julianDate )
{
    // Set Julian date.
    julianDate_ = julianDate;

    // Compute number of centuries past J2000.
    numberOfCenturiesPastJ2000_ = ( julianDate - 2451545.0 ) / 36525.0;

    // Compute and set semi-major axis of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_
            .setSemiMajorAxis( approximatePlanetPositionsDataContainer_.semiMajorAxis_
                               + ( approximatePlanetPositionsDataContainer_
                                   .rateOfChangeOfSemiMajorAxis_ * numberOfCenturiesPastJ2000_ ) );

    // Compute and set eccentricity of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_
            .setEccentricity( approximatePlanetPositionsDataContainer_.eccentricity_
                              + ( approximatePlanetPositionsDataContainer_
                                  .rateOfChangeOfEccentricity_ * numberOfCenturiesPastJ2000_ ) );

    // Compute and set inclination of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_
            .setInclination( approximatePlanetPositionsDataContainer_.inclination_
                             + ( approximatePlanetPositionsDataContainer_
                                 .rateOfChangeOfInclination_ * numberOfCenturiesPastJ2000_ ) );

    // Compute and set longitude of ascending node of planet at given
    // Julian date.
    planetKeplerianElementsAtGivenJulianDate_
            .setLongitudeOfAscendingNode(
                approximatePlanetPositionsDataContainer_.longitudeOfAscendingNode_
                + ( approximatePlanetPositionsDataContainer_
                    .rateOfChangeOfLongitudeOfAscendingNode_ * numberOfCenturiesPastJ2000_ ) );

    // Compute longitude of perihelion of planet at given Julian date.
    longitudeOfPerihelionAtGivenJulianDate_
            = approximatePlanetPositionsDataContainer_.longitudeOfPerihelion_
            + ( approximatePlanetPositionsDataContainer_.rateOfChangeOfLongitudeOfPerihelion_
                * numberOfCenturiesPastJ2000_ );

    // Compute and set argument of periapsis of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_.setArgumentOfPeriapsis(
                longitudeOfPerihelionAtGivenJulianDate_
                - planetKeplerianElementsAtGivenJulianDate_.getLongitudeOfAscendingNode( ) );

    // Compute mean longitude of planet at given Julian date.
    meanLongitudeAtGivenJulianDate_ = approximatePlanetPositionsDataContainer_.meanLongitude_
            + ( approximatePlanetPositionsDataContainer_.rateOfChangeOfMeanLongitude_
                * numberOfCenturiesPastJ2000_ );

    // Compute mean anomaly of planet at given Julian date.
    meanAnomalyAtGivenJulianDate_ = meanLongitudeAtGivenJulianDate_
            - longitudeOfPerihelionAtGivenJulianDate_
            + ( approximatePlanetPositionsDataContainer_.additionalTermB_
                * pow( julianDate, 2.0 ) )
            + ( approximatePlanetPositionsDataContainer_.additionalTermC_
                * cos( approximatePlanetPositionsDataContainer_.additionalTermF_ * julianDate ) )
            + ( approximatePlanetPositionsDataContainer_.additionalTermS_
                * sin( approximatePlanetPositionsDataContainer_.additionalTermF_ * julianDate ) );

    // Compute modulo of mean anomaly for interval :
    // 0 <= meanAnomalyAtGivenJulianDate_ < 360.
    meanAnomalyAtGivenJulianDate_ = mathematics::computeModulo(
                meanAnomalyAtGivenJulianDate_, 360.0 );

    // Translate mean anomaly to:
    // -180 < meanAnomalyAtGivenJulianDate_ <= 180 bounds.
    if ( meanAnomalyAtGivenJulianDate_ > 180.0 )
    {
        meanAnomalyAtGivenJulianDate_ -= 360.0;
    }

    // Set eccentricty for mean anomaly to eccentric anomaly conversion.
    convertMeanAnomalyToEccentricAnomaly_.setEccentricity(
                planetKeplerianElementsAtGivenJulianDate_.getEccentricity( ) );

    // Set mean anomaly for conversion to eccentric anomaly.
    convertMeanAnomalyToEccentricAnomaly_.setMeanAnomaly(
                unit_conversions::convertDegreesToRadians( meanAnomalyAtGivenJulianDate_ ) );

    // Set Newton-Raphson method to use for mean anomaly to eccentric
    // anomaly conversion.
    convertMeanAnomalyToEccentricAnomaly_.setNewtonRaphson( &newtonRaphson_ );

    eccentricAnomalyAtGivenJulianDate_ = convertMeanAnomalyToEccentricAnomaly_.convert( );

    // Convert eccentric anomaly to true anomaly and set in planet elements.
    trueAnomalyAtGivenJulianData_
            = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                eccentricAnomalyAtGivenJulianDate_,
                planetKeplerianElementsAtGivenJulianDate_.getEccentricity( ) );

    planetKeplerianElementsAtGivenJulianDate_.setTrueAnomaly(
                trueAnomalyAtGivenJulianData_ );

    // Convert Keplerian elements to standard units.
    // Convert semi-major axis from AU to meters.
    planetKeplerianElementsAtGivenJulianDate_.setSemiMajorAxis(
                unit_conversions::convertAstronomicalUnitsToMeters(
                    planetKeplerianElementsAtGivenJulianDate_.getSemiMajorAxis( ) ) );

    // Convert inclination from degrees to radians.
    planetKeplerianElementsAtGivenJulianDate_.setInclination(
                unit_conversions::convertDegreesToRadians(
                    planetKeplerianElementsAtGivenJulianDate_.getInclination( ) ) );

    // Convert longitude of ascending node from degrees to radians.
    planetKeplerianElementsAtGivenJulianDate_
            .setLongitudeOfAscendingNode(
                unit_conversions::convertDegreesToRadians(
                    planetKeplerianElementsAtGivenJulianDate_.getLongitudeOfAscendingNode( ) ) );

    // Convert argument of periapsis from degrees to radians.
    planetKeplerianElementsAtGivenJulianDate_.setArgumentOfPeriapsis(
                unit_conversions::convertDegreesToRadians(
                    planetKeplerianElementsAtGivenJulianDate_.getArgumentOfPeriapsis( ) ) );

    // Create predefined Sun.
    Planet predefinedSun_;
    predefinedSun_.setPredefinedPlanetSettings( Planet::sun );

    // Convert planet Elements in Keplerian elements to Cartesian elements.
    planetCartesianElementsAtGivenJulianDate_.state
            = orbital_element_conversions::convertKeplerianToCartesianElements(
                planetKeplerianElementsAtGivenJulianDate_.state,
                predefinedSun_.getGravitationalParameter( ) );

    // Return Cartesian elements of planet at given Julian date.
    return &planetCartesianElementsAtGivenJulianDate_;
}

} // namespace tudat
