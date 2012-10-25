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
 *      110221    K. Kumar          Creation of code.
 *      110224    K. Kumar          Renamed class and file.
 *      120217    K. Kumar          Updated computeModuloForSignedValues() to computeModulo()
 *                                  from Tudat Core.
 *      120322    D. Dirkx          Modified to new Ephemeris interfaces.
 *      120522    P. Musegaas       Fixed bug for coordinates of outer planets.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the
 *          Major Planets, http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf,
 *          last accessed: 24 February, 2011.
 *
 */

#include <cmath>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"

namespace tudat
{
namespace ephemerides
{

//! Get state from ephemeris.
Eigen::VectorXd ApproximatePlanetPositions::getCartesianStateFromEphemeris(
        const double julianDate )
{
    // Convert planet elements in Keplerian elements to Cartesian elements.
    return orbital_element_conversions::convertKeplerianToCartesianElements(
                getKeplerianStateFromEphemeris( julianDate ),
                sunGravitationalParameter );
}

//! Get keplerian state from ephemeris.
Eigen::VectorXd ApproximatePlanetPositions::getKeplerianStateFromEphemeris(
        const double julianDate )
{
    using std::pow;
    using std::sin;
    using std::cos;
    using namespace basic_astrodynamics;

    // Set Julian date.
    julianDate_ = julianDate;

    // Compute number of centuries past J2000.
    numberOfCenturiesPastJ2000_ = ( julianDate - 2451545.0 ) / 36525.0;

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
                * cos( approximatePlanetPositionsDataContainer_.additionalTermF_ *
                       numberOfCenturiesPastJ2000_ ) )
            + ( approximatePlanetPositionsDataContainer_.additionalTermS_
                * sin( approximatePlanetPositionsDataContainer_.additionalTermF_ *
                       numberOfCenturiesPastJ2000_ ) );

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

    // Set Newton-Raphson method to use for mean anomaly to eccentric anomaly conversion.
    boost::shared_ptr< NewtonRaphson > newtonRaphson_ = boost::make_shared< NewtonRaphson >( );

    // Set eccentricty and mean anomaly for mean anomaly to eccentric anomaly conversion.
    orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly
            convertMeanAnomalyToEccentricAnomaly_(
                planetKeplerianElementsAtGivenJulianDate_( eccentricityIndex ),
                unit_conversions::convertDegreesToRadians( meanAnomalyAtGivenJulianDate_ ),
                newtonRaphson_ );

    // Convert mean anomaly to eccentric anomaly.
    eccentricAnomalyAtGivenJulianDate_ = convertMeanAnomalyToEccentricAnomaly_.convert( );

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
