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
 *      http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 */

#include <cmath>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"

namespace tudat
{
namespace ephemerides
{

//! Get state from ephemeris.
Eigen::Vector6d ApproximateJplEphemeris::getCartesianState(
        const double secondsSinceEpoch)
{
    // Convert planet elements in Keplerian elements to Cartesian elements.
    return orbital_element_conversions::
            convertKeplerianToCartesianElements(
                getKeplerianStateFromEphemeris( secondsSinceEpoch ),
                sunGravitationalParameter_ + planetGravitationalParameter_ );
}

//! Get keplerian state from ephemeris.
Eigen::Vector6d ApproximateJplEphemeris::getKeplerianStateFromEphemeris(
        const double secondsSinceEpoch )
{
    using std::pow;
    using std::sin;
    using std::cos;
    using namespace orbital_element_conversions;

    // Set Julian date.
    julianDate_ = basic_astrodynamics::convertSecondsSinceEpochToJulianDay(
                secondsSinceEpoch, basic_astrodynamics::JULIAN_DAY_ON_J2000 );

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
    // Convert semi-major axis from AU to meters
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

ApproximateGtopEphemeris::ApproximateGtopEphemeris( const std::string& bodyName ):
    Ephemeris( "Sun", "ECLIPJ2000" )
{
    bodyIndex_ = getPlanetIndex( bodyName );
}


Eigen::Vector6d getGtopCartesianElements ( const double daysSinceMjd2000,
                                           const int planet )
{
    const double gtopAstronomicalUnit = 149597870660.0; // Astronomical Unit
    const double  gtopSunGravitationalParameter = 1.327124280000000e+020;  //Gravitational constant of Sun;

    double T = daysSinceMjd2000 / ( 100.0 * physical_constants::JULIAN_YEAR_IN_DAYS ) + 1.0;
    Eigen::Vector6d keplerElements;
    double meanMotion;
    switch (planet)
    {
        case(1):// Mercury
            keplerElements( 0 ) = (0.38709860);
            keplerElements( 1 ) = (0.205614210 + 0.000020460 * T - 0.000000030 * T * T);
            keplerElements( 2 ) = (7.002880555555555560 + 1.86083333333333333e-3 * T - 1.83333333333333333e-5 * T * T);
            keplerElements( 4 ) = (4.71459444444444444e+1 + 1.185208333333333330 * T + 1.73888888888888889e-4 * T * T);
            keplerElements( 3 ) = (2.87537527777777778e+1 + 3.70280555555555556e-1 * T +1.20833333333333333e-4 * T * T);
            meanMotion = 1.49472515288888889e+5 + 6.38888888888888889e-6 * T;
            keplerElements( 5 ) = (1.02279380555555556e2 + meanMotion * T);
        break;
        case(2):// Venus
            keplerElements( 0 ) = (0.72333160);
            keplerElements( 1 ) = (0.006820690 - 0.000047740 * T + 0.0000000910 * T * T);
            keplerElements( 2 ) = (3.393630555555555560 + 1.00583333333333333e-3 * T - 9.72222222222222222e-7 * T * T);
            keplerElements( 4 ) = (7.57796472222222222e+1 + 8.9985e-1 * T + 4.1e-4 * T * T);
            keplerElements( 3 ) = (5.43841861111111111e+1 + 5.08186111111111111e-1 * T -1.38638888888888889e-3 * T * T);
            meanMotion = 5.8517803875e+4 + 1.28605555555555556e-3 * T;
            keplerElements( 5 ) = (2.12603219444444444e2 + meanMotion * T);
        break;
        case(3):// Earth
            keplerElements( 0 ) = (1.000000230);
            keplerElements( 1 ) = (0.016751040 - 0.000041800 * T - 0.0000001260 * T * T);
            keplerElements( 2 ) = (0.00);
            keplerElements( 4 ) = (0.00);
            keplerElements( 3 ) = (1.01220833333333333e+2 + 1.7191750 * T + 4.52777777777777778e-4 * T * T + 3.33333333333333333e-6 * T * T * T);
             meanMotion = 3.599904975e+4 - 1.50277777777777778e-4 * T - 3.33333333333333333e-6 * T * T;
            keplerElements( 5 ) = (3.58475844444444444e2 + meanMotion * T);
        break;
        case(4):// Mars
            keplerElements( 0 ) = (1.5236883990);
            keplerElements( 1 ) = (0.093312900 + 0.0000920640 * T - 0.0000000770 * T * T);
            keplerElements( 2 ) = (1.850333333333333330 - 6.75e-4 * T + 1.26111111111111111e-5 * T * T);
            keplerElements( 4 ) = (4.87864416666666667e+1 + 7.70991666666666667e-1 * T - 1.38888888888888889e-6 * T * T - 5.33333333333333333e-6 * T * T * T);
            keplerElements( 3 ) = (2.85431761111111111e+2 + 1.069766666666666670 * T +  1.3125e-4 * T * T + 4.13888888888888889e-6 * T * T * T);
            meanMotion = 1.91398585e+4 + 1.80805555555555556e-4 * T + 1.19444444444444444e-6 * T * T;
            keplerElements( 5 ) = (3.19529425e2 + meanMotion * T);
        break;
        case(5):// Jupiter
            keplerElements( 0 ) = (5.2025610);
            keplerElements( 1 ) = (0.048334750 + 0.000164180 * T  - 0.00000046760 * T * T -0.00000000170 * T * T * T);
            keplerElements( 2 ) = (1.308736111111111110 - 5.69611111111111111e-3 * T +  3.88888888888888889e-6 * T * T);
            keplerElements( 4 ) = (9.94433861111111111e+1 + 1.010530 * T + 3.52222222222222222e-4 * T * T - 8.51111111111111111e-6 * T * T * T);
            keplerElements( 3 ) = (2.73277541666666667e+2 + 5.99431666666666667e-1 * T + 7.0405e-4 * T * T + 5.07777777777777778e-6 * T * T * T);
            meanMotion = 3.03469202388888889e+3 - 7.21588888888888889e-4 * T + 1.78444444444444444e-6 * T * T;
            keplerElements( 5 ) = (2.25328327777777778e2 + meanMotion * T);
        break;
        case(6):// Saturn
            keplerElements( 0 ) = (9.5547470);
            keplerElements( 1 ) = (0.055892320 - 0.00034550 * T - 0.0000007280 * T * T + 0.000000000740 * T * T * T);
            keplerElements( 2 ) = (2.492519444444444440 - 3.91888888888888889e-3 * T - 1.54888888888888889e-5 * T * T + 4.44444444444444444e-8 * T * T * T);
            keplerElements( 4 ) = (1.12790388888888889e+2 + 8.73195138888888889e-1 * T -1.52180555555555556e-4 * T * T - 5.30555555555555556e-6 * T * T * T);
            keplerElements( 3 ) = (3.38307772222222222e+2 + 1.085220694444444440 * T + 9.78541666666666667e-4 * T * T + 9.91666666666666667e-6 * T * T * T);
            meanMotion = 1.22155146777777778e+3 - 5.01819444444444444e-4 * T - 5.19444444444444444e-6 * T * T;
            keplerElements( 5 ) = (1.75466216666666667e2 + meanMotion * T);
        break;
        case(7):// Uranus
            keplerElements( 0 ) = (19.218140);
            keplerElements( 1 ) = (0.04634440 - 0.000026580 * T + 0.0000000770 * T * T);
            keplerElements( 2 ) = (7.72463888888888889e-1 + 6.25277777777777778e-4 * T + 3.95e-5 * T * T);
            keplerElements( 4 ) = (7.34770972222222222e+1 + 4.98667777777777778e-1 * T + 1.31166666666666667e-3 * T * T);
            keplerElements( 3 ) = (9.80715527777777778e+1 + 9.85765e-1 * T - 1.07447222222222222e-3 * T * T - 6.05555555555555556e-7 * T * T * T);
            meanMotion = 4.28379113055555556e+2 + 7.88444444444444444e-5 * T + 1.11111111111111111e-9 * T * T;
            keplerElements( 5 ) = (7.26488194444444444e1 + meanMotion * T);
        break;
        case(8)://Neptune
            keplerElements( 0 ) = (30.109570);
            keplerElements( 1 ) = (0.008997040 + 0.0000063300 * T - 0.0000000020 * T * T);
            keplerElements( 2 ) = (1.779241666666666670 - 9.54361111111111111e-3 * T - 9.11111111111111111e-6 * T * T);
            keplerElements( 4 ) = (1.30681358333333333e+2 + 1.0989350 * T + 2.49866666666666667e-4 * T * T - 4.71777777777777778e-6 * T * T * T);
            keplerElements( 3 ) = (2.76045966666666667e+2 + 3.25639444444444444e-1 * T + 1.4095e-4 * T * T + 4.11333333333333333e-6 * T * T * T);
            meanMotion = 2.18461339722222222e+2 - 7.03333333333333333e-5 * T;
            keplerElements( 5 ) = (3.77306694444444444e1 + meanMotion * T);
        break;
        case(9):// Pluto
            //Fifth order polynomial least square fit generated by Dario Izzo
            //(ESA ACT). JPL405 ephemerides (Charon-Pluto barycenter) have been used to produce the coefficients.
            //This approximation should not be used outside the range 2000-2100;
            T = daysSinceMjd2000 / 36525.00;
            keplerElements( 0 ) = (39.34041961252520 + 4.33305138120726 * T - 22.93749932403733 * T * T + 48.76336720791873 * T * T * T - 45.52494862462379 * T * T * T * T + 15.55134951783384 * T * T * T * T * T);
            keplerElements( 1 ) = (0.24617365396517 + 0.09198001742190 * T - 0.57262288991447 * T * T + 1.39163022881098 * T * T * T - 1.46948451587683 * T * T * T * T + 0.56164158721620 * T * T * T * T * T);
            keplerElements( 2 ) = (17.16690003784702 - 0.49770248790479 * T + 2.73751901890829 * T * T - 6.26973695197547 * T * T * T + 6.36276927397430 * T * T * T * T - 2.37006911673031 * T * T * T * T * T);
            keplerElements( 4 ) = (110.222019291707 + 1.551579150048 * T - 9.701771291171 * T * T + 25.730756810615 * T * T * T - 30.140401383522 * T * T * T * T + 12.796598193159 * T * T * T * T * T);
            keplerElements( 3 ) = (113.368933916592 + 9.436835192183 * T - 35.762300003726 * T * T + 48.966118351549 * T * T * T - 19.384576636609 * T * T * T * T - 3.362714022614 * T * T * T * T * T);
            keplerElements( 5 ) = (15.17008631634665 + 137.023166578486 * T + 28.362805871736 * T * T - 29.677368415909 * T * T * T - 3.585159909117 * T * T * T * T + 13.406844652829 * T * T * T * T * T);
        break;

    }


    keplerElements( 0 ) *= gtopAstronomicalUnit;

    // conversion of DEG into RAD
    double degreeToRadians = unit_conversions::convertDegreesToRadians( 1.0 );
    keplerElements( 2 ) *= degreeToRadians;
    keplerElements( 3 ) *= degreeToRadians;
    keplerElements( 4 ) *= degreeToRadians;
    keplerElements( 5 ) *= degreeToRadians;
    keplerElements( 5 ) =  fmod( keplerElements( 5 ), 2.0 * mathematical_constants::PI );

    // Conversion from Mean Anomaly to Eccentric Anomaly via Kepler's equation
    keplerElements( 5 ) = tudat::orbital_element_conversions::convertMeanAnomalyToTrueAnomaly(keplerElements( 1 ),keplerElements( 5 ) );
    return tudat::orbital_element_conversions::convertKeplerianToCartesianElements(
                keplerElements,  gtopSunGravitationalParameter );


}
} // namespace ephemerides
} // namespace tudat
