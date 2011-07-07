/*! \file approximatePlanetPositions.cpp
 *    This source file contains the definition of an ephemeris class that makes
 *    use of the JPL "Approximate Positions of Major Planets"
 *    ( http://ssd.jpl.nasa.gov/?planet_pos ) to retrieve ephemeris data for a
 *    specific planet. The ephemeris file used is for the period 3000 BC to
 *    3000 AD.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 24 February, 2011
 *    Last modified     : 24 February, 2011
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the
 *          Major Planets, http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf,
 *          last accessed: 24 February, 2011.
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 */

// Include statements.
#include "approximatePlanetPositions.h"
#include "planet.h"

// Using declarations.
using std::cerr;
using std::endl;
using mathematics::raiseToIntegerPower;

//! Default constructor.
ApproximatePlanetPositions::ApproximatePlanetPositions( )
{
    // Set relative path to ephemeris file in file reader.
    ephemerisTextFileReader_.setRelativePath( "External/EphemerisData/" );

    // Set file name of ephemeris data file.
    ephemerisTextFileReader_.setFileName( "p_elem_t2.txt" );

    // Open ephemeris file.
    ephemerisTextFileReader_.openFile( );

    // Skip first 17 lines of file.
    ephemerisTextFileReader_.skipLines( 17 );

    // Read and store next 18 lines of file.
    ephemerisTextFileReader_.readAndStoreData( 18 );

    // Skip next 12 lines of file.
    ephemerisTextFileReader_.skipLines( 12 );

    // Read and store next 5 lines of file.
    ephemerisTextFileReader_.readAndStoreData( 5 );

    // Close file.
    ephemerisTextFileReader_.closeFile( );

    // Get container of data from file.
    containerOfDataFromEphemerisFile_ = ephemerisTextFileReader_.getContainerOfData( );
}

//! Default destructor.
ApproximatePlanetPositions::
        ~ApproximatePlanetPositions( )
{
}

//! Set planet.
void ApproximatePlanetPositions::
        setPlanet( BodiesWithEphemerisData bodyWithEphemerisData )
{
    bodyWithEphemerisData_ = bodyWithEphemerisData;

    switch ( bodyWithEphemerisData )
    {
    case mercury:

        parseEphemerisLineData_( 18 );
        break;

    case venus:

        parseEphemerisLineData_( 20 );
        break;

    case earthMoonBarycenter:

        parseEphemerisLineData_( 22 );
        break;

    case mars:

        parseEphemerisLineData_( 24 );
        break;

    case jupiter:

        parseEphemerisLineData_( 26 );
        parseExtraTermsEphemerisLineData_( 48 );
        break;

    case saturn:

        parseEphemerisLineData_( 28 );
        parseExtraTermsEphemerisLineData_( 49 );
        break;

    case uranus:

        parseEphemerisLineData_( 30 );
        parseExtraTermsEphemerisLineData_( 50 );
        break;

    case neptune:

        parseEphemerisLineData_( 32 );
        parseExtraTermsEphemerisLineData_( 51 );
        break;

    case pluto:

        parseEphemerisLineData_( 34 );
        parseExtraTermsEphemerisLineData_( 52 );
        break;

    default:

        break;
    }
}

//! Get state from ephemeris.
CartesianElements* ApproximatePlanetPositions::
        getStateFromEphemeris( const double& julianDate )
{
    // Set Julian date.
    julianDate_ = julianDate;

    // Compute number of centuries past J2000.
    numberOfCenturiesPastJ2000_ = ( julianDate - 2451545.0 ) / 36525.0;

    // Compute and set semi-major axis of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_
            .setSemiMajorAxis( approximatePlanetPositionsDataContainer_
                               .semiMajorAxis_
                               + ( approximatePlanetPositionsDataContainer_
                                   .rateOfChangeOfSemiMajorAxis_
                                   * numberOfCenturiesPastJ2000_ ) );

    // Compute and set eccentricity of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_
            .setEccentricity( approximatePlanetPositionsDataContainer_
                              .eccentricity_
                              + ( approximatePlanetPositionsDataContainer_
                                  .rateOfChangeOfEccentricity_
                                  * numberOfCenturiesPastJ2000_ ) );

    // Compute and set inclination of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_
            .setInclination( approximatePlanetPositionsDataContainer_
                             .inclination_
                             + ( approximatePlanetPositionsDataContainer_
                                 .rateOfChangeOfInclination_
                                 * numberOfCenturiesPastJ2000_ ) );

    // Compute and set longitude of ascending node of planet at given
    // Julian date.
    planetKeplerianElementsAtGivenJulianDate_
            .setLongitudeOfAscendingNode(
                    approximatePlanetPositionsDataContainer_
                    .longitudeOfAscendingNode_
                    + ( approximatePlanetPositionsDataContainer_
                        .rateOfChangeOfLongitudeOfAscendingNode_
                        * numberOfCenturiesPastJ2000_ ) );

    // Compute longitude of perihelion of planet at given Julian date.
    longitudeOfPerihelionAtGivenJulianDate_
            = approximatePlanetPositionsDataContainer_.longitudeOfPerihelion_
              + ( approximatePlanetPositionsDataContainer_
                  .rateOfChangeOfLongitudeOfPerihelion_
                  * numberOfCenturiesPastJ2000_ );

    // Compute and set argument of periapsis of planet at given Julian date.
    planetKeplerianElementsAtGivenJulianDate_
            .setArgumentOfPeriapsis(
                    longitudeOfPerihelionAtGivenJulianDate_
                    - planetKeplerianElementsAtGivenJulianDate_
                    .getLongitudeOfAscendingNode( ) );

    // Compute mean longitude of planet at given Julian date.
    meanLongitudeAtGivenJulianDate_
            = approximatePlanetPositionsDataContainer_.meanLongitude_
              + ( approximatePlanetPositionsDataContainer_
                  .rateOfChangeOfMeanLongitude_
                  * numberOfCenturiesPastJ2000_ );

    // Compute mean anomaly of planet at given Julian date.
    meanAnomalyAtGivenJulianDate_
            = meanLongitudeAtGivenJulianDate_
              - longitudeOfPerihelionAtGivenJulianDate_
              + ( approximatePlanetPositionsDataContainer_.additionalTermB_
                  * raiseToIntegerPower( julianDate, 2 ) )
              + ( approximatePlanetPositionsDataContainer_.additionalTermC_
                  * cos( approximatePlanetPositionsDataContainer_
                         .additionalTermF_ * julianDate ) )
              + ( approximatePlanetPositionsDataContainer_.additionalTermS_
                  * sin( approximatePlanetPositionsDataContainer_
                         .additionalTermF_ * julianDate ) );

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
            unit_conversions::convertDegreesToRadians(
                    meanAnomalyAtGivenJulianDate_ ) );

    // Set Newton-Raphson method to use for mean anomaly to eccentric
    // anomaly conversion.
    convertMeanAnomalyToEccentricAnomaly_.setNewtonRaphson( &newtonRaphson_ );

    eccentricAnomalyAtGivenJulianDate_ = convertMeanAnomalyToEccentricAnomaly_
                                         .convert( );

    // Convert eccentric anomaly to true anomaly and set in planet elements.
    trueAnomalyAtGivenJulianData_
            = orbital_element_conversions::
              convertEccentricAnomalyToTrueAnomaly(
                      eccentricAnomalyAtGivenJulianDate_,
                      planetKeplerianElementsAtGivenJulianDate_
                      .getEccentricity( ) );

    planetKeplerianElementsAtGivenJulianDate_.setTrueAnomaly(
            trueAnomalyAtGivenJulianData_ );

    // Convert Keplerian elements to standard units.
    // Convert semi-major axis from AU to meters.
    planetKeplerianElementsAtGivenJulianDate_.setSemiMajorAxis(
            unit_conversions::convertAstronomicalUnitsToMeters(
                    planetKeplerianElementsAtGivenJulianDate_
                    .getSemiMajorAxis( ) ) );

    // Convert inclination from degrees to radians.
    planetKeplerianElementsAtGivenJulianDate_.setInclination(
            unit_conversions::convertDegreesToRadians(
                    planetKeplerianElementsAtGivenJulianDate_
                    .getInclination( ) ) );

    // Convert longitude of ascending node from degrees to radians.
    planetKeplerianElementsAtGivenJulianDate_
            .setLongitudeOfAscendingNode(
                    unit_conversions::convertDegreesToRadians(
                            planetKeplerianElementsAtGivenJulianDate_
                            .getLongitudeOfAscendingNode( ) ) );

    // Convert argument of periapsis from degrees to radians.
    planetKeplerianElementsAtGivenJulianDate_.setArgumentOfPeriapsis(
            unit_conversions::convertDegreesToRadians(
                    planetKeplerianElementsAtGivenJulianDate_
                    .getArgumentOfPeriapsis( ) ) );

    // Create predefined Sun.
    Planet predefinedSun_;
    predefinedSun_.setPredefinedPlanetSettings( Planet::sun );

    // Convert planet Elements in Keplerian elements to Cartesian elements.
    planetCartesianElementsAtGivenJulianDate_
            = orbital_element_conversions::
              convertKeplerianToCartesianElements(
                      &planetKeplerianElementsAtGivenJulianDate_,
                      &predefinedSun_ );

    // Return Cartesian elements of planet at given Julian date.
    return &planetCartesianElementsAtGivenJulianDate_;
}

//! Parse ephemeris line data.
void ApproximatePlanetPositions::parseEphemerisLineData_(
        const unsigned int& firstLineNumber )
{
    // Parse data from container of strings using a string stream.

    // Clear stringstream.
    ephemerisLineData_.clear( );

    // Read first line of data.
    ephemerisLineData_
            << containerOfDataFromEphemerisFile_[ firstLineNumber ];

    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_.planetName_;

    // Check if the line number corresponds to that for "EM Bary".
    if ( firstLineNumber == 22 )
    {
        string earthMoonBarycenter_;

        ephemerisLineData_ >> earthMoonBarycenter_;
        approximatePlanetPositionsDataContainer_.planetName_
                += " " + earthMoonBarycenter_;
    }

    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_.semiMajorAxis_;
    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_.eccentricity_;
    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_.inclination_;
    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_.meanLongitude_;
    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_
            .longitudeOfPerihelion_;
    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_
            .longitudeOfAscendingNode_;

    // Clear stringstream.
    ephemerisLineData_.clear( );

    // Read second line of data.
    ephemerisLineData_
            << containerOfDataFromEphemerisFile_[ firstLineNumber + 1 ];

    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_
            .rateOfChangeOfSemiMajorAxis_;
    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_
            .rateOfChangeOfEccentricity_;
    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_
            .rateOfChangeOfInclination_;
    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_
            .rateOfChangeOfMeanLongitude_;
    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_
            .rateOfChangeOfLongitudeOfPerihelion_;
    ephemerisLineData_ >>
            approximatePlanetPositionsDataContainer_
            .rateOfChangeOfLongitudeOfAscendingNode_;
}

//! Parse line data for extra terms for ephemeris.
void ApproximatePlanetPositions::parseExtraTermsEphemerisLineData_(
        const unsigned int& lineNumber )
 {
    // Clear stringstream.
    ephemerisLineData_.clear( );

    // Read second line of data.
    ephemerisLineData_
            << containerOfDataFromEphemerisFile_[ lineNumber ];

    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_
                            .planetName_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_
                            .additionalTermB_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_
                            .additionalTermC_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_
                            .additionalTermS_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_
                            .additionalTermF_;
}

// End of file.
