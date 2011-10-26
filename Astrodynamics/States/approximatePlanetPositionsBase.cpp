/*! \file approximatePlanetPositionsBase.cpp
 *    This source file contains the definition of an ephemeris class that makes use of the JPL
 *    "Approximate Positions of Major Planets" ( http://ssd.jpl.nasa.gov/?planet_pos ) to retrieve
 *    ephemeris data for a specific planet. The ephemeris file used is for the period 3000 BC to
 *    3000 AD.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : L. van der Ham
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.vanderHam@student.tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
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
 *      110221    K. Kumar          First creation of code.
 *      110224    K. Kumar          Renamed class and file.
 *      110629    L. van der Ham    Added circular coplanar case.
 *      110803    L. van der Ham    Created base class and seperated approximatePlanetPositions
 *                                  from approximatePlanetPositionsCircularCoplanar.
 */

// Include statements.
#include "Astrodynamics/Bodies/planet.h"
#include "Astrodynamics/States/approximatePlanetPositionsBase.h"

//! Tudat library namespace.
namespace tudat
{

//! Default constructor.
ApproximatePlanetPositionsBase::ApproximatePlanetPositionsBase( ) : julianDate_( -0.0 ),
    meanLongitudeAtGivenJulianDate_( -0.0 ), numberOfCenturiesPastJ2000_( -0.0 )
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

//! Set planet.
void ApproximatePlanetPositionsBase::setPlanet( BodiesWithEphemerisData bodyWithEphemerisData )
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

//! Parse ephemeris line data.
void ApproximatePlanetPositionsBase::parseEphemerisLineData_( const unsigned int& firstLineNumber )
{
    // Parse data from container of strings using a string stream.

    // Clear stringstream.
    ephemerisLineData_.clear( );

    // Read first line of data.
    ephemerisLineData_ << containerOfDataFromEphemerisFile_[ firstLineNumber ];

    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.planetName_;

    // Check if the line number corresponds to that for "EM Bary".
    if ( firstLineNumber == 22 )
    {
        string earthMoonBarycenter_;

        ephemerisLineData_ >> earthMoonBarycenter_;
        approximatePlanetPositionsDataContainer_.planetName_ += " " + earthMoonBarycenter_;
    }

    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.semiMajorAxis_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.eccentricity_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.inclination_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.meanLongitude_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.longitudeOfPerihelion_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.longitudeOfAscendingNode_;

    // Clear stringstream.
    ephemerisLineData_.clear( );

    // Read second line of data.
    ephemerisLineData_ << containerOfDataFromEphemerisFile_[ firstLineNumber + 1 ];

    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.rateOfChangeOfSemiMajorAxis_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.rateOfChangeOfEccentricity_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.rateOfChangeOfInclination_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.rateOfChangeOfMeanLongitude_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_
                          .rateOfChangeOfLongitudeOfPerihelion_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_
                          .rateOfChangeOfLongitudeOfAscendingNode_;
}

//! Parse line data for extra terms for ephemeris.
void ApproximatePlanetPositionsBase::parseExtraTermsEphemerisLineData_(
    const unsigned int& lineNumber )
{
    // Clear stringstream.
    ephemerisLineData_.clear( );

    // Read second line of data.
    ephemerisLineData_<< containerOfDataFromEphemerisFile_[ lineNumber ];

    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.planetName_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.additionalTermB_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.additionalTermC_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.additionalTermS_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.additionalTermF_;
}

}

// End of file.
