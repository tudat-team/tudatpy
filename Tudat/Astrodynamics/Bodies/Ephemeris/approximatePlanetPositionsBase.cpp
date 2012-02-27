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
 *      110629    L. van der Ham    Added circular coplanar case.
 *      110803    L. van der Ham    Created base class and seperated approximatePlanetPositions
 *                                  from approximatePlanetPositionsCircularCoplanar.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the
 *          Major Planets, http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf,
 *          last accessed: 24 February, 2011.
 *
 */

#include <iostream>
#include <boost/format.hpp>
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>
#include "Tudat/Astrodynamics/Bodies/Ephemeris/approximatePlanetPositionsBase.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{

//! Set planet.
void ApproximatePlanetPositionsBase::setPlanet( BodiesWithEphemerisData bodyWithEphemerisData )
{
    // Check if ephemeris data has been loaded, and reload data if not.
    if ( containerOfDataFromEphemerisFile_.size( ) == 0 )
    {
        reloadData( );
    }

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
        std::string earthMoonBarycenter_;

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

//! Load in ephemeris data for planets.
void ApproximatePlanetPositionsBase::reloadData( )
{

    // Set  path to ephemeris file in file reader.
    std::string filePath_ = input_output::getTudatRootPath( ) +
            "External/EphemerisData/p_elem_t2.txt";

    // Open ephemeris file.
    std::ifstream ephemerisFile_( filePath_.c_str( ) );
    if ( ephemerisFile_.fail( ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            boost::str( boost::format( "Data file '%s' could not be opened." )
                                 % filePath_.c_str( ) ) ) )
            << boost::errinfo_file_name( filePath_.c_str( ) )
            << boost::errinfo_file_open_mode( "std::ios::binary" )
            << boost::errinfo_api_function( "std::ifstream::open" ) );
    }

    // Read the file into a container
    for ( int line = 1; line < 53; line++ )
    {
        std::string lineData;
        getline( ephemerisFile_, lineData );
        containerOfDataFromEphemerisFile_[line] = lineData;
        if ( ephemerisFile_.fail( ) )
        {
            break;
        }
    }

    // Close file
    ephemerisFile_.close( );

}

} // namespace tudat.
