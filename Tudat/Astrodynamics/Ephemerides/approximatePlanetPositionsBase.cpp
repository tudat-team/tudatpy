/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      110629    L. van der Ham    Added circular coplanar case.
 *      110803    L. van der Ham    Created base class and seperated approximatePlanetPositions
 *                                  from approximatePlanetPositionsCircularCoplanar.
 *      120322    D. Dirkx          Modified to new Ephemeris interfaces.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 *    Notes
 *
 */

#include <iostream>

#include <boost/exception/all.hpp>
#include <boost/format.hpp>
#include <boost/throw_exception.hpp>

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositionsBase.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace ephemerides
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

    // Read the file into a container.
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

    // Close file.
    ephemerisFile_.close( );
}

} // namespace ephemerides
} // namespace tudat
