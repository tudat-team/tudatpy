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
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 */

#include "tudat/astro/ephemerides/approximatePlanetPositionsBase.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{
namespace ephemerides
{

int getPlanetIndex( const std::string& bodyName )
{
    int bodyIndex;
    if( bodyName == "Mercury" )
    {
        bodyIndex = 1;
    }
    else if( bodyName == "Venus" )
    {
        bodyIndex = 2;
    }
    else if( bodyName == "Earth" )
    {
        bodyIndex = 3;
    }
    else if( bodyName == "EMB" )
    {
        bodyIndex = 3;
    }
    else if( bodyName == "EarthMoonBarycenter" )
    {
        bodyIndex = 3;
    }
    else if( bodyName == "Mars" )
    {
        bodyIndex = 4;
    }
    else if( bodyName == "Jupiter" )
    {
        bodyIndex = 5;
    }
    else if( bodyName == "Saturn" )
    {
        bodyIndex = 6;
    }
    else if( bodyName == "Uranus" )
    {
        bodyIndex = 7;
    }
    else if( bodyName == "Neptune" )
    {
        bodyIndex = 8;
    }
    else if( bodyName == "Pluto" )
    {
        bodyIndex = 9;
    }
    else if( bodyName == "" )
    {
        bodyIndex = -1;
    }
    else
    {
        throw std::runtime_error( "Error, could find body " + bodyName + " when getting planet index." );
    }
    return bodyIndex;
}

//! Set planet.
void ApproximateJplSolarSystemEphemerisBase::setPlanet( const std::string& bodyName )
{
    // Check if ephemeris data has been loaded, and reload data if not.
    if ( containerOfDataFromEphemerisFile_.size( ) == 0 )
    {
        reloadData( );
    }

    int bodyIndex = getPlanetIndex( bodyName );
    switch ( bodyIndex )
    {
    case 1:

        parseEphemerisLineData_( 18 );
        planetGravitationalParameter_ = celestial_body_constants::MERCURY_GRAVITATIONAL_PARAMETER;
        break;

    case 2:

        parseEphemerisLineData_( 20 );
        planetGravitationalParameter_ = celestial_body_constants::VENUS_GRAVITATIONAL_PARAMETER;
        break;

    case 3:

        parseEphemerisLineData_( 22 );
        planetGravitationalParameter_ = celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER +
                celestial_body_constants::MOON_GRAVITATIONAL_PARAMETER;
        break;

    case 4:

        parseEphemerisLineData_( 24 );
        planetGravitationalParameter_ = celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER;
        break;

    case 5:

        parseEphemerisLineData_( 26 );
        parseExtraTermsEphemerisLineData_( 48 );
        planetGravitationalParameter_ = celestial_body_constants::JUPITER_GRAVITATIONAL_PARAMETER;

        break;

    case 6:

        parseEphemerisLineData_( 28 );
        parseExtraTermsEphemerisLineData_( 49 );
        planetGravitationalParameter_ = celestial_body_constants::SATURN_GRAVITATIONAL_PARAMETER;

        break;

    case 7:

        parseEphemerisLineData_( 30 );
        parseExtraTermsEphemerisLineData_( 50 );
        planetGravitationalParameter_ = celestial_body_constants::URANUS_GRAVITATIONAL_PARAMETER;

        break;

    case 8:

        parseEphemerisLineData_( 32 );
        parseExtraTermsEphemerisLineData_( 51 );
        planetGravitationalParameter_ = celestial_body_constants::NEPTUNE_GRAVITATIONAL_PARAMETER;

        break;

    case 9:

        parseEphemerisLineData_( 34 );
        parseExtraTermsEphemerisLineData_( 52 );
        planetGravitationalParameter_ = celestial_body_constants::PLUTO_GRAVITATIONAL_PARAMETER;

        break;

    default:
        throw std::runtime_error( "Error, did not recognize planet index " + std::to_string( bodyIndex ) );
        break;
    }
}

//! Parse ephemeris line data.
void ApproximateJplSolarSystemEphemerisBase::parseEphemerisLineData_( const unsigned int& firstLineNumber )
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
void ApproximateJplSolarSystemEphemerisBase::parseExtraTermsEphemerisLineData_(
        const unsigned int& lineNumber )
{
    // Clear stringstream.
    ephemerisLineData_.clear( );

    // Read second line of data.
    ephemerisLineData_ << containerOfDataFromEphemerisFile_[ lineNumber ];

    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.planetName_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.additionalTermB_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.additionalTermC_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.additionalTermS_;
    ephemerisLineData_ >> approximatePlanetPositionsDataContainer_.additionalTermF_;
}

//! Load in ephemeris data for planets.
void ApproximateJplSolarSystemEphemerisBase::reloadData( )
{
    // Set  path to ephemeris file in file reader.
    std::string filePath_ = paths::getEphemerisDataFilesPath( ) + "/p_elem_t2.txt";

    // Open ephemeris file.
    std::ifstream ephemerisFile_( filePath_.c_str( ) );
    if ( ephemerisFile_.fail( ) )
    {
        throw std::runtime_error(
                    "Data file could not be opened:" + filePath_ );
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
