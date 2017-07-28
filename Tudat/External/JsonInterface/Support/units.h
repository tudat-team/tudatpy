/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_UNITS_H
#define TUDAT_JSONINTERFACE_UNITS_H

// #include <stdio.h>
#include <map>

#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>

namespace tudat
{

namespace json_interface
{


//! Equivalent SI value of units supported by json_interface.
//! IMPORTANT: keys are all defined in lowercase here, but user can also use uppercase in input files.
static std::map< std::string, double > SIUnits =
{
    // Time [ s ]
    { "s", 1 },
    { "min", 60 },
    { "h", 3600 },
    { "sd", physical_constants::SIDEREAL_DAY  },
    { "d", physical_constants::JULIAN_DAY },
    { "y", physical_constants::JULIAN_YEAR },
    { "sy", physical_constants::SIDEREAL_YEAR },

    // Distance [ m ]
    { "m", 1 },
    { "km", 1e3 },
    { "au", physical_constants::ASTRONOMICAL_UNIT },
};


//! -DOC
template< typename T >
void split( const std::string& string, char delimiter, T result ) {
    std::stringstream stream;
    stream.str( string );
    std::string item;
    while ( std::getline( stream, item, delimiter ) )
    {
        *(result++) = item;
    }
}

//! -DOC
static inline std::vector< std::string > split( const std::string& string, char delimiter )
{
    std::vector< std::string > parts;
    split( string, delimiter, std::back_inserter( parts ) );
    return parts;
}


//! -DOC
template< typename T >
T convertToSecondsSinceJ2000( const std::string& date )
{
    using namespace boost::posix_time;
    ptime pt;
    try
    {
        pt = time_from_string( date );
    }
    catch ( ... )
    {
        try
        {
            pt = from_iso_string( date );
        }
        catch ( ... )
        {
            try
            {
                pt = from_iso_extended_string( date );
            }
            catch ( ... )
            {
                std::cerr << "Unrecognized date format. Supported formats: \"1992-02-14 07:30:00\", \"19920214T0730\"."
                          << std::endl;
                throw;
            }
        }
    }
    time_duration duration = pt - time_from_string( "2000-01-01 12:00:00" );
    return duration.total_seconds( );
}

//! -DOC
template< typename T >
T convertToSIUnits( T number, const std::string& units )
{
    std::string lowercaseUnits = boost::algorithm::to_lower_copy( units );
    if ( SIUnits.count( lowercaseUnits ) == 0 )
    {
        std::cerr << "Unrecognized units \"" << units << "\"." << std::endl;
        std::cerr << "Supported units:";
        for ( auto entry : SIUnits )
        {
            std::cerr << " " << entry.first;
        }
        std::cerr << "." << std::endl;
        throw;
    }
    return number * static_cast< T >( SIUnits.at( lowercaseUnits ) );
}

//! -DOC
template< typename T >
T parseMagnitudeWithUnits( const std::string& text )
{
    std::vector< std::string > parts = split( text, ' ' );
    T number = boost::lexical_cast< T >( parts.at( 0 ) );
    std::string originalUnits = parts.at( 1 );
    return convertToSIUnits< T >( number, originalUnits );
}


} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_UNITS_H
