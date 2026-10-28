/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <cctype>
#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "tudat/astro/earth_orientation/readEarthOrientationData.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace earth_orientation
{
namespace
{
const double arcsecondToRadian = mathematical_constants::PI / ( 180.0 * 3600.0 );

std::string trimWhitespace( const std::string& text )
{
    const std::string whitespace = " \t\r\n";
    const auto first = text.find_first_not_of( whitespace );
    return first == std::string::npos ? "" : text.substr( first, text.find_last_not_of( whitespace ) - first + 1 );
}

std::vector< std::string > splitIntoTokens( const std::string& text )
{
    std::istringstream stream( text );
    stream.imbue( std::locale::classic( ) );
    std::vector< std::string > tokens;
    std::string token;
    while( stream >> token )
    {
        tokens.push_back( token );
    }
    return tokens;
}

double parseEarthOrientationNumber( const std::string& text )
{
    std::size_t consumed = 0;
    double value = TUDAT_NAN;
    try
    {
        value = std::stod( text, &consumed );
    }
    catch( const std::exception& )
    {
        consumed = 0;
    }
    if( consumed != text.size( ) || !std::isfinite( value ) )
    {
        throw std::runtime_error( "invalid finite number: " + text );
    }
    return value;
}

EarthOrientationRecord c04Record( const std::vector< std::string >& tokens )
{
    const bool modern = tokens.size( ) == 21;
    if( !modern && tokens.size( ) != 16 )
    {
        throw std::runtime_error( "C04 requires 16 or 21 columns" );
    }
    std::vector< double > values;
    for( const auto& token : tokens )
    {
        values.push_back( parseEarthOrientationNumber( token ) );
    }

    // The 21-column format has an additional hour column, and the excess length of day after the pole offset rates.
    const int offset = modern ? 1 : 0;
    EarthOrientationRecord record;
    record.utcModifiedJulianDate = values[ 3 + offset ];
    record.polarMotion = { values[ 4 + offset ] * arcsecondToRadian, values[ 5 + offset ] * arcsecondToRadian };
    record.ut1MinusUtc = values[ 6 + offset ];
    record.celestialPoleOffsets = std::array< double, 2 >{ values[ 8 ] * arcsecondToRadian, values[ 9 ] * arcsecondToRadian };
    record.lengthOfDay = modern ? values[ 12 ] : values[ 7 ];
    return record;
}

std::optional< double > optionalNumber( const std::string& line, const int first, const int count, const double scale = 1.0 )
{
    // USNO specification columns are one-based.
    const std::string text = trimWhitespace( line.substr( first - 1, count ) );
    return text.empty( ) ? std::nullopt : std::optional< double >( parseEarthOrientationNumber( text ) * scale );
}

std::optional< EarthOrientationRecord > finalsRecord( const std::string& line )
{
    if( line.size( ) < 134 )
    {
        throw std::runtime_error( "truncated finals2000A record (requires columns 1-134)" );
    }
    EarthOrientationRecord record;
    record.utcModifiedJulianDate = *optionalNumber( line, 8, 8 );
    const auto x = optionalNumber( line, 19, 9, arcsecondToRadian );
    const auto y = optionalNumber( line, 38, 9, arcsecondToRadian );
    const auto ut1 = optionalNumber( line, 59, 10 );
    if( !x && !y && !ut1 )
    {
        return std::nullopt;  // date-only future placeholders
    }
    if( !x || !y || !ut1 )
    {
        throw std::runtime_error( "missing mandatory Bulletin A polar motion or UT1" );
    }
    record.polarMotion = { *x, *y };
    record.ut1MinusUtc = *ut1;
    const auto dX = optionalNumber( line, 98, 9, arcsecondToRadian / 1000.0 );
    const auto dY = optionalNumber( line, 117, 9, arcsecondToRadian / 1000.0 );
    if( dX && dY )
    {
        record.celestialPoleOffsets = std::array< double, 2 >{ *dX, *dY };
    }
    record.lengthOfDay = optionalNumber( line, 80, 7, 0.001 );
    return record;  // Columns 135-185 are Bulletin B and are intentionally not used.
}
}  // namespace

std::map< double, EarthOrientationRecord > readEarthOrientationData( std::istream& stream,
                                                                     const EarthOrientationFormat format,
                                                                     const std::string& source )
{
    std::map< double, EarthOrientationRecord > records;
    std::string line;
    std::size_t lineNumber = 0;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        const std::string content = trimWhitespace( line );

        // C04 header lines do not start with a digit.
        if( content.empty( ) ||
            ( format == EarthOrientationFormat::c04 && !std::isdigit( static_cast< unsigned char >( content.front( ) ) ) ) )
        {
            continue;
        }
        try
        {
            const auto record = format == EarthOrientationFormat::c04
                    ? std::optional< EarthOrientationRecord >( c04Record( splitIntoTokens( content ) ) )
                    : finalsRecord( line );
            if( !record )
            {
                continue;
            }
            if( !records.empty( ) && record->utcModifiedJulianDate <= records.rbegin( )->first )
            {
                throw std::runtime_error( "duplicate or unsorted UTC epoch" );
            }
            records.emplace( record->utcModifiedJulianDate, *record );
        }
        catch( const std::exception& exception )
        {
            throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": " + exception.what( ) );
        }
    }
    if( records.empty( ) )
    {
        throw std::runtime_error( source + ": no Earth orientation records" );
    }
    return records;
}

std::map< double, EarthOrientationRecord > readEarthOrientationData( const std::string& fileName, const EarthOrientationFormat format )
{
    std::ifstream stream( fileName );
    if( !stream )
    {
        throw std::runtime_error( "Cannot open Earth orientation file: " + fileName );
    }
    return readEarthOrientationData( stream, format, fileName );
}

}  // namespace earth_orientation
}  // namespace tudat
