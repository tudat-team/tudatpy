/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/io/readAod1bGravityField.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <boost/date_time/gregorian/gregorian.hpp>

namespace tudat
{
namespace input_output
{
namespace
{

//! Remove header padding without changing retained source lines.
std::string trim( const std::string& value )
{
    const auto first = value.find_first_not_of( " \t\r" );

    return first == std::string::npos ? "" : value.substr( first, value.find_last_not_of( " \t\r" ) - first + 1 );
}

//! Report the product and one-based record location with each malformed input.
[[noreturn]] void fail( const std::string& source, std::size_t line, const std::string& message )
{
    throw std::runtime_error( source + ":" + std::to_string( line ) + ": " + message );
}

//! Parse one finite number, accepting Fortran exponents used by gravity products.
double number( std::string token )
{
    std::replace( token.begin( ), token.end( ), 'D', 'E' );
    std::replace( token.begin( ), token.end( ), 'd', 'e' );

    // Require a complete finite field after translating the optional Fortran exponent.
    std::size_t used = 0;

    const double value = std::stod( token, &used );

    if( used != token.size( ) || !std::isfinite( value ) )
    {
        throw std::runtime_error( "Expected a finite number: " + token );
    }
    return value;
}

//! Parse a complete header integer rather than accepting a numeric prefix with trailing text.
int integer( const std::string& token )
{
    // Retain the consumed length so malformed dimensional headers cannot be truncated silently.
    std::size_t used = 0;

    const int value = std::stoi( token, &used );

    if( used != token.size( ) )
    {
        throw std::runtime_error( "Expected an integer: " + token );
    }
    return value;
}

//! Decode the four format-defined component identifiers without aliases.
Aod1bComponent component( const std::string& token )
{
    if( token == "atm" )
    {
        return Aod1bComponent::atmosphere;
    }
    if( token == "ocn" )
    {
        return Aod1bComponent::ocean;
    }
    if( token == "glo" )
    {
        return Aod1bComponent::global;
    }
    if( token == "oba" )
    {
        return Aod1bComponent::oceanBottomPressure;
    }
    throw std::runtime_error( "Unknown AOD1B component: " + token );
}

//! Use Boost's existing Gregorian calendar to retain exact integral-day GPS epochs.
double gpsSeconds( const std::string& text )
{
    // RL07 samples occur at integral three-hour GPS epochs, without leap-second labels.
    int year = 0, month = 0, day = 0, hour = 0, minute = 0, second = 0;

    char separator = 0;

    std::istringstream input( text );

    if( !( input >> year >> separator ) || separator != '-' || !( input >> month >> separator ) || separator != '-' ||
        !( input >> day >> hour >> separator ) || separator != ':' || !( input >> minute >> separator ) || separator != ':' ||
        !( input >> second ) || !( input >> std::ws ).eof( ) || year < 1975 || year > 9999 || month < 1 || month > 12 || day < 1 ||
        day > 31 || hour < 0 || hour >= 24 || hour % 3 != 0 || minute != 0 || second != 0 )
    {
        throw std::runtime_error( "Expected an RL07 three-hour GPS date: " + text );
    }

    // Calendar subtraction supplies integral days; the GPS reference epoch is at noon.
    const boost::gregorian::date date( year, month, day );

    const long days = ( date - boost::gregorian::date( 2000, 1, 1 ) ).days( );

    return static_cast< double >( days ) * 86400.0 - 43200.0 + hour * 3600.0;
}

}  // namespace

//! Validate mutable input records before building any time-dependent model.
void validateAod1bGravityFieldData( const Aod1bGravityFieldData& data )
{
    if( data.release != 7 || data.softwareVersion.find( "atm_ocean_dealise.07_" ) != 0 || !std::isfinite( data.gravitationalParameter ) ||
        data.gravitationalParameter <= 0.0 || !std::isfinite( data.referenceRadius ) || data.referenceRadius <= 0.0 ||
        data.maximumDegree < 0 || data.maximumDegree > 180 || data.samples.empty( ) )
    {
        fail( data.source, 0, "Expected RL07 metadata, positive GM/radius and a nonempty expansion through degree 180" );
    }

    // Each supplied triangle must be complete, including explicit low-degree zero coefficients.
    std::set< std::pair< Aod1bComponent, double > > epochs;

    const std::size_t expected = ( data.maximumDegree + 1 ) * ( data.maximumDegree + 2 ) / 2;

    for( const auto& sample : data.samples )
    {
        if( sample.component < Aod1bComponent::atmosphere || sample.component > Aod1bComponent::oceanBottomPressure ||
            gpsSeconds( sample.gpsDate ) != sample.gpsSecondsSinceJ2000 ||
            !epochs.emplace( sample.component, sample.gpsSecondsSinceJ2000 ).second || sample.coefficients.size( ) != expected )
        {
            fail( data.source, 0, "Invalid/duplicate component epoch or incomplete coefficient triangle" );
        }
        for( const auto& entry : sample.coefficients )
        {
            if( entry.first.first < 0 || entry.first.first > data.maximumDegree || entry.first.second < 0 ||
                entry.first.second > entry.first.first || !std::isfinite( entry.second.values[ 0 ] ) ||
                !std::isfinite( entry.second.values[ 1 ] ) || ( entry.first.second == 0 && entry.second.values[ 1 ] != 0.0 ) )
            {
                fail( data.source, entry.second.sourceLine, "Invalid degree/order or potential coefficient" );
            }
        }
    }
}

//! Read the header and counted component blocks without inferring their physical metadata.
Aod1bGravityFieldData readAod1bGravityField( std::istream& stream, const std::string& source )
{
    // Keep the original header while interpreting only the format fields required below.
    Aod1bGravityFieldData result;

    result.source = source;

    std::map< std::string, std::string > header;

    std::string line;

    std::size_t lineNumber = 0;

    try
    {
        while( std::getline( stream, line ) )
        {
            ++lineNumber;
            result.headerLines.push_back( line );
            if( trim( line ) == "END OF HEADER" )
            {
                break;
            }

            const auto colon = line.find( ':' );

            if( colon != std::string::npos && !header.emplace( trim( line.substr( 0, colon ) ), trim( line.substr( colon + 1 ) ) ).second )
            {
                throw std::runtime_error( "Duplicate header field" );
            }
        }
        if( trim( line ) != "END OF HEADER" || header.at( "TIME EPOCH (GPS TIME)" ) != "2000-01-01 12:00:00" ||
            header.at( "FILE FORMAT 0=BINARY 1=ASCII" ) != "1" || header.at( "COEFFICIENTS ERRORS (YES/NO)" ) != "NO" ||
            header.at( "SOFTWARE VERSION" ).find( "atm_ocean_dealise.07_" ) != 0 ||
            header.at( "FILENAME" ).find( "_07.asc" ) == std::string::npos )
        {
            throw std::runtime_error( "Expected the explicit ASCII RL07 GPS potential format without error columns" );
        }
        result.release = 7;
        result.softwareVersion = header.at( "SOFTWARE VERSION" );
        result.gravitationalParameter = number( header.at( "CONSTANT GM [M^3/S^2]" ) );
        result.referenceRadius = number( header.at( "CONSTANT A [M]" ) );
        result.maximumDegree = integer( header.at( "MAXIMUM DEGREE" ) );

        // Normalization and record counts are explicit declarations, not coefficient guesses.
        const auto normalization = header.at( "COEFF. NORMALIZED (YES/NO)" );

        if( normalization != "YES" && normalization != "NO" )
        {
            throw std::runtime_error( "Unspecified coefficient normalization" );
        }
        result.fullyNormalized = normalization == "YES";

        const int sets = integer( header.at( "NUMBER OF DATA SETS" ) );

        const int count = integer( header.at( "NUMBER OF DATA RECORDS" ) );

        std::size_t records = 0;

        for( int index = 0; index < sets; ++index )
        {
            if( !std::getline( stream, line ) )
            {
                throw std::runtime_error( "Missing DATA SET" );
            }
            ++lineNumber;

            // The block label supplies the physical component and printed GPS timestamp.
            std::istringstream block( line );

            std::string data, set, id, coefficients, forToken, date, time, of, type, kind;

            int size = 0;

            if( !( block >> data >> set >> id >> size >> coefficients >> forToken >> date >> time >> of >> type >> kind ) ||
                data != "DATA" || set != "SET" || id != ( index + 1 < 10 ? "0" : "" ) + std::to_string( index + 1 ) + ":" ||
                coefficients != "COEFFICIENTS" || forToken != "FOR" || of != "OF" || type != "TYPE" || size <= 0 )
            {
                throw std::runtime_error( "Malformed DATA SET header" );
            }

            Aod1bGravitySample sample;

            sample.component = component( kind );
            sample.gpsDate = date + " " + time;
            sample.gpsSecondsSinceJ2000 = gpsSeconds( sample.gpsDate );
            for( int record = 0; record < size; ++record )
            {
                if( !std::getline( stream, line ) )
                {
                    throw std::runtime_error( "Truncated coefficient block" );
                }
                ++lineNumber;

                // Store values and source position only; keeping every raw line would duplicate the file.
                std::istringstream row( line );

                int degree = 0, order = 0;

                std::string cosine, sine, extra;

                if( !( row >> degree >> order >> cosine >> sine ) || row >> extra )
                {
                    throw std::runtime_error( "Expected degree, order, Cnm and Snm" );
                }
                if( !sample.coefficients
                             .emplace( std::make_pair( degree, order ),
                                       Aod1bGravityCoefficient{ { number( cosine ), number( sine ) }, lineNumber } )
                             .second )
                {
                    throw std::runtime_error( "Duplicate coefficient" );
                }
                ++records;
            }
            result.samples.push_back( std::move( sample ) );
        }
        if( count < 0 || records != static_cast< std::size_t >( count ) )
        {
            throw std::runtime_error( "Header coefficient count disagrees with DATA SET records" );
        }
        while( std::getline( stream, line ) )
        {
            ++lineNumber;
            if( !trim( line ).empty( ) )
            {
                throw std::runtime_error( "Unexpected record after the declared DATA SET blocks" );
            }
        }
        validateAod1bGravityFieldData( result );

        // The header uses the same continuous GPS origin as the component timestamps.
        const auto range = std::minmax_element( result.samples.begin( ), result.samples.end( ), []( const auto& a, const auto& b ) {
            return a.gpsSecondsSinceJ2000 < b.gpsSecondsSinceJ2000;
        } );

        if( number( header.at( "TIME FIRST OBS(SEC PAST EPOCH)" )
                            .substr( 0, header.at( "TIME FIRST OBS(SEC PAST EPOCH)" ).find( ' ' ) ) ) !=
                    range.first->gpsSecondsSinceJ2000 ||
            number( header.at( "TIME LAST OBS(SEC PAST EPOCH)" ).substr( 0, header.at( "TIME LAST OBS(SEC PAST EPOCH)" ).find( ' ' ) ) ) !=
                    range.second->gpsSecondsSinceJ2000 )
        {
            throw std::runtime_error( "GPS observation range disagrees with DATA SET dates" );
        }
    }
    catch( const std::exception& error )
    {
        fail( source, lineNumber, error.what( ) );
    }
    return result;
}

//! Open one uncompressed RL07 product and report missing files explicitly.
Aod1bGravityFieldData readAod1bGravityField( const std::string& fileName )
{
    // Delegate all format validation to the same stream reader.
    std::ifstream stream( fileName );

    if( !stream )
    {
        fail( fileName, 0, "Cannot open AOD1B product" );
    }
    return readAod1bGravityField( stream, fileName );
}

//! Preserve daily source boundaries while moving each parsed product into the returned collection.
std::vector< Aod1bGravityFieldData > readAod1bGravityFields( const std::vector< std::string >& fileNames )
{
    // Reserve the collection before moving each parsed daily product into it.
    std::vector< Aod1bGravityFieldData > result;

    result.reserve( fileNames.size( ) );
    for( const auto& file : fileNames )
    {
        result.push_back( readAod1bGravityField( file ) );
    }
    return result;
}

}  // namespace input_output
}  // namespace tudat
