/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/io/readSp3File.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>
#include <stdexcept>

#include "tudat/io/textParsing.h"
#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{
namespace input_output
{
namespace
{

//! Return a trimmed fixed-width field; omitted trailing columns are blank in SP3.
std::string field( const std::string& line, const std::size_t first, const std::size_t length )
{
    return first < line.size( ) ? text_parsing::trim( line.substr( first, length ) ) : "";
}

//! Parse a finite mandatory SP3 value using the existing text-reader convention.
double number( const std::string& text )
{
    const auto value = text_parsing::parseFiniteNumber( text );

    if( !value )
    {
        throw std::runtime_error( "invalid SP3 numeric field: " + text );
    }

    return *value;
}

//! Preserve fractional seconds separately from the large integral calendar epoch.
Time epoch( const std::string& text )
{
    std::istringstream stream( text );

    stream.imbue( std::locale::classic( ) );

    int year, month, day, hour, minute;

    long double second;

    if( !( stream >> year >> month >> day >> hour >> minute >> second ) || hour < 0 || hour > 23 || minute < 0 || minute > 59 ||
        !std::isfinite( second ) || second < 0.0L || second >= 61.0L )
    {
        throw std::runtime_error( "invalid SP3 calendar epoch: " + text );
    }

    const long double days = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch< long double >(
            year, month, day, 0, 0, 0.0L, basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    return Time( static_cast< int >( days * 24.0L ) + hour, minute * 60.0L + second );
}

}  // namespace

//! Read the declared frame/time metadata and raw satellite records; no clock or frame conversion is inferred.
Sp3Data readSp3File( std::istream& stream, const std::string& source )
{
    Sp3Data output;

    output.source = source;

    std::string line;

    std::size_t lineNumber = 0;

    std::size_t epochCount = 0;

    Time currentEpoch;

    std::map< std::string, bool > velocityRead;

    bool timeSystemRead = false;

    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        try
        {
            if( lineNumber == 1 )
            {
                if( line.size( ) < 60 || line[ 0 ] != '#' || ( line[ 1 ] != 'c' && line[ 1 ] != 'd' ) ||
                    ( line[ 2 ] != 'P' && line[ 2 ] != 'V' ) )
                {
                    throw std::runtime_error( "expected an SP3c/d position or velocity header" );
                }

                output.version = line[ 1 ];
                output.hasVelocityRecords = line[ 2 ] == 'V';
                output.declaredEpochCount = std::stoi( field( line, 32, 7 ) );
                output.dataUsed = field( line, 40, 5 );
                output.coordinateSystem = field( line, 46, 5 );
                output.orbitType = field( line, 52, 3 );
                output.agency = field( line, 56, 4 );
            }
            else if( line.rfind( "##", 0 ) == 0 )
            {
                output.epochInterval = number( field( line, 24, 14 ) );
            }
            else if( line.rfind( "+ ", 0 ) == 0 )
            {
                for( std::size_t index = 9; index < line.size( ); index += 3 )
                {
                    const std::string satellite = field( line, index, 3 );

                    if( !satellite.empty( ) && satellite != "0" )
                    {
                        if( output.records.count( satellite ) )
                        {
                            throw std::runtime_error( "duplicate satellite in SP3 header: " + satellite );
                        }

                        output.satelliteIds.push_back( satellite );
                        output.records[ satellite ] = {};
                    }
                }
            }
            else if( line.rfind( "%c", 0 ) == 0 && !timeSystemRead )
            {
                output.timeSystem = field( line, 9, 3 );
                timeSystemRead = true;
            }
            else if( line.rfind( "/*", 0 ) == 0 || line.rfind( "%/*", 0 ) == 0 )
            {
                output.comments.push_back( line );
            }
            else if( !line.empty( ) && line[ 0 ] == '*' )
            {

                // Complete the previous epoch before accepting another one.
                if( epochCount &&
                    ( velocityRead.size( ) != output.satelliteIds.size( ) ||
                      ( output.hasVelocityRecords &&
                        std::any_of( velocityRead.begin( ), velocityRead.end( ), []( const auto& entry ) { return !entry.second; } ) ) ) )
                {
                    throw std::runtime_error( "incomplete SP3 epoch" );
                }

                const Time nextEpoch = epoch( line.substr( 1 ) );

                if( epochCount && !( nextEpoch > currentEpoch ) )
                {
                    throw std::runtime_error( "SP3 epochs are not strictly increasing" );
                }

                currentEpoch = nextEpoch;
                velocityRead.clear( );
                ++epochCount;
            }
            else if( !line.empty( ) && ( line[ 0 ] == 'P' || line[ 0 ] == 'V' ) )
            {
                const std::string satellite = field( line, 1, 3 );

                if( !epochCount || !output.records.count( satellite ) )
                {
                    throw std::runtime_error( "SP3 record has no declared satellite or epoch" );
                }

                Eigen::Vector3d coordinates;

                for( int axis = 0; axis < 3; ++axis )
                {
                    coordinates[ axis ] = number( field( line, 4 + 14 * axis, 14 ) );
                }

                const double clock = number( field( line, 46, 14 ) );

                const bool available = !coordinates.isZero( 0.0 );

                std::vector< Sp3Record >& records = output.records.at( satellite );

                if( line[ 0 ] == 'P' )
                {
                    if( velocityRead.count( satellite ) )
                    {
                        throw std::runtime_error( "duplicate SP3 position at one epoch" );
                    }

                    Sp3Record record;

                    record.epoch = currentEpoch;
                    record.sourceLine = lineNumber;
                    record.positionAvailable = available;
                    if( available )
                    {
                        record.position = coordinates * 1000.0;
                    }

                    record.clockCorrection = clock >= 999999.0 ? TUDAT_NAN : clock * 1.0e-6;
                    record.clockEvent = field( line, 74, 1 ) == "E";
                    record.clockPredicted = field( line, 75, 1 ) == "P";
                    record.maneuver = field( line, 78, 1 ) == "M";
                    record.orbitPredicted = field( line, 79, 1 ) == "P";
                    records.push_back( record );
                    velocityRead[ satellite ] = false;
                }
                else
                {
                    if( !output.hasVelocityRecords || !velocityRead.count( satellite ) || velocityRead.at( satellite ) )
                    {
                        throw std::runtime_error( "unexpected or duplicate SP3 velocity record" );
                    }

                    Sp3Record& record = records.back( );

                    record.velocityAvailable = available;
                    if( available )
                    {
                        record.velocity = coordinates * 0.1;
                    }

                    record.clockRate = clock >= 999999.0 ? TUDAT_NAN : clock * 1.0e-10;
                    velocityRead[ satellite ] = true;
                }
            }
            else if( line.rfind( "EOF", 0 ) == 0 )
            {
                break;
            }

            // EP/EV covariance records and unused format control lines do not change the samples.
        }
        catch( const std::exception& error )
        {
            throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": " + error.what( ) );
        }
    }

    // ILRSA combined products omit the literal EOF line. Physical EOF is accepted only with the declared epoch count
    // and a complete final set of satellite records, so an incomplete final epoch still fails.
    if( stream.bad( ) || output.records.empty( ) || !timeSystemRead || output.timeSystem.empty( ) ||
        epochCount != static_cast< std::size_t >( output.declaredEpochCount ) || velocityRead.size( ) != output.satelliteIds.size( ) ||
        ( output.hasVelocityRecords &&
          std::any_of( velocityRead.begin( ), velocityRead.end( ), []( const auto& entry ) { return !entry.second; } ) ) )
    {
        throw std::runtime_error( source + ": incomplete SP3 header or samples, or stream read failure" );
    }

    return output;
}

//! Open an uncompressed SP3 product and retain its filename as the source identifier.
Sp3Data readSp3File( const std::string& fileName )
{
    std::ifstream stream( fileName );

    if( !stream )
    {
        throw std::runtime_error( "Cannot open SP3 file: " + fileName );
    }

    return readSp3File( stream, fileName );
}

}  // namespace input_output
}  // namespace tudat
