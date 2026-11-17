/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readSp3File.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{

namespace input_output
{

//! Track which records have been read while assembling one satellite state at one epoch.
struct EpochRecordFlags {
    bool hasPosition = false;
    bool hasVelocity = false;
};

// These helpers implement the SP3 fixed-width grammar. They remain translation-unit local
// because callers outside the reader should receive parsed data instead of depending on file-format details.

//! Extract and trim one fixed-width field from an SP3 record.
static std::string getFixedWidthField( const std::string& line, const unsigned int startIndex, const unsigned int length )
{
    // Accept short records here so that the typed parsers can report a regular invalid-field error.
    if( startIndex >= line.size( ) )
    {
        return "";
    }
    std::string field = line.substr( startIndex, std::min( length, static_cast< unsigned int >( line.size( ) - startIndex ) ) );
    boost::algorithm::trim( field );
    return field;
}

//! Attempt to parse an integer from a fixed-width SP3 field.
static bool tryParseFixedWidthInt( const std::string& line, const unsigned int startIndex, const unsigned int length, int& output )
{
    const std::string value = getFixedWidthField( line, startIndex, length );
    try
    {
        // Require the complete field to be numeric; std::stoi otherwise accepts prefixes such as "2X".
        std::size_t parsedCharacterCount = 0;
        output = std::stoi( value, &parsedCharacterCount );
        return parsedCharacterCount == value.size( );
    }
    catch( ... )
    {
        return false;
    }
}

//! Attempt to parse a floating-point value from a fixed-width SP3 field.
static bool tryParseFixedWidthDouble( const std::string& line, const unsigned int startIndex, const unsigned int length, double& output )
{
    const std::string value = getFixedWidthField( line, startIndex, length );
    try
    {
        // Reject trailing characters and non-finite values before they enter a satellite state.
        std::size_t parsedCharacterCount = 0;
        output = std::stod( value, &parsedCharacterCount );
        return parsedCharacterCount == value.size( ) && std::isfinite( output );
    }
    catch( ... )
    {
        return false;
    }
}

//! Register a non-placeholder satellite identifier encountered in the file.
static void addSatelliteIdToStateMap( const std::string& rawSatelliteId,
                                      Sp3FileContents& fileContents,
                                      std::vector< std::string >& knownSatelliteIds )
{
    // Placeholder fields fill unused slots on '+' records and must not become satellites.
    const std::string trimmedSatelliteId = boost::algorithm::trim_copy( rawSatelliteId );
    if( trimmedSatelliteId.empty( ) || trimmedSatelliteId == "0" || trimmedSatelliteId == "00" || trimmedSatelliteId == "000" )
    {
        return;
    }

    if( fileContents.satelliteStates.count( trimmedSatelliteId ) == 0 )
    {
        fileContents.satelliteStates[ trimmedSatelliteId ] = std::map< double, Eigen::Vector6d >( );
        knownSatelliteIds.push_back( trimmedSatelliteId );
    }
}

//! Validate the records for one epoch and append them to the parsed state histories.
static void validateAndFlushCurrentEpoch( const double currentTime,
                                          const std::string& fileName,
                                          const std::vector< std::string >& knownSatelliteIds,
                                          const bool velocityRecordsExpected,
                                          const std::map< std::string, Eigen::Vector6d >& currentStates,
                                          const std::map< std::string, EpochRecordFlags >& currentRecordFlags,
                                          Sp3FileContents& fileContents )
{
    // An epoch is complete only when every declared satellite has the record types promised by the header.
    for( const std::string& satelliteId : knownSatelliteIds )
    {
        const auto flagIt = currentRecordFlags.find( satelliteId );
        if( flagIt == currentRecordFlags.end( ) || !flagIt->second.hasPosition )
        {
            throw std::runtime_error( "Error when reading SP3 file: missing position record for satellite " + satelliteId + " at epoch " +
                                      std::to_string( currentTime ) + " in " + fileName );
        }

        if( velocityRecordsExpected && !flagIt->second.hasVelocity )
        {
            throw std::runtime_error( "Error when reading SP3 file: missing velocity record for satellite " + satelliteId + " at epoch " +
                                      std::to_string( currentTime ) + " in " + fileName );
        }

        const auto stateIt = currentStates.find( satelliteId );
        if( stateIt == currentStates.end( ) )
        {
            throw std::runtime_error( "Error when reading SP3 file: missing state container for satellite " + satelliteId + " at epoch " +
                                      std::to_string( currentTime ) + " in " + fileName );
        }

        fileContents.satelliteStates[ satelliteId ][ currentTime ] = stateIt->second;
    }
}

//! Reconstruct velocities for a position-only SP3 file using three-point finite differences.
static void deriveVelocitiesFromPositions( Sp3FileContents& fileContents, const std::string& fileName )
{
    for( auto& satellite : fileContents.satelliteStates )
    {
        // A quadratic derivative requires three distinct position samples.
        if( satellite.second.size( ) < 3 )
        {
            throw std::runtime_error(
                    "Error when reading position-only SP3 file: at least three epochs are required to derive "
                    "second-order velocities for satellite " +
                    satellite.first + " in " + fileName );
        }

        // Index the ordered map once so adjacent samples can be selected by position.
        std::vector< std::map< double, Eigen::Vector6d >::iterator > states;
        states.reserve( satellite.second.size( ) );
        for( auto stateIterator = satellite.second.begin( ); stateIterator != satellite.second.end( ); ++stateIterator )
        {
            states.push_back( stateIterator );
        }

        for( unsigned int i = 0; i < states.size( ); ++i )
        {
            // Use a forward stencil first, a backward stencil last, and a centered stencil elsewhere.
            const unsigned int firstStencilIndex = ( i == 0 ) ? 0 : ( i == states.size( ) - 1 ? i - 2 : i - 1 );
            const std::array< double, 3 > times = { states[ firstStencilIndex ]->first,
                                                    states[ firstStencilIndex + 1 ]->first,
                                                    states[ firstStencilIndex + 2 ]->first };
            const std::array< Eigen::Vector3d, 3 > positions = { states[ firstStencilIndex ]->second.segment< 3 >( 0 ),
                                                                 states[ firstStencilIndex + 1 ]->second.segment< 3 >( 0 ),
                                                                 states[ firstStencilIndex + 2 ]->second.segment< 3 >( 0 ) };

            if( positions[ 0 ].allFinite( ) && positions[ 1 ].allFinite( ) && positions[ 2 ].allFinite( ) )
            {
                // Differentiate the quadratic Lagrange interpolant directly at the current epoch.
                Eigen::Vector3d derivative = Eigen::Vector3d::Zero( );
                for( unsigned int stencilIndex = 0; stencilIndex < 3; ++stencilIndex )
                {
                    const unsigned int secondIndex = ( stencilIndex + 1 ) % 3;
                    const unsigned int thirdIndex = ( stencilIndex + 2 ) % 3;
                    const double denominator =
                            ( times[ stencilIndex ] - times[ secondIndex ] ) * ( times[ stencilIndex ] - times[ thirdIndex ] );
                    if( denominator == 0.0 )
                    {
                        throw std::runtime_error( "Cannot derive SP3 velocities from duplicate epochs." );
                    }
                    derivative += positions[ stencilIndex ] * ( 2.0 * states[ i ]->first - times[ secondIndex ] - times[ thirdIndex ] ) /
                            denominator;
                }
                states[ i ]->second.segment< 3 >( 3 ) = derivative;
            }
        }
    }
    fileContents.velocitiesWereDerived = true;
}

//! Retrieve and validate one satellite state history from parsed SP3 contents.
std::map< double, Eigen::Vector6d > Sp3FileContents::getSatelliteStateHistory( const std::string& satelliteId ) const
{
    // Validate at the data boundary so ephemeris creation receives a complete Cartesian state history.
    const auto satelliteIterator = satelliteStates.find( satelliteId );
    if( satelliteIterator == satelliteStates.end( ) )
    {
        throw std::invalid_argument( "Satellite '" + satelliteId + "' is not present in the SP3 file." );
    }
    if( satelliteIterator->second.empty( ) )
    {
        throw std::invalid_argument( "Satellite '" + satelliteId + "' has no states in the SP3 file." );
    }
    for( const auto& state : satelliteIterator->second )
    {
        if( !state.second.allFinite( ) )
        {
            throw std::runtime_error( "SP3 state history requires finite position and velocity records for satellite '" + satelliteId +
                                      "'. The selected SP3 file contains a missing value." );
        }
    }
    return satelliteIterator->second;
}

//! Read and validate an SP3 file, converting its position and velocity records to SI units.
std::shared_ptr< Sp3FileContents > readSp3File( const std::string& fileName, const double referenceJulianDay )
{
    // Store raw-file metadata and SI-unit state histories in the data-input object.
    std::shared_ptr< Sp3FileContents > fileContents = std::make_shared< Sp3FileContents >( );
    fileContents->referenceJulianDay = referenceJulianDay;

    // Open the file before initializing per-epoch parser state.
    std::ifstream stream( fileName, std::ios_base::in );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening SP3 file: " + fileName );
    }

    std::string currentLine;

    // These containers accumulate the P/V records belonging to the current '*' epoch.
    std::vector< std::string > knownSatelliteIds;
    std::map< std::string, Eigen::Vector6d > currentStates;
    std::map< std::string, EpochRecordFlags > currentRecordFlags;

    // These flags retain file-level declarations needed to validate records as they are encountered.
    double currentTime = TUDAT_NAN;
    bool hasCurrentEpoch = false;
    bool velocityRecordsExpected = false;
    bool headerSeen = false;
    bool satelliteCountSeen = false;
    int parsedEpochCount = 0;

    while( std::getline( stream, currentLine ) )
    {
        // Dispatch records by prefix while preserving the original line for fixed-width extraction.
        const std::string trimmedLine = boost::algorithm::trim_copy( currentLine );
        if( trimmedLine.empty( ) )
        {
            continue;
        }

        if( trimmedLine.rfind( "EOF", 0 ) == 0 )
        {
            break;
        }

        const char firstCharacter = currentLine.at( 0 );

        // Header line '#': fixed-width parsing only.
        if( firstCharacter == '#' && ( currentLine.size( ) < 2 || currentLine.at( 1 ) != '#' ) )
        {
            if( headerSeen )
            {
                throw std::runtime_error( "Error when reading SP3 file: duplicate first header line in " + fileName );
            }
            if( currentLine.size( ) < 3 || currentLine.at( 1 ) < 'a' || currentLine.at( 1 ) > 'd' )
            {
                throw std::runtime_error( "Error when reading SP3 file: only SP3 versions a, b, c, and d are supported in " + fileName );
            }
            if( currentLine.at( 2 ) != 'P' && currentLine.at( 2 ) != 'V' )
            {
                throw std::runtime_error( "Error when reading SP3 file: header must declare position ('P') or velocity ('V') mode in " +
                                          fileName );
            }

            headerSeen = true;
            fileContents->formatVersion = currentLine.at( 1 );
            velocityRecordsExpected = ( currentLine.at( 2 ) == 'V' );
            fileContents->hasVelocityRecords = velocityRecordsExpected;

            // Convert the calendar fields to seconds from the caller-selected Julian-day reference.
            int year, month, day, hour, minute;
            double second;
            if( !( tryParseFixedWidthInt( currentLine, 3, 4, year ) && tryParseFixedWidthInt( currentLine, 8, 2, month ) &&
                   tryParseFixedWidthInt( currentLine, 11, 2, day ) && tryParseFixedWidthInt( currentLine, 14, 2, hour ) &&
                   tryParseFixedWidthInt( currentLine, 17, 2, minute ) && tryParseFixedWidthDouble( currentLine, 20, 11, second ) ) )
            {
                throw std::runtime_error( "Error when reading SP3 file: invalid epoch in first header line in " + fileName );
            }
            fileContents->startEpoch = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch< double >(
                                               year, month, day, 0, 0, 0.0, referenceJulianDay ) *
                            physical_constants::JULIAN_DAY +
                    static_cast< double >( hour ) * 3600.0 + static_cast< double >( minute ) * 60.0 + second;

            int declaredEpochCount;
            if( !tryParseFixedWidthInt( currentLine, 32, 7, declaredEpochCount ) )
            {
                throw std::runtime_error( "Error when reading SP3 file: invalid declared epoch count in " + fileName );
            }
            fileContents->declaredNumberOfEpochs = declaredEpochCount;

            const std::string frameName = getFixedWidthField( currentLine, 46, 5 );
            if( !frameName.empty( ) )
            {
                fileContents->frameName = frameName;
            }

            const std::string analysisCenter = getFixedWidthField( currentLine, 56, 4 );
            if( !analysisCenter.empty( ) )
            {
                fileContents->analysisCenter = analysisCenter;
            }

            continue;
        }

        // Header line '##': fixed-width parsing only.
        if( currentLine.rfind( "##", 0 ) == 0 )
        {
            double declaredEpochInterval;
            if( !tryParseFixedWidthDouble( currentLine, 24, 14, declaredEpochInterval ) )
            {
                throw std::runtime_error( "Error when reading SP3 file: invalid declared epoch interval in " + fileName );
            }
            fileContents->declaredEpochInterval = declaredEpochInterval;
            continue;
        }

        // Satellite ID list line '+': fixed-width parsing only.
        if( firstCharacter == '+' && ( currentLine.size( ) < 2 || currentLine.at( 1 ) != '+' ) )
        {
            if( !satelliteCountSeen )
            {
                int declaredSatelliteCount;
                if( !tryParseFixedWidthInt( currentLine, 3, 3, declaredSatelliteCount ) )
                {
                    throw std::runtime_error( "Error when reading SP3 file: invalid declared satellite count in " + fileName );
                }
                fileContents->declaredNumberOfSatellites = declaredSatelliteCount;
                satelliteCountSeen = true;
            }
            for( unsigned int i = 0; i < 17; i++ )
            {
                addSatelliteIdToStateMap( getFixedWidthField( currentLine, 9 + 3 * i, 3 ), *fileContents, knownSatelliteIds );
            }
            continue;
        }

        // Time-system line '%c': fixed-width parsing only.
        if( currentLine.rfind( "%c", 0 ) == 0 )
        {
            if( fileContents->timeScale.empty( ) )
            {
                const std::string timeScale = getFixedWidthField( currentLine, 9, 3 );
                if( !timeScale.empty( ) )
                {
                    fileContents->timeScale = boost::iequals( timeScale, "ccc" ) ? "GPS" : boost::to_upper_copy( timeScale );
                }
            }
            continue;
        }

        // Epoch line '*': fixed-width parsing only.
        if( firstCharacter == '*' )
        {
            int year, month, day, hour, minute;
            double second;
            if( !( tryParseFixedWidthInt( currentLine, 3, 4, year ) && tryParseFixedWidthInt( currentLine, 8, 2, month ) &&
                   tryParseFixedWidthInt( currentLine, 11, 2, day ) && tryParseFixedWidthInt( currentLine, 14, 2, hour ) &&
                   tryParseFixedWidthInt( currentLine, 17, 2, minute ) && tryParseFixedWidthDouble( currentLine, 20, 11, second ) ) )
            {
                throw std::runtime_error( "Error when reading SP3 file: invalid epoch line in " + fileName );
            }

            const double nextEpochTime = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch< double >(
                                                 year, month, day, 0, 0, 0.0, referenceJulianDay ) *
                            physical_constants::JULIAN_DAY +
                    static_cast< double >( hour ) * 3600.0 + static_cast< double >( minute ) * 60.0 + second;

            if( hasCurrentEpoch )
            {
                // Store the previous epoch before resetting the temporary state containers.
                if( nextEpochTime <= currentTime )
                {
                    throw std::runtime_error( "Error when reading SP3 file: epochs are not strictly increasing in " + fileName );
                }
                validateAndFlushCurrentEpoch( currentTime,
                                              fileName,
                                              knownSatelliteIds,
                                              velocityRecordsExpected,
                                              currentStates,
                                              currentRecordFlags,
                                              *fileContents );
            }

            currentTime = nextEpochTime;
            hasCurrentEpoch = true;
            parsedEpochCount++;
            currentStates.clear( );
            currentRecordFlags.clear( );

            for( const std::string& satelliteId : knownSatelliteIds )
            {
                currentStates[ satelliteId ] = Eigen::Vector6d::Constant( TUDAT_NAN );
                currentRecordFlags[ satelliteId ] = EpochRecordFlags( );
            }

            continue;
        }

        // State lines 'P' and 'V': fixed-width parsing only.
        if( firstCharacter == 'P' || firstCharacter == 'V' )
        {
            if( !hasCurrentEpoch )
            {
                throw std::runtime_error( "Error when reading SP3 file: encountered state record before first epoch in " + fileName );
            }

            const std::string satelliteId = getFixedWidthField( currentLine, 1, 3 );
            if( satelliteId.empty( ) )
            {
                throw std::runtime_error( "Error when reading SP3 file: invalid state record satellite id in " + fileName );
            }

            addSatelliteIdToStateMap( satelliteId, *fileContents, knownSatelliteIds );
            if( currentStates.count( satelliteId ) == 0 )
            {
                currentStates[ satelliteId ] = Eigen::Vector6d::Constant( TUDAT_NAN );
                currentRecordFlags[ satelliteId ] = EpochRecordFlags( );
            }

            double x, y, z;
            const bool readSucceeded = tryParseFixedWidthDouble( currentLine, 4, 14, x ) &&
                    tryParseFixedWidthDouble( currentLine, 18, 14, y ) && tryParseFixedWidthDouble( currentLine, 32, 14, z );

            if( !readSucceeded )
            {
                throw std::runtime_error( "Error when reading SP3 file: invalid numeric state record in " + fileName );
            }

            // SP3 uses an all-zero position or velocity vector for a bad or absent state.
            if( x == 0.0 && y == 0.0 && z == 0.0 )
            {
                x = TUDAT_NAN;
                y = TUDAT_NAN;
                z = TUDAT_NAN;
            }

            // Convert raw SP3 positions from km and velocities from dm/s to SI units.
            const int indexOffset = ( firstCharacter == 'P' ) ? 0 : 3;
            const double unitScaling = ( firstCharacter == 'P' ) ? 1000.0 : 0.1;

            if( firstCharacter == 'P' )
            {
                if( currentRecordFlags[ satelliteId ].hasPosition )
                {
                    throw std::runtime_error( "Error when reading SP3 file: duplicate position record for satellite " + satelliteId +
                                              " in " + fileName );
                }
                currentRecordFlags[ satelliteId ].hasPosition = true;
            }
            else
            {
                if( !velocityRecordsExpected )
                {
                    throw std::runtime_error( "Error when reading SP3 file: velocity record encountered in a position-only file in " +
                                              fileName );
                }
                if( currentRecordFlags[ satelliteId ].hasVelocity )
                {
                    throw std::runtime_error( "Error when reading SP3 file: duplicate velocity record for satellite " + satelliteId +
                                              " in " + fileName );
                }
                currentRecordFlags[ satelliteId ].hasVelocity = true;
            }

            currentStates[ satelliteId ]( indexOffset + 0 ) = x * unitScaling;
            currentStates[ satelliteId ]( indexOffset + 1 ) = y * unitScaling;
            currentStates[ satelliteId ]( indexOffset + 2 ) = z * unitScaling;
            continue;
        }

        // Ignore other SP3 records that are not currently used.
        if( firstCharacter == '%' || firstCharacter == '/' || firstCharacter == 'E' || firstCharacter == '+' )
        {
            continue;
        }
    }

    if( !headerSeen )
    {
        throw std::runtime_error( "Error when reading SP3 file: first header line is missing in " + fileName );
    }
    if( !hasCurrentEpoch )
    {
        throw std::runtime_error( "Error when reading SP3 file: no epochs were found in " + fileName );
    }

    if( hasCurrentEpoch )
    {
        // EOF does not start another epoch, so explicitly store the final accumulated one.
        validateAndFlushCurrentEpoch(
                currentTime, fileName, knownSatelliteIds, velocityRecordsExpected, currentStates, currentRecordFlags, *fileContents );
    }

    // Cross-check the complete file against its header declarations after continuation records are known.
    if( !satelliteCountSeen )
    {
        throw std::runtime_error( "Error when reading SP3 file: satellite-list header is missing in " + fileName );
    }
    if( fileContents->declaredNumberOfSatellites != static_cast< int >( knownSatelliteIds.size( ) ) )
    {
        throw std::runtime_error( "Error when reading SP3 file: parsed satellite count (" + std::to_string( knownSatelliteIds.size( ) ) +
                                  ") does not match declared satellite count (" +
                                  std::to_string( fileContents->declaredNumberOfSatellites ) + ") in " + fileName );
    }
    if( fileContents->timeScale.empty( ) )
    {
        throw std::runtime_error( "Error when reading SP3 file: time-system code is missing in " + fileName );
    }
    static const std::set< std::string > supportedTimeSystems = { "GPS", "GLO", "GAL", "BDT", "TAI", "UTC", "IRN", "QZS" };
    if( supportedTimeSystems.count( fileContents->timeScale ) == 0 )
    {
        throw std::runtime_error( "Error when reading SP3 file: unsupported time-system code '" + fileContents->timeScale + "' in " +
                                  fileName );
    }

    if( fileContents->declaredNumberOfEpochs > 0 && parsedEpochCount != fileContents->declaredNumberOfEpochs )
    {
        throw std::runtime_error( "Error when reading SP3 file: parsed epoch count (" + std::to_string( parsedEpochCount ) +
                                  ") does not match declared epoch count (" + std::to_string( fileContents->declaredNumberOfEpochs ) +
                                  ") in " + fileName );
    }

    if( !fileContents->hasVelocityRecords )
    {
        // Position-only products become complete state histories at the data-input boundary.
        deriveVelocitiesFromPositions( *fileContents, fileName );
    }

    return fileContents;
}

}  // namespace input_output

}  // namespace tudat
