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

#include <array>
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

namespace
{

struct EpochRecordFlags {
    bool hasPosition = false;
    bool hasVelocity = false;
};

std::string getFixedWidthField( const std::string& line, const unsigned int startIndex, const unsigned int length )
{
    if( startIndex >= line.size( ) )
    {
        return "";
    }
    std::string field = line.substr( startIndex, std::min( length, static_cast< unsigned int >( line.size( ) - startIndex ) ) );
    boost::algorithm::trim( field );
    return field;
}

bool tryParseInt( const std::string& value, int& output )
{
    try
    {
        output = std::stoi( value );
        return true;
    }
    catch( ... )
    {
        return false;
    }
}

bool tryParseDouble( const std::string& value, double& output )
{
    try
    {
        output = std::stod( value );
        return true;
    }
    catch( ... )
    {
        return false;
    }
}

bool tryParseFixedWidthInt( const std::string& line, const unsigned int startIndex, const unsigned int length, int& output )
{
    return tryParseInt( getFixedWidthField( line, startIndex, length ), output );
}

bool tryParseFixedWidthDouble( const std::string& line, const unsigned int startIndex, const unsigned int length, double& output )
{
    return tryParseDouble( getFixedWidthField( line, startIndex, length ), output );
}

bool isMissingSp3State( const double x, const double y, const double z )
{
    // SP3 uses an all-zero position or velocity vector for a bad or absent state.
    return x == 0.0 && y == 0.0 && z == 0.0;
}

bool isSupportedSp3TimeSystem( const std::string& timeSystem )
{
    static const std::set< std::string > supportedTimeSystems = { "GPS", "GLO", "GAL", "BDT", "TAI", "UTC", "IRN", "QZS" };
    return supportedTimeSystems.count( timeSystem ) > 0;
}

void addSatelliteIdToStateMap( const std::string& rawSatelliteId,
                               Sp3FileContents& fileContents,
                               std::vector< std::string >& knownSatelliteIds )
{
    const std::string trimmedSatelliteId = boost::algorithm::trim_copy( rawSatelliteId );
    if( trimmedSatelliteId.empty( ) || trimmedSatelliteId == "0" || trimmedSatelliteId == "00" || trimmedSatelliteId == "000" )
    {
        return;
    }

    if( fileContents.satelliteStates.count( trimmedSatelliteId ) == 0 )
    {
        fileContents.satelliteStates[ trimmedSatelliteId ] = std::map< double, Eigen::VectorXd >( );
        knownSatelliteIds.push_back( trimmedSatelliteId );
    }
}

void validateAndFlushCurrentEpoch( const double currentTime,
                                   const std::string& fileName,
                                   const std::vector< std::string >& knownSatelliteIds,
                                   const bool velocityRecordsExpected,
                                   const std::map< std::string, Eigen::VectorXd >& currentStates,
                                   const std::map< std::string, EpochRecordFlags >& currentRecordFlags,
                                   Sp3FileContents& fileContents )
{
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

Eigen::Vector3d calculateThreePointDerivative( const std::array< double, 3 >& times,
                                               const std::array< Eigen::Vector3d, 3 >& positions,
                                               const double evaluationTime )
{
    Eigen::Vector3d derivative = Eigen::Vector3d::Zero( );
    for( unsigned int i = 0; i < 3; ++i )
    {
        const unsigned int j = ( i + 1 ) % 3;
        const unsigned int k = ( i + 2 ) % 3;
        const double denominator = ( times[ i ] - times[ j ] ) * ( times[ i ] - times[ k ] );
        if( denominator == 0.0 )
        {
            throw std::runtime_error( "Cannot derive SP3 velocities from duplicate epochs." );
        }
        derivative += positions[ i ] * ( 2.0 * evaluationTime - times[ j ] - times[ k ] ) / denominator;
    }
    return derivative;
}

void deriveVelocitiesFromPositions( Sp3FileContents& fileContents, const std::string& fileName )
{
    for( auto& satellite : fileContents.satelliteStates )
    {
        if( satellite.second.size( ) < 3 )
        {
            throw std::runtime_error(
                    "Error when reading position-only SP3 file: at least three epochs are required to derive "
                    "second-order velocities for satellite " +
                    satellite.first + " in " + fileName );
        }

        std::vector< std::map< double, Eigen::VectorXd >::iterator > states;
        states.reserve( satellite.second.size( ) );
        for( auto stateIterator = satellite.second.begin( ); stateIterator != satellite.second.end( ); ++stateIterator )
        {
            states.push_back( stateIterator );
        }

        for( unsigned int i = 0; i < states.size( ); ++i )
        {
            const unsigned int firstStencilIndex = ( i == 0 ) ? 0 : ( i == states.size( ) - 1 ? i - 2 : i - 1 );
            const std::array< double, 3 > times = { states[ firstStencilIndex ]->first,
                                                    states[ firstStencilIndex + 1 ]->first,
                                                    states[ firstStencilIndex + 2 ]->first };
            const std::array< Eigen::Vector3d, 3 > positions = { states[ firstStencilIndex ]->second.segment< 3 >( 0 ),
                                                                 states[ firstStencilIndex + 1 ]->second.segment< 3 >( 0 ),
                                                                 states[ firstStencilIndex + 2 ]->second.segment< 3 >( 0 ) };

            if( positions[ 0 ].allFinite( ) && positions[ 1 ].allFinite( ) && positions[ 2 ].allFinite( ) )
            {
                states[ i ]->second.segment< 3 >( 3 ) = calculateThreePointDerivative( times, positions, states[ i ]->first );
            }
        }
    }
    fileContents.velocitiesWereDerived = true;
}

}  // namespace

std::shared_ptr< Sp3FileContents > readSp3File( const std::string& fileName, const double referenceJulianDay )
{
    std::shared_ptr< Sp3FileContents > fileContents = std::make_shared< Sp3FileContents >( );
    fileContents->referenceJulianDay = referenceJulianDay;

    std::ifstream stream( fileName, std::ios_base::in );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening SP3 file: " + fileName );
    }

    std::string currentLine;

    std::vector< std::string > knownSatelliteIds;
    std::map< std::string, Eigen::VectorXd > currentStates;
    std::map< std::string, EpochRecordFlags > currentRecordFlags;

    double currentTime = TUDAT_NAN;
    bool hasCurrentEpoch = false;
    bool velocityRecordsExpected = false;
    bool headerSeen = false;
    bool satelliteCountSeen = false;
    int parsedEpochCount = 0;

    while( std::getline( stream, currentLine ) )
    {
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
            if( tryParseFixedWidthInt( currentLine, 32, 7, declaredEpochCount ) )
            {
                fileContents->declaredNumberOfEpochs = declaredEpochCount;
            }

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
            if( tryParseFixedWidthDouble( currentLine, 24, 14, declaredEpochInterval ) )
            {
                fileContents->declaredEpochInterval = declaredEpochInterval;
            }
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
                currentStates[ satelliteId ] = Eigen::VectorXd::Constant( 6, TUDAT_NAN );
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
                currentStates[ satelliteId ] = Eigen::VectorXd::Constant( 6, TUDAT_NAN );
                currentRecordFlags[ satelliteId ] = EpochRecordFlags( );
            }

            double x, y, z;
            const bool readSucceeded = tryParseFixedWidthDouble( currentLine, 4, 14, x ) &&
                    tryParseFixedWidthDouble( currentLine, 18, 14, y ) && tryParseFixedWidthDouble( currentLine, 32, 14, z );

            if( !readSucceeded )
            {
                throw std::runtime_error( "Error when reading SP3 file: invalid numeric state record in " + fileName );
            }

            if( isMissingSp3State( x, y, z ) )
            {
                x = TUDAT_NAN;
                y = TUDAT_NAN;
                z = TUDAT_NAN;
            }

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
        validateAndFlushCurrentEpoch(
                currentTime, fileName, knownSatelliteIds, velocityRecordsExpected, currentStates, currentRecordFlags, *fileContents );
    }

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
    if( !isSupportedSp3TimeSystem( fileContents->timeScale ) )
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
        deriveVelocitiesFromPositions( *fileContents, fileName );
    }

    return fileContents;
}

}  // namespace input_output

}  // namespace tudat
