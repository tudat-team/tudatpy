/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readSinexFile.h"

#include <cmath>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{

namespace input_output
{

namespace
{

std::string trimCopy( const std::string& value )
{
    std::string trimmed = value;
    boost::algorithm::trim( trimmed );
    return trimmed;
}

bool tryParseInt( const std::string& value, int& output )
{
    try
    {
        output = boost::lexical_cast< int >( value );
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
        output = boost::lexical_cast< double >( value );
        return true;
    }
    catch( ... )
    {
        return false;
    }
}

std::vector< std::string > splitTokens( const std::string& line )
{
    std::vector< std::string > tokens;
    boost::algorithm::split( tokens, trimCopy( line ), boost::algorithm::is_any_of( " \t" ), boost::algorithm::token_compress_on );
    return tokens;
}

bool isStateEntry( const std::string& entry )
{
    return ( entry == "STAX" || entry == "STAY" || entry == "STAZ" || entry == "VELX" || entry == "VELY" || entry == "VELZ" );
}

bool parseSiteIdStyleEntry( const std::vector< std::string >& tokens, int& stationCode, std::string& domesId )
{
    if( tokens.size( ) < 3 )
    {
        return false;
    }

    if( !tryParseInt( tokens.at( 0 ), stationCode ) )
    {
        return false;
    }

    domesId = tokens.at( 2 );
    if( domesId.empty( ) || domesId == "---" )
    {
        return false;
    }

    return true;
}

int getStateIndex( const std::string& stateType )
{
    if( stateType == "STAX" )
    {
        return 0;
    }
    if( stateType == "STAY" )
    {
        return 1;
    }
    if( stateType == "STAZ" )
    {
        return 2;
    }
    if( stateType == "VELX" )
    {
        return 3;
    }
    if( stateType == "VELY" )
    {
        return 4;
    }
    if( stateType == "VELZ" )
    {
        return 5;
    }
    return -1;
}

std::string extractDomesIdFromSiteIdLine( const std::string& line )
{
    if( line.size( ) >= 16 )
    {
        // SINEX SITE/ID DOMES id starts at column 8 (1-based), i.e. index 7.
        std::string domesId = trimCopy( line.substr( 7, 9 ) );
        if( !domesId.empty( ) )
        {
            return domesId;
        }
    }
    return "";
}

void parseSiteIdBlockLine( const std::string& line, std::map< std::string, std::string >& siteCodeToDomes, std::map< int, std::string >& siteIndexToDomes )
{
    std::string siteCode;
    std::string domesId;

    if( line.size( ) >= 17 )
    {
        siteCode = trimCopy( line.substr( 0, 4 ) );
        domesId = extractDomesIdFromSiteIdLine( line );
    }

    const std::vector< std::string > tokens = splitTokens( line );
    if( tokens.size( ) >= 3 && !tokens.at( 2 ).empty( ) && tokens.at( 2 ) != "---" )
    {
        domesId = tokens.at( 2 );
    }

    if( siteCode.empty( ) && !tokens.empty( ) )
    {
        siteCode = tokens.at( 0 );
    }

    if( !siteCode.empty( ) && !domesId.empty( ) )
    {
        siteCodeToDomes[ siteCode ] = domesId;
    }

    if( !tokens.empty( ) )
    {
        int siteIndex = -1;
        if( tryParseInt( tokens.at( 0 ), siteIndex ) && !domesId.empty( ) )
        {
            siteIndexToDomes[ siteIndex ] = domesId;
        }
    }
}

double parseSinexEpochWithOpenEnd( const std::string& epochString, const double referenceJulianDay, bool& hasOpenEnd )
{
    hasOpenEnd = false;
    const std::string trimmedEpoch = trimCopy( epochString );
    if( trimmedEpoch.empty( ) )
    {
        return TUDAT_NAN;
    }
    if( trimmedEpoch == "00:000:00000" )
    {
        hasOpenEnd = true;
        return TUDAT_NAN;
    }
    return convertSinexDateTimeToSecondsSinceEpoch( trimmedEpoch, referenceJulianDay );
}

double convertDmsToRadians( const double degrees, const double arcMinutes, const double arcSeconds )
{
    const double sign = ( degrees < 0.0 ) ? -1.0 : 1.0;
    return sign * ( std::fabs( degrees ) + arcMinutes / 60.0 + arcSeconds / 3600.0 ) *
            mathematical_constants::PI / 180.0;
}

double normalizeLongitude( const double longitude )
{
    double normalizedLongitude = std::fmod( longitude + mathematical_constants::PI, 2.0 * mathematical_constants::PI );
    if( normalizedLongitude < 0.0 )
    {
        normalizedLongitude += 2.0 * mathematical_constants::PI;
    }
    return normalizedLongitude - mathematical_constants::PI;
}

bool parseDmsTriplet(
        const std::vector< std::string >& tokens, const unsigned int firstIndex, double& degrees, double& arcMinutes, double& arcSeconds )
{
    if( tokens.size( ) <= firstIndex + 2 )
    {
        return false;
    }
    return tryParseDouble( tokens.at( firstIndex ), degrees ) && tryParseDouble( tokens.at( firstIndex + 1 ), arcMinutes ) &&
            tryParseDouble( tokens.at( firstIndex + 2 ), arcSeconds );
}

bool parseApproximateCoordinatesFromFixedWidth(
        const std::string& line, double& longitude, double& latitude, double& height )
{
    if( line.size( ) < 77 )
    {
        return false;
    }

    const std::vector< std::string > longitudeTokens = splitTokens( line.substr( 43, 12 ) );
    const std::vector< std::string > latitudeTokens = splitTokens( line.substr( 56, 12 ) );
    const std::string heightToken = trimCopy( line.substr( 69, 8 ) );

    double longitudeDegrees = TUDAT_NAN;
    double longitudeMinutes = TUDAT_NAN;
    double longitudeSeconds = TUDAT_NAN;
    double latitudeDegrees = TUDAT_NAN;
    double latitudeMinutes = TUDAT_NAN;
    double latitudeSeconds = TUDAT_NAN;
    if( !parseDmsTriplet( longitudeTokens, 0, longitudeDegrees, longitudeMinutes, longitudeSeconds ) ||
        !parseDmsTriplet( latitudeTokens, 0, latitudeDegrees, latitudeMinutes, latitudeSeconds ) ||
        !tryParseDouble( heightToken, height ) )
    {
        return false;
    }

    longitude = normalizeLongitude( convertDmsToRadians( longitudeDegrees, longitudeMinutes, longitudeSeconds ) );
    latitude = convertDmsToRadians( latitudeDegrees, latitudeMinutes, latitudeSeconds );
    return true;
}

bool parseApproximateCoordinatesFromTokens(
        const std::vector< std::string >& tokens, double& longitude, double& latitude, double& height )
{
    if( tokens.size( ) < 11 )
    {
        return false;
    }

    int tailShift = 0;
    double maybeSod = TUDAT_NAN;
    if( tryParseDouble( tokens.back( ), maybeSod ) && std::fabs( maybeSod ) > 1.0E4 )
    {
        tailShift = 1;
    }

    const int heightIndex = static_cast< int >( tokens.size( ) ) - 1 - tailShift;
    if( heightIndex < 6 )
    {
        return false;
    }

    double longitudeDegrees = TUDAT_NAN;
    double longitudeMinutes = TUDAT_NAN;
    double longitudeSeconds = TUDAT_NAN;
    double latitudeDegrees = TUDAT_NAN;
    double latitudeMinutes = TUDAT_NAN;
    double latitudeSeconds = TUDAT_NAN;
    if( !tryParseDouble( tokens.at( static_cast< unsigned int >( heightIndex ) ), height ) ||
        !tryParseDouble( tokens.at( static_cast< unsigned int >( heightIndex - 1 ) ), latitudeSeconds ) ||
        !tryParseDouble( tokens.at( static_cast< unsigned int >( heightIndex - 2 ) ), latitudeMinutes ) ||
        !tryParseDouble( tokens.at( static_cast< unsigned int >( heightIndex - 3 ) ), latitudeDegrees ) ||
        !tryParseDouble( tokens.at( static_cast< unsigned int >( heightIndex - 4 ) ), longitudeSeconds ) ||
        !tryParseDouble( tokens.at( static_cast< unsigned int >( heightIndex - 5 ) ), longitudeMinutes ) ||
        !tryParseDouble( tokens.at( static_cast< unsigned int >( heightIndex - 6 ) ), longitudeDegrees ) )
    {
        return false;
    }

    longitude = normalizeLongitude( convertDmsToRadians( longitudeDegrees, longitudeMinutes, longitudeSeconds ) );
    latitude = convertDmsToRadians( latitudeDegrees, latitudeMinutes, latitudeSeconds );
    return true;
}

std::string extractStationNameFromSiteIdLine( const std::string& line, const std::vector< std::string >& tokens )
{
    std::string stationName;
    if( line.size( ) >= 42 )
    {
        stationName = trimCopy( line.substr( 20, 22 ) );
    }
    if( stationName.empty( ) && tokens.size( ) > 4 )
    {
        stationName = tokens.at( 4 );
    }

    std::vector< std::string > nameTokens = splitTokens( stationName );
    return ( nameTokens.empty( ) ? stationName : nameTokens.at( 0 ) );
}

}  // namespace

double convertSinexDateTimeToSecondsSinceEpoch( const std::string& dateTime, const double referenceJulianDay )
{
    if( dateTime.size( ) < 11 )
    {
        throw std::runtime_error( "Error when converting SINEX epoch: expected format YY:DDD:SSSSS, got " + dateTime + "." );
    }

    int year = boost::lexical_cast< int >( dateTime.substr( 0, 2 ) );
    if( year < 50 )
    {
        year += 2000;
    }
    else
    {
        year += 1900;
    }

    int dayOfYear = boost::lexical_cast< int >( dateTime.substr( 3, 3 ) );
    double secondsInDay = boost::lexical_cast< double >( dateTime.substr( 7, 5 ) );
    boost::gregorian::date date = basic_astrodynamics::convertYearAndDaysInYearToDate( year, dayOfYear - 1 );
    return basic_astrodynamics::calculateJulianDaySinceEpoch( date, secondsInDay / physical_constants::JULIAN_DAY, referenceJulianDay ) *
            physical_constants::JULIAN_DAY;
}

std::map< std::string, SinexStationState > readSinexStationData( const std::string& fileName, const double referenceJulianDay )
{
    std::ifstream stream( fileName, std::ios_base::in );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening SINEX file: " + fileName );
    }

    std::map< std::string, std::string > siteCodeToDomes;
    std::map< int, std::string > siteIndexToDomes;
    std::map< std::string, Eigen::Matrix< double, 6, 1 > > siteStates;
    std::map< std::string, double > siteReferenceEpochs;

    bool isInSiteIdBlock = false;
    bool isInSolutionEstimateBlock = false;

    std::string line;
    while( std::getline( stream, line ) )
    {
        std::string trimmedLine = trimCopy( line );
        if( trimmedLine.empty( ) )
        {
            continue;
        }

        if( trimmedLine == "+SITE/ID" )
        {
            isInSiteIdBlock = true;
            continue;
        }
        if( trimmedLine == "-SITE/ID" )
        {
            isInSiteIdBlock = false;
            continue;
        }
        if( trimmedLine == "+SOLUTION/ESTIMATE" )
        {
            isInSolutionEstimateBlock = true;
            continue;
        }
        if( trimmedLine == "-SOLUTION/ESTIMATE" )
        {
            isInSolutionEstimateBlock = false;
            continue;
        }

        if( trimmedLine.at( 0 ) == '*' )
        {
            continue;
        }

        if( isInSiteIdBlock )
        {
            parseSiteIdBlockLine( line, siteCodeToDomes, siteIndexToDomes );
            continue;
        }

        if( isInSolutionEstimateBlock )
        {
            std::string stateType;
            std::string siteIdentifier;
            std::string referenceEpoch;
            double stateValue = TUDAT_NAN;
            bool hasParsedValue = false;

            if( line.size( ) >= 68 )
            {
                stateType = trimCopy( line.substr( 6, 6 ) );
                siteIdentifier = trimCopy( line.substr( 13, 4 ) );
                referenceEpoch = trimCopy( line.substr( 26, 12 ) );
                hasParsedValue = tryParseDouble( trimCopy( line.substr( 47, 21 ) ), stateValue );
            }

            if( !isStateEntry( stateType ) || !hasParsedValue )
            {
                std::vector< std::string > tokens;
                boost::algorithm::split(
                        tokens, trimmedLine, boost::algorithm::is_any_of( " \t" ), boost::algorithm::token_compress_on );

                int stateTokenIndex = -1;
                for( unsigned int i = 0; i < tokens.size( ); i++ )
                {
                    if( isStateEntry( tokens.at( i ) ) )
                    {
                        stateTokenIndex = static_cast< int >( i );
                        break;
                    }
                }

                if( stateTokenIndex < 0 || tokens.size( ) <= static_cast< unsigned int >( stateTokenIndex + 1 ) )
                {
                    continue;
                }

                stateType = tokens.at( stateTokenIndex );
                siteIdentifier = tokens.at( stateTokenIndex + 1 );

                if( tokens.size( ) > static_cast< unsigned int >( stateTokenIndex + 4 ) )
                {
                    referenceEpoch = tokens.at( stateTokenIndex + 4 );
                }

                int parsedFromEnd = 0;
                for( int i = static_cast< int >( tokens.size( ) ) - 1; i >= 0; i-- )
                {
                    double parsedValue = TUDAT_NAN;
                    if( tryParseDouble( tokens.at( static_cast< unsigned int >( i ) ), parsedValue ) )
                    {
                        parsedFromEnd++;
                        if( parsedFromEnd == 2 )
                        {
                            stateValue = parsedValue;
                            hasParsedValue = true;
                            break;
                        }
                    }
                }
            }

            if( !isStateEntry( stateType ) || !hasParsedValue )
            {
                continue;
            }

            std::string domesId;
            if( siteCodeToDomes.count( siteIdentifier ) > 0 )
            {
                domesId = siteCodeToDomes.at( siteIdentifier );
            }
            else
            {
                int siteIndex = -1;
                if( tryParseInt( siteIdentifier, siteIndex ) && siteIndexToDomes.count( siteIndex ) > 0 )
                {
                    domesId = siteIndexToDomes.at( siteIndex );
                }
                else
                {
                    domesId = siteIdentifier;
                }
            }

            if( siteStates.count( domesId ) == 0 )
            {
                siteStates[ domesId ] = Eigen::Matrix< double, 6, 1 >::Constant( TUDAT_NAN );
            }

            int stateIndex = getStateIndex( stateType );
            if( stateIndex < 0 )
            {
                continue;
            }
            siteStates[ domesId ]( stateIndex ) = stateValue;

            if( !referenceEpoch.empty( ) && referenceEpoch.find( ':' ) != std::string::npos )
            {
                siteReferenceEpochs[ domesId ] = convertSinexDateTimeToSecondsSinceEpoch( referenceEpoch, referenceJulianDay );
            }
        }
    }

    std::map< std::string, SinexStationState > stationStates;
    for( const auto& stateEntry: siteStates )
    {
        SinexStationState currentState;
        currentState.domesId_ = stateEntry.first;
        currentState.position_ = stateEntry.second.segment( 0, 3 );
        currentState.velocity_ = stateEntry.second.segment( 3, 3 );
        currentState.referenceEpoch_ =
                ( siteReferenceEpochs.count( stateEntry.first ) > 0 ) ? siteReferenceEpochs.at( stateEntry.first ) : TUDAT_NAN;

        if( currentState.velocity_( 0 ) == currentState.velocity_( 0 ) )
        {
            currentState.velocity_ /= physical_constants::JULIAN_YEAR;
        }

        stationStates[ stateEntry.first ] = currentState;
    }

    return stationStates;
}

std::map< std::string, std::vector< SinexStationEccentricity > > readSinexStationEccentricities(
        const std::string& fileName,
        const double referenceJulianDay )
{
    std::ifstream stream( fileName, std::ios_base::in );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening SINEX eccentricity file: " + fileName );
    }

    std::map< std::string, std::string > siteCodeToDomes;
    std::map< int, std::string > siteIndexToDomes;
    std::map< std::string, std::vector< SinexStationEccentricity > > eccentricityData;

    bool isInSiteIdBlock = false;
    bool isInSiteEccentricityBlock = false;

    std::string line;
    while( std::getline( stream, line ) )
    {
        const std::string trimmedLine = trimCopy( line );
        if( trimmedLine.empty( ) )
        {
            continue;
        }

        if( trimmedLine == "+SITE/ID" )
        {
            isInSiteIdBlock = true;
            continue;
        }
        if( trimmedLine == "-SITE/ID" )
        {
            isInSiteIdBlock = false;
            continue;
        }
        if( trimmedLine == "+SITE/ECCENTRICITY" )
        {
            isInSiteEccentricityBlock = true;
            continue;
        }
        if( trimmedLine == "-SITE/ECCENTRICITY" )
        {
            isInSiteEccentricityBlock = false;
            continue;
        }
        if( trimmedLine.at( 0 ) == '*' )
        {
            continue;
        }

        if( isInSiteIdBlock )
        {
            parseSiteIdBlockLine( line, siteCodeToDomes, siteIndexToDomes );
            continue;
        }

        if( !isInSiteEccentricityBlock )
        {
            continue;
        }

        const std::vector< std::string > tokens = splitTokens( line );
        if( tokens.size( ) < 10 )
        {
            continue;
        }

        int stationCode = -1;
        if( !tryParseInt( tokens.at( 0 ), stationCode ) )
        {
            continue;
        }

        const std::string coordinateType = tokens.at( 6 );
        if( coordinateType != "XYZ" )
        {
            continue;
        }

        double xEccentricity = 0.0;
        double yEccentricity = 0.0;
        double zEccentricity = 0.0;
        if( !tryParseDouble( tokens.at( 7 ), xEccentricity ) || !tryParseDouble( tokens.at( 8 ), yEccentricity ) ||
            !tryParseDouble( tokens.at( 9 ), zEccentricity ) )
        {
            continue;
        }

        std::string domesId;
        if( siteIndexToDomes.count( stationCode ) > 0 )
        {
            domesId = siteIndexToDomes.at( stationCode );
        }
        else if( siteCodeToDomes.count( std::to_string( stationCode ) ) > 0 )
        {
            domesId = siteCodeToDomes.at( std::to_string( stationCode ) );
        }
        else
        {
            continue;
        }

        SinexStationEccentricity eccentricityEntry;
        eccentricityEntry.domesId_ = domesId;
        eccentricityEntry.stationCode_ = stationCode;
        eccentricityEntry.eccentricity_ = Eigen::Vector3d( xEccentricity, yEccentricity, zEccentricity );

        bool hasOpenEnd = false;
        bool unusedStartOpenFlag = false;
        eccentricityEntry.startEpoch_ = parseSinexEpochWithOpenEnd( tokens.at( 4 ), referenceJulianDay, unusedStartOpenFlag );
        eccentricityEntry.endEpoch_ = parseSinexEpochWithOpenEnd( tokens.at( 5 ), referenceJulianDay, hasOpenEnd );
        eccentricityEntry.hasOpenEnd_ = hasOpenEnd;

        eccentricityData[ domesId ].push_back( eccentricityEntry );
    }

    return eccentricityData;
}

std::map< int, IlrsStationRegistryEntry > readIlrsStationRegistryFromSinexSiteId( const std::string& fileName )
{
    std::ifstream stream( fileName, std::ios_base::in );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening SINEX SITE/ID file: " + fileName );
    }

    std::map< int, IlrsStationRegistryEntry > stationRegistry;
    bool isInSiteIdBlock = false;

    std::string line;
    while( std::getline( stream, line ) )
    {
        const std::string trimmedLine = trimCopy( line );
        if( trimmedLine.empty( ) )
        {
            continue;
        }
        if( trimmedLine == "+SITE/ID" )
        {
            isInSiteIdBlock = true;
            continue;
        }
        if( trimmedLine == "-SITE/ID" )
        {
            isInSiteIdBlock = false;
            continue;
        }
        if( !isInSiteIdBlock || trimmedLine.at( 0 ) == '*' || trimmedLine.at( 0 ) == '#' || trimmedLine.at( 0 ) == '%' )
        {
            continue;
        }

        const std::vector< std::string > tokens = splitTokens( line );
        if( tokens.size( ) < 3 )
        {
            continue;
        }

        int stationCode = -1;
        if( !tryParseInt( tokens.at( 0 ), stationCode ) )
        {
            continue;
        }

        IlrsStationRegistryEntry registryEntry;
        registryEntry.stationCode_ = stationCode;
        registryEntry.stationName_ = extractStationNameFromSiteIdLine( line, tokens );
        registryEntry.domesId_ = tokens.at( 2 );

        double longitude = TUDAT_NAN;
        double latitude = TUDAT_NAN;
        double height = TUDAT_NAN;
        const bool hasParsedCoordinates = parseApproximateCoordinatesFromTokens( tokens, longitude, latitude, height ) ||
                                          parseApproximateCoordinatesFromFixedWidth( line, longitude, latitude, height );
        if( hasParsedCoordinates )
        {
            registryEntry.approximateLongitude_ = longitude;
            registryEntry.approximateLatitude_ = latitude;
            registryEntry.approximateHeight_ = height;
        }

        stationRegistry[ stationCode ] = registryEntry;
    }

    return stationRegistry;
}

std::map< std::string, std::string > readDomesIdNumbers( const std::string& fileName )
{
    std::ifstream stream( fileName, std::ios_base::in );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening station DOMES file: " + fileName );
    }

    std::map< std::string, std::string > stationNameMap;
    bool hasSiteIdBlock = false;
    bool isInSiteIdBlock = false;
    std::string line;
    while( std::getline( stream, line ) )
    {
        std::string trimmed = trimCopy( line );
        if( trimmed.empty( ) )
        {
            continue;
        }

        if( trimmed == "+SITE/ID" )
        {
            hasSiteIdBlock = true;
            isInSiteIdBlock = true;
            continue;
        }
        if( trimmed == "-SITE/ID" )
        {
            isInSiteIdBlock = false;
            continue;
        }

        if( hasSiteIdBlock && !isInSiteIdBlock )
        {
            continue;
        }
        if( trimmed.at( 0 ) == '*' || trimmed.at( 0 ) == '#' || trimmed.at( 0 ) == '%' )
        {
            continue;
        }

        std::string stationName;
        std::string domesId;

        std::vector< std::string > tokens;
        boost::algorithm::split( tokens, trimmed, boost::algorithm::is_any_of( " \t" ), boost::algorithm::token_compress_on );
        int stationCode = -1;
        if( parseSiteIdStyleEntry( tokens, stationCode, domesId ) )
        {
            stationNameMap[ std::to_string( stationCode ) ] = domesId;
            continue;
        }

        if( line.size( ) >= 42 )
        {
            stationName = trimCopy( line.substr( 0, 24 ) );
            domesId = trimCopy( line.substr( 33, 9 ) );
        }

        if( stationName.empty( ) && !tokens.empty( ) )
        {
            stationName = tokens.at( 0 );
        }
        if( domesId.empty( ) && tokens.size( ) > 1 )
        {
            domesId = tokens.at( 1 );
        }

        if( !stationName.empty( ) && !domesId.empty( ) )
        {
            stationNameMap[ stationName ] = domesId;
        }
    }

    return stationNameMap;
}

std::map< int, std::string > readMonumentNumbers( const std::string& fileName )
{
    std::ifstream stream( fileName, std::ios_base::in );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening station monument file: " + fileName );
    }

    std::map< int, std::string > stationNameMap;
    bool hasSiteIdBlock = false;
    bool isInSiteIdBlock = false;
    std::string line;
    while( std::getline( stream, line ) )
    {
        std::string trimmed = trimCopy( line );
        if( trimmed.empty( ) )
        {
            continue;
        }

        if( trimmed == "+SITE/ID" )
        {
            hasSiteIdBlock = true;
            isInSiteIdBlock = true;
            continue;
        }
        if( trimmed == "-SITE/ID" )
        {
            isInSiteIdBlock = false;
            continue;
        }

        if( hasSiteIdBlock && !isInSiteIdBlock )
        {
            continue;
        }
        if( trimmed.at( 0 ) == '*' || trimmed.at( 0 ) == '#' || trimmed.at( 0 ) == '%' )
        {
            continue;
        }

        std::vector< std::string > tokens;
        boost::algorithm::split(
                tokens, trimmed, boost::algorithm::is_any_of( " \t" ), boost::algorithm::token_compress_on );

        int stationCode = -1;
        std::string domesId;
        if( parseSiteIdStyleEntry( tokens, stationCode, domesId ) )
        {
            stationNameMap[ stationCode ] = std::to_string( stationCode );
            continue;
        }

        std::string stationName;
        int monumentId = -1;
        bool hasMonumentId = false;

        if( line.size( ) >= 52 )
        {
            stationName = trimCopy( line.substr( 0, 24 ) );
            hasMonumentId = tryParseInt( trimCopy( line.substr( 48, 4 ) ), monumentId );
        }

        if( stationName.empty( ) || !hasMonumentId )
        {
            if( tokens.size( ) >= 2 )
            {
                stationName = tokens.at( 0 );
                hasMonumentId = tryParseInt( tokens.back( ), monumentId );
            }
        }

        if( !stationName.empty( ) && hasMonumentId )
        {
            stationNameMap[ monumentId ] = stationName;
        }
    }

    return stationNameMap;
}

std::map< std::string, std::string > readGroundStationNames( const std::string& fileName )
{
    std::map< std::string, std::string > stationNameToDomes = readDomesIdNumbers( fileName );
    std::map< std::string, std::string > domesToStationName;

    for( const auto& stationEntry: stationNameToDomes )
    {
        domesToStationName[ stationEntry.second ] = stationEntry.first;
    }

    return domesToStationName;
}

}  // namespace input_output

}  // namespace tudat
