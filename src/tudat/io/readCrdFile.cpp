/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readCrdFile.h"

#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"

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

std::vector< std::string > splitLine( const std::string& line )
{
    std::vector< std::string > tokens;
    std::string trimmedLine = trimCopy( line );
    boost::algorithm::split( tokens, trimmedLine, boost::algorithm::is_any_of( " \t" ), boost::algorithm::token_compress_on );
    return tokens;
}

bool passHasAnyData( const CrdPass& pass )
{
    return ( !pass.data_.fullRateData_.empty( ) || !pass.data_.normalPointData_.empty( ) || !pass.data_.meteorologicalData_.empty( ) ||
             pass.configuration_.cdpPadId_ >= 0 );
}

}  // namespace

double convertCrdTwoWayTimeOfFlightToSlrRange( const double twoWayTimeOfFlight )
{
    return 0.5 * physical_constants::SPEED_OF_LIGHT * twoWayTimeOfFlight;
}

std::vector< CrdPass > readCrdFile( const std::string& fileName )
{
    std::ifstream stream( fileName, std::ios_base::in );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening CRD file: " + fileName );
    }

    std::vector< CrdPass > fileContents;
    CrdPass currentPass;

    std::string line;
    while( std::getline( stream, line ) )
    {
        std::vector< std::string > tokens = splitLine( line );
        if( tokens.empty( ) )
        {
            continue;
        }

        const std::string lineId = tokens.at( 0 );
        if( lineId == "H8" || lineId == "h8" )
        {
            if( passHasAnyData( currentPass ) )
            {
                fileContents.push_back( currentPass );
            }
            currentPass = CrdPass( );
            continue;
        }
        if( lineId == "H9" || lineId == "h9" )
        {
            if( passHasAnyData( currentPass ) )
            {
                fileContents.push_back( currentPass );
            }
            currentPass = CrdPass( );
            break;
        }

        if( lineId == "H2" || lineId == "h2" )
        {
            if( tokens.size( ) > 1 )
            {
                currentPass.configuration_.stationName_ = tokens.at( 1 );
            }
            if( tokens.size( ) > 2 )
            {
                tryParseInt( tokens.at( 2 ), currentPass.configuration_.cdpPadId_ );
            }
            continue;
        }
        if( lineId == "H3" || lineId == "h3" )
        {
            if( tokens.size( ) > 1 )
            {
                currentPass.configuration_.targetName_ = tokens.at( 1 );
            }
            continue;
        }
        if( lineId == "H4" || lineId == "h4" )
        {
            if( tokens.size( ) > 7 )
            {
                tryParseInt( tokens.at( 2 ), currentPass.configuration_.startYear_ );
                tryParseInt( tokens.at( 3 ), currentPass.configuration_.startMonth_ );
                tryParseInt( tokens.at( 4 ), currentPass.configuration_.startDay_ );
            }
            if( tokens.size( ) > 13 )
            {
                tryParseInt( tokens.at( 8 ), currentPass.configuration_.endYear_ );
                tryParseInt( tokens.at( 9 ), currentPass.configuration_.endMonth_ );
                tryParseInt( tokens.at( 10 ), currentPass.configuration_.endDay_ );
            }
            continue;
        }
        if( lineId == "C0" || lineId == "c0" )
        {
            if( tokens.size( ) > 2 )
            {
                tryParseDouble( tokens.at( 2 ), currentPass.configuration_.transmitWavelengthNm_ );
            }
            continue;
        }
        if( lineId == "11" )
        {
            if( tokens.size( ) < 3 )
            {
                continue;
            }

            CrdNormalPointRecord currentRecord;
            tryParseDouble( tokens.at( 1 ), currentRecord.secondOfDay_ );
            tryParseDouble( tokens.at( 2 ), currentRecord.twoWayTimeOfFlight_ );
            currentRecord.oneWayRange_ = convertCrdTwoWayTimeOfFlightToSlrRange( currentRecord.twoWayTimeOfFlight_ );
            if( tokens.size( ) > 3 )
            {
                currentRecord.systemConfigurationId_ = tokens.at( 3 );
            }
            if( tokens.size( ) > 4 )
            {
                tryParseInt( tokens.at( 4 ), currentRecord.epochEvent_ );
            }
            if( tokens.size( ) > 5 )
            {
                tryParseDouble( tokens.at( 5 ), currentRecord.normalPointWindowLength_ );
            }
            if( tokens.size( ) > 6 )
            {
                tryParseInt( tokens.at( 6 ), currentRecord.numberOfReturns_ );
            }
            if( tokens.size( ) > 7 )
            {
                tryParseDouble( tokens.at( 7 ), currentRecord.binRms_ );
            }

            currentPass.data_.normalPointData_.push_back( currentRecord );
            continue;
        }
        if( lineId == "10" )
        {
            if( tokens.size( ) < 3 )
            {
                continue;
            }

            CrdFullRateRecord currentRecord;
            tryParseDouble( tokens.at( 1 ), currentRecord.secondOfDay_ );
            tryParseDouble( tokens.at( 2 ), currentRecord.twoWayTimeOfFlight_ );
            currentRecord.oneWayRange_ = convertCrdTwoWayTimeOfFlightToSlrRange( currentRecord.twoWayTimeOfFlight_ );
            if( tokens.size( ) > 3 )
            {
                currentRecord.systemConfigurationId_ = tokens.at( 3 );
            }
            if( tokens.size( ) > 4 )
            {
                tryParseInt( tokens.at( 4 ), currentRecord.epochEvent_ );
            }
            if( tokens.size( ) > 5 )
            {
                tryParseInt( tokens.at( 5 ), currentRecord.filterFlag_ );
            }
            if( tokens.size( ) > 6 )
            {
                tryParseInt( tokens.at( 6 ), currentRecord.detectorChannel_ );
            }
            if( tokens.size( ) > 7 )
            {
                tryParseInt( tokens.at( 7 ), currentRecord.stopNumber_ );
            }

            currentPass.data_.fullRateData_.push_back( currentRecord );
            continue;
        }
        if( lineId == "20" )
        {
            if( tokens.size( ) < 5 )
            {
                continue;
            }

            CrdMeteoRecord currentRecord;
            tryParseDouble( tokens.at( 1 ), currentRecord.secondOfDay_ );
            tryParseDouble( tokens.at( 2 ), currentRecord.pressure_ );
            tryParseDouble( tokens.at( 3 ), currentRecord.temperature_ );
            tryParseDouble( tokens.at( 4 ), currentRecord.humidity_ );
            currentPass.data_.meteorologicalData_.push_back( currentRecord );
            continue;
        }
    }

    if( passHasAnyData( currentPass ) )
    {
        fileContents.push_back( currentPass );
    }

    return fileContents;
}

std::vector< CrdPass > readCrdFiles( const std::vector< std::string >& fileNames )
{
    std::vector< CrdPass > allData;
    for( const std::string& fileName: fileNames )
    {
        std::vector< CrdPass > currentFileData = readCrdFile( fileName );
        allData.insert( allData.end( ), currentFileData.begin( ), currentFileData.end( ) );
    }
    return allData;
}

std::map< std::string, std::vector< CrdPass > > groupCrdDataPerTarget( const std::vector< CrdPass >& crdPasses )
{
    std::map< std::string, std::vector< CrdPass > > groupedData;
    for( const CrdPass& pass: crdPasses )
    {
        if( !pass.configuration_.targetName_.empty( ) )
        {
            groupedData[ pass.configuration_.targetName_ ].push_back( pass );
        }
    }
    return groupedData;
}

std::map< std::string, std::vector< CrdPass > > groupCrdDataPerStation(
        const std::vector< CrdPass >& crdPasses,
        const std::map< int, std::string >& monumentIdToGroundStationNameMap )
{
    std::map< std::string, std::vector< CrdPass > > groupedData;
    for( const CrdPass& pass: crdPasses )
    {
        if( monumentIdToGroundStationNameMap.count( pass.configuration_.cdpPadId_ ) > 0 )
        {
            groupedData[ monumentIdToGroundStationNameMap.at( pass.configuration_.cdpPadId_ ) ].push_back( pass );
        }
    }
    return groupedData;
}

std::map< double, double > extractNormalPointMeasurements( const CrdPass& passData )
{
    if( passData.configuration_.startYear_ < 0 || passData.configuration_.startMonth_ < 0 || passData.configuration_.startDay_ < 0 )
    {
        throw std::runtime_error(
                "Error when extracting CRD normal-point measurements: pass start calendar date is not available." );
    }

    std::map< double, double > normalPointData;
    const double startOfDayInSeconds =
            basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch( passData.configuration_.startYear_,
                                                                             passData.configuration_.startMonth_,
                                                                             passData.configuration_.startDay_,
                                                                             0,
                                                                             0,
                                                                             0.0,
                                                                             basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    double currentDayStart = startOfDayInSeconds;
    double previousSecondOfDay = TUDAT_NAN;

    for( const CrdNormalPointRecord& point: passData.data_.normalPointData_ )
    {
        if( previousSecondOfDay == previousSecondOfDay && point.secondOfDay_ + 1.0E4 < previousSecondOfDay )
        {
            currentDayStart += physical_constants::JULIAN_DAY;
        }

        normalPointData[ currentDayStart + point.secondOfDay_ ] = point.oneWayRange_;
        previousSecondOfDay = point.secondOfDay_;
    }

    return normalPointData;
}

std::map< double, double > extractNormalPointMeasurements( const std::vector< CrdPass >& passData )
{
    std::map< double, double > normalPointDataAtUtcSinceJ2000;
    for( const CrdPass& pass: passData )
    {
        std::map< double, double > singlePassData = extractNormalPointMeasurements( pass );
        normalPointDataAtUtcSinceJ2000.insert( singlePassData.begin( ), singlePassData.end( ) );
    }
    return normalPointDataAtUtcSinceJ2000;
}

std::map< double, double > extractFullRateMeasurements( const CrdPass& passData )
{
    if( passData.configuration_.startYear_ < 0 || passData.configuration_.startMonth_ < 0 || passData.configuration_.startDay_ < 0 )
    {
        throw std::runtime_error(
                "Error when extracting CRD full-rate measurements: pass start calendar date is not available." );
    }

    std::map< double, double > fullRateData;
    const double startOfDayInSeconds =
            basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch( passData.configuration_.startYear_,
                                                                             passData.configuration_.startMonth_,
                                                                             passData.configuration_.startDay_,
                                                                             0,
                                                                             0,
                                                                             0.0,
                                                                             basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    double currentDayStart = startOfDayInSeconds;
    double previousSecondOfDay = TUDAT_NAN;

    for( const CrdFullRateRecord& point: passData.data_.fullRateData_ )
    {
        if( previousSecondOfDay == previousSecondOfDay && point.secondOfDay_ + 1.0E4 < previousSecondOfDay )
        {
            currentDayStart += physical_constants::JULIAN_DAY;
        }

        fullRateData[ currentDayStart + point.secondOfDay_ ] = point.oneWayRange_;
        previousSecondOfDay = point.secondOfDay_;
    }

    return fullRateData;
}

std::map< double, double > extractFullRateMeasurements( const std::vector< CrdPass >& passData )
{
    std::map< double, double > fullRateDataAtUtcSinceJ2000;
    for( const CrdPass& pass: passData )
    {
        std::map< double, double > singlePassData = extractFullRateMeasurements( pass );
        fullRateDataAtUtcSinceJ2000.insert( singlePassData.begin( ), singlePassData.end( ) );
    }
    return fullRateDataAtUtcSinceJ2000;
}

std::map< std::string, double > getStationWavelengths( const std::map< std::string, std::vector< CrdPass > >& groupedData )
{
    std::map< std::string, double > wavelengthMap;

    for( const auto& stationData: groupedData )
    {
        if( stationData.second.empty( ) )
        {
            continue;
        }

        double firstWavelengthNm = stationData.second.at( 0 ).configuration_.transmitWavelengthNm_;
        for( unsigned int i = 1; i < stationData.second.size( ); i++ )
        {
            if( stationData.second.at( i ).configuration_.transmitWavelengthNm_ != firstWavelengthNm )
            {
                throw std::runtime_error( "Error when extracting station wavelength from CRD data for station " +
                                          stationData.first + ": wavelength is not consistent over passes." );
            }
        }
        wavelengthMap[ stationData.first ] = firstWavelengthNm * 1.0E-9;
    }

    return wavelengthMap;
}

}  // namespace input_output

}  // namespace tudat
