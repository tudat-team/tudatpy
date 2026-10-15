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

#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>
#include <stdexcept>

#include "tudat/io/textParsing.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{
namespace input_output
{
namespace
{

//! Read a numeric field; absent optional fields and CRDv2's na retain the record default.
template< typename Number >
void readNumber( const std::vector< std::string >& fields, const std::size_t index, Number& value )
{
    if( index >= fields.size( ) || fields[ index ] == "na" || fields[ index ] == "NA" )
    {
        return;
    }

    std::istringstream stream( fields[ index ] );

    stream.imbue( std::locale::classic( ) );
    if( !( stream >> value ) || !stream.eof( ) || !std::isfinite( value ) )
    {
        throw std::runtime_error( "invalid numeric field: " + fields[ index ] );
    }
}

//! Reject truncated data records instead of returning a partially populated measurement.
void requireFields( const std::vector< std::string >& fields, const std::size_t count )
{
    if( fields.size( ) < count )
    {
        throw std::runtime_error( "incomplete " + fields.front( ) + " record" );
    }
}

//! Retain only sessions with actual measurement or ancillary records.
bool passHasAnyData( const CrdPass& pass )
{
    return !pass.data_.fullRateData_.empty( ) || !pass.data_.normalPointData_.empty( ) || !pass.data_.meteorologicalData_.empty( ) ||
            !pass.data_.rangeSupplementData_.empty( ) || !pass.data_.calibrationData_.empty( );
}

//! Form precise UTC-coordinate epochs while retaining input order, duplicate events and midnight rollover.
template< typename Record >
std::vector< std::pair< Time, double > > extractMeasurements( const CrdPassConfigurationData& configuration,
                                                              const std::vector< Record >& records )
{
    if( configuration.rangeType_ != 2 ||
        ( configuration.stationEpochTimeScale_ != 3 && configuration.stationEpochTimeScale_ != 4 &&
          configuration.stationEpochTimeScale_ != 7 ) )
    {
        throw std::runtime_error( "CRD extraction requires two-way ranging and a declared station UTC reference (3, 4 or 7)." );
    }

    if( configuration.startYear_ < 0 || configuration.startMonth_ < 1 || configuration.startDay_ < 1 )
    {
        throw std::runtime_error( "CRD extraction requires the session start calendar date." );
    }

    const long double day =
            basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch< long double >( configuration.startYear_,
                                                                                           configuration.startMonth_,
                                                                                           configuration.startDay_,
                                                                                           0,
                                                                                           0,
                                                                                           0.0L,
                                                                                           basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    Time dayStart( static_cast< int >( day * 24.0L ), 0.0L );

    long double previousSeconds = TUDAT_NAN;

    std::vector< std::pair< Time, double > > output;

    output.reserve( records.size( ) );
    for( const Record& record : records )
    {
        if( std::isfinite( previousSeconds ) && record.secondOfDay_ + 1.0e4L < previousSeconds )
        {
            dayStart += physical_constants::JULIAN_DAY;
        }

        output.emplace_back( dayStart + record.secondOfDay_, record.oneWayRange_ );
        previousSeconds = record.secondOfDay_;
    }

    return output;
}

}  // namespace

//! Convert a two-way measured interval in seconds to a one-way range in metres.
double convertCrdTwoWayTimeOfFlightToSlrRange( const double twoWayTimeOfFlight )
{
    return 0.5 * physical_constants::SPEED_OF_LIGHT * twoWayTimeOfFlight;
}

//! Read CRD records without applying the corrections described by their flags.
std::vector< CrdPass > readCrdFile( std::istream& stream, const std::string& source )
{
    std::vector< CrdPass > output;

    CrdPass current;

    CrdPassConfigurationData sharedConfiguration;

    sharedConfiguration.source_ = source;
    current.configuration_ = sharedConfiguration;

    std::string line;

    std::size_t lineNumber = 0;

    while( std::getline( stream, line ) )
    {
        ++lineNumber;

        const std::vector< std::string > fields = text_parsing::splitIntoTokens( line );

        if( fields.empty( ) )
        {
            continue;
        }

        std::string id = fields.front( );

        if( id.front( ) == 'h' || id.front( ) == 'c' )
        {
            id.front( ) -= 'a' - 'A';
        }

        try
        {
            CrdPassConfigurationData& configuration = current.configuration_;

            if( id == "H8" || id == "H9" )
            {
                if( passHasAnyData( current ) )
                {
                    output.push_back( current );
                }

                current = CrdPass( );
                current.configuration_ = sharedConfiguration;
                if( id == "H9" )
                {
                    break;
                }
            }
            else if( id == "H1" )
            {
                requireFields( fields, 3 );
                sharedConfiguration = CrdPassConfigurationData( );
                sharedConfiguration.source_ = source;
                readNumber( fields, 2, sharedConfiguration.formatVersion_ );
                if( sharedConfiguration.formatVersion_ != 1 && sharedConfiguration.formatVersion_ != 2 )
                {
                    throw std::runtime_error( "only CRD versions 1 and 2 are supported" );
                }

                configuration = sharedConfiguration;
            }
            else if( id == "H2" )
            {
                requireFields( fields, 6 );
                configuration.stationName_ = fields[ 1 ];
                readNumber( fields, 2, configuration.cdpPadId_ );
                readNumber( fields, 3, configuration.systemNumber_ );
                readNumber( fields, 4, configuration.occupancy_ );
                readNumber( fields, 5, configuration.stationEpochTimeScale_ );
                configuration.stationNetwork_ = fields.size( ) > 6 ? fields[ 6 ] : "";
                sharedConfiguration = configuration;
            }
            else if( id == "H3" )
            {
                requireFields( fields, 3 );
                configuration.targetName_ = fields[ 1 ];
                configuration.ilrsSatelliteId_ = fields[ 2 ];
                sharedConfiguration = configuration;
            }
            else if( id == "H4" )
            {
                requireFields( fields, 22 );
                readNumber( fields, 1, configuration.dataType_ );
                readNumber( fields, 2, configuration.startYear_ );
                readNumber( fields, 3, configuration.startMonth_ );
                readNumber( fields, 4, configuration.startDay_ );
                readNumber( fields, 5, configuration.startHour_ );
                readNumber( fields, 6, configuration.startMinute_ );
                readNumber( fields, 7, configuration.startSecond_ );
                readNumber( fields, 8, configuration.endYear_ );
                readNumber( fields, 9, configuration.endMonth_ );
                readNumber( fields, 10, configuration.endDay_ );
                readNumber( fields, 11, configuration.endHour_ );
                readNumber( fields, 12, configuration.endMinute_ );
                readNumber( fields, 13, configuration.endSecond_ );
                readNumber( fields, 14, configuration.dataRelease_ );
                readNumber( fields, 15, configuration.troposphericCorrectionApplied_ );
                readNumber( fields, 16, configuration.centerOfMassCorrectionApplied_ );
                readNumber( fields, 17, configuration.receiveAmplitudeCorrectionApplied_ );
                readNumber( fields, 18, configuration.stationSystemDelayApplied_ );
                readNumber( fields, 19, configuration.spacecraftSystemDelayApplied_ );
                readNumber( fields, 20, configuration.rangeType_ );
                readNumber( fields, 21, configuration.dataQualityAlert_ );
            }
            else if( id == "C0" )
            {
                requireFields( fields, 4 );

                CrdSystemConfigurationRecord record;

                readNumber( fields, 2, record.transmitWavelengthNm_ );
                record.systemConfigurationId_ = fields[ 3 ];
                record.componentConfigurationIds_.assign( fields.begin( ) + 4, fields.end( ) );
                configuration.systemConfigurations_[ fields[ 3 ] ] = record;
                sharedConfiguration.systemConfigurations_ = configuration.systemConfigurations_;
            }
            else if( id == "10" || id == "11" )
            {
                requireFields( fields, id == "11" ? 8 : 5 );

                CrdNormalPointRecord normal;

                CrdFullRateRecord full;

                long double seconds = TUDAT_NAN;

                double flightTime = TUDAT_NAN;

                int event = -1;

                readNumber( fields, 1, seconds );
                readNumber( fields, 2, flightTime );
                readNumber( fields, 4, event );
                if( !std::isfinite( seconds ) || !std::isfinite( flightTime ) )
                {
                    throw std::runtime_error( "missing observation epoch or time of flight" );
                }

                // A half-range is meaningful only when H4 explicitly identifies a two-way observation.
                const double halfRange = configuration.rangeType_ == 2 ? convertCrdTwoWayTimeOfFlightToSlrRange( flightTime ) : TUDAT_NAN;

                if( id == "11" )
                {
                    normal.secondOfDay_ = seconds;
                    normal.twoWayTimeOfFlight_ = flightTime;
                    normal.oneWayRange_ = halfRange;
                    normal.systemConfigurationId_ = fields[ 3 ];
                    normal.epochEvent_ = event;
                    readNumber( fields, 5, normal.normalPointWindowLength_ );
                    readNumber( fields, 6, normal.numberOfReturns_ );
                    readNumber( fields, 7, normal.binRms_ );
                    readNumber( fields, 8, normal.binSkew_ );
                    readNumber( fields, 9, normal.binKurtosis_ );
                    readNumber( fields, 10, normal.binPeakMinusMean_ );
                    readNumber( fields, 11, normal.returnRate_ );
                    readNumber( fields, 12, normal.detectorChannel_ );
                    readNumber( fields, 13, normal.signalToNoiseRatio_ );
                    current.data_.normalPointData_.push_back( normal );
                }
                else
                {
                    full.secondOfDay_ = seconds;
                    full.twoWayTimeOfFlight_ = flightTime;
                    full.oneWayRange_ = halfRange;
                    full.systemConfigurationId_ = fields[ 3 ];
                    full.epochEvent_ = event;
                    readNumber( fields, 5, full.filterFlag_ );
                    readNumber( fields, 6, full.detectorChannel_ );
                    readNumber( fields, 7, full.stopNumber_ );
                    current.data_.fullRateData_.push_back( full );
                }
            }
            else if( id == "20" )
            {
                requireFields( fields, 5 );

                CrdMeteoRecord record;

                readNumber( fields, 1, record.secondOfDay_ );
                readNumber( fields, 2, record.pressure_ );
                readNumber( fields, 3, record.temperature_ );
                readNumber( fields, 4, record.humidity_ );
                readNumber( fields, 5, record.originOfValues_ );
                current.data_.meteorologicalData_.push_back( record );
            }
            else if( id == "12" )
            {
                requireFields( fields, 7 );

                CrdRangeSupplementRecord record;

                readNumber( fields, 1, record.secondOfDay_ );
                record.systemConfigurationId_ = fields[ 2 ];
                readNumber( fields, 3, record.troposphericCorrectionPs_ );
                readNumber( fields, 4, record.centerOfMassCorrection_ );
                readNumber( fields, 5, record.neutralDensityFilter_ );
                readNumber( fields, 6, record.timeBiasApplied_ );
                readNumber( fields, 7, record.rangeRate_ );
                current.data_.rangeSupplementData_.push_back( record );
            }
            else if( id == "40" || id == "41" )
            {
                requireFields( fields, 16 );

                CrdCalibrationRecord record;

                record.recordType_ = id == "40" ? 40 : 41;
                readNumber( fields, 1, record.secondOfDay_ );
                readNumber( fields, 2, record.dataType_ );
                record.systemConfigurationId_ = fields[ 3 ];
                readNumber( fields, 4, record.numberRecorded_ );
                readNumber( fields, 5, record.numberUsed_ );
                readNumber( fields, 6, record.targetDistance_ );
                readNumber( fields, 7, record.systemDelayPs_ );
                readNumber( fields, 8, record.delayShiftPs_ );
                readNumber( fields, 9, record.rmsPs_ );
                readNumber( fields, 10, record.skew_ );
                readNumber( fields, 11, record.kurtosis_ );
                readNumber( fields, 12, record.peakMinusMeanPs_ );
                readNumber( fields, 13, record.calibrationType_ );
                readNumber( fields, 14, record.shiftType_ );
                readNumber( fields, 15, record.detectorChannel_ );
                readNumber( fields, 16, record.span_ );
                readNumber( fields, 17, record.returnRate_ );
                current.data_.calibrationData_.push_back( record );
            }
        }
        catch( const std::exception& error )
        {
            throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": " + error.what( ) );
        }
    }

    if( passHasAnyData( current ) )
    {
        output.push_back( current );
    }

    return output;
}

//! Open a CRD file and retain its source name in the returned sessions.
std::vector< CrdPass > readCrdFile( const std::string& fileName )
{
    std::ifstream stream( fileName );

    if( !stream )
    {
        throw std::runtime_error( "Error when opening CRD file: " + fileName );
    }

    return readCrdFile( stream, fileName );
}

//! Concatenate parsed sessions without changing their input order or correction metadata.
std::vector< CrdPass > readCrdFiles( const std::vector< std::string >& fileNames )
{
    std::vector< CrdPass > output;

    for( const std::string& file : fileNames )
    {
        const std::vector< CrdPass > sessions = readCrdFile( file );

        output.insert( output.end( ), sessions.begin( ), sessions.end( ) );
    }

    return output;
}

//! Group sessions by the target name supplied in H3.
std::map< std::string, std::vector< CrdPass > > groupCrdDataPerTarget( const std::vector< CrdPass >& crdPasses )
{
    std::map< std::string, std::vector< CrdPass > > grouped;

    for( const CrdPass& pass : crdPasses )
    {
        grouped[ pass.configuration_.targetName_ ].push_back( pass );
    }

    return grouped;
}

//! Group sessions using the caller's station-pad identifiers.
std::map< std::string, std::vector< CrdPass > > groupCrdDataPerStation( const std::vector< CrdPass >& crdPasses,
                                                                        const std::map< int, std::string >& names )
{
    std::map< std::string, std::vector< CrdPass > > grouped;

    for( const CrdPass& pass : crdPasses )
    {
        if( names.count( pass.configuration_.cdpPadId_ ) )
        {
            grouped[ names.at( pass.configuration_.cdpPadId_ ) ].push_back( pass );
        }
    }

    return grouped;
}

//! Extract normal points at their original event epochs without merging equal timestamps.
std::vector< std::pair< Time, double > > extractNormalPointMeasurements( const CrdPass& pass )
{
    return extractMeasurements( pass.configuration_, pass.data_.normalPointData_ );
}

//! Concatenate normal-point pairs from multiple sessions without changing input order.
std::vector< std::pair< Time, double > > extractNormalPointMeasurements( const std::vector< CrdPass >& passes )
{
    std::vector< std::pair< Time, double > > output;

    for( const CrdPass& pass : passes )
    {
        const auto measurements = extractNormalPointMeasurements( pass );

        output.insert( output.end( ), measurements.begin( ), measurements.end( ) );
    }

    return output;
}

//! Extract full-rate observations at their original event epochs without merging equal timestamps.
std::vector< std::pair< Time, double > > extractFullRateMeasurements( const CrdPass& pass )
{
    return extractMeasurements( pass.configuration_, pass.data_.fullRateData_ );
}

//! Concatenate full-rate pairs from multiple sessions without changing input order.
std::vector< std::pair< Time, double > > extractFullRateMeasurements( const std::vector< CrdPass >& passes )
{
    std::vector< std::pair< Time, double > > output;

    for( const CrdPass& pass : passes )
    {
        const auto measurements = extractFullRateMeasurements( pass );

        output.insert( output.end( ), measurements.begin( ), measurements.end( ) );
    }

    return output;
}

//! Return a wavelength only when all recorded C0 systems at a station use the same wavelength.
std::map< std::string, double > getStationWavelengths( const std::map< std::string, std::vector< CrdPass > >& groupedData )
{
    std::map< std::string, double > output;

    for( const auto& station : groupedData )
    {
        double wavelength = TUDAT_NAN;

        for( const CrdPass& pass : station.second )
        {
            for( const auto& system : pass.configuration_.systemConfigurations_ )
            {
                const double current = system.second.transmitWavelengthNm_;

                if( std::isfinite( wavelength ) && current != wavelength )
                {
                    throw std::runtime_error( "CRD station " + station.first + " uses more than one transmit wavelength." );
                }

                wavelength = current;
            }
        }

        output[ station.first ] = wavelength * 1.0e-9;
    }

    return output;
}

}  // namespace input_output
}  // namespace tudat
