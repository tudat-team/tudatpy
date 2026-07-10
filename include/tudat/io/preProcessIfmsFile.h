/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PREPROCESS_IFMS_FILE_H
#define TUDAT_PREPROCESS_IFMS_FILE_H

#include <Eigen/Core>

#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/basics/timeType.h"
#include "tudat/basics/utilities.h"
#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingSupplementaryData.h"

namespace tudat
{

namespace input_output
{

data::PlainLinkDefinition getIfmsLinkEnds( const std::string& spacecraftName,
                                           const std::string& earthName,
                                           const std::string& groundStationName );

double getIfmsFrequencyRampRate( const double frequencyRampRate );

template< typename ObservationScalarType >
std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > makeScalarIfmsObservations(
        const std::vector< double >& rawObservationValues )
{
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
    observations.reserve( rawObservationValues.size( ) );
    for( double rawObservationValue : rawObservationValues )
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observation( 1 );
        observation( 0 ) = static_cast< ObservationScalarType >( rawObservationValue );
        observations.push_back( observation );
    }
    return observations;
}

void addIfmsFrequencyRamps( const std::shared_ptr< TrackingTxtFileContents >& rawIfmsData,
                            const std::shared_ptr< data::RampedFrequencySupplementaryData >& frequencySupplementaryData );

template< typename ObservationScalarType = double, typename TimeType = Time >
std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > extractTrackingDataFromRawIfms(
        const std::vector< std::shared_ptr< TrackingTxtFileContents > >& rawIfmsDataVector,
        const std::string& spacecraftName,
        const std::vector< std::string >& groundStationNames,
        const std::string& earthName = "Earth" )
{
    if( rawIfmsDataVector.size( ) != groundStationNames.size( ) )
    {
        throw std::runtime_error( "Error when reading IFMS files, list of files has different size than list of stations." );
    }

    std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > trackingDataSets;
    for( unsigned int i = 0; i < rawIfmsDataVector.size( ); ++i )
    {
        const std::vector< double > utcReceptionTimes =
                rawIfmsDataVector.at( i )->getDoubleDataColumn( TrackingDataType::utc_reception_time_j2000 );
        const std::vector< double > dopplerObservations =
                rawIfmsDataVector.at( i )->getDoubleDataColumn( TrackingDataType::doppler_averaged_frequency );

        auto trackingData = std::make_shared< data::TrackingData< ObservationScalarType, TimeType > >(
                "DsnNWayAveragedDoppler",
                getIfmsLinkEnds( spacecraftName, earthName, groundStationNames.at( i ) ),
                makeScalarIfmsObservations< ObservationScalarType >( dopplerObservations ),
                utilities::staticCastVector< TimeType, double >( utcReceptionTimes ),
                "receiver",
                "UTC" );

        const auto& metaDataDoubleMap = rawIfmsDataVector.at( i )->getMetaDataDoubleMap( );
        auto integrationTimeIterator = metaDataDoubleMap.find( TrackingDataType::doppler_integration_time );
        if( integrationTimeIterator != metaDataDoubleMap.end( ) && std::isfinite( integrationTimeIterator->second ) )
        {
            trackingData->addAncillarySettings( "Doppler observable integration time", integrationTimeIterator->second );
        }

        trackingDataSets.push_back( trackingData );
    }

    return trackingDataSets;
}

inline std::vector< std::shared_ptr< data::TrackingSupplementaryData > > extractTrackingSupplementaryDataFromRawIfms(
        const std::vector< std::shared_ptr< TrackingTxtFileContents > >& rawIfmsDataVector,
        const std::vector< std::string >& groundStationNames,
        const std::string& earthName = "Earth" )
{
    if( rawIfmsDataVector.size( ) != groundStationNames.size( ) )
    {
        throw std::runtime_error( "Error when reading IFMS files, list of files has different size than list of stations." );
    }

    std::map< std::string, std::shared_ptr< data::RampedFrequencySupplementaryData > > frequencyDataPerStation;
    for( unsigned int i = 0; i < rawIfmsDataVector.size( ); ++i )
    {
        std::shared_ptr< data::RampedFrequencySupplementaryData >& frequencyData = frequencyDataPerStation[ groundStationNames.at( i ) ];
        if( frequencyData == nullptr )
        {
            frequencyData = std::make_shared< data::RampedFrequencySupplementaryData >( );
        }

        addIfmsFrequencyRamps( rawIfmsDataVector.at( i ), frequencyData );
    }

    std::vector< std::shared_ptr< data::TrackingSupplementaryData > > trackingSupplementaryDataSets;
    for( const auto& stationFrequencyData : frequencyDataPerStation )
    {
        auto supplementaryData = std::make_shared< data::TrackingSupplementaryData >( earthName, stationFrequencyData.first );
        supplementaryData->setFrequencySupplementaryData(
                std::vector< std::shared_ptr< data::FrequencySupplementaryData > >( { stationFrequencyData.second } ) );
        trackingSupplementaryDataSets.push_back( supplementaryData );
    }

    return trackingSupplementaryDataSets;
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
convertRawIfmsFiles( const std::vector< std::shared_ptr< TrackingTxtFileContents > >& rawIfmsDataVector,
                     const std::string& spacecraftName,
                     const std::vector< std::string >& groundStationNames,
                     const std::string& earthName = "Earth" )
{
    return std::make_pair( extractTrackingDataFromRawIfms< ObservationScalarType, TimeType >(
                                   rawIfmsDataVector, spacecraftName, groundStationNames, earthName ),
                           extractTrackingSupplementaryDataFromRawIfms( rawIfmsDataVector, groundStationNames, earthName ) );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
readIfmsFiles( const std::vector< std::string >& ifmsFileNames,
               const std::string& spacecraftName,
               const std::vector< std::string >& groundStationNames,
               const std::string& earthName = "Earth",
               const bool applyTroposphereCorrection = true,
               const bool filterInvalidLines = true )
{
    if( ifmsFileNames.size( ) != groundStationNames.size( ) )
    {
        throw std::runtime_error( "Error when reading IFMS files, list of files has different size than list of stations." );
    }

    std::vector< std::shared_ptr< TrackingTxtFileContents > > rawIfmsDataVector;
    rawIfmsDataVector.reserve( ifmsFileNames.size( ) );
    for( const std::string& ifmsFileName : ifmsFileNames )
    {
        rawIfmsDataVector.push_back( readIfmsFile( ifmsFileName, applyTroposphereCorrection, filterInvalidLines ) );
    }

    return convertRawIfmsFiles< ObservationScalarType, TimeType >( rawIfmsDataVector, spacecraftName, groundStationNames, earthName );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
readIfmsFiles( const std::vector< std::string >& ifmsFileNames,
               const std::string& spacecraftName,
               const std::string& groundStationName,
               const std::string& earthName = "Earth",
               const bool applyTroposphereCorrection = true,
               const bool filterInvalidLines = true )
{
    return readIfmsFiles< ObservationScalarType, TimeType >( ifmsFileNames,
                                                             spacecraftName,
                                                             std::vector< std::string >( ifmsFileNames.size( ), groundStationName ),
                                                             earthName,
                                                             applyTroposphereCorrection,
                                                             filterInvalidLines );
}

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_PREPROCESS_IFMS_FILE_H
