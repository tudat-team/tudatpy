/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PREPROCESS_FDETS_FILE_H
#define TUDAT_PREPROCESS_FDETS_FILE_H

#include <Eigen/Core>

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

data::PlainLinkDefinition getFdetsLinkEnds( const std::string& spacecraftName,
                                            const std::string& earthName,
                                            const std::string& transmittingStationName,
                                            const std::string& receivingStationName );

template< typename ObservationScalarType >
std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > makeScalarFdetsObservations(
        const std::vector< double >& rawObservationValues,
        const double baseFrequency )
{
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
    observations.reserve( rawObservationValues.size( ) );
    for( double rawObservationValue : rawObservationValues )
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observation( 1 );
        observation( 0 ) = static_cast< ObservationScalarType >( baseFrequency + rawObservationValue );
        observations.push_back( observation );
    }
    return observations;
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > extractTrackingDataFromRawFdets(
        const std::vector< std::shared_ptr< TrackingTxtFileContents > >& rawFdetsDataVector,
        const std::vector< double >& baseFrequencies,
        const std::string& spacecraftName,
        const std::vector< std::string >& transmittingStationNames,
        const std::vector< std::string >& receivingStationNames,
        const std::string& earthName = "Earth" )
{
    if( rawFdetsDataVector.size( ) != baseFrequencies.size( ) || rawFdetsDataVector.size( ) != transmittingStationNames.size( ) ||
        rawFdetsDataVector.size( ) != receivingStationNames.size( ) )
    {
        throw std::runtime_error(
                "Error when reading FDETS files, lists of files, base frequencies and stations do not have the same size." );
    }

    std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > trackingDataSets;
    for( unsigned int i = 0; i < rawFdetsDataVector.size( ); ++i )
    {
        const std::vector< double > utcReceptionTimes =
                rawFdetsDataVector.at( i )->getDoubleDataColumn( TrackingDataType::utc_reception_time_j2000 );
        const std::vector< double > dopplerObservations =
                rawFdetsDataVector.at( i )->getDoubleDataColumn( TrackingDataType::doppler_measured_frequency );

        auto trackingData = std::make_shared< data::TrackingData< ObservationScalarType, TimeType > >(
                "DopplerMeasuredFrequency",
                getFdetsLinkEnds( spacecraftName, earthName, transmittingStationNames.at( i ), receivingStationNames.at( i ) ),
                makeScalarFdetsObservations< ObservationScalarType >( dopplerObservations, baseFrequencies.at( i ) ),
                utilities::staticCastVector< TimeType, double >( utcReceptionTimes ),
                "receiver",
                "UTC" );

        trackingData->addAncillarySettings( "Doppler base frequency", baseFrequencies.at( i ) );

        trackingDataSets.push_back( trackingData );
    }

    return trackingDataSets;
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
convertRawFdetsFiles( const std::vector< std::shared_ptr< TrackingTxtFileContents > >& rawFdetsDataVector,
                      const std::vector< double >& baseFrequencies,
                      const std::string& spacecraftName,
                      const std::vector< std::string >& transmittingStationNames,
                      const std::vector< std::string >& receivingStationNames,
                      const std::string& earthName = "Earth" )
{
    return std::make_pair(
            extractTrackingDataFromRawFdets< ObservationScalarType, TimeType >(
                    rawFdetsDataVector, baseFrequencies, spacecraftName, transmittingStationNames, receivingStationNames, earthName ),
            std::vector< std::shared_ptr< data::TrackingSupplementaryData > >( ) );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
readFdetsFiles( const std::vector< std::string >& fdetsFileNames,
                const std::vector< double >& baseFrequencies,
                const FdetDateFormat dateFormat,
                const std::string& spacecraftName,
                const std::vector< std::string >& transmittingStationNames,
                const std::vector< std::string >& receivingStationNames,
                const std::string& earthName = "Earth" )
{
    if( fdetsFileNames.size( ) != baseFrequencies.size( ) )
    {
        throw std::runtime_error( "Error when reading FDETS files, list of files has different size than list of base frequencies." );
    }

    std::vector< std::shared_ptr< TrackingTxtFileContents > > rawFdetsDataVector;
    rawFdetsDataVector.reserve( fdetsFileNames.size( ) );
    for( const std::string& fdetsFileName : fdetsFileNames )
    {
        rawFdetsDataVector.push_back( readFdetsFile( fdetsFileName, dateFormat ) );
    }

    return convertRawFdetsFiles< ObservationScalarType, TimeType >(
            rawFdetsDataVector, baseFrequencies, spacecraftName, transmittingStationNames, receivingStationNames, earthName );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
readFdetsFiles( const std::vector< std::string >& fdetsFileNames,
                const double baseFrequency,
                const FdetDateFormat dateFormat,
                const std::string& spacecraftName,
                const std::string& transmittingStationName,
                const std::string& receivingStationName,
                const std::string& earthName = "Earth" )
{
    return readFdetsFiles< ObservationScalarType, TimeType >( fdetsFileNames,
                                                              std::vector< double >( fdetsFileNames.size( ), baseFrequency ),
                                                              dateFormat,
                                                              spacecraftName,
                                                              std::vector< std::string >( fdetsFileNames.size( ), transmittingStationName ),
                                                              std::vector< std::string >( fdetsFileNames.size( ), receivingStationName ),
                                                              earthName );
}

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_PREPROCESS_FDETS_FILE_H
