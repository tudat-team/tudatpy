/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_COMPRESSDOPPLEROBSERVATIONCOLLECTION_H
#define TUDAT_COMPRESSDOPPLEROBSERVATIONCOLLECTION_H

#include <memory>

#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"

namespace tudat
{
namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > compressDopplerData(
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > originalDopplerData,
        const unsigned int setId,
        const unsigned int compressionRatio,
        const std::map< std::string, Eigen::Vector3d >& approximateGroundStationPositions =
                simulation_setup::getCombinedApproximateGroundStationPositions( ) )
{
    if( compressionRatio == 0 )
    {
        throw std::runtime_error( "Error in Doppler data compression, compression ratio must be positive." );
    }

    const double cadenceTolerance = 0.01;
    ObservationScalarType floatingCompressionRatio =
            mathematical_constants::getFloatingInteger< ObservationScalarType >( compressionRatio );

    const observation_models::ObservationSetMetadata< ObservationScalarType, TimeType >& metadata =
            originalDopplerData->getObservationSetMetadata( setId );
    const observation_models::LinkDefinition& linkDefinition = originalDopplerData->getLinkDefinition( metadata.linkDefinitionId_ );
    const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >& originalAncillarySettings =
            originalDopplerData->getAncillarySettings( metadata.ancillarySettingsId_ );

    double currentCompressionTime = originalAncillarySettings->getAncillaryDoubleData(
            observation_models::ObservationAncillarySimulationVariable::doppler_integration_time );

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > originalObservations =
            originalDopplerData->getObservationsForSet( setId );
    std::vector< TimeType > originalObservationTimesTdb = originalDopplerData->getObservationTimesForSet( setId );

    earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter = earth_orientation::TerrestrialTimeScaleConverter( );
    std::string stationName = linkDefinition.linkEnds_.at( observation_models::LinkEndType::receiver ).getReferencePointName( );
    auto stationPositionIterator = approximateGroundStationPositions.find( stationName );
    if( stationPositionIterator == approximateGroundStationPositions.end( ) )
    {
        throw std::runtime_error( "Error in Doppler data compression, could not retrieve approximate station position for " + stationName );
    }
    Eigen::Vector3d stationPosition = stationPositionIterator->second;

    std::vector< TimeType > originalObservationTimesUtc =
            timeScaleConverter.getCurrentTimesFromSinglePosition< TimeType >( basic_astrodynamics::TimeScales::tdb_scale,
                                                                              basic_astrodynamics::TimeScales::utc_scale,
                                                                              originalObservationTimesTdb,
                                                                              stationPosition );

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > compressedObservations;
    std::vector< TimeType > compressedObservationTimesUtc;

    for( unsigned int runStart = 0; runStart < originalObservations.size( ); )
    {
        unsigned int runEnd = runStart + 1;
        while( runEnd < originalObservations.size( ) )
        {
            double currentTimeStep =
                    static_cast< double >( originalObservationTimesUtc.at( runEnd ) - originalObservationTimesUtc.at( runEnd - 1 ) );
            if( std::fabs( currentTimeStep - currentCompressionTime ) <= cadenceTolerance )
            {
                runEnd++;
            }
            else
            {
                break;
            }
        }

        for( unsigned int blockStart = runStart; blockStart + compressionRatio <= runEnd; blockStart += compressionRatio )
        {
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > newObservable = originalObservations.at( blockStart );
            TimeType newTime = originalObservationTimesUtc.at( blockStart );

            for( unsigned int j = 1; j < compressionRatio; j++ )
            {
                newObservable += originalObservations.at( blockStart + j );
                newTime += originalObservationTimesUtc.at( blockStart + j );
            }

            newObservable /= floatingCompressionRatio;
            newTime = newTime / floatingCompressionRatio;

            compressedObservations.push_back( newObservable );
            compressedObservationTimesUtc.push_back( newTime );
        }

        runStart = runEnd;
    }

    std::vector< Eigen::Vector3d > compressedEarthFixedPositions;
    for( unsigned int i = 0; i < compressedObservationTimesUtc.size( ); ++i )
    {
        compressedEarthFixedPositions.push_back( stationPosition );
    }
    std::vector< TimeType > compressedObservationTimesTdb =
            timeScaleConverter.getCurrentTimes< TimeType >( basic_astrodynamics::TimeScales::utc_scale,
                                                            basic_astrodynamics::TimeScales::tdb_scale,
                                                            compressedObservationTimesUtc,
                                                            compressedEarthFixedPositions );

    std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySimulationSettings =
            std::make_shared< observation_models::ObservationAncillarySimulationSettings >( *originalAncillarySettings );
    double originalIntegrationTime = ancillarySimulationSettings->getAncillaryDoubleData(
            observation_models::ObservationAncillarySimulationVariable::doppler_integration_time );
    ancillarySimulationSettings->setAncillaryDoubleData(
            observation_models::ObservationAncillarySimulationVariable::doppler_integration_time,
            originalIntegrationTime * static_cast< double >( floatingCompressionRatio ) );

    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > compressedDataset =
            std::make_shared< observation_models::ObservationDataset< ObservationScalarType, TimeType > >( );
    compressedDataset->addObservationSet( metadata.observableType_,
                                          linkDefinition,
                                          compressedObservations,
                                          compressedObservationTimesTdb,
                                          metadata.referenceLinkEnd_,
                                          std::vector< Eigen::VectorXd >( ),
                                          originalDopplerData->getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ ),
                                          ancillarySimulationSettings,
                                          std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                          std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                                          true );
    return compressedDataset;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > createCompressedDopplerDataset(
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > originalDopplerData,
        const unsigned int compressionRatio,
        const unsigned int minNumberObservations = 10,
        const double maxArcGap = 300.0,
        const std::map< std::string, Eigen::Vector3d >& approximateGroundStationPositions =
                simulation_setup::getCombinedApproximateGroundStationPositions( ) )
{
    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > compressedData =
            std::make_shared< observation_models::ObservationDataset< ObservationScalarType, TimeType > >( );

    for( unsigned int setId = 0; setId < originalDopplerData->getNumberOfObservationSets( ); ++setId )
    {
        const observation_models::ObservationSetMetadata< ObservationScalarType, TimeType >& metadata =
                originalDopplerData->getObservationSetMetadata( setId );
        if( metadata.observableType_ != observation_models::ObservableType::dsn_n_way_averaged_doppler )
        {
            compressedData->addObservationSetFromDataset( *originalDopplerData, setId );
            continue;
        }

        const observation_models::LinkDefinition& linkDefinition = originalDopplerData->getLinkDefinition( metadata.linkDefinitionId_ );
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations =
                originalDopplerData->getObservationsForSet( setId );
        const std::vector< TimeType > observationTimes = originalDopplerData->getObservationTimesForSet( setId );

        for( unsigned int arcStart = 0; arcStart < observationTimes.size( ); )
        {
            unsigned int arcEnd = arcStart + 1;
            while( arcEnd < observationTimes.size( ) &&
                   static_cast< double >( observationTimes.at( arcEnd ) - observationTimes.at( arcEnd - 1 ) ) <= maxArcGap )
            {
                ++arcEnd;
            }

            if( arcEnd - arcStart >= minNumberObservations )
            {
                std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > arcObservations( observations.begin( ) + arcStart,
                                                                                                          observations.begin( ) + arcEnd );
                std::vector< TimeType > arcObservationTimes( observationTimes.begin( ) + arcStart, observationTimes.begin( ) + arcEnd );

                std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > arcDataset =
                        std::make_shared< observation_models::ObservationDataset< ObservationScalarType, TimeType > >( );
                arcDataset->addObservationSet( metadata.observableType_,
                                               linkDefinition,
                                               arcObservations,
                                               arcObservationTimes,
                                               metadata.referenceLinkEnd_,
                                               std::vector< Eigen::VectorXd >( ),
                                               originalDopplerData->getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ ),
                                               originalDopplerData->getAncillarySettings( metadata.ancillarySettingsId_ ),
                                               std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                               std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                                               true );

                std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > compressedArcDataset =
                        compressDopplerData< ObservationScalarType, TimeType >(
                                arcDataset, 0, compressionRatio, approximateGroundStationPositions );
                if( compressedArcDataset->getNumberOfObservationSets( ) > 0 &&
                    compressedArcDataset->getNumberOfObservationsForSet( 0 ) > 0 )
                {
                    compressedData->addObservationSetFromDataset( *compressedArcDataset, 0 );
                }
            }

            arcStart = arcEnd;
        }
    }

    return compressedData;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > compressDopplerData(
        const std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > originalDopplerData,
        const unsigned int compressionRatio,
        const std::map< std::string, Eigen::Vector3d >& approximateGroundStationPositions =
                simulation_setup::getCombinedApproximateGroundStationPositions( ) )
{
    return observation_models::createSingleObservationSet< ObservationScalarType, TimeType >(
            compressDopplerData< ObservationScalarType, TimeType >(
                    originalDopplerData->getObservationDataset( ),
                    originalDopplerData->getObservationSetId( ),
                    compressionRatio,
                    approximateGroundStationPositions ) );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createCompressedDopplerCollection(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > originalDopplerData,
        const unsigned int compressionRatio,
        const unsigned int minNumberObservations = 10,
        const double maxArcGap = 300.0,
        const std::map< std::string, Eigen::Vector3d >& approximateGroundStationPositions =
                simulation_setup::getCombinedApproximateGroundStationPositions( ) )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createCompressedDopplerDataset< ObservationScalarType, TimeType >(
                    originalDopplerData->getObservationDataset( ),
                    compressionRatio,
                    minNumberObservations,
                    maxArcGap,
                    approximateGroundStationPositions ) );
}

}  // namespace observation_models
}  // namespace tudat

#endif  // TUDAT_COMPRESSDOPPLEROBSERVATIONCOLLECTION_H
