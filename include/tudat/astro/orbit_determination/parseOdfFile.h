/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PARSEODFFILE_H
#define TUDAT_PARSEODFFILE_H

#include "tudat/basics/utilities.h"
#include "tudat/io/readOdfFile.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationFrequencies.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/simulation/estimation_setup/observations.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/math/quadrature/trapezoidQuadrature.h"

namespace tudat
{

namespace orbit_determination
{

observation_models::ObservableType getObservableTypeForOdfId( const int odfId );

observation_models::FrequencyBands getFrequencyBandForOdfId ( const int odfId );

std::string getStationNameFromStationId ( const int networkId, const int stationId );

class ProcessedOdfFileSingleLinkData
{
public:

    ProcessedOdfFileSingleLinkData( ){ }

    virtual ~ProcessedOdfFileSingleLinkData( ){ }

    std::vector< double > observationTimes_;
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observableValues_;
    std::vector< double > receiverDownlinkDelays_;

    std::vector< int > downlinkBandIds_;
    std::vector< int > uplinkBandIds_;
    std::vector< int > referenceBandIds_;

    std::vector< std::string > originFiles_;

    observation_models::ObservableType observableType_;

    std::string transmittingStation_;
    std::string receivingStation_;

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > getUnprocessedObservables( )
    {
        return utilities::createMapFromVectors( observationTimes_, observableValues_ );
    }

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getUnprocessedObservablesVector( )
    {
        return observableValues_;
    }

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > getProcessedObservables( )
    {
        return utilities::createMapFromVectors( observationTimes_, getProcessedObservablesVector( ) );
    }

    virtual std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getProcessedObservablesVector( )
    {
        return observableValues_;
    }

    std::vector< double > getObservationTimesUtc (  )
    {
        return observationTimes_;
    }

    std::vector< double > getObservationTimesTdb (
            const simulation_setup::SystemOfBodies& bodies )
    {
        earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter =
                earth_orientation::TerrestrialTimeScaleConverter( );

        std::vector< Eigen::Vector3d > earthFixedPositions;
        for ( unsigned int i = 0; i < observationTimes_.size( ); ++i )
        {
            earthFixedPositions.push_back(
                    bodies.getBody( "Earth" )->getGroundStation( receivingStation_ )->getStateInPlanetFixedFrame< double, double >(
                            observationTimes_.at( i ) ).segment( 0, 3 )
                            );
        }

        return timeScaleConverter.getCurrentTimes(
                basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, observationTimes_,
                earthFixedPositions );
    }

};


class ProcessedOdfFileDopplerData: public ProcessedOdfFileSingleLinkData
{
public:

    ~ProcessedOdfFileDopplerData( ){ }

    std::vector< int > receiverChannels_;
    std::vector< double > referenceFrequencies_;
    std::vector< double > countInterval_;
    std::vector< double > transmitterUplinkDelays_;
    std::vector< bool > receiverRampingFlags_;

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getProcessedObservablesVector( )
    {
        throw std::runtime_error("Getting range rate observables not implemented");
        return observableValues_;
    }

    std::map< double, bool > getReceiverRampingFlags( )
    {
        return utilities::createMapFromVectors( observationTimes_, receiverRampingFlags_ );
    }

    std::map< double, double > getReferenceFrequencies( )
    {
        return utilities::createMapFromVectors( observationTimes_, referenceFrequencies_ );
    }

    std::map< double, double > getCountInterval( )
    {
        return utilities::createMapFromVectors( observationTimes_, countInterval_ );
    }
};


class ProcessedOdfFileContents
{
public:

    std::string spacecraftName_;

    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > > > processedDataBlocks_;

    std::map< std::string, std::shared_ptr< observation_models::PiecewiseLinearFrequencyInterpolator > > rampInterpolators_;
};


std::shared_ptr< observation_models::PiecewiseLinearFrequencyInterpolator > mergeRampDataInterpolators(
        const std::vector< std::shared_ptr< observation_models::PiecewiseLinearFrequencyInterpolator > >& interpolatorList );

void addOdfFileContentsToMergedContents(
        const observation_models::ObservableType observableType,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > mergedOdfFileContents,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > blockToAdd );

std::shared_ptr< ProcessedOdfFileContents > mergeOdfFileContents(
        const std::vector< std::shared_ptr< ProcessedOdfFileContents > > odfFileContents );

void addOdfDataBlockToProcessedData(
        const observation_models::ObservableType currentObservableType,
        const std::shared_ptr< input_output::OdfDataBlock > rawDataBlock,
        const std::shared_ptr< ProcessedOdfFileSingleLinkData > processedDataBlock );

observation_models::LinkEnds getLinkEndsFromOdfBlock (
        const std::shared_ptr< input_output::OdfDataBlock > dataBlock,
        std::string spacecraftName );

std::shared_ptr< ProcessedOdfFileContents > processOdfFileContents(
        const std::shared_ptr< input_output::OdfRawFileContents > rawOdfData,
        bool verbose = true );

template< typename TimeType = double >
observation_models::ObservationAncilliarySimulationSettings< TimeType > createOdfAncillarySettings
        ( std::shared_ptr< ProcessedOdfFileSingleLinkData > odfDataContents,
          unsigned int dataIndex )
{
    if ( dataIndex >= odfDataContents->observationTimes_.size( ) )
    {
        throw std::runtime_error("Error when creating ODF data ancillary settings: specified data index is larger than data size.");
    }

    observation_models::ObservationAncilliarySimulationSettings< TimeType > ancillarySettings =
            observation_models::ObservationAncilliarySimulationSettings< TimeType >( );

    observation_models::ObservableType currentObservableType = odfDataContents->observableType_;

    // Set common ancillary settings
    ancillarySettings.setAncilliaryDoubleData(
            observation_models::uplink_band, getFrequencyBandForOdfId(
                    odfDataContents->uplinkBandIds_.at( dataIndex ) ) );
    ancillarySettings.setAncilliaryDoubleData(
            observation_models::downlink_band, getFrequencyBandForOdfId(
                    odfDataContents->downlinkBandIds_.at( dataIndex ) ) );

    if ( std::dynamic_pointer_cast< ProcessedOdfFileDopplerData >( odfDataContents ) != nullptr )
    {
        std::shared_ptr< ProcessedOdfFileDopplerData > dopplerDataBlock
                = std::dynamic_pointer_cast< ProcessedOdfFileDopplerData >( odfDataContents );

        ancillarySettings.setAncilliaryDoubleData(
                observation_models::doppler_integration_time, dopplerDataBlock->getCountInterval( ).at( dataIndex ) );
        ancillarySettings.setAncilliaryDoubleData(
                observation_models::doppler_reference_frequency, dopplerDataBlock->referenceFrequencies_.at( dataIndex ) );

        if ( currentObservableType == observation_models::n_way_differenced_range )
        {
            ancillarySettings.setAncilliaryDoubleVectorData(
                observation_models::retransmission_delays,  std::vector< double >{
                    dopplerDataBlock->transmitterUplinkDelays_.at( dataIndex ), 0.0,
                    dopplerDataBlock->receiverDownlinkDelays_.at( dataIndex ) } );
        }
        else
        {
            ancillarySettings.setAncilliaryDoubleVectorData(
                observation_models::retransmission_delays, std::vector< double >{
                    dopplerDataBlock->transmitterUplinkDelays_.at( dataIndex ),
                    dopplerDataBlock->receiverDownlinkDelays_.at( dataIndex ) } );
        }

    }
    else
    {
        throw std::runtime_error("Error when casting ODF processed data: data type not identified.");
    }
}

template< typename ObservationScalarType = double, typename TimeType = double >
void separateSingleLinkOdfData(
        observation_models::ObservableType currentObservableType,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > odfSingleLinkData,
        std::vector< std::vector< TimeType > >& observationTimes,
        std::vector< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > >& observables,
        std::vector< observation_models::ObservationAncilliarySimulationSettings< TimeType > >& ancillarySettings,
        const simulation_setup::SystemOfBodies& bodies )
{
    // Initialize vectors
    observationTimes.clear( );
    observables.clear( );
    ancillarySettings.clear( );

    // Get time and observables vectors
    std::vector< double > observationTimesTdb = odfSingleLinkData->getObservationTimesTdb( bodies );
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observablesVector =
            odfSingleLinkData->getProcessedObservablesVector( );

    for ( unsigned int i = 0; i < odfSingleLinkData->observationTimes_.size( ); ++i )
    {
        observation_models::ObservationAncilliarySimulationSettings< TimeType > currentAncillarySettings =
                createOdfAncillarySettings( odfSingleLinkData, i );

        bool newAncillarySettings = true;

        for ( unsigned int j = 0; j < ancillarySettings.size( ); ++j )
        {
            if ( ancillarySettings.at( j ) == currentAncillarySettings )
            {
                newAncillarySettings = false;
                observationTimes.at( j ).push_back( observationTimesTdb.at( i ) );
                observables.at( j ).push_back( observablesVector.at( i ) );
                break;
            }
        }

        if ( newAncillarySettings )
        {
            observationTimes.push_back ( std::vector< TimeType >{ observationTimesTdb.at( i ) } );
            observables.push_back( std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >{
                observablesVector.at( i ) } );
            ancillarySettings.push_back( currentAncillarySettings );
        }
    }
}

void setGroundStationsTransmittingFrequencies(
        std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents,
        const simulation_setup::SystemOfBodies& bodies )
{
    for( auto it = processedOdfFileContents->rampInterpolators_.begin( ); it != processedOdfFileContents->rampInterpolators_.end( ); it++ )
    {
       bodies.getBody( "Earth" )->getGroundStation( it->first )->setTransmittingFrequencyCalculator(
               it->second );
    }
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createOdfObservationCollection(
        std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr )
{

    // Create and fill single observation sets
    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, std::vector< std::shared_ptr<
            observation_models::SingleObservationSet< ObservationScalarType, TimeType > > > > > sortedObservationSets;

    for ( auto observableTypeIterator = processedOdfFileContents->processedDataBlocks_.begin( );
            observableTypeIterator != processedOdfFileContents->processedDataBlocks_.end( ); ++observableTypeIterator )
    {
        observation_models::ObservableType currentObservableType = observableTypeIterator->first;

        for ( auto linkEndsIterator = observableTypeIterator->second.begin( );
                linkEndsIterator != observableTypeIterator->second.end( ); ++linkEndsIterator )
        {
            observation_models::LinkEnds currentLinkEnds = linkEndsIterator->first;
            std::shared_ptr< ProcessedOdfFileSingleLinkData > currentOdfSingleLinkData = linkEndsIterator->second;

            // Reset vector of observation sets
            sortedObservationSets.at( currentObservableType ).at( currentLinkEnds ).clear( );

            // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
            std::vector< std::vector< TimeType > > observationTimes;
            std::vector< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > observables;
            std::vector< observation_models::ObservationAncilliarySimulationSettings< TimeType > > ancillarySettings;

            // Fill vectors
            separateSingleLinkOdfData(
                    currentObservableType, currentOdfSingleLinkData, observationTimes, observables, ancillarySettings,
                    bodies );

            // Create the single observation sets and save them
            for ( unsigned int i = 0; i < observationTimes.size( ); ++i )
            {
                if ( dependentVariableCalculator != nullptr )
                {
//                    sortedObservationSets.at( currentObservableType ).at( currentLinkEnds ).push_back(
//                        currentObservableType, currentLinkEnds, observables.at( i ), observationTimes.at( i ),
//                        observation_models::receiver,
//                        dependentVariableCalculator->calculateDependentVariables( ),
//                        dependentVariableCalculator, ancillarySettings.at( i ) );
                    throw std::runtime_error( "Computation of dependent variables for ODF observables is not implemented." );
                }
                else
                {
                    sortedObservationSets.at( currentObservableType ).at( currentLinkEnds ).push_back(
                        currentObservableType, currentLinkEnds, observables.at( i ), observationTimes.at( i ),
                        observation_models::receiver,
                        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                        nullptr, ancillarySettings.at( i ) );
                }

            }
        }
    }

    // Add transmitting stations to ground stations
    for( auto it = processedOdfFileContents->rampInterpolators_.begin( );
            it != processedOdfFileContents->rampInterpolators_.end( ); it++ )
    {
       bodies.getBody( "Earth" )->getGroundStation( it->first )->setTransmittingFrequencyCalculator( it->second );
    }

}

} // namespace orbit_determination

} // namespace tudat

#endif // TUDAT_PARSEODFFILE_H
