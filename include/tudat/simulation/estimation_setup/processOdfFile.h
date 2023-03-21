/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROCESSODFFILE_H
#define TUDAT_PROCESSODFFILE_H

#include "tudat/basics/utilities.h"
#include "tudat/io/readOdfFile.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/simulation/estimation_setup/observations.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/math/quadrature/trapezoidQuadrature.h"

namespace tudat
{

namespace observation_models
{

observation_models::ObservableType getObservableTypeForOdfId( const int odfId );

observation_models::FrequencyBands getFrequencyBandForOdfId ( const int odfId );

std::string getStationNameFromStationId ( const int networkId, const int stationId );

class ProcessedOdfFileSingleLinkData
{
public:

    ProcessedOdfFileSingleLinkData( ){ }

    virtual ~ProcessedOdfFileSingleLinkData( ){ }

    std::vector< double > unprocessedObservationTimes_;
    std::vector< double > processedObservationTimes_;

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observableValues_;
    std::vector< double > receiverDownlinkDelays_;

    std::vector< int > downlinkBandIds_;
    std::vector< int > uplinkBandIds_;
    std::vector< int > referenceBandIds_;

    std::vector< std::string > originFiles_;

    observation_models::ObservableType observableType_;

    std::string transmittingStation_;
    std::string receivingStation_;

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > getObservables( )
    {
        return utilities::createMapFromVectors( processedObservationTimes_, observableValues_ );
    }

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getObservablesVector( )
    {
        return observableValues_;
    }

    std::vector< double > getObservationTimesVector( )
    {
        return processedObservationTimes_;
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

    std::map< double, bool > getReceiverRampingFlags( )
    {
        return utilities::createMapFromVectors( processedObservationTimes_, receiverRampingFlags_ );
    }

    std::map< double, double > getReferenceFrequencies( )
    {
        return utilities::createMapFromVectors( processedObservationTimes_, referenceFrequencies_ );
    }

    std::vector< double > getReferenceFrequenciesVector( )
    {
        return referenceFrequencies_;
    }

    std::map< double, double > getCountInterval( )
    {
        return utilities::createMapFromVectors( processedObservationTimes_, countInterval_ );
    }
};


class ProcessedOdfFileContents
{
public:

    // Constructor for single raw ODF data object
    ProcessedOdfFileContents(
            const std::shared_ptr< const input_output::OdfRawFileContents > rawOdfData,
            const std::shared_ptr< const simulation_setup::Body > bodyWithGroundStations,
            bool verbose = true,
            std::string spacecraftName = "" ):
        spacecraftName_( spacecraftName ),
        bodyWithGroundStations_ ( bodyWithGroundStations ),
        verbose_( verbose )
    {
        if ( spacecraftName != "" )
        {
            spacecraftName_ = spacecraftName;
        }
        else
        {
            // Spacecraft name selected to be the "NAIF Id", which is equal to -"JPL Id" (for a spacecraft)
            spacecraftName_ = std::to_string( - rawOdfData->spacecraftId_ );
        }

        extractRawOdfOrbitData( rawOdfData );
        updateProcessedObservationTimes( );
        extractRawOdfRampData( rawOdfData );
    }

    // Constructor for multiple ODF data objects
    ProcessedOdfFileContents(
            std::vector< std::shared_ptr< const input_output::OdfRawFileContents > > rawOdfDataVector,
            const std::shared_ptr< const simulation_setup::Body > bodyWithGroundStations,
            bool verbose = true,
            std::string spacecraftName = "" ):
        spacecraftName_( spacecraftName ),
        bodyWithGroundStations_ ( bodyWithGroundStations ),
        verbose_( verbose )
    {
        unsigned int spacecraftId = rawOdfDataVector.front( )->spacecraftId_;

        if ( spacecraftName != "" )
        {
            spacecraftName_ = spacecraftName;
        }
        else
        {
            // Spacecraft name selected to be the "NAIF Id", which is equal to -"JPL Id" (for a spacecraft)
            spacecraftName_ = std::to_string( - spacecraftId );
        }

        for ( unsigned int i = 0; i < rawOdfDataVector.size( ); ++i )
        {
            // Check if spacecraft ID is valid
            if ( rawOdfDataVector.at( i )->spacecraftId_ != spacecraftId )
            {
                throw std::runtime_error("Error when creating processed ODF object from raw data: multiple spacecraft IDs"
                                         "found (" + std::to_string( spacecraftId ) + " and " +
                                         std::to_string( rawOdfDataVector.at( i )->spacecraftId_ ) + ").");
            }

            extractRawOdfOrbitData( rawOdfDataVector.at( i ) );
        }

        updateProcessedObservationTimes( );

        extractMultipleRawOdfRampData( rawOdfDataVector );
    }

    std::string getSpacecraftName( )
    {
        return spacecraftName_;
    }

    std::vector< std::string > getGroundStationsNames( );

    std::vector< observation_models::ObservableType > getProcessedObservableTypes( );

    std::pair< double, double > getStartAndEndTime( const simulation_setup::SystemOfBodies& bodies );

    std::vector< int > getNotReadRawOdfObservableTypes( )
    {
        return notReadRawOdfObservableTypes_;
    }

    const std::map< std::string, std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > >& getRampInterpolators( )
    {
        return rampInterpolators_;
    }

    const std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > > >& getProcessedDataBlocks( )
    {
        return processedDataBlocks_;
    }

private:

    void extractRawOdfOrbitData( std::shared_ptr< const input_output::OdfRawFileContents > rawOdfData );

    void extractRawOdfRampData( std::shared_ptr< const input_output::OdfRawFileContents > rawOdfData );

    void extractMultipleRawOdfRampData(
            std::vector< std::shared_ptr< const input_output::OdfRawFileContents > > rawOdfDataVector );

    void addOdfRawDataBlockToProcessedData(
            const observation_models::ObservableType currentObservableType,
            const std::shared_ptr< const input_output::OdfDataBlock > rawDataBlock,
            const std::shared_ptr< ProcessedOdfFileSingleLinkData > singleLinkProcessedData,
            const std::string rawDataFileName );

    void updateProcessedObservationTimes( );

    std::vector< double > computeObservationTimesUtcFromJ2000( std::vector< double > observationTimesUtcFromEME1950 );

    std::vector< double > computeObservationTimesTdbFromJ2000(
            std::string groundStation, std::vector< double > observationTimesUtcFromEME1950 );

    std::string spacecraftName_;

    std::shared_ptr< const simulation_setup::Body > bodyWithGroundStations_;

    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > > > processedDataBlocks_;

    std::map< std::string, std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > > rampInterpolators_;

    // Vector keeping the invalid observable types that were found in the raw ODF data
    std::vector< int > notReadRawOdfObservableTypes_;

    bool verbose_;
};


//std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > mergeRampDataInterpolators(
//        const std::vector< std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > >& interpolatorList );
//
//void addOdfFileContentsToMergedContents(
//        const observation_models::ObservableType observableType,
//        std::shared_ptr< ProcessedOdfFileSingleLinkData > mergedOdfFileContents,
//        std::shared_ptr< ProcessedOdfFileSingleLinkData > blockToAdd );
//
//std::shared_ptr< ProcessedOdfFileContents > mergeOdfFileContents(
//        const std::vector< std::shared_ptr< ProcessedOdfFileContents > > odfFileContents );

observation_models::LinkEnds getLinkEndsFromOdfBlock (
        const std::shared_ptr< input_output::OdfDataBlock > dataBlock,
        std::string spacecraftName );

std::shared_ptr< ProcessedOdfFileContents > processOdfFileContents(
        const std::shared_ptr< input_output::OdfRawFileContents > rawOdfData,
        std::string spacecraftName = "",
        bool verbose = true );

template< typename TimeType = double >
observation_models::ObservationAncilliarySimulationSettings< TimeType > createOdfAncillarySettings
        ( std::shared_ptr< ProcessedOdfFileSingleLinkData > odfDataContents,
          unsigned int dataIndex )
{
    if ( dataIndex >= odfDataContents->unprocessedObservationTimes_.size( ) )
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
                observation_models::doppler_integration_time, dopplerDataBlock->countInterval_.at( dataIndex ) );
        ancillarySettings.setAncilliaryDoubleData(
                 observation_models::doppler_reference_frequency, dopplerDataBlock->referenceFrequencies_.at( dataIndex ) );

        if ( currentObservableType == observation_models::dsn_n_way_averaged_doppler )
        {
//            std::cout << "Delays (" << odfDataContents->processedObservationTimes_.at( dataIndex ) << "): " <<
//                dopplerDataBlock->transmitterUplinkDelays_.at( dataIndex ) << " " << 0.0 << " " <<
//                dopplerDataBlock->receiverDownlinkDelays_.at( dataIndex ) << std::endl;
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

    return ancillarySettings;
}

template< typename ObservationScalarType = double, typename TimeType = double >
void separateSingleLinkOdfData(
        observation_models::ObservableType currentObservableType,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > odfSingleLinkData,
        std::vector< std::vector< TimeType > >& observationTimes,
        std::vector< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > >& observables,
        std::vector< observation_models::ObservationAncilliarySimulationSettings< TimeType > >& ancillarySettings )
{
    // Initialize vectors
    observationTimes.clear( );
    observables.clear( );
    ancillarySettings.clear( );

    // Get time and observables vectors
    std::vector< double > observationTimesTdb = odfSingleLinkData->getObservationTimesVector( );
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observablesVector =
            odfSingleLinkData->getObservablesVector( );

    for ( unsigned int i = 0; i < odfSingleLinkData->unprocessedObservationTimes_.size( ); ++i )
    {
        observation_models::ObservationAncilliarySimulationSettings< TimeType > currentAncillarySettings =
                createOdfAncillarySettings< TimeType >( odfSingleLinkData, i );

        bool newAncillarySettings = true;

        for ( unsigned int j = 0; j < ancillarySettings.size( ); ++j )
        {
            if ( ancillarySettings.at( j ) == currentAncillarySettings )
            {
                newAncillarySettings = false;
                observationTimes.at( j ).push_back( static_cast< TimeType >( observationTimesTdb.at( i ) ) );
                observables.at( j ).push_back( observablesVector.at( i ).template cast< ObservationScalarType >( ) );
                break;
            }
        }

        if ( newAncillarySettings )
        {
            observationTimes.push_back ( std::vector< TimeType >{ static_cast< TimeType >( observationTimesTdb.at( i ) ) } );
            observables.push_back( std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >{
                observablesVector.at( i ).template cast< ObservationScalarType >( ) } );
            ancillarySettings.push_back( currentAncillarySettings );
        }
    }
}

// Add transmitting stations to ground stations
inline void setGroundStationsTransmittingFrequencies(
        std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents,
        std::shared_ptr< simulation_setup::Body > bodyWithGroundStations )
{
    for( auto it = processedOdfFileContents->getRampInterpolators( ).begin( );
            it != processedOdfFileContents->getRampInterpolators( ).end( ); it++ )
    {
       bodyWithGroundStations->getGroundStation( it->first )->setTransmittingFrequencyCalculator( it->second );
    }
}


template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createOdfObservedObservationCollection(
        std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents,
        simulation_setup::SystemOfBodies& bodies,
        bool setGroundStationsFrequencies = true,
        std::string bodyWithGroundStations = "Earth")
{

    // Create and fill single observation sets
    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, std::vector< std::shared_ptr<
            observation_models::SingleObservationSet< ObservationScalarType, TimeType > > > > > sortedObservationSets;

    for ( auto observableTypeIterator = processedOdfFileContents->getProcessedDataBlocks( ).begin( );
            observableTypeIterator != processedOdfFileContents->getProcessedDataBlocks( ).end( ); ++observableTypeIterator )
    {
        observation_models::ObservableType currentObservableType = observableTypeIterator->first;

        for ( auto linkEndsIterator = observableTypeIterator->second.begin( );
                linkEndsIterator != observableTypeIterator->second.end( ); ++linkEndsIterator )
        {
            observation_models::LinkEnds currentLinkEnds = linkEndsIterator->first;
            std::shared_ptr< ProcessedOdfFileSingleLinkData > currentOdfSingleLinkData = linkEndsIterator->second;

            // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
            std::vector< std::vector< TimeType > > observationTimes;
            std::vector< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > observables;
            std::vector< observation_models::ObservationAncilliarySimulationSettings< TimeType > > ancillarySettings;

            // Fill vectors
            separateSingleLinkOdfData(
                    currentObservableType, currentOdfSingleLinkData, observationTimes, observables, ancillarySettings );

            // Create the single observation sets and save them
            for ( unsigned int i = 0; i < observationTimes.size( ); ++i )
            {
                sortedObservationSets[ currentObservableType ][ currentLinkEnds ].push_back(
                    std::make_shared< observation_models::SingleObservationSet< ObservationScalarType, TimeType > >(
                        currentObservableType, currentLinkEnds, observables.at( i ), observationTimes.at( i ),
                        observation_models::receiver,
                        std::vector< Eigen::VectorXd >( ),
                        nullptr, std::make_shared< observation_models::ObservationAncilliarySimulationSettings< TimeType > >(
                            ancillarySettings.at( i ) ) ) );
            }
        }
    }

    if ( setGroundStationsFrequencies )
    {
        setGroundStationsTransmittingFrequencies( processedOdfFileContents, bodies.getBody( bodyWithGroundStations ) );
    }

    return std::make_shared< observation_models::ObservationCollection< ObservationScalarType, TimeType > >(
            sortedObservationSets );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< TimeType > > > createOdfObservationSimulationSettingsList(
        std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observedObservationCollection )
{
    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< TimeType > > > observationSimulationSettings;

    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, std::vector< std::shared_ptr<
            observation_models::SingleObservationSet< ObservationScalarType, TimeType > > > > > observationSetList =
                    observedObservationCollection->getObservations( );

    for ( auto observableTypeIterator = observationSetList.begin( ); observableTypeIterator != observationSetList.end( );
            ++observableTypeIterator )
    {
        observation_models::ObservableType currentObservableType = observableTypeIterator->first;

        for ( auto linkEndsIterator = observableTypeIterator->second.begin( );
              linkEndsIterator != observableTypeIterator->second.end( ); ++linkEndsIterator )
        {
            observation_models::LinkEnds currentLinkEnds = linkEndsIterator->first;
            std::vector< std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObservationSets = linkEndsIterator->second;

            for ( unsigned int i = 0; i < singleObservationSets.size( ); ++i )
            {
                observationSimulationSettings.push_back(
                        std::make_shared< simulation_setup::TabulatedObservationSimulationSettings< TimeType > >(
                                currentObservableType, currentLinkEnds, singleObservationSets.at( i )->getObservationTimes( ),
                                singleObservationSets.at( i )->getReferenceLinkEnd( ),
                                std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ),
                                nullptr,
                                singleObservationSets.at( i )->getAncilliarySettings( )
                                )
                        );
            }

        }
    }

    return observationSimulationSettings;
}

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_PROCESSODFFILE_H
