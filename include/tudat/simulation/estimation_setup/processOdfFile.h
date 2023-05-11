/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References: 820-013, TRK-2-18 Tracking System Interfaces Orbit Data File Interface, Revision E, 2008, JPL/DSN
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

/*! Get the observable type associated with a given ODF data ID.
 *
 * @param odfId ODF data ID
 * @return Observable type
 */
observation_models::ObservableType getObservableTypeForOdfId( const int odfId );

/*! Get the frequency band associated with a given ODF frequency band ID.
 *
 * @param odfId ODF frequency band ID
 * @return Frequency band
 */
observation_models::FrequencyBands getFrequencyBandForOdfId ( const int odfId );

/*! Get the ground station name associated with a given ODF network and station ID.
 *
 * @param networkId ODF network ID
 * @param stationId ODF ground station ID
 * @return Ground station name
 */
std::string getStationNameFromStationId ( const int networkId, const int stationId );

class ProcessedOdfFileSingleLinkData
{
public:

    ProcessedOdfFileSingleLinkData( observation_models::ObservableType observableType,
                                    std::string receivingStation ):
        receivingStation_( receivingStation ),
        observableType_( observableType )
    { }

    virtual ~ProcessedOdfFileSingleLinkData( ){ }

    std::vector< double > unprocessedObservationTimes_;
    std::vector< double > processedObservationTimes_;

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observableValues_;
    std::vector< double > receiverDownlinkDelays_;

    std::vector< int > downlinkBandIds_;
    std::vector< int > uplinkBandIds_;
    std::vector< int > referenceBandIds_;

    std::vector< std::string > originFiles_;

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

    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

private:

    observation_models::ObservableType observableType_;

};


class ProcessedOdfFileDopplerData: public ProcessedOdfFileSingleLinkData
{
public:

    ProcessedOdfFileDopplerData( observation_models::ObservableType observableType,
                                 std::string receivingStation,
                                 std::string transmittingStation ):
        ProcessedOdfFileSingleLinkData( observableType, receivingStation ),
        transmittingStation_( transmittingStation )
    { }

    ~ProcessedOdfFileDopplerData( ){ }

    std::string transmittingStation_;

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


bool compareRawOdfDataByStartDate( std::shared_ptr< input_output::OdfRawFileContents > rawOdfData1,
                                   std::shared_ptr< input_output::OdfRawFileContents > rawOdfData2 );

class ProcessedOdfFileContents
{
public:

    // Constructor for single raw ODF data object
    ProcessedOdfFileContents(
            const std::shared_ptr< input_output::OdfRawFileContents > rawOdfData,
            const std::shared_ptr< const simulation_setup::Body > bodyWithGroundStations,
            bool verbose = true,
            std::string spacecraftName = "" ):
        ProcessedOdfFileContents(
                std::vector< std::shared_ptr< input_output::OdfRawFileContents > >{ rawOdfData },
                bodyWithGroundStations, verbose, spacecraftName )
    { }

    // Constructor for multiple ODF data objects
    ProcessedOdfFileContents(
            std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector,
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
            // Spacecraft name selected to be "NAIF Id", which is equal to -"JPL Id" (for a spacecraft)
            spacecraftName_ = std::to_string( - rawOdfDataVector.front( )->spacecraftId_ );
        }

        // Sort ODF data files
        sortAndValidateOdfDataVector( rawOdfDataVector );

        // Extract and process raw ODF data
        extractMultipleRawOdfRampData( rawOdfDataVector );
        for ( unsigned int i = 0; i < rawOdfDataVector.size( ); ++i )
        {
            extractRawOdfOrbitData( rawOdfDataVector.at( i ) );
        }
        updateProcessedObservationTimes( );
    }

    std::string getSpacecraftName( )
    {
        return spacecraftName_;
    }

    std::vector< std::string > getGroundStationsNames( );

    std::vector< observation_models::ObservableType > getProcessedObservableTypes( );

    std::pair< double, double > getStartAndEndTime( );

    std::vector< int > getIgnoredRawOdfObservableTypes( )
    {
        return ignoredRawOdfObservableTypes_;
    }

    std::vector< std::string > getIgnoredGroundStations( )
    {
        return ignoredGroundStations_;
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

    void sortAndValidateOdfDataVector( std::vector< std::shared_ptr< input_output::OdfRawFileContents > >& rawOdfDataVector );

    bool isObservationValid( std::shared_ptr< input_output::OdfDataBlock > rawDataBlock,
                             observation_models::LinkEnds linkEnds,
                             observation_models::ObservableType currentObservableType );

    void extractRawOdfOrbitData( std::shared_ptr< input_output::OdfRawFileContents > rawOdfData );

    void extractMultipleRawOdfRampData(
            std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector );

    void addOdfRawDataBlockToProcessedData(
            const std::shared_ptr< input_output::OdfDataBlock > rawDataBlock,
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

    // Unprocessed ramp start and end times. Used for pre-processing observations
    std::map< std::string, std::vector< double > > unprocessedRampStartTimesPerStation_;
    std::map< std::string, std::vector< double > > unprocessedRampEndTimesPerStation_;

    // Vector keeping the invalid observable types that were found in the raw ODF data
    std::vector< int > ignoredRawOdfObservableTypes_;

    // Vector keeping the invalid ground stations that were found in the raw ODF data
    std::vector< std::string > ignoredGroundStations_;

    // Vector keeping the invalid ground stations that were found in the raw ODF data
    std::vector< std::shared_ptr< input_output::OdfDataBlock > > ignoredOdfRawDataBlocks_;

    // Flag indicating whether to print warnings
    bool verbose_;
};

observation_models::LinkEnds getLinkEndsFromOdfBlock (
        const std::shared_ptr< input_output::OdfDataBlock > dataBlock,
        std::string spacecraftName );

template< typename TimeType = double >
observation_models::ObservationAncilliarySimulationSettings createOdfAncillarySettings(
        std::shared_ptr< ProcessedOdfFileSingleLinkData > odfDataContents,
        unsigned int dataIndex )
{
    if ( dataIndex >= odfDataContents->unprocessedObservationTimes_.size( ) )
    {
        throw std::runtime_error("Error when creating ODF data ancillary settings: specified data index is larger than data size.");
    }

    observation_models::ObservationAncilliarySimulationSettings ancillarySettings =
            observation_models::ObservationAncilliarySimulationSettings( );

    observation_models::ObservableType currentObservableType = odfDataContents->getObservableType( );

    // Set common ancillary settings
    observation_models::FrequencyBands uplinkBand = getFrequencyBandForOdfId( odfDataContents->uplinkBandIds_.at( dataIndex ) );
    observation_models::FrequencyBands downlinkBand = getFrequencyBandForOdfId( odfDataContents->downlinkBandIds_.at( dataIndex ) );

    ancillarySettings.setAncilliaryDoubleData( observation_models::uplink_band, uplinkBand );
    ancillarySettings.setAncilliaryDoubleData( observation_models::downlink_band, downlinkBand);

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
            ancillarySettings.setAncilliaryDoubleVectorData(
                    observation_models::transmission_reception_delays, std::vector< double >{
                    dopplerDataBlock->transmitterUplinkDelays_.at( dataIndex ), 0.0,
                    dopplerDataBlock->receiverDownlinkDelays_.at( dataIndex ) } );
        }
        else
        {
            ancillarySettings.setAncilliaryDoubleVectorData(
                    observation_models::transmission_reception_delays, std::vector< double >{
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
        std::vector< observation_models::ObservationAncilliarySimulationSettings >& ancillarySettings )
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
        observation_models::ObservationAncilliarySimulationSettings currentAncillarySettings =
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
        std::vector< observation_models::ObservableType > observableTypesToProcess = std::vector< observation_models::ObservableType >( ) )
{
    // Set observables to process
    if ( observableTypesToProcess.empty( ) )
    {
        observableTypesToProcess = processedOdfFileContents->getProcessedObservableTypes( );
    }

    // Create and fill single observation sets
    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, std::vector< std::shared_ptr<
            observation_models::SingleObservationSet< ObservationScalarType, TimeType > > > > > sortedObservationSets;

    for ( auto observableTypeIterator = processedOdfFileContents->getProcessedDataBlocks( ).begin( );
            observableTypeIterator != processedOdfFileContents->getProcessedDataBlocks( ).end( ); ++observableTypeIterator )
    {
        observation_models::ObservableType currentObservableType = observableTypeIterator->first;

        // Check if an observation set should be created for the current observable type
        if ( std::count( observableTypesToProcess.begin( ), observableTypesToProcess.end( ), currentObservableType ) == 0 )
        {
            continue;
        }

        for ( auto linkEndsIterator = observableTypeIterator->second.begin( );
                linkEndsIterator != observableTypeIterator->second.end( ); ++linkEndsIterator )
        {
            observation_models::LinkEnds currentLinkEnds = linkEndsIterator->first;
            std::shared_ptr< ProcessedOdfFileSingleLinkData > currentOdfSingleLinkData = linkEndsIterator->second;

            // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
            std::vector< std::vector< TimeType > > observationTimes;
            std::vector< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > observables;
            std::vector< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings;

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
                        nullptr, std::make_shared< observation_models::ObservationAncilliarySimulationSettings >(
                            ancillarySettings.at( i ) ) ) );
            }
        }
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

void setOdfInformationInBodies(
        const std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents,
        simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyWithGroundStations = "Earth",
        const std::function< double ( FrequencyBands uplinkBand, FrequencyBands downlinkBand ) > getTurnaroundRatio =
                &getDsnDefaultTurnaroundRatios );

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_PROCESSODFFILE_H
