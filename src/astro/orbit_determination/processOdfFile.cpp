/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>

#include "tudat/astro/orbit_determination/processOdfFile.h"

namespace tudat
{

namespace orbit_determination
{

observation_models::ObservableType getObservableTypeForOdfId(
        const int odfId )
{
    observation_models::ObservableType observableType;

    switch( odfId )
    {
//    case 11:
//        observableType = observation_models::dsn_one_way_averaged_doppler;
//        break;
    case 12:
        observableType = observation_models::dsn_n_way_averaged_doppler;
        break;
    case 13:
        observableType = observation_models::dsn_n_way_averaged_doppler;
        break;
    default:
        throw std::runtime_error( "Error when getting observable type for ODF ID, ID: " +
                                  std::to_string( odfId ) + " not recognized." );
    }

    return observableType;
}

observation_models::FrequencyBands getFrequencyBandForOdfId ( const int odfId )
{
    observation_models::FrequencyBands frequencyBand;

    switch( odfId )
    {
    case 0:
        frequencyBand = observation_models::ku_band;
        break;
    case 1:
        frequencyBand = observation_models::s_band;
        break;
    case 2:
        frequencyBand = observation_models::x_band;
        break;
    case 3:
        frequencyBand = observation_models::ka_band;
        break;
    default:
        throw std::runtime_error( "Error when getting observable type for ODF ID, ID: " +
                                  std::to_string( odfId ) + " not recognized." );
    }

    return frequencyBand;
}

std::string getStationNameFromStationId ( const int networkId, const int stationId )
{
    std::string stationName;

    if ( networkId == 0 )
    {
        stationName = "DSS-" + std::to_string( stationId );
    }
    else if ( networkId == 3 )
    {
        stationName = "UPL-" + std::to_string( stationId );
    }
    else
    {
        stationName = "Station-" + std::to_string( stationId );
    }

    return stationName;
}


std::vector< std::string > ProcessedOdfFileContents::getGroundStationsNames( )
{
    std::vector< std::string > groundStations;

    for ( auto observableIt = processedDataBlocks_.begin( ); observableIt != processedDataBlocks_.end( );
            ++observableIt )
    {
        for ( auto linkEndIt = observableIt->second.begin( ); linkEndIt != observableIt->second.end( );
                ++linkEndIt )
        {
            observation_models::LinkEnds linkEnd = linkEndIt->first;

            for ( auto linkEndTypeIt = linkEnd.begin( ); linkEndTypeIt != linkEnd.end( ); ++linkEndTypeIt )
            {
                // Check if linkEndId is a ground station
                if ( linkEndTypeIt->second.stationName_ != "" )
                {
                    if ( !std::count( groundStations.begin( ), groundStations.end( ), linkEndTypeIt->second.stationName_)  )
                    {
                        groundStations.push_back( linkEndTypeIt->second.stationName_ );
                    }
                }
            }
        }
    }

    return groundStations;
}

std::vector< observation_models::ObservableType > ProcessedOdfFileContents::getProcessedObservableTypes( )
{
    std::vector< observation_models::ObservableType > observableTypes;

    for ( auto observableIt = processedDataBlocks_.begin( ); observableIt != processedDataBlocks_.end( );
            ++observableIt )
    {
        observableTypes.push_back( observableIt->first );
    }

    return observableTypes;
}

std::pair< double, double > ProcessedOdfFileContents::getStartAndEndTime( const simulation_setup::SystemOfBodies& bodies )
    {
        // Reset variables
        double startTimeTdbSinceJ2000 = TUDAT_NAN;
        double endTimeTdbSinceJ2000 = TUDAT_NAN;

        // Loop over data
        for ( auto observableIt = processedDataBlocks_.begin( ); observableIt != processedDataBlocks_.end( );
                ++observableIt )
        {
            for ( auto linkEndIt = observableIt->second.begin( ); linkEndIt != observableIt->second.end( );
                  ++linkEndIt )
            {
                std::shared_ptr< ProcessedOdfFileSingleLinkData > processedSingleLinkData = linkEndIt->second;

                // Extract the start and end times
                std::vector< double > timeVector = processedSingleLinkData->processedObservationTimes_;

                if ( timeVector.front( ) < startTimeTdbSinceJ2000 || std::isnan( startTimeTdbSinceJ2000 ) )
                {
                    startTimeTdbSinceJ2000 = timeVector.front( );
                }

                if ( timeVector.back( ) > endTimeTdbSinceJ2000 || std::isnan( endTimeTdbSinceJ2000 ) )
                {
                    endTimeTdbSinceJ2000 = timeVector.back( );
                }
            }
        }

        return std::make_pair( startTimeTdbSinceJ2000, endTimeTdbSinceJ2000 );
    }

//std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > mergeRampDataInterpolators(
//        const std::vector< std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > >& interpolatorList )
//{
//    std::map< double, double > rampEndTimesPerStation;
//    std::map< double, double > rampRatesPerStation;
//    std::map< double, double > rampStartFrequenciesPerStation;
//
//    for( unsigned int i = 0; i < interpolatorList.size( ); i++ )
//    {
//        for( unsigned int j = 0; j < interpolatorList.at( i )->getStartTimes().size( ); j++ )
//        {
//            rampEndTimesPerStation[ interpolatorList.at( i )->getStartTimes().at( j ) ] =
//                    interpolatorList.at( i )->getEndTimes().at( j );
//            rampRatesPerStation[ interpolatorList.at( i )->getStartTimes().at( j ) ] =
//                    interpolatorList.at( i )->getRampRates().at( j );
//            rampStartFrequenciesPerStation[ interpolatorList.at( i )->getStartTimes().at( j ) ] =
//                    interpolatorList.at( i )->getStartFrequencies().at( j );
//        }
//    }
//    return std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >(
//                utilities::createVectorFromMapKeys( rampEndTimesPerStation ),
//                utilities::createVectorFromMapValues( rampEndTimesPerStation ),
//                utilities::createVectorFromMapValues( rampRatesPerStation ),
//                utilities::createVectorFromMapValues( rampStartFrequenciesPerStation ) );
//}
//
//void addOdfFileContentsToMergedContents(
//        const observation_models::ObservableType observableType,
//        std::shared_ptr< ProcessedOdfFileSingleLinkData > mergedOdfFileContents,
//        std::shared_ptr< ProcessedOdfFileSingleLinkData > blockToAdd )
//{
//    mergedOdfFileContents->downlinkBandIds_.insert(
//            mergedOdfFileContents->downlinkBandIds_.end( ),
//            blockToAdd->downlinkBandIds_.begin( ), blockToAdd->downlinkBandIds_.end( ) );
//
//    mergedOdfFileContents->observableValues_.insert(
//            mergedOdfFileContents->observableValues_.end( ),
//            blockToAdd->observableValues_.begin( ), blockToAdd->observableValues_.end( ) );
//
//    mergedOdfFileContents->unprocessedObservationTimes_.insert(
//            mergedOdfFileContents->unprocessedObservationTimes_.end( ),
//            blockToAdd->unprocessedObservationTimes_.begin( ), blockToAdd->unprocessedObservationTimes_.end( ) );
//
//    mergedOdfFileContents->originFiles_.insert(
//            mergedOdfFileContents->originFiles_.end( ),
//            blockToAdd->originFiles_.begin( ), blockToAdd->originFiles_.end( ) );
//
//    if( observableType == observation_models::dsn_one_way_averaged_doppler ||
//            observableType == observation_models::dsn_n_way_averaged_doppler )
//    {
//        std::shared_ptr< ProcessedOdfFileDopplerData > dopplerBlockToAdd
//                = std::dynamic_pointer_cast< ProcessedOdfFileDopplerData >(
//                    blockToAdd );
//        std::shared_ptr< ProcessedOdfFileDopplerData > currentDopplerObservableMergedData =
//                std::dynamic_pointer_cast< ProcessedOdfFileDopplerData >(
//                    mergedOdfFileContents );
//
//
//        currentDopplerObservableMergedData->referenceFrequencies_.insert(
//                currentDopplerObservableMergedData->referenceFrequencies_.end( ),
//                dopplerBlockToAdd->referenceFrequencies_.begin( ), dopplerBlockToAdd->referenceFrequencies_.end( ) );
//
//        currentDopplerObservableMergedData->receiverRampingFlags_.insert(
//                currentDopplerObservableMergedData->receiverRampingFlags_.end( ),
//                dopplerBlockToAdd->receiverRampingFlags_.begin( ), dopplerBlockToAdd->receiverRampingFlags_.end( ) );
//
//        currentDopplerObservableMergedData->countInterval_.insert(
//                currentDopplerObservableMergedData->countInterval_.end( ),
//                dopplerBlockToAdd->countInterval_.begin( ), dopplerBlockToAdd->countInterval_.end( ) );
//
//        currentDopplerObservableMergedData->receiverChannels_.insert(
//                currentDopplerObservableMergedData->receiverChannels_.end( ),
//                dopplerBlockToAdd->receiverChannels_.begin( ), dopplerBlockToAdd->receiverChannels_.end( ) );
//    }
//}
//
//std::shared_ptr< ProcessedOdfFileContents > mergeOdfFileContents(
//        const std::vector< std::shared_ptr< ProcessedOdfFileContents > > odfFileContents )
//{
//    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds,
//            std::shared_ptr< ProcessedOdfFileSingleLinkData > > > mergedOdfFileContents;
//    std::map< std::string, std::vector< std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > > > rampInterpolatorList;
//
//    // Iterate over all ODF files
//    for( unsigned int i = 0; i < odfFileContents.size( ); i++ )
//    {
//        // Retrieve contents of current file.
//        std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds,
//                std::shared_ptr< ProcessedOdfFileSingleLinkData > > >  dataBlocks =
//                odfFileContents.at( i )->processedDataBlocks_;
//
//        for( auto it = dataBlocks.begin( ); it != dataBlocks.end( ); it++ )
//        {
//            for( auto linkIt = it->second.begin( ); linkIt != it->second.end( ); linkIt++ )
//            {
//                bool setNewObject = false;
//                if( mergedOdfFileContents.count( it->first ) == 0 )
//                {
//                    setNewObject = true;
//                }
//                else if( mergedOdfFileContents.at( it->first ).count( linkIt->first ) == 0 )
//                {
//                    setNewObject = true;
//                }
//
//                if( setNewObject )
//                {
//                    mergedOdfFileContents[ it->first ][ linkIt->first ] = linkIt->second;
//                }
//                else
//                {
//                    addOdfFileContentsToMergedContents(
//                                it->first, mergedOdfFileContents[ it->first ][ linkIt->first ], linkIt->second );
//                }
//            }
//        }
//
//        std::map< std::string, std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > > currentRampInterpolators =
//                odfFileContents.at( i )->rampInterpolators_;
//
//        for( auto it = currentRampInterpolators.begin( ); it != currentRampInterpolators.end( ); it++ )
//        {
//            rampInterpolatorList[ it->first ].push_back( it->second );
//        }
//    }
//
//    std::map< std::string, std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > > mergedRampInterpolators;
//    for( auto it = rampInterpolatorList.begin( ); it != rampInterpolatorList.end( ); it++ )
//    {
//        mergedRampInterpolators[ it->first ] = mergeRampDataInterpolators(
//                    it->second );
//    }
//
//
//    std::shared_ptr< ProcessedOdfFileContents > processedOdfFile =
//            std::make_shared< ProcessedOdfFileContents >( );
//
//    processedOdfFile->processedDataBlocks_ = mergedOdfFileContents;
//    processedOdfFile->rampInterpolators_ = mergedRampInterpolators;
//    processedOdfFile->spacecraftName_ = "AAA";
//    return processedOdfFile;
//}

observation_models::LinkEnds getLinkEndsFromOdfBlock (
        const std::shared_ptr< input_output::OdfDataBlock > dataBlock,
        std::string spacecraftName )
{
    int currentObservableId = dataBlock->observableSpecificDataBlock_->dataType_;

    observation_models::LinkEnds linkEnds;

    if ( currentObservableId == 11 )
    {
        linkEnds[ observation_models::transmitter ] = observation_models::LinkEndId ( spacecraftName );
        linkEnds[ observation_models::receiver ] = observation_models::LinkEndId ( "Earth", getStationNameFromStationId(
                0, dataBlock->commonDataBlock_->receivingStationId_ ) );
    }
    else if ( currentObservableId == 12 )
    {
        linkEnds[ observation_models::transmitter ] = observation_models::LinkEndId (
                "Earth", getStationNameFromStationId( dataBlock->commonDataBlock_->transmittingStationNetworkId_,
                                                      dataBlock->commonDataBlock_->transmittingStationId_ ) );
        linkEnds[ observation_models::reflector1 ] = observation_models::LinkEndId ( spacecraftName );
        linkEnds[ observation_models::receiver ] = observation_models::LinkEndId (
                "Earth", getStationNameFromStationId( 0, dataBlock->commonDataBlock_->receivingStationId_ ) );
    }
    else if ( currentObservableId == 13 )
    {
        linkEnds[ observation_models::transmitter ] = observation_models::LinkEndId (
                "Earth", getStationNameFromStationId( dataBlock->commonDataBlock_->transmittingStationNetworkId_,
                                                      dataBlock->commonDataBlock_->transmittingStationId_ ) );
        linkEnds[ observation_models::reflector1 ] = observation_models::LinkEndId ( spacecraftName );
        linkEnds[ observation_models::receiver ] = observation_models::LinkEndId (
                "Earth", getStationNameFromStationId( 0, dataBlock->commonDataBlock_->receivingStationId_ ) );
    }
    else
    {
        throw std::runtime_error(
                "Error when getting link ends from ODF data blocks, data type " +
                std::to_string( currentObservableId ) + " not recognized." );
    }

    return linkEnds;
}

void ProcessedOdfFileContents::addOdfRawDataBlockToProcessedData(
        const observation_models::ObservableType currentObservableType,
        const std::shared_ptr< const input_output::OdfDataBlock > rawDataBlock,
        const std::shared_ptr< ProcessedOdfFileSingleLinkData > singleLinkProcessedData,
        const std::string rawDataFileName )
{
    // Add properties to data block if data is valid
    if ( rawDataBlock->commonDataBlock_->validity_ == 0 )
    {

        // Add common properties to data object
        singleLinkProcessedData->downlinkBandIds_.push_back( rawDataBlock->commonDataBlock_->downlinkBandId_ );
        singleLinkProcessedData->uplinkBandIds_.push_back( rawDataBlock->commonDataBlock_->uplinkBandId_ );
        singleLinkProcessedData->referenceBandIds_.push_back( rawDataBlock->commonDataBlock_->referenceBandId_ );
        singleLinkProcessedData->unprocessedObservationTimes_.push_back( rawDataBlock->commonDataBlock_->getObservableTime( ) );
        singleLinkProcessedData->receiverDownlinkDelays_.push_back( rawDataBlock->commonDataBlock_->getReceivingStationDownlinkDelay( ) );
        singleLinkProcessedData->originFiles_.push_back( rawDataFileName );

        // Add properties to data object for Doppler data
        if ( currentObservableType == observation_models::dsn_one_way_averaged_doppler ||
             currentObservableType == observation_models::dsn_n_way_averaged_doppler )
        {
            std::shared_ptr< input_output::OdfDopplerDataBlock > odfDopplerDataBlock =
                    std::dynamic_pointer_cast< input_output::OdfDopplerDataBlock >(
                            rawDataBlock->observableSpecificDataBlock_ );
            std::shared_ptr< ProcessedOdfFileDopplerData > odfParsedDopplerDataBlock =
                    std::dynamic_pointer_cast< ProcessedOdfFileDopplerData >(
                            singleLinkProcessedData );

            singleLinkProcessedData->observableValues_.push_back(
                    ( Eigen::Matrix< double, 1, 1 >( ) << rawDataBlock->commonDataBlock_->getObservableValue( )
                    ).finished( ) );

            odfParsedDopplerDataBlock->countInterval_.push_back( odfDopplerDataBlock->getCompressionTime( ) );
            odfParsedDopplerDataBlock->receiverChannels_.push_back( odfDopplerDataBlock->receiverChannel_ );
            odfParsedDopplerDataBlock->receiverRampingFlags_.push_back( odfDopplerDataBlock->receiverExciterFlag_ );
            odfParsedDopplerDataBlock->referenceFrequencies_.push_back( odfDopplerDataBlock->getReferenceFrequency( ) );
            odfParsedDopplerDataBlock->transmitterUplinkDelays_.push_back( odfDopplerDataBlock->transmittingStationUplinkDelay_ );
        }
    }
}

void ProcessedOdfFileContents::updateProcessedObservationTimes( )
{
    // Loop over saved data and convert time to TDB wrt J2000
    for ( auto observableTypeIterator = processedDataBlocks_.begin( );
            observableTypeIterator != processedDataBlocks_.end( ); ++observableTypeIterator )
    {
        observation_models::ObservableType currentObservableType = observableTypeIterator->first;

        for ( auto linkEndsIterator = observableTypeIterator->second.begin( );
              linkEndsIterator != observableTypeIterator->second.end( ); ++linkEndsIterator )
        {
            linkEndsIterator->second->processedObservationTimes_ = computeObservationTimesTdbFromJ2000(
                    linkEndsIterator->second->transmittingStation_,
                    linkEndsIterator->second->unprocessedObservationTimes_ );
        }
    }
}

std::vector< double > ProcessedOdfFileContents::computeObservationTimesUtcFromJ2000(
            std::vector< double > observationTimesUtcFromEME1950 )
    {
        std::vector < double > observationTimesUtcFromJ2000;

        double EME1950ToJ2000Offset = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch< double >(
                1950, 1, 1, 0, 0, 0, basic_astrodynamics::JULIAN_DAY_ON_J2000 )
                        * physical_constants::JULIAN_DAY;

        for ( unsigned int i = 0; i < observationTimesUtcFromEME1950.size( ); ++i )
        {
            observationTimesUtcFromJ2000.push_back( observationTimesUtcFromEME1950.at( i ) + EME1950ToJ2000Offset );
        }

        return observationTimesUtcFromJ2000;
    }

std::vector< double > ProcessedOdfFileContents::computeObservationTimesTdbFromJ2000(
            std::string groundStation,
            std::vector< double > observationTimesUtcFromEME1950 )
    {
        earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter =
                earth_orientation::TerrestrialTimeScaleConverter( );

        std::vector< double > observationTimesUtcFromJ2000 = computeObservationTimesUtcFromJ2000( observationTimesUtcFromEME1950 );

        std::vector< Eigen::Vector3d > earthFixedPositions;
        for ( unsigned int i = 0; i < observationTimesUtcFromJ2000.size( ); ++i )
        {
            // Approximation: UTC time used to retrieve the ground station's position
            earthFixedPositions.push_back(
                    bodyWithGroundStations_->getGroundStation( groundStation )->getStateInPlanetFixedFrame< double, double >(
                            observationTimesUtcFromJ2000.at( i ) ).segment( 0, 3 ) );
        }

        std::vector< double > observationTimesTdbFromJ2000 = timeScaleConverter.getCurrentTimes(
                basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, observationTimesUtcFromJ2000,
                earthFixedPositions );

        return observationTimesTdbFromJ2000;
    }

void ProcessedOdfFileContents::extractRawOdfOrbitData(
        std::shared_ptr< const input_output::OdfRawFileContents > rawOdfData )
{
    // Retrieve data blocks from ODF file raw contents
    std::vector< std::shared_ptr< input_output::OdfDataBlock > > rawDataBlocks = rawOdfData->dataBlocks_;

    // Iterate over all block of ODF file.
    for( unsigned int i = 0; i < rawDataBlocks.size( ); i++ )
    {
        // Retrieve observable type and link ends
        int currentObservableId = rawDataBlocks.at( i )->observableSpecificDataBlock_->dataType_;

        // Get current observable type and throw warning if not implemented
        observation_models::ObservableType currentObservableType;
        try
        {
            currentObservableType = getObservableTypeForOdfId( currentObservableId );
        }
        catch( const std::runtime_error& )
        {
            if ( verbose_ )
            {
                if ( std::find( notReadRawOdfObservableTypes_.begin( ), notReadRawOdfObservableTypes_.end( ),
                                currentObservableId ) == notReadRawOdfObservableTypes_.end( ) )
                {
                    notReadRawOdfObservableTypes_.push_back( currentObservableId );
                    std::cerr << "Warning: processing of ODF data type " << currentObservableId <<
                        " is not implemented, ignoring the corresponding data." << std::endl;
                }
            }
            continue;
        }

        observation_models::LinkEnds linkEnds = getLinkEndsFromOdfBlock(
                rawDataBlocks.at( i ), spacecraftName_ );

        // Check if data object already exists for current observable/link ends
        bool createNewObject = false;
        if( processedDataBlocks_.count( currentObservableType ) == 0 )
        {
            createNewObject = true;
        }
        else if( processedDataBlocks_.at( currentObservableType ).count( linkEnds ) == 0 )
        {
            createNewObject = true;
        }

        // Create new data object, if required
        if( createNewObject )
        {
            if( currentObservableType == observation_models::dsn_one_way_averaged_doppler ||
                    currentObservableType == observation_models::dsn_n_way_averaged_doppler )
            {
                processedDataBlocks_[ currentObservableType ][ linkEnds ] = std::make_shared< ProcessedOdfFileDopplerData >( );
            }
            else
            {
                throw std::runtime_error( "Error when processing ODF file contents: processing of the selected observable type (" +
                    std::to_string( currentObservableType ) + ") is not implemented.");
            }

            processedDataBlocks_[ currentObservableType ][ linkEnds ]->transmittingStation_ =
                    getStationNameFromStationId( rawDataBlocks.at( i )->commonDataBlock_->transmittingStationNetworkId_,
                                                 rawDataBlocks.at( i )->commonDataBlock_->transmittingStationId_ );
            processedDataBlocks_[ currentObservableType ][ linkEnds ]->receivingStation_ =
                    getStationNameFromStationId( 0, rawDataBlocks.at( i )->commonDataBlock_->receivingStationId_ );
            processedDataBlocks_[ currentObservableType ][ linkEnds ]->observableType_ = currentObservableType;
        }

        addOdfRawDataBlockToProcessedData(
                currentObservableType, rawDataBlocks.at( i ),
                processedDataBlocks_[ currentObservableType ][ linkEnds ],
                rawOdfData->fileName_ );
    }

}

void ProcessedOdfFileContents::extractRawOdfRampData( std::shared_ptr< const input_output::OdfRawFileContents > rawOdfData )
{

    for( auto it = rawOdfData->rampBlocks_.begin( ); it != rawOdfData->rampBlocks_.end( ); it++ )
    {
        std::string stationName = getStationNameFromStationId( 0, it->first );

        std::vector< std::shared_ptr< input_output::OdfRampBlock > > rampBlocks = it->second;

        std::vector< double > startTimes, endTimes, rampRates, startFrequencies;
        for( unsigned int i = 0; i < it->second.size( ); i++ )
        {
            startTimes.push_back( rampBlocks.at( i )->getRampStartTime( ) );
            endTimes.push_back( rampBlocks.at( i )->getRampEndTime( ) );
            rampRates.push_back( rampBlocks.at( i )->getRampRate( ) );
            startFrequencies.push_back( rampBlocks.at( i )->getRampStartFrequency( ) );
        }

        rampInterpolators_[ stationName ] = std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >(
                computeObservationTimesTdbFromJ2000( stationName, startTimes ),
                computeObservationTimesTdbFromJ2000( stationName, endTimes ),
                rampRates, startFrequencies );
    }

}

void ProcessedOdfFileContents::extractMultipleRawOdfRampData(
        std::vector< std::shared_ptr< const input_output::OdfRawFileContents > > rawOdfDataVector )
{

    std::map< std::string, std::vector< double > > startTimesPerStation, endTimesPerStation,
        rampRatesPerStation, startFrequenciesPerStation;

    for( unsigned int i = 0; i < rawOdfDataVector.size( ); ++i )
    {
        for( auto it = rawOdfDataVector.at( i )->rampBlocks_.begin( );
                it != rawOdfDataVector.at( i )->rampBlocks_.end( ); it++ )
        {
            std::string stationName = getStationNameFromStationId( 0, it->first );

            std::vector< std::shared_ptr< input_output::OdfRampBlock > > rampBlocks = it->second;

            for( unsigned int j = 0; i < it->second.size( ); j++ )
            {
                startTimesPerStation[ stationName ].push_back( rampBlocks.at( j )->getRampStartTime( ) );
                endTimesPerStation[ stationName ].push_back( rampBlocks.at( j )->getRampEndTime( ) );
                rampRatesPerStation[ stationName ].push_back( rampBlocks.at( j )->getRampRate( ) );
                startFrequenciesPerStation[ stationName ].push_back( rampBlocks.at( j )->getRampStartFrequency( ) );
            }
        }
    }

    for( auto it = startTimesPerStation.begin( ); it != startTimesPerStation.end( ); ++it )
    {
        std::string stationName = it->first;
        rampInterpolators_[ stationName ] = std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >(
                computeObservationTimesTdbFromJ2000( stationName, startTimesPerStation[ stationName ] ),
                computeObservationTimesTdbFromJ2000( stationName, endTimesPerStation[ stationName ] ),
                rampRatesPerStation[ stationName ], startFrequenciesPerStation[ stationName ] );
    }

}

} // namespace orbit_determination

} // namespace tudat
