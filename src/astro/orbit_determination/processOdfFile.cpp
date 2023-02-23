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
    case 11:
        observableType = observation_models::dsn_one_way_averaged_doppler;
        break;
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

std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > mergeRampDataInterpolators(
        const std::vector< std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > >& interpolatorList )
{
    std::map< double, double > rampEndTimesPerStation;
    std::map< double, double > rampRatesPerStation;
    std::map< double, double > rampStartFrequenciesPerStation;

    for( unsigned int i = 0; i < interpolatorList.size( ); i++ )
    {
        for( unsigned int j = 0; j < interpolatorList.at( i )->getStartTimes().size( ); j++ )
        {
            rampEndTimesPerStation[ interpolatorList.at( i )->getStartTimes().at( j ) ] =
                    interpolatorList.at( i )->getEndTimes().at( j );
            rampRatesPerStation[ interpolatorList.at( i )->getStartTimes().at( j ) ] =
                    interpolatorList.at( i )->getRampRates().at( j );
            rampStartFrequenciesPerStation[ interpolatorList.at( i )->getStartTimes().at( j ) ] =
                    interpolatorList.at( i )->getStartFrequencies().at( j );
        }
    }
    return std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >(
                utilities::createVectorFromMapKeys( rampEndTimesPerStation ),
                utilities::createVectorFromMapValues( rampEndTimesPerStation ),
                utilities::createVectorFromMapValues( rampRatesPerStation ),
                utilities::createVectorFromMapValues( rampStartFrequenciesPerStation ) );
}

void addOdfFileContentsToMergedContents(
        const observation_models::ObservableType observableType,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > mergedOdfFileContents,
        std::shared_ptr< ProcessedOdfFileSingleLinkData > blockToAdd )
{
    mergedOdfFileContents->downlinkBandIds_.insert(
            mergedOdfFileContents->downlinkBandIds_.end( ),
            blockToAdd->downlinkBandIds_.begin( ), blockToAdd->downlinkBandIds_.end( ) );

    mergedOdfFileContents->observableValues_.insert(
            mergedOdfFileContents->observableValues_.end( ),
            blockToAdd->observableValues_.begin( ), blockToAdd->observableValues_.end( ) );

    mergedOdfFileContents->observationTimes_.insert(
            mergedOdfFileContents->observationTimes_.end( ),
            blockToAdd->observationTimes_.begin( ), blockToAdd->observationTimes_.end( ) );

    mergedOdfFileContents->originFiles_.insert(
            mergedOdfFileContents->originFiles_.end( ),
            blockToAdd->originFiles_.begin( ), blockToAdd->originFiles_.end( ) );

    if( observableType == observation_models::dsn_one_way_averaged_doppler ||
            observableType == observation_models::dsn_n_way_averaged_doppler )
    {
        std::shared_ptr< ProcessedOdfFileDopplerData > dopplerBlockToAdd
                = std::dynamic_pointer_cast< ProcessedOdfFileDopplerData >(
                    blockToAdd );
        std::shared_ptr< ProcessedOdfFileDopplerData > currentDopplerObservableMergedData =
                std::dynamic_pointer_cast< ProcessedOdfFileDopplerData >(
                    mergedOdfFileContents );


        currentDopplerObservableMergedData->referenceFrequencies_.insert(
                currentDopplerObservableMergedData->referenceFrequencies_.end( ),
                dopplerBlockToAdd->referenceFrequencies_.begin( ), dopplerBlockToAdd->referenceFrequencies_.end( ) );

        currentDopplerObservableMergedData->receiverRampingFlags_.insert(
                currentDopplerObservableMergedData->receiverRampingFlags_.end( ),
                dopplerBlockToAdd->receiverRampingFlags_.begin( ), dopplerBlockToAdd->receiverRampingFlags_.end( ) );

        currentDopplerObservableMergedData->countInterval_.insert(
                currentDopplerObservableMergedData->countInterval_.end( ),
                dopplerBlockToAdd->countInterval_.begin( ), dopplerBlockToAdd->countInterval_.end( ) );

        currentDopplerObservableMergedData->receiverChannels_.insert(
                currentDopplerObservableMergedData->receiverChannels_.end( ),
                dopplerBlockToAdd->receiverChannels_.begin( ), dopplerBlockToAdd->receiverChannels_.end( ) );
    }
}

std::shared_ptr< ProcessedOdfFileContents > mergeOdfFileContents(
        const std::vector< std::shared_ptr< ProcessedOdfFileContents > > odfFileContents )
{
    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds,
            std::shared_ptr< ProcessedOdfFileSingleLinkData > > > mergedOdfFileContents;
    std::map< std::string, std::vector< std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > > > rampInterpolatorList;

    // Iterate over all ODF files
    for( unsigned int i = 0; i < odfFileContents.size( ); i++ )
    {
        // Retrieve contents of current file.
        std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds,
                std::shared_ptr< ProcessedOdfFileSingleLinkData > > >  dataBlocks =
                odfFileContents.at( i )->processedDataBlocks_;

        for( auto it = dataBlocks.begin( ); it != dataBlocks.end( ); it++ )
        {
            for( auto linkIt = it->second.begin( ); linkIt != it->second.end( ); linkIt++ )
            {
                bool setNewObject = false;
                if( mergedOdfFileContents.count( it->first ) == 0 )
                {
                    setNewObject = true;
                }
                else if( mergedOdfFileContents.at( it->first ).count( linkIt->first ) == 0 )
                {
                    setNewObject = true;
                }

                if( setNewObject )
                {
                    mergedOdfFileContents[ it->first ][ linkIt->first ] = linkIt->second;
                }
                else
                {
                    addOdfFileContentsToMergedContents(
                                it->first, mergedOdfFileContents[ it->first ][ linkIt->first ], linkIt->second );
                }
            }
        }

        std::map< std::string, std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > > currentRampInterpolators =
                odfFileContents.at( i )->rampInterpolators_;

        for( auto it = currentRampInterpolators.begin( ); it != currentRampInterpolators.end( ); it++ )
        {
            rampInterpolatorList[ it->first ].push_back( it->second );
        }
    }

    std::map< std::string, std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > > mergedRampInterpolators;
    for( auto it = rampInterpolatorList.begin( ); it != rampInterpolatorList.end( ); it++ )
    {
        mergedRampInterpolators[ it->first ] = mergeRampDataInterpolators(
                    it->second );
    }


    std::shared_ptr< ProcessedOdfFileContents > processedOdfFile =
            std::make_shared< ProcessedOdfFileContents >( );

    processedOdfFile->processedDataBlocks_ = mergedOdfFileContents;
    processedOdfFile->rampInterpolators_ = mergedRampInterpolators;
    processedOdfFile->spacecraftName_ = "AAA";
    return processedOdfFile;
}

void addOdfDataBlockToProcessedData(
        const observation_models::ObservableType currentObservableType,
        const std::shared_ptr< input_output::OdfDataBlock > rawDataBlock,
        const std::shared_ptr< ProcessedOdfFileSingleLinkData > processedDataBlock )
{
    // Add properties to data block if data is valid
    if ( rawDataBlock->commonDataBlock_->validity_ == 0 )
    {

        // Add common properties to data object
        processedDataBlock->downlinkBandIds_.push_back( rawDataBlock->commonDataBlock_->downlinkBandId_ );
        processedDataBlock->uplinkBandIds_.push_back( rawDataBlock->commonDataBlock_->uplinkBandId_ );
        processedDataBlock->referenceBandIds_.push_back( rawDataBlock->commonDataBlock_->referenceBandId_ );
        processedDataBlock->observationTimes_.push_back( rawDataBlock->commonDataBlock_->getObservableTime( ) );
        processedDataBlock->receiverDownlinkDelays_.push_back( rawDataBlock->commonDataBlock_->getReceivingStationDownlinkDelay( ) );

        // Add properties to data object for Doppler data
        if ( currentObservableType == observation_models::dsn_one_way_averaged_doppler ||
             currentObservableType == observation_models::dsn_n_way_averaged_doppler )
        {
            std::shared_ptr< input_output::OdfDopplerDataBlock > odfDopplerDataBlock =
                    std::dynamic_pointer_cast< input_output::OdfDopplerDataBlock >(
                            rawDataBlock->observableSpecificDataBlock_ );
            std::shared_ptr< ProcessedOdfFileDopplerData > odfParsedDopplerDataBlock =
                    std::dynamic_pointer_cast< ProcessedOdfFileDopplerData >(
                            processedDataBlock );

            processedDataBlock->observableValues_.push_back(
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

std::shared_ptr< ProcessedOdfFileContents > processOdfFileContents(
        const std::shared_ptr< input_output::OdfRawFileContents > rawOdfData,
        bool verbose )
{
    // Create output object
    std::shared_ptr< ProcessedOdfFileContents > processedOdfFile = std::make_shared< ProcessedOdfFileContents >( );
    std::string spacecraftName = std::to_string( rawOdfData->spacecraftId_ );

    // Retrieve data blocks from ODF file raw contents
    std::vector< std::shared_ptr< input_output::OdfDataBlock > > rawDataBlocks = rawOdfData->dataBlocks_;

    // Create list of data, sorted by observable type and link ends; single object per combination of the two
    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds,
            std::shared_ptr< ProcessedOdfFileSingleLinkData > > > processedDataBlocks;

    // Vector keeping the invalid observable types that were found
    std::vector< int > invalidObservableTypes;

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
            if ( verbose )
            {
                if ( std::find( invalidObservableTypes.begin( ), invalidObservableTypes.end( ), currentObservableId )
                    == invalidObservableTypes.end( ) )
                {
                    invalidObservableTypes.push_back( currentObservableId );
                    std::cerr << "Warning: processing of ODF data type " << currentObservableId << " is not implemented, ignoring the corresponding data." << std::endl;
                }
            }
            continue;
        }

        observation_models::LinkEnds linkEnds = getLinkEndsFromOdfBlock(
                rawDataBlocks.at( i ), spacecraftName );

        // Check if data object already exists for current observable/link ends
        bool createNewObject = false;
        if( processedDataBlocks.count( currentObservableType ) == 0 )
        {
            createNewObject = true;
        }
        else if( processedDataBlocks.at( currentObservableType ).count( linkEnds ) == 0 )
        {
            createNewObject = true;
        }

        // Create new data object, if required
        if( createNewObject )
        {
            if( currentObservableType == observation_models::dsn_one_way_averaged_doppler ||
                    currentObservableType == observation_models::dsn_n_way_averaged_doppler )
            {
                processedDataBlocks[ currentObservableType ][ linkEnds ] = std::make_shared< ProcessedOdfFileDopplerData >( );
            }
            else
            {
                throw std::runtime_error( "Error when processing ODF file contents: processing of the selected observable type (" +
                    std::to_string( currentObservableType ) + ") is not implemented.");
            }

            processedDataBlocks[ currentObservableType ][ linkEnds ]->transmittingStation_ =
                    getStationNameFromStationId( rawDataBlocks.at( i )->commonDataBlock_->transmittingStationNetworkId_,
                                                 rawDataBlocks.at( i )->commonDataBlock_->transmittingStationId_ );
            processedDataBlocks[ currentObservableType ][ linkEnds ]->receivingStation_ =
                    getStationNameFromStationId( 0, rawDataBlocks.at( i )->commonDataBlock_->receivingStationId_ );
            processedDataBlocks[ currentObservableType ][ linkEnds ]->observableType_ = currentObservableType;
        }

        addOdfDataBlockToProcessedData(
                currentObservableType, rawDataBlocks.at( i ),
                processedDataBlocks[ currentObservableType ][ linkEnds ] );
    }

    // Save output
    processedOdfFile->processedDataBlocks_ = processedDataBlocks;
    processedOdfFile->spacecraftName_ = spacecraftName;

    // Process ramp blocks
    std::map< std::string, std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > > rampInterpolators;

    for( auto it = rawOdfData->rampBlocks_.begin( ); it != rawOdfData->rampBlocks_.end( ); it++ )
    {
        std::string stationName = getStationNameFromStationId( 0, it->first );
        rampInterpolators[ stationName ] = std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >( it->second );
    }
    processedOdfFile->rampInterpolators_ = rampInterpolators;

    return processedOdfFile;
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

} // namespace orbit_determination

} // namespace tudat
