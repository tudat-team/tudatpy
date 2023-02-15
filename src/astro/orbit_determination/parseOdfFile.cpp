/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/parseOdfFile.h"
#include "tudat/simulation/estimation_setup/observations.h"

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
        observableType = observation_models::one_way_differenced_range;
        break;
    case 12:
        observableType = observation_models::n_way_differenced_range;
        break;
    case 13:
        observableType = observation_models::n_way_differenced_range;
        break;
    default:
        throw std::runtime_error( "Error when getting observable type for ODF ID, ID: " +
                                  std::to_string( odfId ) + " not recognized." );
    }
    return observableType;
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

std::shared_ptr< RampedReferenceFrequencyInterpolator > mergeRampDataInterpolators(
        const std::vector< std::shared_ptr< RampedReferenceFrequencyInterpolator > >& interpolatorList )
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
    return std::make_shared< RampedReferenceFrequencyInterpolator >(
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
    mergedOdfFileContents->downlinkBand.insert(
                mergedOdfFileContents->downlinkBand.end( ),
                blockToAdd->downlinkBand.begin( ), blockToAdd->downlinkBand.end( ) );

    mergedOdfFileContents->observableValues.insert(
                mergedOdfFileContents->observableValues.end( ),
                blockToAdd->observableValues.begin( ), blockToAdd->observableValues.end( ) );

    mergedOdfFileContents->observationTimes.insert(
                mergedOdfFileContents->observationTimes.end( ),
                blockToAdd->observationTimes.begin( ), blockToAdd->observationTimes.end( ) );

    mergedOdfFileContents->originFile.insert(
                mergedOdfFileContents->originFile.end( ),
                blockToAdd->originFile.begin( ), blockToAdd->originFile.end( ) );

    if( observableType == observation_models::one_way_differenced_range ||
            observableType == observation_models::n_way_differenced_range )
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

        currentDopplerObservableMergedData->compressionTimes_.insert(
                currentDopplerObservableMergedData->compressionTimes_.end( ),
                dopplerBlockToAdd->compressionTimes_.begin( ), dopplerBlockToAdd->compressionTimes_.end( ) );

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
    std::map< int, std::vector< std::shared_ptr< RampedReferenceFrequencyInterpolator > > > rampInterpolatorList;

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

        std::map< int, std::shared_ptr< RampedReferenceFrequencyInterpolator > > currentRampInterpolators =
                odfFileContents.at( i )->rampInterpolators_;

        for( auto it = currentRampInterpolators.begin( ); it != currentRampInterpolators.end( ); it++ )
        {
            rampInterpolatorList[ it->first ].push_back( it->second );
        }
    }

    std::map< int, std::shared_ptr< RampedReferenceFrequencyInterpolator > > mergedRampInterpolators;
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
        processedDataBlock->downlinkBand.push_back( rawDataBlock->commonDataBlock_->downlinkBandId_ );
        processedDataBlock->uplinkBand.push_back( rawDataBlock->commonDataBlock_->uplinkBandId_ );
        processedDataBlock->referenceBand.push_back( rawDataBlock->commonDataBlock_->referenceBandId_ );
        processedDataBlock->observableValues.push_back( rawDataBlock->commonDataBlock_->getObservableValue( ) );
        processedDataBlock->observationTimes.push_back( rawDataBlock->commonDataBlock_->getObservableTime( ) );
        processedDataBlock->receiverDownlinkDelay.push_back( rawDataBlock->commonDataBlock_->getReceivingStationDownlinkDelay( ) );

        // Add properties to data object for Doppler data
        if ( currentObservableType == observation_models::one_way_differenced_range ||
             currentObservableType == observation_models::n_way_differenced_range )
        {
            std::shared_ptr< input_output::OdfDopplerDataBlock > odfDopplerDataBlock =
                    std::dynamic_pointer_cast< input_output::OdfDopplerDataBlock >(
                            rawDataBlock->observableSpecificDataBlock_ );
            std::shared_ptr< ProcessedOdfFileDopplerData > odfParsedDopplerDataBlock =
                    std::dynamic_pointer_cast< ProcessedOdfFileDopplerData >(
                            processedDataBlock );

            odfParsedDopplerDataBlock->compressionTimes_.push_back( odfDopplerDataBlock->getCompressionTime( ) );
            odfParsedDopplerDataBlock->receiverChannels_.push_back( odfDopplerDataBlock->receiverChannel_ );
            odfParsedDopplerDataBlock->receiverRampingFlags_.push_back( odfDopplerDataBlock->receiverExciterFlag_ );
            odfParsedDopplerDataBlock->referenceFrequencies_.push_back( odfDopplerDataBlock->getReferenceFrequency( ) );
            odfParsedDopplerDataBlock->uplinkDelays_.push_back( odfDopplerDataBlock->transmittingStationUplinkDelay_ );
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
                "Earth", getStationNameFromStationId( dataBlock->commonDataBlock_->transmittingStationNetworkId_,
                                                      dataBlock->commonDataBlock_->receivingStationId_ ) );
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
                "Error when getting link definition from ODF data blocks, data type " +
                std::to_string( currentObservableId ) + " not recognized." );
    }

    return linkEnds;
}

std::shared_ptr< ProcessedOdfFileContents > processOdfFileContents(
        const std::shared_ptr< input_output::OdfRawFileContents > rawOdfData )
{
    // Create output object
    std::shared_ptr< ProcessedOdfFileContents > processedOdfFile = std::make_shared< ProcessedOdfFileContents >( );
    std::string spacecraftName = std::to_string( rawOdfData->spacecraftId_ );

    // Retrieve data blocks from ODF file raw contents
    std::vector< std::shared_ptr< input_output::OdfDataBlock > > rawDataBlocks = rawOdfData->dataBlocks_;

    // Create list of data, sorted by observable type and link ends; single object per combination of the two
    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds,
            std::shared_ptr< ProcessedOdfFileSingleLinkData > > > processedDataBlocks;

    int currentObservableId;
    observation_models::ObservableType currentObservableType;

    // Iterate over all block of ODF file.
    for( unsigned int i = 0; i < rawDataBlocks.size( ); i++ )
    {
        // Retrieve observable type and link ends
        currentObservableId = rawDataBlocks.at( i )->observableSpecificDataBlock_->dataType_;
        currentObservableType = getObservableTypeForOdfId( currentObservableId );
        int appendedTransmittingStationId =
                rawDataBlocks.at( i )->commonDataBlock_->transmittingStationId_ +
                1000 * rawDataBlocks.at( i )->commonDataBlock_->transmittingStationNetworkId_;

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
            if( currentObservableType == observation_models::one_way_differenced_range ||
                    currentObservableType == observation_models::n_way_differenced_range )
            {
                processedDataBlocks[ currentObservableType ][ linkEnds ] = std::make_shared< ProcessedOdfFileDopplerData >( );
                processedDataBlocks[ currentObservableType ][ linkEnds ]->transmittingStation = appendedTransmittingStationId;
                processedDataBlocks[ currentObservableType ][ linkEnds ]->receivingStation =
                        std::to_string( rawDataBlocks.at( i )->commonDataBlock_->receivingStationId_ );
                processedDataBlocks[ currentObservableType ][ linkEnds ]->observableType = currentObservableType;
            }
            else
            {
                throw std::runtime_error( "Error when parsing ODF file contents, can currently only handle Doppler data" );
            }
        }

        addOdfDataBlockToProcessedData(
                currentObservableType, rawDataBlocks.at( i ),
                processedDataBlocks[ currentObservableType ][ linkEnds ] );
    }

    // Save output and return
    processedOdfFile->processedDataBlocks_ = processedDataBlocks;
    processedOdfFile->spacecraftName_ = spacecraftName;

    std::map< int, std::vector< std::shared_ptr< input_output::OdfRampBlock > > > rampDataBlocks = rawOdfData->rampBlocks_;
    std::map< int, std::shared_ptr< RampedReferenceFrequencyInterpolator > > rampInterpolators;

    for( auto it = rampDataBlocks.begin( ); it != rampDataBlocks.end( ); it++ )
    {
        rampInterpolators[ it->first ] =
                std::make_shared< RampedReferenceFrequencyInterpolator >( it->second );
    }
    processedOdfFile->rampInterpolators_ = rampInterpolators;

    return processedOdfFile;
}

} // namespace orbit_determination

} // namespace tudat
