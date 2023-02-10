#include "tudat/astro/orbit_determination/parseOdfFile.h"

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
        throw std::runtime_error( "Error wen getting observable type for ODF ID, ID: " +
                                  std::to_string( odfId ) + " not recognized." );
    }
    return observableType;
}


std::shared_ptr< RampedReferenceFrequencyInterpolator > mergeRampDataInterpolators(
        const std::vector< std::shared_ptr< RampedReferenceFrequencyInterpolator > >& interpolatorList )
{
    std::map< double, double > rampEndTimesPerStation;
    std::map< double, double > rampRatesPerStation;
    std::map< double, double > rampStartFrequenciesPerStation;

    for( unsigned int i = 0; i < interpolatorList.size( ); i++ )
    {
        for( unsigned int j = 0; j < interpolatorList.at( i )->startTimes.size( ); j++ )
        {
            rampEndTimesPerStation[ interpolatorList.at( i )->startTimes.at( j ) ] =
                    interpolatorList.at( i )->endTimes.at( j );
            rampRatesPerStation[ interpolatorList.at( i )->startTimes.at( j ) ] =
                    interpolatorList.at( i )->rampRates.at( j );
            rampStartFrequenciesPerStation[ interpolatorList.at( i )->startTimes.at( j ) ] =
                    interpolatorList.at( i )->startFrequency.at( j );
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
        std::shared_ptr< ProcessdOdfFileSingleLinkData > mergedOdfFileContents,
        std::shared_ptr< ProcessdOdfFileSingleLinkData > blockToAdd )
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
        std::shared_ptr< ProcessdOdfFileDopplerData > dopplerBlockToAdd
                = std::dynamic_pointer_cast< ProcessdOdfFileDopplerData >(
                    blockToAdd );
        std::shared_ptr< ProcessdOdfFileDopplerData > currentDopplerObservableMergedData =
                std::dynamic_pointer_cast< ProcessdOdfFileDopplerData >(
                    mergedOdfFileContents );


        currentDopplerObservableMergedData->referenceFrequency.insert(
                    currentDopplerObservableMergedData->referenceFrequency.end( ),
                    dopplerBlockToAdd->referenceFrequency.begin( ), dopplerBlockToAdd->referenceFrequency.end( ) );

        currentDopplerObservableMergedData->rampingFlag.insert(
                    currentDopplerObservableMergedData->rampingFlag.end( ),
                    dopplerBlockToAdd->rampingFlag.begin( ), dopplerBlockToAdd->rampingFlag.end( ) );

        currentDopplerObservableMergedData->compressionTimes.insert(
                    currentDopplerObservableMergedData->compressionTimes.end( ),
                    dopplerBlockToAdd->compressionTimes.begin( ), dopplerBlockToAdd->compressionTimes.end( ) );

        currentDopplerObservableMergedData->receiverChannels.insert(
                    currentDopplerObservableMergedData->receiverChannels.end( ),
                    dopplerBlockToAdd->receiverChannels.begin( ), dopplerBlockToAdd->receiverChannels.end( ) );
    }
}

std::shared_ptr< ProcessedOdfFileContents > mergeOdfFileContents(
        const std::vector< std::shared_ptr< ProcessedOdfFileContents > > odfFileContents )
{
    std::map< observation_models::ObservableType, std::map< std::pair< std::string, std::string >,
            std::shared_ptr< ProcessdOdfFileSingleLinkData > > > mergedOdfFileContents;
    std::map< int, std::vector< std::shared_ptr< RampedReferenceFrequencyInterpolator > > > rampInterpolatorList;

    // Iterate over all ODF files
    for( unsigned int i = 0; i < odfFileContents.size( ); i++ )
    {
        // Retrieve contents of current file.
        std::map< observation_models::ObservableType, std::map< std::pair< std::string, std::string >,
                std::shared_ptr< ProcessdOdfFileSingleLinkData > > >  dataBlocks =
                odfFileContents.at( i )->processedDataBlocks;

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
                odfFileContents.at( i )->rampInterpolators;

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

    processedOdfFile->processedDataBlocks = mergedOdfFileContents;
    processedOdfFile->rampInterpolators = mergedRampInterpolators;
    processedOdfFile->spacecraftName = "AAA";
    return processedOdfFile;
}

void addOdfDataBlockToParsedData(
        const observation_models::ObservableType currentObservableType,
        const std::shared_ptr< input_output::OdfDataBlock > rawDataBlock,
        const std::shared_ptr< ProcessdOdfFileSingleLinkData > processedDataBlock )
{
    // Add common properties to data object
    processedDataBlock->downlinkBand.push_back( rawDataBlock->commonDataBlock->downlinkBand );
    processedDataBlock->uplinkBand.push_back( rawDataBlock->commonDataBlock->uplinkBand );
    processedDataBlock->referenceBand.push_back( rawDataBlock->commonDataBlock->referenceBand );
    processedDataBlock->observableValues.push_back( rawDataBlock->commonDataBlock->getObservableValue( ) );
    processedDataBlock->observationTimes.push_back( rawDataBlock->commonDataBlock->getObservableTime( ) );
    processedDataBlock->receiverDownlinkDelay.push_back( rawDataBlock->commonDataBlock->receivingStationDownlinkDelay );

    // Add properties to data object for Doppler data
    if( currentObservableType == observation_models::one_way_differenced_range ||
            currentObservableType == observation_models::n_way_differenced_range )
    {
        std::shared_ptr< input_output::OdfDopplerDataBlock > odfDopplerDataBlock =
                std::dynamic_pointer_cast< input_output::OdfDopplerDataBlock >(
                    rawDataBlock->observableSpecificDataBlock );
        std::shared_ptr< ProcessdOdfFileDopplerData > odfParsedDopplerDataBlock =
                std::dynamic_pointer_cast< ProcessdOdfFileDopplerData >(
                    processedDataBlock );

        odfParsedDopplerDataBlock->compressionTimes.push_back( odfDopplerDataBlock->compressionTime );
        odfParsedDopplerDataBlock->receiverChannels.push_back( odfDopplerDataBlock->receiverChannel );
        odfParsedDopplerDataBlock->rampingFlag.push_back( odfDopplerDataBlock->receiverExciterFlag );
        odfParsedDopplerDataBlock->referenceFrequency.push_back( odfDopplerDataBlock->getReferenceFrequency( ) );
        odfParsedDopplerDataBlock->reservedData.push_back( odfDopplerDataBlock->reservedSegment );
        odfParsedDopplerDataBlock->uplinkDelays.push_back( odfDopplerDataBlock->transmittingStationUplinkDelay );
    }
}

std::shared_ptr< ProcessedOdfFileContents > parseOdfFileContents(
        const std::shared_ptr< input_output::OdfRawFileContents > rawOdfData )
{
    // Create output object
    std::shared_ptr< ProcessedOdfFileContents > processedOdfFile =
            std::make_shared< ProcessedOdfFileContents >( );
    std::string spacecraftName = std::to_string( rawOdfData->spacecraftId );

    // Retrieve data blocks from ODF file raw contents
    std::vector< std::shared_ptr< input_output::OdfDataBlock > > rawDataBlocks = rawOdfData->dataBlocks;

    // Create list of data, sorted by observable type and link ends; single object per combination of the two
    std::map< observation_models::ObservableType, std::map< std::pair< std::string, std::string >,
            std::shared_ptr< ProcessdOdfFileSingleLinkData > > > processedDataBlocks;

    bool createNewObject = false;
    int currentObservableId;
    observation_models::ObservableType currentObservableType;
    std::pair< std::string, std::string > stationIds;

    // Iterate over all block of ODF file.
    for( unsigned int i = 0; i < rawDataBlocks.size( ); i++ )
    {
        // Retrieve observable type and link end names
        currentObservableId = rawDataBlocks.at( i )->observableSpecificDataBlock->dataType_;
        currentObservableType = getObservableTypeForOdfId( currentObservableId );
        int appendedTransmittingStationId =
                rawDataBlocks.at( i )->commonDataBlock->transmittingStation + 100 *
                rawDataBlocks.at( i )->commonDataBlock->transmittingStationNetworkId;
        stationIds = std::make_pair( std::to_string( appendedTransmittingStationId ),
                                     std::to_string( rawDataBlocks.at( i )->commonDataBlock->receivingStation ) );

        // Check if data object already exists for current observable/link ends
        createNewObject = false;
        if( processedDataBlocks.count( currentObservableType ) == 0 )
        {
            createNewObject = true;
        }
        else if( processedDataBlocks.at( currentObservableType ).count( stationIds ) == 0 )
        {
            createNewObject = true;
        }

        // Create new data object, if required
        if( createNewObject )
        {
            if( currentObservableType == observation_models::one_way_differenced_range ||
                    currentObservableType == observation_models::n_way_differenced_range )
            {
                processedDataBlocks[ currentObservableType ][ stationIds ] = std::make_shared< ProcessdOdfFileDopplerData >( );
                processedDataBlocks[ currentObservableType ][ stationIds ]->transmittingStation = appendedTransmittingStationId;
                processedDataBlocks[ currentObservableType ][ stationIds ]->receivingStation =
                        std::to_string( rawDataBlocks.at( i )->commonDataBlock->receivingStation );
                processedDataBlocks[ currentObservableType ][ stationIds ]->observableType = currentObservableType;
            }
            else
            {
                throw std::runtime_error( "Error when parsing ODF file contents, can currently only handle Doppler data" );
            }
        }

        addOdfDataBlockToParsedData(
                    currentObservableType, rawDataBlocks.at( i ), processedDataBlocks[ currentObservableType ][ stationIds ] );
    }

    // Save output and return
    processedOdfFile->processedDataBlocks = processedDataBlocks;
    processedOdfFile->spacecraftName = spacecraftName;

    std::map< int, std::vector< std::shared_ptr< input_output::OdfRampBlock > > > rampDataBlocks = rawOdfData->odfRampBlocks;
    std::map< int, std::shared_ptr< RampedReferenceFrequencyInterpolator > > rampInterpolators;

    for( auto it = rampDataBlocks.begin( ); it != rampDataBlocks.end( ); it++ )
    {
        rampInterpolators[ it->first ] =
                std::make_shared< RampedReferenceFrequencyInterpolator >( it->second );
    }
    processedOdfFile->rampInterpolators = rampInterpolators;

    return processedOdfFile;
}

} // namespace orbit_determination

} // namespace tudat
