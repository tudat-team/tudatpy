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

#include "tudat/simulation/estimation_setup/processOdfFile.h"

namespace tudat
{

namespace observation_models
{

observation_models::ObservableType getObservableTypeForOdfId( const int odfId )
{
    observation_models::ObservableType observableType;

    switch( odfId )
    {
    // TODO: don't forget to remove ProcessedOdfFileContentsPrivateFunctionTest class after implementing processing of data type 11 (1-way Doppler)
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
                if ( linkEndTypeIt->second.stationName_ != "" && linkEndTypeIt->second.bodyName_ != spacecraftName_ )
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

std::pair< double, double > ProcessedOdfFileContents::getStartAndEndTime( )
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

observation_models::LinkEnds getLinkEndsFromOdfBlock (
        const std::shared_ptr< input_output::OdfDataBlock > dataBlock,
        std::string spacecraftName )
{
    int currentObservableId = dataBlock->getObservableSpecificDataBlock( )->dataType_;

    observation_models::LinkEnds linkEnds;

    if ( currentObservableId == 11 )
    {
        linkEnds[ observation_models::transmitter ] = observation_models::LinkEndId ( spacecraftName, "Antenna" );
        linkEnds[ observation_models::receiver ] = observation_models::LinkEndId ( "Earth", getStationNameFromStationId(
                0, dataBlock->getCommonDataBlock( )->receivingStationId_ ) );
    }
    else if ( currentObservableId == 12 )
    {
        linkEnds[ observation_models::transmitter ] = observation_models::LinkEndId (
                "Earth", getStationNameFromStationId( dataBlock->getCommonDataBlock( )->transmittingStationNetworkId_,
                                                      dataBlock->getCommonDataBlock( )->transmittingStationId_ ) );
        linkEnds[ observation_models::reflector1 ] = observation_models::LinkEndId ( spacecraftName, "Antenna" );
        linkEnds[ observation_models::receiver ] = observation_models::LinkEndId (
                "Earth", getStationNameFromStationId( 0, dataBlock->getCommonDataBlock( )->receivingStationId_ ) );
    }
    else if ( currentObservableId == 13 )
    {
        linkEnds[ observation_models::transmitter ] = observation_models::LinkEndId (
                "Earth", getStationNameFromStationId( dataBlock->getCommonDataBlock( )->transmittingStationNetworkId_,
                                                      dataBlock->getCommonDataBlock( )->transmittingStationId_ ) );
        linkEnds[ observation_models::reflector1 ] = observation_models::LinkEndId ( spacecraftName, "Antenna" );
        linkEnds[ observation_models::receiver ] = observation_models::LinkEndId (
                "Earth", getStationNameFromStationId( 0, dataBlock->getCommonDataBlock( )->receivingStationId_ ) );
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
        const std::shared_ptr< input_output::OdfDataBlock > rawDataBlock,
        const std::shared_ptr< ProcessedOdfFileSingleLinkData > singleLinkProcessedData,
        const std::string rawDataFileName )
{
    // Add properties to data block if data is valid
    if ( rawDataBlock->getCommonDataBlock( )->validity_ == 0 )
    {

        // Add common properties to data object
        singleLinkProcessedData->downlinkBandIds_.push_back( rawDataBlock->getCommonDataBlock( )->downlinkBandId_ );
        singleLinkProcessedData->uplinkBandIds_.push_back( rawDataBlock->getCommonDataBlock( )->uplinkBandId_ );
        singleLinkProcessedData->referenceBandIds_.push_back( rawDataBlock->getCommonDataBlock( )->referenceBandId_ );
        singleLinkProcessedData->unprocessedObservationTimes_.push_back( rawDataBlock->getCommonDataBlock( )->getObservableTime( ) );
        singleLinkProcessedData->receiverDownlinkDelays_.push_back( rawDataBlock->getCommonDataBlock( )->getReceivingStationDownlinkDelay( ) );
        singleLinkProcessedData->originFiles_.push_back( rawDataFileName );

        // Add properties to data object for Doppler data
        if ( singleLinkProcessedData->getObservableType( ) == observation_models::dsn_n_way_averaged_doppler )
        {
            std::shared_ptr< input_output::OdfDopplerDataBlock > odfDopplerDataBlock =
                    std::dynamic_pointer_cast< input_output::OdfDopplerDataBlock >(
                            rawDataBlock->getObservableSpecificDataBlock( ) );
            std::shared_ptr< ProcessedOdfFileDopplerData > odfParsedDopplerDataBlock =
                    std::dynamic_pointer_cast< ProcessedOdfFileDopplerData >(
                            singleLinkProcessedData );

            singleLinkProcessedData->observableValues_.push_back(
                    ( Eigen::Matrix< double, 1, 1 >( ) << rawDataBlock->getCommonDataBlock( )->getObservableValue( )
                    ).finished( ) );

            odfParsedDopplerDataBlock->countInterval_.push_back( odfDopplerDataBlock->getCompressionTime( ) );
            odfParsedDopplerDataBlock->receiverChannels_.push_back( odfDopplerDataBlock->getReceiverChannel( ) );
            odfParsedDopplerDataBlock->receiverRampingFlags_.push_back( odfDopplerDataBlock->getReceiverExciterFlag( ) );
            odfParsedDopplerDataBlock->referenceFrequencies_.push_back( odfDopplerDataBlock->getReferenceFrequency( ) );
            odfParsedDopplerDataBlock->transmitterUplinkDelays_.push_back( odfDopplerDataBlock->getTransmittingStationUplinkDelay( ) );
        }
    }
}

void ProcessedOdfFileContents::updateProcessedObservationTimes( )
{
    // Loop over saved data and convert time to TDB wrt J2000
    for ( auto observableTypeIterator = processedDataBlocks_.begin( );
            observableTypeIterator != processedDataBlocks_.end( ); ++observableTypeIterator )
    {
        for ( auto linkEndsIterator = observableTypeIterator->second.begin( );
              linkEndsIterator != observableTypeIterator->second.end( ); ++linkEndsIterator )
        {
            linkEndsIterator->second->processedObservationTimes_ = computeObservationTimesTdbFromJ2000(
                    linkEndsIterator->second->receivingStation_,
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

        if ( approximateEarthFixedGroundStationPositions_.count( groundStation ) == 0 )
        {
            throw std::runtime_error( "Error when processing ODF file, converting time from UTC to TDB: the position of "
                                      "the ground station " + groundStation + " was not specified." );
        }

        std::vector< Eigen::Vector3d > earthFixedPositions;
        for ( unsigned int i = 0; i < observationTimesUtcFromJ2000.size( ); ++i )
        {
            // TDB time wrt Earth center used to retrieve the ground station's position
            earthFixedPositions.push_back( approximateEarthFixedGroundStationPositions_.at( groundStation ) );
        }

        std::vector< double > observationTimesTdbFromJ2000 = timeScaleConverter.getCurrentTimes(
                basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, observationTimesUtcFromJ2000,
                earthFixedPositions );

        return observationTimesTdbFromJ2000;
    }

bool compareRawOdfDataByStartDate( std::shared_ptr< input_output::OdfRawFileContents > rawOdfData1,
                                   std::shared_ptr< input_output::OdfRawFileContents > rawOdfData2 )
{
    if ( rawOdfData1->getDataBlocks( ).at( 0 )->getCommonDataBlock( )->getObservableTime( ) <
        rawOdfData2->getDataBlocks( ).at( 0 )->getCommonDataBlock( )->getObservableTime( ) )
    {
        return true;
    }
    else
    {
        return false;
    }
}

void ProcessedOdfFileContents::sortAndValidateOdfDataVector(
        std::vector< std::shared_ptr< input_output::OdfRawFileContents > >& rawOdfDataVector )
{
    unsigned int spacecraftId = rawOdfDataVector.front( )->spacecraftId_;

    for ( unsigned int i = 0; i < rawOdfDataVector.size( ); ++i )
    {
        // Check if spacecraft ID is valid
        if ( rawOdfDataVector.at( i )->spacecraftId_ != spacecraftId )
        {
            throw std::runtime_error( "Error when creating processed ODF object from raw data: multiple spacecraft IDs"
                                      "found (" + std::to_string( spacecraftId ) + " and " +
                                      std::to_string( rawOdfDataVector.at( i )->spacecraftId_ ) + ")." );
        }
    }

    std::stable_sort( rawOdfDataVector.begin( ), rawOdfDataVector.end( ), &compareRawOdfDataByStartDate );
}

bool ProcessedOdfFileContents::isObservationValid(
        std::shared_ptr< input_output::OdfDataBlock > rawDataBlock,
        observation_models::LinkEnds linkEnds,
        observation_models::ObservableType currentObservableType )
{
    int currentObservableId = rawDataBlock->getObservableSpecificDataBlock( )->dataType_;

    std::string transmittingStation, receivingStation;

    if ( requiresTransmittingStation( currentObservableType ) )
    {
        transmittingStation = linkEnds.at( transmitter ).stationName_;

        // Check if transmitting station is in ramp tables
        if ( rampInterpolators_.count( transmittingStation ) == 0 )
        {
            if ( std::count( ignoredGroundStations_.begin( ), ignoredGroundStations_.end( ), transmittingStation ) == 0 )
            {
                ignoredGroundStations_.push_back( transmittingStation );
                if ( verbose_ )
                {
                    std::cerr << "Warning: ground station " << transmittingStation << " not available in ramp tables," <<
                        " ignoring corresponding data." << std::endl;
                }
            }
            ignoredOdfRawDataBlocks_.push_back( rawDataBlock );
            return false;
        }

        // Check if observation time is covered by ramp tables
        if ( rawDataBlock->getCommonDataBlock( )->getObservableTime( ) <
            unprocessedRampStartTimesPerStation_[ transmittingStation ].front( ) )
        {
            if ( verbose_ )
            {
                std::cerr << "Warning: observation of ODF type " << currentObservableId << " not covered by ramp table of station " <<
                    transmittingStation << ", ignoring it." << std::endl;
            }
            ignoredOdfRawDataBlocks_.push_back( rawDataBlock );
            return false;
        }
    }
    if ( requiresFirstReceivingStation( currentObservableType ) )
    {
        receivingStation = linkEnds.at( receiver ).stationName_;

        // Check if receiving station is in ramp tables
        if ( rampInterpolators_.count( receivingStation ) == 0 )
        {
            if ( std::count( ignoredGroundStations_.begin( ), ignoredGroundStations_.end( ), receivingStation ) == 0 )
            {
                ignoredGroundStations_.push_back( receivingStation );
                if ( verbose_ )
                {
                    std::cerr << "Warning: observation of ODF type " << currentObservableId << " not covered by ramp table of station " <<
                        receivingStation << ", ignoring it." << std::endl;
                }
            }
            ignoredOdfRawDataBlocks_.push_back( rawDataBlock );
            return false;
        }

        // Check if observation time is covered by ramp tables
        if ( rawDataBlock->getCommonDataBlock( )->getObservableTime( ) < unprocessedRampStartTimesPerStation_[ receivingStation ].front( ) ||
            rawDataBlock->getCommonDataBlock( )->getObservableTime( ) > unprocessedRampStartTimesPerStation_[ receivingStation ].back( ) )
        {
            if ( verbose_ )
            {
                std::cerr << "Warning: observation of ODF type " << currentObservableId << " not covered by ramp tables," <<
                    " ignoring it." << std::endl;
            }
            ignoredOdfRawDataBlocks_.push_back( rawDataBlock );
            return false;
        }
    }

    return true;
}

void ProcessedOdfFileContents::extractRawOdfOrbitData(
        std::shared_ptr< input_output::OdfRawFileContents > rawOdfData )
{
    // Retrieve data blocks from ODF file raw contents
    std::vector< std::shared_ptr< input_output::OdfDataBlock > > rawDataBlocks = rawOdfData->getDataBlocks( );

    // Iterate over all block of ODF file.
    for( unsigned int i = 0; i < rawDataBlocks.size( ); i++ )
    {
        // Retrieve observable type and link ends
        int currentObservableId = rawDataBlocks.at( i )->getObservableSpecificDataBlock( )->dataType_;

        // Get current observable type and throw warning if not implemented
        observation_models::ObservableType currentObservableType;
        try
        {
            currentObservableType = getObservableTypeForOdfId( currentObservableId );
        }
        catch( const std::runtime_error& )
        {
            if ( std::find( ignoredRawOdfObservableTypes_.begin( ), ignoredRawOdfObservableTypes_.end( ),
                            currentObservableId ) == ignoredRawOdfObservableTypes_.end( ) )
            {
                ignoredRawOdfObservableTypes_.push_back( currentObservableId );
                if ( verbose_ )
                {
                    std::cerr << "Warning: processing of ODF data type " << currentObservableId <<
                        " is not implemented, ignoring the corresponding data." << std::endl;
                }
            }
            ignoredOdfRawDataBlocks_.push_back( rawDataBlocks.at( i ) );
            continue;
        }

        observation_models::LinkEnds linkEnds = getLinkEndsFromOdfBlock(
                rawDataBlocks.at( i ), spacecraftName_ );

        // Check if observation is valid and should be processed
        if ( isObservationValid( rawDataBlocks.at( i ), linkEnds, currentObservableType ) )
        {
            // Check if data object already exists for current observable/link ends
            bool createNewObject = false;
            if ( processedDataBlocks_.count( currentObservableType ) == 0 )
            {
                createNewObject = true;
            }
            else if ( processedDataBlocks_.at( currentObservableType ).count( linkEnds ) == 0 )
            {
                createNewObject = true;
            }

            // Create new data object, if required
            if ( createNewObject )
            {
                processedDataBlocks_[ currentObservableType ][ linkEnds ] = std::make_shared< ProcessedOdfFileDopplerData >(
                        currentObservableType, linkEnds.at( receiver ).stationName_, linkEnds.at( transmitter ).stationName_ );
            }

            addOdfRawDataBlockToProcessedData(
                    rawDataBlocks.at( i ), processedDataBlocks_[ currentObservableType ][ linkEnds ],
                    rawOdfData->fileName_ );
        }
    }

}

void ProcessedOdfFileContents::extractMultipleRawOdfRampData(
        std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector )
{

    std::map< std::string, std::vector< double > > rampRatesPerStation, startFrequenciesPerStation;

    for( unsigned int i = 0; i < rawOdfDataVector.size( ); ++i )
    {
        std::map< int, std::vector< std::shared_ptr< input_output::OdfRampBlock > > >
                rampBlocksPerStation = rawOdfDataVector.at( i )->getRampBlocks( );
        for( auto it = rampBlocksPerStation.begin( ); it != rampBlocksPerStation.end( ); it++ )
        {
            std::string stationName = getStationNameFromStationId( 0, it->first );

            std::vector< std::shared_ptr< input_output::OdfRampBlock > > rampBlocks = it->second;

            for( unsigned int j = 0; j < it->second.size( ); j++ )
            {
                // Check if zero time ramp
                if ( rampBlocks.at( j )->getRampStartTime( ) == rampBlocks.at( j )->getRampEndTime( ) )
                {
                    continue;
                }

                // Check if adding ramp block vector to previously existing vector: add connection point
                if ( j == 0 && !unprocessedRampStartTimesPerStation_[ stationName ].empty( ) )
                {
                    unprocessedRampStartTimesPerStation_[ stationName ].push_back( unprocessedRampEndTimesPerStation_[ stationName ].back( ) );
                    unprocessedRampEndTimesPerStation_[ stationName ].push_back( rampBlocks.at( j )->getRampStartTime( ) );
                    rampRatesPerStation[ stationName ].push_back( TUDAT_NAN );
                    startFrequenciesPerStation[ stationName ].push_back( TUDAT_NAN );
                }

                unprocessedRampStartTimesPerStation_[ stationName ].push_back( rampBlocks.at( j )->getRampStartTime( ) );
                unprocessedRampEndTimesPerStation_[ stationName ].push_back( rampBlocks.at( j )->getRampEndTime( ) );
                rampRatesPerStation[ stationName ].push_back( rampBlocks.at( j )->getRampRate( ) );
                startFrequenciesPerStation[ stationName ].push_back( rampBlocks.at( j )->getRampStartFrequency( ) );
            }
        }
    }

    for( auto it = unprocessedRampStartTimesPerStation_.begin( ); it != unprocessedRampStartTimesPerStation_.end( ); ++it )
    {
        std::string stationName = it->first;
        rampInterpolators_[ stationName ] = std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >(
                computeObservationTimesTdbFromJ2000( stationName, unprocessedRampStartTimesPerStation_[ stationName ] ),
                computeObservationTimesTdbFromJ2000( stationName, unprocessedRampEndTimesPerStation_[ stationName ] ),
                rampRatesPerStation[ stationName ], startFrequenciesPerStation[ stationName ] );
    }

}

void setOdfInformationInBodies(
        const std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents,
        simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyWithGroundStations,
        const std::function< double ( FrequencyBands uplinkBand, FrequencyBands downlinkBand ) > getTurnaroundRatio )
{

    // Set transmitting frequency objects in ground stations
    setTransmittingFrequenciesInGroundStations( processedOdfFileContents, bodies.getBody( bodyWithGroundStations ) );

    // Set turnaround ratios in spacecraft body
    std::shared_ptr< system_models::VehicleSystems > vehicleSystems =  bodies.getBody( processedOdfFileContents->getSpacecraftName( ) )->getVehicleSystems( );
    vehicleSystems->setTransponderTurnaroundRatio( getTurnaroundRatio );
}

} // namespace observation_models

} // namespace tudat
