/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/processOdfFile.h"

#include <algorithm>

namespace tio = tudat::input_output;

namespace tudat
{

namespace observation_models
{

observation_models::ObservableType getObservableTypeForOdfId( const input_output::OdfDataType observable_type_id )
{
    observation_models::ObservableType observableType;

    switch( observable_type_id )
    {
            // TODO: don't forget to remove
            // ProcessedOdfFileContentsPrivateFunctionTest class after implementing
            // processing of data type 11 (1-way Doppler)
            //    case 11:
            //        observableType =
            //        observation_models::dsn_one_way_averaged_doppler; break;
        case tio::OdfDataType::two_way_doppler:
        case tio::OdfDataType::three_way_doppler:
            observableType = observation_models::dsn_n_way_averaged_doppler;
            break;
        case tio::OdfDataType::sra_planetary_operational_discrete_spectrum_range:
            observableType = observation_models::dsn_n_way_range;
            break;
        default:
            throw std::runtime_error( "Error when getting observable type from ODF ID: ID " +
                                      std::to_string( static_cast< int >( observable_type_id ) ) + " not recognized" );
    }

    return observableType;
}

observation_models::FrequencyBands getFrequencyBandForOdfId( const int odfId )
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
            throw std::runtime_error( "Error when getting observable type for ODF ID, ID: " + std::to_string( odfId ) +
                                      " not recognized." );
    }

    return frequencyBand;
}

std::string getStationNameFromStationId( const int networkId, const int stationId )
{
    std::string stationName;

    if( networkId == 0 )
    {
        stationName = "DSS-" + std::to_string( stationId );
    }
    else if( networkId == 3 )
    {
        stationName = "UPL-" + std::to_string( stationId );
    }
    else
    {
        stationName = "Station-" + std::to_string( stationId );
    }

    return stationName;
}

observation_models::LinkEnds getLinkEndsFromOdfBlock( const std::shared_ptr< input_output::OdfDataBlock > dataBlock,
                                                      std::string spacecraftName )
{
    std::shared_ptr< input_output::OdfCommonDataBlock > commonDataBlock = dataBlock->getCommonDataBlock( );

    observation_models::LinkEnds linkEnds;

    switch( commonDataBlock->dataType_ )
    {
        case input_output::OdfDataType::one_way_doppler:
            linkEnds[ observation_models::LinkEndType::transmitter ] = observation_models::LinkEndId( spacecraftName );
            linkEnds[ observation_models::LinkEndType::receiver ] =
                    observation_models::LinkEndId( "Earth", getStationNameFromStationId( 0, commonDataBlock->receivingStationId_ ) );
            break;
        case input_output::OdfDataType::two_way_doppler:
        case input_output::OdfDataType::three_way_doppler:
        case input_output::OdfDataType::sra_planetary_operational_discrete_spectrum_range:
            linkEnds[ observation_models::LinkEndType::transmitter ] =
                    observation_models::LinkEndId( "Earth",
                                                   getStationNameFromStationId( commonDataBlock->transmittingStationNetworkId_,
                                                                                commonDataBlock->transmittingStationId_ ) );
            linkEnds[ observation_models::LinkEndType::reflector1 ] = observation_models::LinkEndId( spacecraftName );
            linkEnds[ observation_models::LinkEndType::receiver ] =
                    observation_models::LinkEndId( "Earth", getStationNameFromStationId( 0, commonDataBlock->receivingStationId_ ) );
            break;
        default:
            throw std::runtime_error( "Error when getting link ends from ODF data blocks: Data type " +
                                      std::to_string( static_cast< int >( commonDataBlock->dataType_ ) ) + " not recognized." );
            break;
    }

    return linkEnds;
}

bool compareRawOdfDataByStartDate( std::shared_ptr< input_output::OdfRawFileContents > rawOdfData1,
                                   std::shared_ptr< input_output::OdfRawFileContents > rawOdfData2 )
{
    if( rawOdfData1->getDataBlocks( ).at( 0 )->getCommonDataBlock( )->getObservableTime( ) <
        rawOdfData2->getDataBlocks( ).at( 0 )->getCommonDataBlock( )->getObservableTime( ) )
    {
        return true;
    }
    else
    {
        return false;
    }
}

}  // namespace observation_models

}  // namespace tudat
