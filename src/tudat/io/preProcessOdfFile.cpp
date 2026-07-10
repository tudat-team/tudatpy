/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/preProcessOdfFile.h"
#include <string>
#include <utility>
#include <vector>

namespace tio = tudat::input_output;

namespace tudat
{

namespace input_output
{

std::string getObservableNameForOdfId( const input_output::OdfDataType observable_type_id )
{
    std::string observableName;

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
            observableName = "DsnNWayAveragedDoppler";
            break;
        case tio::OdfDataType::sra_planetary_operational_discrete_spectrum_range:
            observableName = "DsnNWayRange";
            break;
        default:
            throw std::runtime_error( "Error when getting observable name from ODF ID: ID " +
                                      std::to_string( static_cast< int >( observable_type_id ) ) + " not recognized" );
    }

    return observableName;
}

std::string getFrequencyBandNameForOdfId( const int bandId )
{
    std::string frequencyBandName;

    switch( bandId )
    {
        case 0:
            frequencyBandName = "Ku-band";
            break;
        case 1:
            frequencyBandName = "S-band";
            break;
        case 2:
            frequencyBandName = "X-band";
            break;
        case 3:
            frequencyBandName = "Ka-band";
            break;
        default:
            throw std::runtime_error( "Error when getting frequency band name for ODF band ID, ID: " + std::to_string( bandId ) +
                                      " not recognized." );
    }

    return frequencyBandName;
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

data::PlainLinkDefinition getLinkEndsFromOdfBlock( const std::shared_ptr< input_output::OdfDataBlock > dataBlock,
                                                   const std::string& spacecraftName )
{
    std::shared_ptr< input_output::OdfCommonDataBlock > commonDataBlock = dataBlock->getCommonDataBlock( );

    data::PlainLinkDefinition linkEnds;

    switch( commonDataBlock->dataType_ )
    {
        case input_output::OdfDataType::one_way_doppler:
            linkEnds = { std::make_pair( std::make_pair( spacecraftName, "" ), "transmitter" ),
                         std::make_pair( std::make_pair( "Earth", getStationNameFromStationId( 0, commonDataBlock->receivingStationId_ ) ),
                                         "receiver" ) };
            break;
        case input_output::OdfDataType::two_way_doppler:
        case input_output::OdfDataType::three_way_doppler:
        case input_output::OdfDataType::sra_planetary_operational_discrete_spectrum_range:
            linkEnds = { std::make_pair( std::make_pair( "Earth",
                                                         getStationNameFromStationId( commonDataBlock->transmittingStationNetworkId_,
                                                                                      commonDataBlock->transmittingStationId_ ) ),
                                         "transmitter" ),
                         std::make_pair( std::make_pair( spacecraftName, "" ), "reflector_1" ),
                         std::make_pair( std::make_pair( "Earth", getStationNameFromStationId( 0, commonDataBlock->receivingStationId_ ) ),
                                         "receiver" ) };

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

}  // namespace input_output

}  // namespace tudat
