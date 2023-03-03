/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <vector>

#include "tudat/basics/testMacros.h"
#include "tudat/io/readOdfFile.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"
#include "tudat/simulation/estimation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;

std::string getAbsolutePath ( std::string relativePath )
{
    std::string currFilePath = __FILE__;
    std::string currFileName = __FILE_NAME__;

    std::string currDir = currFilePath;
    currDir.resize( currFilePath.size( ) - currFileName.size( ) - 1 );

    return currDir + relativePath;
}

BOOST_AUTO_TEST_SUITE( test_odf_file_reader )

//! Checks parsed binary file (odf07155.odf) against values in txt file (odf07155.txt)
//! Files available at https://pds-geosciences.wustl.edu/dataserv/radio_science.htm (see radio science documentation bundle)
BOOST_AUTO_TEST_CASE( testSingleOdfFileReader )
{
    spice_interface::loadStandardSpiceKernels( );

    // Create Earth and it ground stations
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );

    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate );

    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    std::string file = getAbsolutePath( "/data_files/odf07155.odf" );
    std::shared_ptr< input_output::OdfRawFileContents > rawOdfContents = input_output::readOdfFile(
            file );

    // Saved file name
    BOOST_CHECK_EQUAL ( rawOdfContents->fileName_, file );

    // File label group
    BOOST_CHECK_EQUAL ( rawOdfContents->systemId_, "TDDS    " );
    BOOST_CHECK_EQUAL ( rawOdfContents->programId_, "AMMOS   " );
    BOOST_CHECK_EQUAL ( rawOdfContents->spacecraftId_, 236 );
    BOOST_CHECK_EQUAL ( rawOdfContents->fileCreationDate_, 1071106 );
    BOOST_CHECK_EQUAL ( rawOdfContents->fileCreationTime_, 230913 );
    BOOST_CHECK_EQUAL ( rawOdfContents->fileReferenceDate_, 19500101 );
    BOOST_CHECK_EQUAL ( rawOdfContents->fileReferenceTime_, 000000 );

    // Identifier group
    BOOST_CHECK_EQUAL ( rawOdfContents->identifierGroupStringA_, "TIMETAG " );
    BOOST_CHECK_EQUAL ( rawOdfContents->identifierGroupStringB_, "OBSRVBL " );
    BOOST_CHECK_EQUAL ( rawOdfContents->identifierGroupStringC_, "FREQ,ANCILLARY-DATA " );

    // Block 0: 1-way doppler
    std::shared_ptr< input_output::OdfCommonDataBlock > commonDataBlock = rawOdfContents->dataBlocks_.at( 0 )->commonDataBlock_;
    BOOST_CHECK_EQUAL ( commonDataBlock->getObservableTime( ), 1812103240.000 );
    BOOST_CHECK_EQUAL ( commonDataBlock->getReceivingStationDownlinkDelay( ), 0.000000000 );
    BOOST_CHECK_EQUAL ( commonDataBlock->getObservableValue( ), -382738.663803100 );
    BOOST_CHECK_EQUAL ( commonDataBlock->formatId_, 2 );
    BOOST_CHECK_EQUAL ( commonDataBlock->receivingStationId_, 63 );
    BOOST_CHECK_EQUAL ( commonDataBlock->transmittingStationId_, 0 );
    BOOST_CHECK_EQUAL ( commonDataBlock->transmittingStationNetworkId_, 0 );
    BOOST_CHECK_EQUAL ( commonDataBlock->downlinkBandId_, 2 );
    BOOST_CHECK_EQUAL ( commonDataBlock->uplinkBandId_, 0 );
    BOOST_CHECK_EQUAL ( commonDataBlock->referenceBandId_, 2 );
    BOOST_CHECK_EQUAL ( commonDataBlock->validity_, 0 );

    std::shared_ptr< input_output::OdfDopplerDataBlock > dopplerDataBlock =
            std::dynamic_pointer_cast< input_output::OdfDopplerDataBlock >(
                    rawOdfContents->dataBlocks_.at( 0 )->observableSpecificDataBlock_ );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->dataType_, 11 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->receiverChannel_, 1 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->spacecraftId_, 236 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->receiverExciterFlag_, 1 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getReferenceFrequency( ), 2299812417.000 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->reservedSegment_, 0 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getCompressionTime( ), 60 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getTransmittingStationUplinkDelay( ), 0 );

    // Block 19: 2-way doppler
    commonDataBlock = rawOdfContents->dataBlocks_.at( 19 )->commonDataBlock_;
    BOOST_CHECK_EQUAL ( commonDataBlock->getObservableTime( ), 1812104598.000 );
    BOOST_CHECK_EQUAL ( commonDataBlock->getReceivingStationDownlinkDelay( ), 0.000000000 );
    BOOST_CHECK_EQUAL ( commonDataBlock->getObservableValue( ), -157.702220916 );
    BOOST_CHECK_EQUAL ( commonDataBlock->formatId_, 2 );
    BOOST_CHECK_EQUAL ( commonDataBlock->receivingStationId_, 63 );
    BOOST_CHECK_EQUAL ( commonDataBlock->transmittingStationId_, 63 );
    BOOST_CHECK_EQUAL ( commonDataBlock->transmittingStationNetworkId_, 0 );
    BOOST_CHECK_EQUAL ( commonDataBlock->downlinkBandId_, 2 );
    BOOST_CHECK_EQUAL ( commonDataBlock->uplinkBandId_, 2 );
    BOOST_CHECK_EQUAL ( commonDataBlock->referenceBandId_, 2 );
    BOOST_CHECK_EQUAL ( commonDataBlock->validity_, 0 );

    dopplerDataBlock = std::dynamic_pointer_cast< input_output::OdfDopplerDataBlock >(
            rawOdfContents->dataBlocks_.at( 19 )->observableSpecificDataBlock_ );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->dataType_, 12 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->receiverChannel_, 1 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->spacecraftId_, 236 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->receiverExciterFlag_, 1 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getReferenceFrequency( ), 7177648275.000 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->reservedSegment_, 0 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getCompressionTime( ), 60 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getTransmittingStationUplinkDelay( ), 0 );

    // Block 23: sequential range
    commonDataBlock = rawOdfContents->dataBlocks_.at( 23 )->commonDataBlock_;
    BOOST_CHECK_EQUAL ( commonDataBlock->getObservableTime( ), 1812104814.000 );
    BOOST_CHECK_EQUAL ( commonDataBlock->getReceivingStationDownlinkDelay( ), 0.000000000 );
    BOOST_CHECK_EQUAL ( commonDataBlock->getObservableValue( ), 587993.568119415 );
    BOOST_CHECK_EQUAL ( commonDataBlock->formatId_, 2 );
    BOOST_CHECK_EQUAL ( commonDataBlock->receivingStationId_, 63 );
    BOOST_CHECK_EQUAL ( commonDataBlock->transmittingStationId_, 63 );
    BOOST_CHECK_EQUAL ( commonDataBlock->transmittingStationNetworkId_, 0 );
    BOOST_CHECK_EQUAL ( commonDataBlock->downlinkBandId_, 2 );
    BOOST_CHECK_EQUAL ( commonDataBlock->uplinkBandId_, 2 );
    BOOST_CHECK_EQUAL ( commonDataBlock->referenceBandId_, 2 );
    BOOST_CHECK_EQUAL ( commonDataBlock->validity_, 0 );

    std::shared_ptr< input_output::OdfSequentialRangeDataBlock > rangeDataBlock =
            std::dynamic_pointer_cast< input_output::OdfSequentialRangeDataBlock >(
                    rawOdfContents->dataBlocks_.at( 23 )->observableSpecificDataBlock_ );
    BOOST_CHECK_EQUAL ( rangeDataBlock->dataType_, 37 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->lowestRangingComponent_, 14 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->spacecraftId_, 236 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->reservedBlock_, 1 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->getReferenceFrequency( ), 7177004669.452 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->coderInPhaseTimeOffset_, 774 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->compositeTwo_, 400000 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->getTransmittingStationUplinkDelay( ), 0 );

    // Check number of data blocks
    BOOST_CHECK_EQUAL ( rawOdfContents->dataBlocks_.size( ), 2233 - 4 - 1 );

    // Check number of ramp groups
    BOOST_CHECK_EQUAL ( rawOdfContents->rampBlocks_.size( ), 3 );

    // Ramp block for DSS-63
    std::shared_ptr< input_output::OdfRampBlock > rampBlock =
            std::dynamic_pointer_cast< input_output::OdfRampBlock >(
                    rawOdfContents->rampBlocks_.at( 63 ).at( 5 ) );
    BOOST_CHECK_EQUAL ( rampBlock->getRampStartTime( ), 1812101116.000000000 );
    BOOST_CHECK_EQUAL ( rampBlock->getRampRate( ), 0.095680000 );
    BOOST_CHECK_EQUAL ( rampBlock->getRampStartFrequency( ), 7177004073.170830727 );
    BOOST_CHECK_EQUAL ( rampBlock->getRampEndTime( ), 1812101962.000000000 );

    // Check number of ramp blocks
    BOOST_CHECK_EQUAL ( rawOdfContents->rampBlocks_.at( 63 ).size( ), 2331 - 2233 - 1 );

}

    //   boost::shared_ptr< orbit_determination::ProcessedOdfFileContents > odfContents =
    //           orbit_determination::parseOdfFileContents(
    //               input_output::readOdfFile( "/home/dominic/Downloads/mromagr2017_117_0745xmmmv1.odf" ) );

    //   std::map< observation_models::ObservableType,
    //           std::vector< boost::shared_ptr< orbit_determination::ProcessdOdfFileSingleLinkData > > > dataBlocks =
    //           odfContents->dataBlocks;

    //   for( auto it = dataBlocks.begin( ); it != dataBlocks.end( ); it++ )
    //   {
    //       int counter = 0;
    //       for( unsigned int i = 0; i < it->second.size( ); i++ )
    //       {
    //           boost::shared_ptr< orbit_determination::ProcessdOdfFileDopplerData > currentDopplerData =
    //                   boost::dynamic_pointer_cast< orbit_determination::ProcessdOdfFileDopplerData >(
    //                       it->second.at( i ) );
    //           std::string fileSuffix = std::to_string( it->first ) + "_" + std::to_string( counter );

    //           input_output::writeDataMapToTextFile(
    //                       currentDopplerData->getCompressionTimes( ),
    //                           "odfTestCompressionTimes_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
    //           input_output::writeDataMapToTextFile(
    //                       currentDopplerData->getObservationData( ),
    //                           "odfTestObservations_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
    //           input_output::writeDataMapToTextFile(
    //                       currentDopplerData->getReferenceFrequencies( ),
    //                           "odfTestReferenceFrequencies_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
    //           input_output::writeDataMapToTextFile(
    //                       currentDopplerData->getRampFlags( ),
    //                           "odfTestRampFlags_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
    //           counter++;
    //       }
    //   }

//    std::vector< boost::filesystem::path > files = input_output::listAllFilesInDirectory(
//                "/Users/pipas/Documents/dsn_trk-2-18/", false ); // "/home/dominic/Software/MercuryData/odf/"

//    for( unsigned int i = 0; i < files.size( ); i++ )
//    {
//        if( i % 100 == 0 )
//        {
//            //std::cout<<i<<std::endl;
//        }
//        std::string fileString = files.at( i ).string( );
//        int stringSize = fileString.size( );
//
//        if( fileString.substr( stringSize - 3, stringSize -1 ) == "dat" )
//        {
////            input_output::readOdfFile( "/Users/pipas/Documents/dsn_trk-2-18/" + fileString );
//
//            odfContentsList.push_back( orbit_determination::ProcessedOdfFileContents(
//                    input_output::readOdfFile( "/Users/pipas/Documents/dsn_trk-2-18/" + fileString ) ) ); // "/home/dominic/Software/MercuryData/odf/"
//        }
//    }

//    observation_models::ProcessedOdfFileContents odfContents = observation_models::ProcessedOdfFileContents(
//            input_output::readOdfFile( "/Users/pipas/Documents/dsn_trk-2-18/odf07155.dat" ),
//            bodies.getBody( "Earth" ) );

//
//    std::shared_ptr< orbit_determination::ProcessedOdfFileContents > mergedData =
//            orbit_determination::mergeOdfFileContents(
//                odfContentsList );
//
//    std::map< observation_models::ObservableType, std::map< std::pair< std::string, std::string >,
//            std::shared_ptr< orbit_determination::ProcessdOdfFileSingleLinkData > > > dataBlocksMerged =
//            mergedData->processedDataBlocks;
//
//    double filterStartTime = 58.0 * 365.25 * 86400.0;
//    double filterEndTime = 62.0 * 365.25 * 86400.0;
//
//    double filterFrequency = 200.0E3;
//
//    for( auto it = dataBlocksMerged.begin( ); it != dataBlocksMerged.end( ); it++ )
//    {
//        int counter = 0;
//
//        std::vector< double > observationTimes, observables, referenceFrequencies;
//        std::vector< std::string > receivingStation, transmittingStation;
//        std::vector< double > receivingRampFrequency, transmittingRampFrequency;
//        std::vector< std::string > originFiles;
//        std::vector< bool > rampFlags;
//
//
//        for( auto linkIt = it->second.begin( ); linkIt != it->second.end( ); linkIt++ )
//        {
//            std::shared_ptr< orbit_determination::ProcessdOdfFileDopplerData > currentDopplerData =
//                    std::dynamic_pointer_cast< orbit_determination::ProcessdOdfFileDopplerData >( linkIt->second );
//
//            std::vector< double > currentObservationTimes = currentDopplerData->observationTimes,
//                    currentObservables = currentDopplerData->observableValues,
//                    currentReferenceFrequencies = currentDopplerData->referenceFrequency;
//            std::vector< std::string > currentOriginFiles = currentDopplerData->originFile;
//            std::vector< bool > currentRampFlags = currentDopplerData->rampingFlag;
//            std::cout<<it->first<<" "<<linkIt->first.first<<" "<<linkIt->first.second<<" "<<currentObservationTimes.size( )<<std::endl;
//
//            std::shared_ptr< orbit_determination::RampedReferenceFrequencyInterpolator > transmitterRampInterpolator;
//            if( mergedData->rampInterpolators.count(
//                        boost::lexical_cast< int >( linkIt->first.first ) ) != 0 )
//            {
//                transmitterRampInterpolator = mergedData->rampInterpolators.at(
//                            boost::lexical_cast< int >( linkIt->first.first ) );
//            }
//            else
//            {
//                transmitterRampInterpolator = NULL;
//            }
//
//            std::shared_ptr< orbit_determination::RampedReferenceFrequencyInterpolator > receiverRampInterpolator;
//            if( mergedData->rampInterpolators.count(
//                        boost::lexical_cast< int >( linkIt->first.second ) ) != 0 )
//            {
//                receiverRampInterpolator = mergedData->rampInterpolators.at(
//                            boost::lexical_cast< int >( linkIt->first.second ) );
//            }
//            else
//            {
//                receiverRampInterpolator = NULL;
//            }
//
//            for( unsigned int j = 0; j < currentObservationTimes.size( ); j++ )
//            {
//                if( currentObservationTimes.at( j ) > filterStartTime &&
//                        currentObservationTimes.at( j ) < filterEndTime &&
//                        std::fabs( currentObservables.at( j ) ) < filterFrequency
//                        && std::stoi( linkIt->first.second ) > 0 && std::stoi( linkIt->first.second ) < 100
//                         && std::stoi( linkIt->first.first ) > 0 && std::stoi( linkIt->first.first ) < 100 )
//                {
//
//                    observationTimes.push_back( currentObservationTimes.at( j ) );
//                    observables.push_back( currentObservables.at( j ) );
//                    //originFiles.push_back( currentOriginFiles.at( j ) );
//                    rampFlags.push_back( currentRampFlags.at( j ) );
//                    referenceFrequencies.push_back( currentReferenceFrequencies.at( j ) );
//
//                    receivingStation.push_back( linkIt->first.second );
//                    transmittingStation.push_back( linkIt->first.first );
//
//                    bool testBoolean;
//
//                    if( receiverRampInterpolator != NULL )
//                    {
//                       receivingRampFrequency.push_back( receiverRampInterpolator->getCurrentReferenceFrequency(
//                                                             currentObservationTimes.at( j ), testBoolean ) );
//                    }
//                    else
//                    {
//                        receivingRampFrequency.push_back( TUDAT_NAN );
//                    }
//
//                    if( transmitterRampInterpolator != NULL )
//                    {
//                       transmittingRampFrequency.push_back( transmitterRampInterpolator->getCurrentReferenceFrequency(
//                                                                currentObservationTimes.at( j ), testBoolean ) );
//                    }
//                    else
//                    {
//                        transmittingRampFrequency.push_back( TUDAT_NAN );
//                    }
//
//
//
//                    //           std::string fileSuffix = std::to_string( it->first ) + "_" + std::to_string( counter );
//
//                    //           input_output::writeDataMapToTextFile(
//                    //                       currentDopplerData->getCompressionTimes( ),
//                    //                           "odfTestCompressionTimes_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
//                    //           input_output::writeDataMapToTextFile(
//                    //                       currentDopplerData->getObservationData( ),
//                    //                           "odfTestObservations_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
//                    //           input_output::writeDataMapToTextFile(
//                    //                       currentDopplerData->getReferenceFrequencies( ),
//                    //                           "odfTestReferenceFrequencies_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
//                    //           input_output::writeDataMapToTextFile(
//                    //                       currentDopplerData->getRampFlags( ),
//                    //                           "odfTestRampFlags_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
//                    counter++;
//                }
//            }
//        }
//
//        Eigen::MatrixXd dataMatrix = Eigen::MatrixXd( rampFlags.size( ), 8 );
//        for( unsigned int j = 0; j < rampFlags.size( ); j++ )
//        {
//            dataMatrix( j, 0 ) = observationTimes.at( j );
//            dataMatrix( j, 1 ) = observables.at( j );
//            dataMatrix( j, 2 ) = boost::lexical_cast< double >( receivingStation.at( j ) );
//            dataMatrix( j, 3 ) = boost::lexical_cast< double >( transmittingStation.at( j ) );
//            dataMatrix( j, 4 ) = referenceFrequencies.at( j );
//            dataMatrix( j, 5 ) = static_cast< double >( rampFlags.at( j ) );
//            dataMatrix( j, 6 ) = static_cast< double >( receivingRampFrequency.at( j ) );
//            dataMatrix( j, 7 ) = static_cast< double >( transmittingRampFrequency.at( j ) );
//        }
//
//        std::cout << dataMatrix << std::endl << std::endl;
//
//
//        input_output::writeMatrixToFile(
//                    dataMatrix, "odfFileSummary_" + std::to_string( it->first ), 16,
//                    "/Users/pipas/Documents/dsn_trk-2-18/" ); // "/home/dominic/Documents/"
//    }


BOOST_AUTO_TEST_SUITE_END( )

}

}