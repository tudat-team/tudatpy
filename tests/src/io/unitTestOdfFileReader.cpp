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

#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_odf_file_reader )

//! Checks parsed binary file (odf07155.odf) against values in txt file (odf07155.txt)
//! Files available at https://pds-geosciences.wustl.edu/dataserv/radio_science.htm (see radio science documentation bundle)
BOOST_AUTO_TEST_CASE( testSingleOdfFileReader )
{
    std::string file = tudat::paths::getTudatTestDataPath( )  + "/odf07155.odf";
    std::shared_ptr< input_output::OdfRawFileContents > rawOdfContents =
            std::make_shared< input_output::OdfRawFileContents >( file );

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
    std::shared_ptr< input_output::OdfCommonDataBlock > commonDataBlock = rawOdfContents->getDataBlocks( ).at( 0 )->getCommonDataBlock( );
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
                    rawOdfContents->getDataBlocks( ).at( 0 )->getObservableSpecificDataBlock( ) );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->dataType_, 11 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getReceiverChannel(), 1 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getSpacecraftId(), 236 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getReceiverExciterFlag(), 1 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getReferenceFrequency( ), 2299812417.000 );
//    BOOST_CHECK_EQUAL ( dopplerDataBlock->reservedSegment_, 0 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getCompressionTime( ), 60 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getTransmittingStationUplinkDelay( ), 0 );

    // Block 19: 2-way doppler
    commonDataBlock = rawOdfContents->getDataBlocks( ).at( 19 )->getCommonDataBlock( );
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
            rawOdfContents->getDataBlocks( ).at( 19 )->getObservableSpecificDataBlock( ) );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->dataType_, 12 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getReceiverChannel( ), 1 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getSpacecraftId( ), 236 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getReceiverExciterFlag( ), 1 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getReferenceFrequency( ), 7177648275.000 );
//    BOOST_CHECK_EQUAL ( dopplerDataBlock->reservedSegment_, 0 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getCompressionTime( ), 60 );
    BOOST_CHECK_EQUAL ( dopplerDataBlock->getTransmittingStationUplinkDelay( ), 0 );

    // Block 23: sequential range
    commonDataBlock = rawOdfContents->getDataBlocks( ).at( 23 )->getCommonDataBlock( );
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
                    rawOdfContents->getDataBlocks( ).at( 23 )->getObservableSpecificDataBlock( ) );
    BOOST_CHECK_EQUAL ( rangeDataBlock->dataType_, 37 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->lowestRangingComponent_, 14 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->getSpacecraftId( ), 236 );
    //BOOST_CHECK_EQUAL ( rangeDataBlock->reservedBlock_, 1 );
    BOOST_CHECK_EQUAL ( int(rangeDataBlock->reservedBlock_), 1 );
    //BOOST_CHECK_EQUAL ( rangeDataBlock->getReferenceFrequency( ), 7177004669.452 );
    BOOST_CHECK_EQUAL ( std::round( 1000.0 * rangeDataBlock->getReferenceFrequency( ) ), 7177004669452 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->uplinkCoderInPhaseTimeOffset_, 774 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->compositeTwo_, 400000 );
    BOOST_CHECK_EQUAL ( rangeDataBlock->getTransmittingStationUplinkDelay( ), 0 );

    // Check number of data blocks
    BOOST_CHECK_EQUAL ( rawOdfContents->getDataBlocks( ).size( ), 2233 - 4 - 1 );

    // Check number of ramp groups
    BOOST_CHECK_EQUAL ( rawOdfContents->getRampBlocks( ).size( ), 3 );

    // Ramp block for DSS-63
    std::shared_ptr< input_output::OdfRampBlock > rampBlock =
            std::dynamic_pointer_cast< input_output::OdfRampBlock >(
                    rawOdfContents->getRampBlocks( ).at( 63 ).at( 5 ) );
    BOOST_CHECK_EQUAL ( rampBlock->getRampStartTime( ), 1812101116.000000000 );
    BOOST_CHECK_EQUAL ( rampBlock->getRampRate( ), 0.095680000 );
    BOOST_CHECK_EQUAL ( rampBlock->getRampStartFrequency( ), 7177004073.170830727 );
    BOOST_CHECK_EQUAL ( rampBlock->getRampEndTime( ), 1812101962.000000000 );

    // Check number of ramp blocks
    BOOST_CHECK_EQUAL ( rawOdfContents->getRampBlocks( ).at( 63 ).size( ), 2331 - 2233 - 1 );

    // Write to txt file
    rawOdfContents->writeOdfToTextFile(tudat::paths::getTudatTestDataPath( ) + "/odf07155_out.txt");
}

BOOST_AUTO_TEST_CASE( testProcessSingleOdfFile )
{
    spice_interface::loadStandardSpiceKernels( );

    // Create system of bodies
    std::string spacecraftName = "Messenger";

    std::vector< std::string > bodiesToCreate = { "Earth" };

    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate );
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Load ODF file
    std::string file = tudat::paths::getTudatTestDataPath( )  + "/mess_rs_09121_121_odf.dat";
    std::shared_ptr< input_output::OdfRawFileContents > rawOdfContents =
            std::make_shared< input_output::OdfRawFileContents >( file );

    // Process ODF file
    std::shared_ptr< observation_models::ProcessedOdfFileContents > processedOdfFileContents =
            std::make_shared< observation_models::ProcessedOdfFileContents >( rawOdfContents, spacecraftName );

    std::pair< double, double > startAndEndTimeTdb = processedOdfFileContents->getStartAndEndTime( );

    // double startTime = startAndEndTimeTdb.first;
    // TODO: delete following line and replace by previous one after processing of ODF data type 11 (1-way Doppler) is implemented
    double startTime = observation_models::ProcessedOdfFileContentsPrivateFunctionTest::computeObservationTimesTdbFromJ2000(
            processedOdfFileContents,
            "DSS-25",
            rawOdfContents->getDataBlocks( ).front( )->getCommonDataBlock( )->getObservableTime( ) );

    double endTime = startAndEndTimeTdb.second;

    // Compare start and end time with values in LBL file
    for ( unsigned int i = 0; i < 2; ++i )
    {
        double time, expectedFractionOfDay;
        int expectedYear, expectedMonth, expectedDay;
        // Expected start time: 2009-05-01T12:41:18.000 UTC
        if ( i == 0 )
        {
            time = startTime;
            expectedYear = 2009;
            expectedMonth = 5;
            expectedDay = 1;
            expectedFractionOfDay = 12.0 / 24.0 + 41.0 / ( 24.0 * 60.0 ) + 18.0 / ( 24.0 * 3600.0 );
        }
        // Expected end time: 2009-05-01T22:44:35.000
        else
        {
            time = endTime;
            expectedYear = 2009;
            expectedMonth = 5;
            expectedDay = 1;
            expectedFractionOfDay = 22.0 / 24.0 + 44.0 / ( 24.0 * 60.0 ) + 35.0 / ( 24.0 * 3600.0 );
        }

        // Get UTC time
        earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter = earth_orientation::TerrestrialTimeScaleConverter( );
        double timeUtc = timeScaleConverter.getCurrentTime< double >( basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, time );

        // Get UTC calendar date
        int year, month, day;
        double fractionOfDay;
        iauJd2cal( basic_astrodynamics::JULIAN_DAY_ON_J2000, timeUtc / physical_constants::JULIAN_DAY, &year, &month, &day, &fractionOfDay );

        BOOST_CHECK_EQUAL ( year, expectedYear );
        BOOST_CHECK_EQUAL ( month, expectedMonth );
        BOOST_CHECK_EQUAL ( day, expectedDay );
        BOOST_CHECK_CLOSE_FRACTION( fractionOfDay, expectedFractionOfDay, 1e-10 );
    }

}

BOOST_AUTO_TEST_CASE( testProcessMultipleOdfFile )
{
    spice_interface::loadStandardSpiceKernels( );

    // Create system of bodies
    std::string spacecraftName = "MRO";

    std::vector< std::string > bodiesToCreate = { "Earth" };

    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate );
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Load ODF files
    std::shared_ptr< input_output::OdfRawFileContents > rawOdfContents1 = std::make_shared< input_output::OdfRawFileContents >(
            tudat::paths::getTudatTestDataPath( )  + "/mromagr2017_097_1335xmmmv1.odf" );
    std::shared_ptr< input_output::OdfRawFileContents > rawOdfContents2 = std::make_shared< input_output::OdfRawFileContents >(
            tudat::paths::getTudatTestDataPath( )  + "/mromagr2017_098_1555xmmmv1.odf" );
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector = { rawOdfContents2, rawOdfContents1 };

    // Process ODF files
    std::shared_ptr< observation_models::ProcessedOdfFileContents > processedOdfFileContents =
            std::make_shared< observation_models::ProcessedOdfFileContents >( rawOdfDataVector, spacecraftName );

    std::pair< double, double > startAndEndTimeTdb = processedOdfFileContents->getStartAndEndTime( );

    // double startTime = startAndEndTimeTdb.first;
    // TODO: delete following line and replace by previous one after processing of ODF data type 11 (1-way Doppler) is implemented
    double startTime = observation_models::ProcessedOdfFileContentsPrivateFunctionTest::computeObservationTimesTdbFromJ2000(
            processedOdfFileContents,
            "DSS-55",
            rawOdfContents1->getDataBlocks( ).front( )->getCommonDataBlock( )->getObservableTime( ) );

    double endTime = startAndEndTimeTdb.second;

    // Compare start and end time with values in LBL file
    for ( unsigned int i = 0; i < 2; ++i )
    {
        double time, expectedFractionOfDay;
        int expectedYear, expectedMonth, expectedDay;
        // Expected start time: 2017-097T13:35:34.5 UTC
        if ( i == 0 )
        {
            time = startTime;
            expectedYear = 2017;
            expectedMonth = 4;
            expectedDay = 7;
            expectedFractionOfDay = 13.0 / 24.0 + 35.0 / ( 24.0 * 60.0 ) + 34.5 / ( 24.0 * 3600.0 );
        }
        // Expected end time: 2017-100T03:50:21.5 UTC
        else
        {
            time = endTime;
            expectedYear = 2017;
            expectedMonth = 4;
            expectedDay = 10;
            expectedFractionOfDay = 3.0 / 24.0 + 50.0 / ( 24.0 * 60.0 ) + 21.5 / ( 24.0 * 3600.0 );
        }

        // Get UTC time
        earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter = earth_orientation::TerrestrialTimeScaleConverter( );
        double timeUtc = timeScaleConverter.getCurrentTime< double >( basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, time );

        // Get UTC calendar date
        int year, month, day;
        double fractionOfDay;
        iauJd2cal( basic_astrodynamics::JULIAN_DAY_ON_J2000, timeUtc / physical_constants::JULIAN_DAY, &year, &month, &day, &fractionOfDay );

        BOOST_CHECK_EQUAL ( year, expectedYear );
        BOOST_CHECK_EQUAL ( month, expectedMonth );
        BOOST_CHECK_EQUAL ( day, expectedDay );
        BOOST_CHECK_CLOSE_FRACTION( fractionOfDay, expectedFractionOfDay, 1e-9 );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}