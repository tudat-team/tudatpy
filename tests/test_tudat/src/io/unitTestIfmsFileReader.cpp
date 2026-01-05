/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>
#include <utility>
#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"

#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"

// Some simplifications for shorter lines
using namespace tudat::input_output;
using namespace tudat::basic_astrodynamics;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::ground_stations;
using namespace tudat;

namespace tudat
{
namespace unit_tests
{

//! Starting the entire test suite
BOOST_AUTO_TEST_SUITE( test_ifms_file_reader );

//! Test reading of mars express IFMS files.
BOOST_AUTO_TEST_CASE( testIfmsFileReader )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies.

    simulation_setup::BodyListSettings bodySettings = simulation_setup::getDefaultBodySettings( { "Earth" }, "SSB", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->groundStationSettings = getRadioTelescopeStationSettings( );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "MeX" );

    std::vector< std::shared_ptr< TrackingTxtFileContents > > rawIfmsFiles;

    rawIfmsFiles.push_back( readIfmsFile(
            paths::getTudatTestDataPath( ) + "/estrack_n_way_doppler_observation_model/M32ICL2L02_D2X_133630120_00.TAB.txt", true ) );
    rawIfmsFiles.push_back( readIfmsFile(
            paths::getTudatTestDataPath( ) + "/estrack_n_way_doppler_observation_model/M32ICL2L02_D2X_133630203_00.TAB.txt", false ) );
    rawIfmsFiles.push_back( readIfmsFile(
            paths::getTudatTestDataPath( ) + "/estrack_n_way_doppler_observation_model/M32ICL2L02_D2X_133631902_00.TAB.txt", true ) );

    std::vector< int > rawIfmsFileSizes = { 2509, 5934, 8460 };
    std::vector< double > tenthObservation = { 8421021939.244300 + 0.010986, 8420857180.689899, 8420892014.600800 - 0.023160 };
    std::vector< std::string > tenthObservationTimeUtc = { "2013-12-29T01:20:55.500",
                                                           "2013-12-29T02:03:50.500",
                                                           "2013-12-29T19:02:57.500" };

    Eigen::Vector3d stationPosition = getCombinedApproximateGroundStationPositions( ).at( "NWNORCIA" );
    auto timeConverter = earth_orientation::createDefaultTimeConverter( );

    std::vector< int > linesToBeSkipped = { 1, 1, 0 };

    for( unsigned int i = 0; i < rawIfmsFiles.size( ) + 1; i++ )
    {
        std::vector< std::shared_ptr< ProcessedTrackingTxtFileContents< double, Time > > > processedIfmsFiles;
        if( i < rawIfmsFiles.size( ) )
        {
            rawIfmsFiles.at( i )->addMetaData( TrackingDataType::receiving_station_name, "NWNORCIA" );
            rawIfmsFiles.at( i )->addMetaData( TrackingDataType::transmitting_station_name, "NWNORCIA" );
            processedIfmsFiles.push_back( std::make_shared< observation_models::ProcessedTrackingTxtFileContents< double, Time > >(
                    rawIfmsFiles.at( i ), "MeX", simulation_setup::getCombinedApproximateGroundStationPositions( ) ) );
        }
        else
        {
            for( unsigned int j = 0; j < rawIfmsFiles.size( ); j++ )
            {
                processedIfmsFiles.push_back( std::make_shared< observation_models::ProcessedTrackingTxtFileContents< double, Time > >(
                        rawIfmsFiles.at( j ), "MeX", simulation_setup::getCombinedApproximateGroundStationPositions( ) ) );
            }
        }

        setTrackingDataInformationInBodies( processedIfmsFiles, bodies, { dsn_n_way_averaged_doppler } );

        auto observationCollection = observation_models::createTrackingTxtFilesObservationCollection< double, Time >(
                processedIfmsFiles, { dsn_n_way_averaged_doppler } );

        std::vector< Time > observationTimes = observationCollection->getConcatenatedTimeVector( );
        Eigen::VectorXd observationValues = observationCollection->getObservationVector( );

        if( i < rawIfmsFiles.size( ) )
        {
            BOOST_CHECK_EQUAL( observationValues.rows( ), rawIfmsFileSizes.at( i ) - linesToBeSkipped.at( i ) );
            BOOST_CHECK_CLOSE_FRACTION( observationValues( 9 - linesToBeSkipped.at( i ) ),
                                        tenthObservation.at( i ),
                                        std::numeric_limits< double >::epsilon( ) );
            Time utcTime = timeConverter->getCurrentTime< Time >(
                    tdb_scale, utc_scale, observationTimes.at( 9 - linesToBeSkipped.at( i ) ), stationPosition );
            Time utcTimeTest = timeFromIsoString< Time >( tenthObservationTimeUtc.at( i ) );
            BOOST_CHECK_SMALL( static_cast< double >( utcTime - utcTimeTest ), 1.0E-12 );
        }
        else
        {
            BOOST_CHECK_EQUAL( observationValues.rows( ),
                               rawIfmsFileSizes.at( 0 ) + rawIfmsFileSizes.at( 1 ) + rawIfmsFileSizes.at( 2 ) - linesToBeSkipped.at( 0 ) -
                                       linesToBeSkipped.at( 1 ) - linesToBeSkipped.at( 2 ) );
            BOOST_CHECK_CLOSE_FRACTION( observationValues( 9 - linesToBeSkipped.at( 0 ) ),
                                        tenthObservation.at( 0 ),
                                        std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION(
                    observationValues( 9 + rawIfmsFileSizes.at( 0 ) - linesToBeSkipped.at( 0 ) - linesToBeSkipped.at( 1 ) ),
                    tenthObservation.at( 1 ),
                    std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( observationValues( 9 + rawIfmsFileSizes.at( 0 ) + rawIfmsFileSizes.at( 1 ) -
                                                           linesToBeSkipped.at( 0 ) - linesToBeSkipped.at( 1 ) - linesToBeSkipped.at( 2 ) ),
                                        tenthObservation.at( 2 ),
                                        std::numeric_limits< double >::epsilon( ) );

            Time utcTime1 = timeConverter->getCurrentTime< Time >(
                    tdb_scale, utc_scale, observationTimes.at( 9 - linesToBeSkipped.at( 0 ) ), stationPosition );
            Time utcTimeTest1 = timeFromIsoString< Time >( tenthObservationTimeUtc.at( 0 ) );
            BOOST_CHECK_SMALL( static_cast< double >( utcTime1 - utcTimeTest1 ), 1.0E-12 );

            Time utcTime2 = timeConverter->getCurrentTime< Time >(
                    tdb_scale,
                    utc_scale,
                    observationTimes.at( 9 + rawIfmsFileSizes.at( 0 ) - linesToBeSkipped.at( 0 ) - linesToBeSkipped.at( 1 ) ),
                    stationPosition );
            Time utcTimeTest2 = timeFromIsoString< Time >( tenthObservationTimeUtc.at( 1 ) );
            BOOST_CHECK_SMALL( static_cast< double >( utcTime2 - utcTimeTest2 ), 1.0E-12 );

            Time utcTime3 = timeConverter->getCurrentTime< Time >(
                    tdb_scale,
                    utc_scale,
                    observationTimes.at( 9 + rawIfmsFileSizes.at( 0 ) + rawIfmsFileSizes.at( 1 ) - linesToBeSkipped.at( 0 ) -
                                         linesToBeSkipped.at( 1 ) - linesToBeSkipped.at( 2 ) ),
                    stationPosition );
            Time utcTimeTest3 = timeFromIsoString< Time >( tenthObservationTimeUtc.at( 2 ) );
            BOOST_CHECK_SMALL( static_cast< double >( utcTime3 - utcTimeTest3 ), 1.0E-12 );
        }
    }
}

std::map< Time, double > loadIfmsFilesCombined( const std::vector< std::string >& ifmsFileNames, SystemOfBodies& bodies )
{
    std::map< Time, double > uplinkFrequencies;

    // Build a single observation collection from all IFMS files
    std::shared_ptr< ObservationCollection< double, Time > > observationCollectionCombined =
            createIfmsObservedObservationCollectionFromFiles< double, Time >(
                    ifmsFileNames, bodies, "spacecraft", "NWNORCIA", FrequencyBands::x_band, FrequencyBands::x_band );

    // Retrieve frequency calculator
    std::shared_ptr< ground_stations::GroundStation > nwnorcia = bodies.at( "Earth" )->getGroundStation( "NWNORCIA" );
    std::shared_ptr< StationFrequencyInterpolator > freqCalc = nwnorcia->getTransmittingFrequencyCalculator( );

    // Extract concatenated observation times and compute frequencies
    std::vector< Time > epochs = observationCollectionCombined->getConcatenatedObservationTimes( );
    for( double epoch : epochs )
    {
        uplinkFrequencies[ epoch ] = freqCalc->getTemplatedCurrentFrequency( epoch );
    }
    return uplinkFrequencies;
}

std::map< Time, double > loadIfmsFilesSeparate( const std::vector< std::string >& ifmsFileNames, SystemOfBodies& bodies )
{
    std::map< Time, double > uplinkFrequencies;

    // Build a separate observation collection from each IFMS file
    std::vector< std::shared_ptr< ObservationCollection< double, Time > > > observationCollections;
    for( const std::string& fileName : ifmsFileNames )
    {
        observationCollections.push_back(
                createIfmsObservedObservationCollectionFromFiles< double, Time >( std::vector< std::string >( { fileName } ),
                                                                                  bodies,
                                                                                  "spacecraft",
                                                                                  "NWNORCIA",
                                                                                  FrequencyBands::x_band,
                                                                                  FrequencyBands::x_band ) );
    }

    // Retrieve frequency calculator
    std::shared_ptr< ground_stations::GroundStation > nwnorcia = bodies.at( "Earth" )->getGroundStation( "NWNORCIA" );
    std::shared_ptr< StationFrequencyInterpolator > freqCalc = nwnorcia->getTransmittingFrequencyCalculator( );

    for( unsigned int i = 0; i < observationCollections.size( ); i++ )
    {
        // Extract concatenated observation times and compute frequencies
        std::vector< Time > epochs = observationCollections.at( i )->getConcatenatedObservationTimes( );
        for( double epoch : epochs )
        {
            uplinkFrequencies[ epoch ] = freqCalc->getTemplatedCurrentFrequency( epoch );
        }
    }
    return uplinkFrequencies;
}

BOOST_AUTO_TEST_CASE( testMultipleIfmsFileReader )
{
    // Adjust this to point to the directory containing your IFMS files
    const std::string SOURCE_DIR = paths::getTudatTestDataPath( ) + "/estrack_n_way_doppler_observation_model/";

    // Define IFMS file paths
    std::vector< std::string > ifmsFilesOrdered;
    ifmsFilesOrdered.push_back( SOURCE_DIR + "M32ICL1L02_D2X_133612305_00.TAB" );
    ifmsFilesOrdered.push_back( SOURCE_DIR + "M32ICL1L02_D2X_133621819_00.TAB" );
    ifmsFilesOrdered.push_back( SOURCE_DIR + "M32ICL1L02_D2X_133630120_00.TAB" );

    std::vector< std::string > ifmsFilesUnordered;
    ifmsFilesUnordered.push_back( ifmsFilesOrdered.at( 2 ) );
    ifmsFilesUnordered.push_back( ifmsFilesOrdered.at( 1 ) );
    ifmsFilesUnordered.push_back( ifmsFilesOrdered.at( 0 ) );

    std::vector< std::map< Time, double > > transmittedFrequencies;
    for( int test = 0; test < 4; test++ )
    {
        spice_interface::loadStandardSpiceKernels( );

        // Create bodies
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Earth" );
        BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, "SSB", "J2000" );
        bodySettings.at( "Earth" )->groundStationSettings = getRadioTelescopeStationSettings( );
        bodySettings.addSettings( "spacecraft" );
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Load files separately or one by one
        if( test == 0 )
        {
            transmittedFrequencies.push_back( loadIfmsFilesCombined( ifmsFilesOrdered, bodies ) );
        }
        else if( test == 1 )
        {
            transmittedFrequencies.push_back( loadIfmsFilesSeparate( ifmsFilesOrdered, bodies ) );
        }
        else if( test == 2 )
        {
            transmittedFrequencies.push_back( loadIfmsFilesCombined( ifmsFilesUnordered, bodies ) );
        }
        else if( test == 3 )
        {
            transmittedFrequencies.push_back( loadIfmsFilesSeparate( ifmsFilesUnordered, bodies ) );
        }
    }

    BOOST_CHECK_EQUAL( transmittedFrequencies.at( 0 ).size( ), transmittedFrequencies.at( 1 ).size( ) );
    BOOST_CHECK_EQUAL( transmittedFrequencies.at( 0 ).size( ), transmittedFrequencies.at( 2 ).size( ) );
    BOOST_CHECK_EQUAL( transmittedFrequencies.at( 0 ).size( ), transmittedFrequencies.at( 3 ).size( ) );

    for( auto it : transmittedFrequencies.at( 0 ) )
    {
        BOOST_CHECK_CLOSE_FRACTION( transmittedFrequencies.at( 0 ).at( it.first ),
                                    transmittedFrequencies.at( 1 ).at( it.first ),
                                    10.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( transmittedFrequencies.at( 0 ).at( it.first ),
                                    transmittedFrequencies.at( 2 ).at( it.first ),
                                    10.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( transmittedFrequencies.at( 0 ).at( it.first ),
                                    transmittedFrequencies.at( 3 ).at( it.first ),
                                    10.0 * std::numeric_limits< double >::epsilon( ) );
    }
}

// End test suite
BOOST_AUTO_TEST_SUITE_END( );

}  // namespace unit_tests

}  // namespace tudat
