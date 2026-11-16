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

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include <iostream>
#include <utility>
#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/io/preProcessIfmsFile.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/createObservationCollection.h"

// Some simplifications for shorter lines
using namespace tudat::input_output;
using namespace tudat::basic_astrodynamics;
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
    spice_interface::loadStandardSpiceKernels( );
    std::vector< std::string > ifmsFiles;
    ifmsFiles.push_back( paths::getTudatTestDataPath( ) + "/estrack_n_way_doppler_observation_model/M32ICL2L02_D2X_133630120_00.TAB.txt" );
    ifmsFiles.push_back( paths::getTudatTestDataPath( ) + "/estrack_n_way_doppler_observation_model/M32ICL2L02_D2X_133630203_00.TAB.txt" );
    ifmsFiles.push_back( paths::getTudatTestDataPath( ) + "/estrack_n_way_doppler_observation_model/M32ICL2L02_D2X_133631902_00.TAB.txt" );

    std::vector< bool > applyTroposphereCorrection = { true, false, true };

    std::vector< int > rawIfmsFileSizes = { 2509, 5934, 8460 };
    std::vector< double > tenthObservation = { 8421021939.244300 + 0.010986, 8420857180.689899, 8420892014.600800 - 0.023160 };
    std::vector< std::string > tenthObservationTimeUtc = { "2013-12-29T01:20:55.500",
                                                           "2013-12-29T02:03:50.500",
                                                           "2013-12-29T19:02:57.500" };

    std::vector< int > linesToBeSkipped = { 1, 1, 0 };

    for( unsigned int i = 0; i < ifmsFiles.size( ) + 1; i++ )
    {
        std::vector< std::shared_ptr< data::TrackingData< double, Time > > > trackingData;
        if( i < ifmsFiles.size( ) )
        {
            trackingData = readIfmsFiles< double, Time >( std::vector< std::string >( { ifmsFiles.at( i ) } ),
                                                          "MeX",
                                                          std::vector< std::string >( { "NWNORCIA" } ),
                                                          "Earth",
                                                          applyTroposphereCorrection.at( i ) )
                                   .first;
        }
        else
        {
            for( unsigned int j = 0; j < ifmsFiles.size( ); j++ )
            {
                std::vector< std::shared_ptr< data::TrackingData< double, Time > > > currentTrackingData =
                        readIfmsFiles< double, Time >( std::vector< std::string >( { ifmsFiles.at( j ) } ),
                                                       "MeX",
                                                       std::vector< std::string >( { "NWNORCIA" } ),
                                                       "Earth",
                                                       applyTroposphereCorrection.at( j ) )
                                .first;
                trackingData.insert( trackingData.end( ), currentTrackingData.begin( ), currentTrackingData.end( ) );
            }
        }

        if( i < ifmsFiles.size( ) )
        {
            BOOST_REQUIRE_EQUAL( trackingData.size( ), 1 );
            BOOST_CHECK_EQUAL( trackingData.at( 0 )->getObservableType( ), "DsnNWayAveragedDoppler" );
            BOOST_CHECK_EQUAL( trackingData.at( 0 )->getReferenceLinkEnd( ), "receiver" );
            BOOST_CHECK_EQUAL( trackingData.at( 0 )->getTimeScale( ), "UTC" );
            BOOST_CHECK_EQUAL( trackingData.at( 0 )->getNumberOfObservations( ), rawIfmsFileSizes.at( i ) - linesToBeSkipped.at( i ) );
            BOOST_CHECK_CLOSE_FRACTION( trackingData.at( 0 )->getObservations( ).at( 9 - linesToBeSkipped.at( i ) )( 0 ),
                                        tenthObservation.at( i ),
                                        std::numeric_limits< double >::epsilon( ) );
            Time utcTimeTest = timeFromIsoString< Time >( tenthObservationTimeUtc.at( i ) );
            BOOST_CHECK_SMALL(
                    static_cast< double >( trackingData.at( 0 )->getObservationEpochs( ).at( 9 - linesToBeSkipped.at( i ) ) - utcTimeTest ),
                    1.0E-12 );

            simulation_setup::BodyListSettings bodySettings = simulation_setup::getDefaultBodySettings( { "Earth" } );
            bodySettings.at( "Earth" )->groundStationSettings = simulation_setup::getRadioTelescopeStationSettings( );
            simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );
            auto observationCollection = observation_models::createObservationCollection< double, Time >( trackingData, bodies );
            std::vector< Time > observationCollectionEpochs = observationCollection->getConcatenatedTimeVector( );
            const Eigen::Vector3d earthFixedPosition =
                    bodies.getBody( "Earth" )->getGroundStation( "NWNORCIA" )->getNominalStationState( )->getNominalCartesianPosition( );
            Time tdbTimeTest = earth_orientation::TerrestrialTimeScaleConverter( ).getCurrentTime< Time >(
                    utc_scale, tdb_scale, utcTimeTest, earthFixedPosition );
            BOOST_CHECK_SMALL( static_cast< double >( observationCollectionEpochs.at( 9 - linesToBeSkipped.at( i ) ) - tdbTimeTest ),
                               1.0E-12 );
        }
        else
        {
            BOOST_REQUIRE_EQUAL( trackingData.size( ), 3 );
            BOOST_CHECK_EQUAL( trackingData.at( 0 )->getNumberOfObservations( ), rawIfmsFileSizes.at( 0 ) - linesToBeSkipped.at( 0 ) );
            BOOST_CHECK_EQUAL( trackingData.at( 1 )->getNumberOfObservations( ), rawIfmsFileSizes.at( 1 ) - linesToBeSkipped.at( 1 ) );
            BOOST_CHECK_EQUAL( trackingData.at( 2 )->getNumberOfObservations( ), rawIfmsFileSizes.at( 2 ) - linesToBeSkipped.at( 2 ) );

            BOOST_CHECK_CLOSE_FRACTION( trackingData.at( 0 )->getObservations( ).at( 9 - linesToBeSkipped.at( 0 ) )( 0 ),
                                        tenthObservation.at( 0 ),
                                        std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( trackingData.at( 1 )->getObservations( ).at( 9 - linesToBeSkipped.at( 1 ) )( 0 ),
                                        tenthObservation.at( 1 ),
                                        std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( trackingData.at( 2 )->getObservations( ).at( 9 - linesToBeSkipped.at( 2 ) )( 0 ),
                                        tenthObservation.at( 2 ),
                                        std::numeric_limits< double >::epsilon( ) );

            Time utcTimeTest1 = timeFromIsoString< Time >( tenthObservationTimeUtc.at( 0 ) );
            Time utcTimeTest2 = timeFromIsoString< Time >( tenthObservationTimeUtc.at( 1 ) );
            Time utcTimeTest3 = timeFromIsoString< Time >( tenthObservationTimeUtc.at( 2 ) );
            BOOST_CHECK_SMALL( static_cast< double >( trackingData.at( 0 )->getObservationEpochs( ).at( 9 - linesToBeSkipped.at( 0 ) ) -
                                                      utcTimeTest1 ),
                               1.0E-12 );
            BOOST_CHECK_SMALL( static_cast< double >( trackingData.at( 1 )->getObservationEpochs( ).at( 9 - linesToBeSkipped.at( 1 ) ) -
                                                      utcTimeTest2 ),
                               1.0E-12 );
            BOOST_CHECK_SMALL( static_cast< double >( trackingData.at( 2 )->getObservationEpochs( ).at( 9 - linesToBeSkipped.at( 2 ) ) -
                                                      utcTimeTest3 ),
                               1.0E-12 );
        }
    }
}

std::map< Time, double > getIfmsRampStartFrequencies(
        const std::vector< std::shared_ptr< data::TrackingSupplementaryData > >& supplementaryData )
{
    std::map< Time, double > uplinkFrequencies;

    for( const std::shared_ptr< data::TrackingSupplementaryData >& currentSupplementaryData : supplementaryData )
    {
        BOOST_CHECK_EQUAL( currentSupplementaryData->getBodyName( ), "Earth" );
        BOOST_CHECK_EQUAL( currentSupplementaryData->getReferencePointName( ), "NWNORCIA" );
        for( const std::shared_ptr< data::FrequencySupplementaryData >& frequencyData :
             currentSupplementaryData->getFrequencySupplementaryData( ) )
        {
            std::shared_ptr< data::RampedFrequencySupplementaryData > rampedFrequencyData =
                    std::dynamic_pointer_cast< data::RampedFrequencySupplementaryData >( frequencyData );
            BOOST_REQUIRE( rampedFrequencyData != nullptr );
            for( const data::RampedFrequencySupplementaryData::FrequencyRamp& ramp : rampedFrequencyData->getFrequencyRamps( ) )
            {
                uplinkFrequencies[ Time( ramp.startTime_ ) ] = ramp.startFrequency_;
            }
        }
    }

    return uplinkFrequencies;
}

std::map< Time, double > loadIfmsFilesCombined( const std::vector< std::string >& ifmsFileNames )
{
    return getIfmsRampStartFrequencies( readIfmsFiles< double, Time >( ifmsFileNames, "spacecraft", "NWNORCIA", "Earth" ).second );
}

std::map< Time, double > loadIfmsFilesSeparate( const std::vector< std::string >& ifmsFileNames )
{
    std::map< Time, double > uplinkFrequencies;
    for( const std::string& fileName : ifmsFileNames )
    {
        std::map< Time, double > currentUplinkFrequencies = getIfmsRampStartFrequencies(
                readIfmsFiles< double, Time >( std::vector< std::string >( { fileName } ), "spacecraft", "NWNORCIA", "Earth" ).second );
        for( const auto& frequency : currentUplinkFrequencies )
        {
            uplinkFrequencies[ frequency.first ] = frequency.second;
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
        // Load files separately or one by one
        if( test == 0 )
        {
            transmittedFrequencies.push_back( loadIfmsFilesCombined( ifmsFilesOrdered ) );
        }
        else if( test == 1 )
        {
            transmittedFrequencies.push_back( loadIfmsFilesSeparate( ifmsFilesOrdered ) );
        }
        else if( test == 2 )
        {
            transmittedFrequencies.push_back( loadIfmsFilesCombined( ifmsFilesUnordered ) );
        }
        else if( test == 3 )
        {
            transmittedFrequencies.push_back( loadIfmsFilesSeparate( ifmsFilesUnordered ) );
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
