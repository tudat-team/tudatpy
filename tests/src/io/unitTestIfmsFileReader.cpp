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
using namespace tudat;

namespace tudat
{
namespace unit_tests
{



//! Starting the entire test suite
BOOST_AUTO_TEST_SUITE(test_ifms_file_reader);

//! Test reading of mars express IFMS files.
BOOST_AUTO_TEST_CASE(testIfmsFileReader)
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies.

    simulation_setup::BodyListSettings bodySettings =
        simulation_setup::getDefaultBodySettings( { "Earth"}, "SSB", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->groundStationSettings = getRadioTelescopeStationSettings( );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "MeX" );

    std::vector< std::shared_ptr< TrackingTxtFileContents > > rawIfmsFiles;

    rawIfmsFiles.push_back( readIfmsFile( paths::getTudatTestDataPath( )  + "/estrack_n_way_doppler_observation_model/M32ICL2L02_D2X_133630120_00.TAB.txt", true ) );
    rawIfmsFiles.push_back( readIfmsFile( paths::getTudatTestDataPath( )  + "/estrack_n_way_doppler_observation_model/M32ICL2L02_D2X_133630203_00.TAB.txt", false ) );
    rawIfmsFiles.push_back( readIfmsFile( paths::getTudatTestDataPath( )  + "/estrack_n_way_doppler_observation_model/M32ICL2L02_D2X_133631902_00.TAB.txt", true ) );

    std::vector< int > rawIfmsFileSizes = { 2509, 5934, 8460 };
    std::vector< double > tenthObservation = { 8421021939.244300 + 0.010986, 8420857180.689899, 8420892014.600800 - 0.023160};
    std::vector< std::string > tenthObservationTimeUtc = { "2013-12-29T01:20:55.500", "2013-12-29T02:03:50.500", "2013-12-29T19:02:57.500" };

    Eigen::Vector3d stationPosition = getCombinedApproximateGroundStationPositions( ).at( "NWNORCIA" );
    auto timeConverter = earth_orientation::defaultTimeConverter;

    for( unsigned int i = 0; i < rawIfmsFiles.size( ) + 1; i++ )
    {
        std::vector< std::shared_ptr< ProcessedTrackingTxtFileContents< double, Time > > > processedIfmsFiles;
        if( i < rawIfmsFiles.size( ) )
        {
            rawIfmsFiles.at( i )->addMetaData( TrackingDataType::receiving_station_name, "NWNORCIA" );
            rawIfmsFiles.at( i )->addMetaData( TrackingDataType::transmitting_station_name, "NWNORCIA" );
            processedIfmsFiles.push_back( std::make_shared<observation_models::ProcessedTrackingTxtFileContents< double, Time > >(
                rawIfmsFiles.at( i ), "MeX", simulation_setup::getCombinedApproximateGroundStationPositions( ) ) );
        }
        else
        {
            for( unsigned int j = 0; j < rawIfmsFiles.size( ); j++ )
            {
                processedIfmsFiles.push_back( std::make_shared<observation_models::ProcessedTrackingTxtFileContents< double, Time  > >(
                    rawIfmsFiles.at( j ), "MeX", simulation_setup::getCombinedApproximateGroundStationPositions( ) ) );
            }
        }

        setTrackingDataInformationInBodies( processedIfmsFiles, bodies, { dsn_n_way_averaged_doppler } );

        auto observationCollection = observation_models::createTrackingTxtFilesObservationCollection<double, Time >(
            processedIfmsFiles, { dsn_n_way_averaged_doppler } );


        std::vector< Time > observationTimes = observationCollection->getConcatenatedTimeVector( );
        Eigen::VectorXd observationValues = observationCollection->getObservationVector( );

        if( i < rawIfmsFiles.size( ) )
        {
            BOOST_CHECK_EQUAL( observationValues.rows( ), rawIfmsFileSizes.at( i ) );
            BOOST_CHECK_CLOSE_FRACTION( observationValues( 9 ), tenthObservation.at( i ), std::numeric_limits< double >::epsilon( ) );
            Time utcTime = timeConverter->getCurrentTime< Time >( tdb_scale, utc_scale, observationTimes.at( 9 ), stationPosition );
            Time utcTimeTest = timeFromIsoString< Time >( tenthObservationTimeUtc.at( i ) );
            BOOST_CHECK_SMALL( static_cast< double >( utcTime - utcTimeTest ), 1.0E-12 );
        }
        else
        {
            BOOST_CHECK_EQUAL( observationValues.rows( ), rawIfmsFileSizes.at( 0 ) + rawIfmsFileSizes.at( 1 ) + rawIfmsFileSizes.at( 2 ) );
            BOOST_CHECK_CLOSE_FRACTION( observationValues( 9 ), tenthObservation.at( 0 ), std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( observationValues( 9 + rawIfmsFileSizes.at( 0 ) ), tenthObservation.at( 1 ), std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( observationValues( 9 + rawIfmsFileSizes.at( 0 )  + rawIfmsFileSizes.at( 1 ) ), tenthObservation.at( 2 ), std::numeric_limits< double >::epsilon( ) );

            Time utcTime1 =  timeConverter->getCurrentTime< Time >( tdb_scale, utc_scale, observationTimes.at( 9 ), stationPosition );
            Time utcTimeTest1 = timeFromIsoString< Time >( tenthObservationTimeUtc.at( 0 ) );
            BOOST_CHECK_SMALL( static_cast< double >( utcTime1 - utcTimeTest1 ), 1.0E-12 );

            Time utcTime2 =  timeConverter->getCurrentTime< Time >( tdb_scale, utc_scale, observationTimes.at( 9 + rawIfmsFileSizes.at( 0 ) ), stationPosition );
            Time utcTimeTest2 = timeFromIsoString< Time >( tenthObservationTimeUtc.at( 1 ) );
            BOOST_CHECK_SMALL( static_cast< double >( utcTime2 - utcTimeTest2 ), 1.0E-12 );

            Time utcTime3 =  timeConverter->getCurrentTime< Time >( tdb_scale, utc_scale, observationTimes.at( 9 + rawIfmsFileSizes.at( 0 ) + rawIfmsFileSizes.at( 1 ) ), stationPosition );
            Time utcTimeTest3 = timeFromIsoString< Time >( tenthObservationTimeUtc.at( 2 ) );
            BOOST_CHECK_SMALL( static_cast< double >( utcTime3 - utcTimeTest3 ), 1.0E-12 );

        }
    }
}

// End test suite
BOOST_AUTO_TEST_SUITE_END();

}// namespace unit_tests

}// namespace tudat
