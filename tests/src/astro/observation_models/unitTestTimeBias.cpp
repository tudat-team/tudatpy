/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>


#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/angularPositionObservationModel.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;


BOOST_AUTO_TEST_SUITE( test_time_bias_model )


BOOST_AUTO_TEST_CASE( testConstantTimeBias )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );

    // Specify initial time
    double initialTime = 0.0;
    double finalTime = initialTime + 7.0 * 86400.0;
    double buffer = 10.0 * 3600.0;

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings =
            getDefaultBodySettings( bodiesToCreate, initialTime - buffer, finalTime + buffer );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Define link ends for observations.
    LinkDefinition linkEnds;
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "" );
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars" , ""  );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimePerturbingBodies ) );

    // Create observation settings
    std::shared_ptr< ObservationModelSettings > idealObservableSettings = std::make_shared< ObservationModelSettings >
            ( angular_position, linkEnds, lightTimeCorrectionSettings,
              std::make_shared< ConstantTimeBiasSettings >( 0.0, receiver ) );

    double timeBias = 2.0;
    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
    biasSettingsList.push_back( std::make_shared< ConstantTimeBiasSettings >( timeBias, receiver ) );
    biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( ( Eigen::Vector2d( ) << 0.0, 0.0 ).finished( ), true ) );
    std::shared_ptr< MultipleObservationBiasSettings > biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );

    std::shared_ptr< ObservationModelSettings > biasedObservableSettings = std::make_shared< ObservationModelSettings >
            ( angular_position, linkEnds, lightTimeCorrectionSettings, biasSettings );

    // Create observation models.
    std::shared_ptr< ObservationModel< 2, double, double > > idealObservationModel =
           ObservationModelCreator< 2, double, double >::createObservationModel( idealObservableSettings, bodies );

    std::shared_ptr< ObservationModel< 2, double, double > > biasedObservationModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( biasedObservableSettings, bodies );

    // Compute observation separately with two functions.
    double obsTime = ( finalTime + initialTime ) / 2.0;
    std::vector< double > unbiasedLinkEndTimes, biasedLinkEndTimes, linkEndTimesBiasedTime;
    std::vector< Eigen::Vector6d > unbiasedLinkEndStates, biasedLinkEndStates, linkEndStatesBiasedTime;

    Eigen::Vector2d biasedObservation = biasedObservationModel->computeObservationsWithLinkEndData( obsTime, receiver, biasedLinkEndTimes, biasedLinkEndStates );
    idealObservationModel->computeObservationsWithLinkEndData( obsTime, receiver, unbiasedLinkEndTimes, unbiasedLinkEndStates );
    Eigen::Vector2d observationAtBiasedTime = idealObservationModel->computeObservationsWithLinkEndData( obsTime - timeBias, receiver, linkEndTimesBiasedTime, linkEndStatesBiasedTime );


    // Check link end times
    std::cout<<unbiasedLinkEndTimes[ 1 ]<<" "<<biasedLinkEndTimes[ 1 ]<<" "<<timeBias<<std::endl;
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( unbiasedLinkEndTimes[ 1 ] - biasedLinkEndTimes[ 1 ] ), timeBias, std::numeric_limits< double >::epsilon( ) );

    // Check link end states
    for ( unsigned int i = 0 ; i < biasedLinkEndStates.size( ) ; i++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                biasedLinkEndStates.at( i ), linkEndStatesBiasedTime.at( i ), std::numeric_limits< double >::epsilon( ) );
    }

    // Check observations
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            biasedObservation, observationAtBiasedTime, std::numeric_limits< double >::epsilon( ) );

}

BOOST_AUTO_TEST_CASE( testArcWiseTimeBias )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );

    // Specify initial time
    double initialTime = 0.0;
    double finalTime = initialTime + 7.0 * 86400.0;
    double buffer = 10.0 * 3600.0;

    std::vector< double > arcs;
    arcs.push_back( initialTime );
    arcs.push_back( initialTime + 2.0 * 86400.0 );
    arcs.push_back( initialTime + 6.0 * 86400.0 );

    // Define arc-wise time biases
    std::vector< double > zeroBiasesPerArc = { 0.0, 0.0, 0.0 };
    std::vector< double > biasesPerArc = { 2.0, 5.0, 10.0 };

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings = getDefaultBodySettings( bodiesToCreate, initialTime - buffer, finalTime + buffer );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Define link ends for observations.
    LinkDefinition linkEnds;
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "" );
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars" , ""  );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimePerturbingBodies ) );

    // Create observation settings
    std::shared_ptr< ObservationModelSettings > idealObservableSettings = std::make_shared< ObservationModelSettings >
            ( angular_position, linkEnds, lightTimeCorrectionSettings,
              std::make_shared< ArcWiseTimeBiasSettings >( arcs, zeroBiasesPerArc, receiver ) );

    std::shared_ptr< ObservationModelSettings > biasedObservableSettings = std::make_shared< ObservationModelSettings >
            ( angular_position, linkEnds, lightTimeCorrectionSettings,
              std::make_shared< ArcWiseTimeBiasSettings >( arcs, biasesPerArc, receiver ) );

    for ( unsigned int k = 0 ; k < biasesPerArc.size( ) ; k++ )
    {
        std::cout << "biases per arc: " << biasesPerArc.at( k ) << "\n\n";
    }

    // Create observation models.
    std::shared_ptr< ObservationModel< 2, double, double > > idealObservationModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( idealObservableSettings, bodies );

    std::shared_ptr< ObservationModel< 2, double, double > > biasedObservationModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( biasedObservableSettings, bodies );


    // Compute observation separately with two functions.
    std::vector< double > obsTimes;
    obsTimes.push_back( ( arcs[ 0 ] + arcs[ 1 ] ) / 2.0 );
    obsTimes.push_back( ( arcs[ 1 ] + arcs[ 2 ] ) / 2.0 );
    obsTimes.push_back( ( arcs[ 2 ] + finalTime ) / 2.0 );

    for ( unsigned int j = 0 ; j < obsTimes.size( ) ; j++ )
    {
        std::vector< double > unbiasedLinkEndTimes, biasedLinkEndTimes, linkEndTimesBiasedTime;
        std::vector< Eigen::Vector6d > unbiasedLinkEndStates, biasedLinkEndStates, linkEndStatesBiasedTime;

        Eigen::Vector2d biasedObservation = biasedObservationModel->computeObservationsWithLinkEndData( obsTimes[ j ], receiver, biasedLinkEndTimes, biasedLinkEndStates );
        idealObservationModel->computeObservationsWithLinkEndData( obsTimes[ j ], receiver, unbiasedLinkEndTimes, unbiasedLinkEndStates );
        Eigen::Vector2d observationAtBiasedTime = idealObservationModel->computeObservationsWithLinkEndData(
                obsTimes[ j ] - biasesPerArc[ j ], receiver, linkEndTimesBiasedTime, linkEndStatesBiasedTime );

        // Check link end times
        BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( unbiasedLinkEndTimes[ 1 ] - biasedLinkEndTimes[ 1 ] ), biasesPerArc[ j ], std::numeric_limits< double >::epsilon( ) );

        // Check link end states
        for ( unsigned int i = 0 ; i < biasedLinkEndStates.size( ) ; i++ )
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    biasedLinkEndStates.at( i ), linkEndStatesBiasedTime.at( i ), std::numeric_limits< double >::epsilon( ) );
        }

        // Check observations
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                biasedObservation, observationAtBiasedTime, std::numeric_limits< double >::epsilon( ) );
    }





}

BOOST_AUTO_TEST_SUITE_END( )

}

}

