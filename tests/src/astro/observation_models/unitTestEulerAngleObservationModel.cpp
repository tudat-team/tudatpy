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
#include <boost/make_shared.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
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


BOOST_AUTO_TEST_SUITE( test_position_obsevable_model )


BOOST_AUTO_TEST_CASE( testPositionObsevableModel )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ observed_body ] = std::make_pair( "Earth" , ""  );


    // Create observation settings
    std::shared_ptr< ObservationModelSettings > observableSettings = std::make_shared< ObservationModelSettings >
            ( euler_angle_313_observable, linkEnds, std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
              std::make_shared< ConstantObservationBiasSettings >(
                  ( Eigen::Vector3d( ) << 0.6, -0.3, 0.2 ).finished( ), true ) );

    // Create observation model.
    std::shared_ptr< ObservationModel< 3, double, double > > observationModel =
           ObservationModelCreator< 3, double, double >::createObservationModel(
                observableSettings, bodies );
    std::shared_ptr< ObservationBias< 3 > > observationBias = observationModel->getObservationBiasCalculator( );


    // Compute observation separately with two functions.
    double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    Eigen::Vector3d observation = observationModel->computeObservations(
                observationTime, observed_body );

    Eigen::Vector3d observation2 = observationModel->computeObservationsWithLinkEndData(
                observationTime, observed_body, linkEndTimes, linkEndStates );

    // Check size of link end state/time.
    BOOST_CHECK_EQUAL( linkEndTimes.size( ), 1 );
    BOOST_CHECK_EQUAL( linkEndStates.size( ), 1 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                observation, observation2, std::numeric_limits< double >::epsilon( ) );

    // Check link end state/time.
    BOOST_CHECK_CLOSE_FRACTION( observationTime, linkEndTimes[ 0 ], std::numeric_limits< double >::epsilon( ) );

    // Check biased observable
    Eigen::Quaterniond currentRotationToBodyFixedFrame = bodies.at( "Earth" )->getRotationalEphemeris( )->
                getRotationToTargetFrame( observationTime );

    Eigen::Vector3d currentBias =
            observationBias->getObservationBias(
                  std::vector< double >( ), std::vector< Eigen::Vector6d>( ) );
    Eigen::Quaterniond computedRotationToBodyFixedFrame =
            Eigen::Quaterniond(
                            Eigen::AngleAxisd( -( observation( 0 ) - currentBias( 0 ) ), Eigen::Vector3d::UnitZ( ) ) *
                            Eigen::AngleAxisd( -( observation( 1 ) - currentBias( 1 ) ), Eigen::Vector3d::UnitX( ) ) *
                            Eigen::AngleAxisd( -( observation( 2 ) - currentBias( 2 ) ), Eigen::Vector3d::UnitZ( ) ) );

    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( currentRotationToBodyFixedFrame.toRotationMatrix( )( i, j ) -
                    computedRotationToBodyFixedFrame.toRotationMatrix( )( i, j ) ),
                    5.0 * std::numeric_limits< double >::epsilon( ) );
        }
    }

    // Check ideal observable
    observation = observationModel->computeIdealObservations(
                observationTime, observed_body );
    observation2 = observationModel->computeIdealObservationsWithLinkEndData(
                observationTime, observed_body, linkEndTimes, linkEndStates );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                observation, observation2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( observationTime, linkEndTimes[ 0 ], std::numeric_limits< double >::epsilon( ) );


    computedRotationToBodyFixedFrame =
                Eigen::Quaterniond(
                                Eigen::AngleAxisd( -observation( 0 ), Eigen::Vector3d::UnitZ( ) ) *
                                Eigen::AngleAxisd( -observation( 1 ), Eigen::Vector3d::UnitX( ) ) *
                                Eigen::AngleAxisd( -observation( 2 ), Eigen::Vector3d::UnitZ( ) ) );
    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( currentRotationToBodyFixedFrame.toRotationMatrix( )( i, j ) -
                    computedRotationToBodyFixedFrame.toRotationMatrix( )( i, j ) ),
                    5.0 * std::numeric_limits< double >::epsilon( ) );
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}

