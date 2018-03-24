/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/ObservationModels/angularPositionObservationModel.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"

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
    std::map< std::string, boost::shared_ptr< BodySettings > > defaultBodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( defaultBodySettings );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ observed_body ] = std::make_pair( "Earth" , ""  );


    // Create observation settings
    boost::shared_ptr< ObservationSettings > observableSettings = boost::make_shared< ObservationSettings >
            ( position_observable, std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
              boost::make_shared< ConstantObservationBiasSettings >(
                  ( Eigen::Vector3d( ) << 543.2454, -34.244, 3431.24345 ).finished( ), true ) );

    // Create observation model.
    boost::shared_ptr< ObservationModel< 3, double, double > > observationModel =
           ObservationModelCreator< 3, double, double >::createObservationModel(
                linkEnds, observableSettings, bodyMap );
    boost::shared_ptr< ObservationBias< 3 > > observationBias = observationModel->getObservationBiasCalculator( );


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

    // Check link end state/time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( observationTime ),
                linkEndStates[ 0 ], std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( observationTime, linkEndTimes[ 0 ], std::numeric_limits< double >::epsilon( ) );

    // Check biased observable
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( observationTime ).segment( 0, 3 ) +
                observationBias->getObservationBias(
                      std::vector< double >( ), std::vector< Eigen::Vector6d>( ) ) ),
                observation, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                observation, observation2, std::numeric_limits< double >::epsilon( ) );

    // Check ideal observable
    observation = observationModel->computeIdealObservations(
                observationTime, observed_body );
    observation2 = observationModel->computeIdealObservationsWithLinkEndData(
                observationTime, observed_body, linkEndTimes, linkEndStates );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( observationTime ),
                linkEndStates[ 0 ], std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( observationTime ).segment( 0, 3 ),
                observation, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                observation, observation2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( observationTime, linkEndTimes[ 0 ], std::numeric_limits< double >::epsilon( ) );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}

