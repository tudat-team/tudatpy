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


BOOST_AUTO_TEST_SUITE( test_angular_position_model )


BOOST_AUTO_TEST_CASE( testAngularPositionModel )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );

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
    linkEnds[ transmitter ] = std::make_pair( "Earth" , ""  );
    linkEnds[ receiver ] = std::make_pair( "Mars" , ""  );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                lightTimePerturbingBodies ) );

    // Create observation settings
    std::shared_ptr< ObservationModelSettings > observableSettings = std::make_shared< ObservationModelSettings >
            ( angular_position, linkEnds, lightTimeCorrectionSettings,
              std::make_shared< ConstantObservationBiasSettings >(
                  ( Eigen::Vector2d( ) << 3.2E-9, -1.5E-8 ).finished( ), true ) );

    // Create observation model.
    std::shared_ptr< ObservationModel< 2, double, double > > observationModel =
           ObservationModelCreator< 2, double, double >::createObservationModel(
                observableSettings, bodies );
    std::shared_ptr< ObservationBias< 2 > > observationBias = observationModel->getObservationBiasCalculator( );


    // Compute observation separately with two functions.
    double receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    Eigen::Vector2d observationFromReceptionTime = observationModel->computeObservations(
                receiverObservationTime, receiver );
    Eigen::Vector2d observationFromReceptionTime2 = observationModel->computeObservationsWithLinkEndData(
                receiverObservationTime, receiver, linkEndTimes, linkEndStates );
    BOOST_CHECK_EQUAL( linkEndTimes.size( ), 2 );
    BOOST_CHECK_EQUAL( linkEndStates.size( ), 2 );

    // Manually create and compute light time corrections
    std::shared_ptr< LightTimeCorrection > lightTimeCorrectionCalculator =
            createLightTimeCorrections(
                lightTimeCorrectionSettings.at( 0 ), bodies, linkEnds[ transmitter ], linkEnds[ receiver ] );
    double lightTimeCorrection = lightTimeCorrectionCalculator->calculateLightTimeCorrection(
                linkEndStates.at( 0 ), linkEndStates.at( 1 ), linkEndTimes.at( 0 ), linkEndTimes.at( 1 ) );

    // Check equality of computed observations.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                observationFromReceptionTime, observationFromReceptionTime2, std::numeric_limits< double >::epsilon( ) );

    // Check consistency of link end states and time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 0 ), bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 0 ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 1 ), bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 1 ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Check that reception time is kept fixed.
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( receiverObservationTime ),
                                linkEndTimes[ 1 ], std::numeric_limits< double >::epsilon( ) );

    // Manually compute light time
    Eigen::Vector3d positionDifference = ( linkEndStates[ 0 ] - linkEndStates[ 1 ] ).segment( 0, 3 );
    BOOST_CHECK_CLOSE_FRACTION(
                positionDifference.norm( ) / physical_constants::SPEED_OF_LIGHT + lightTimeCorrection,
                linkEndTimes[ 1 ] - linkEndTimes[ 0 ],
                std::numeric_limits< double >::epsilon( ) * 1000.0 );
                        // Poor tolerance due to rounding errors when subtracting times

    // Check computed right ascension/declination from link end states
    Eigen::Vector3d sphericalRelativeCoordinates = coordinate_conversions::convertCartesianToSpherical(
                positionDifference );
    BOOST_CHECK_CLOSE_FRACTION(
                sphericalRelativeCoordinates.z( ) + observationBias->getObservationBias(
                    std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ),
                observationFromReceptionTime( 0 ),
                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
                mathematical_constants::PI / 2.0 - sphericalRelativeCoordinates.y( ) +
                observationBias->getObservationBias(
                    std::vector< double >( ), std::vector< Eigen::Vector6d>( ) ).y( ),
                observationFromReceptionTime( 1 ),
                std::numeric_limits< double >::epsilon( ) );

    // Compute transmission time from light time.
    double transmitterObservationTime = receiverObservationTime - ( linkEndTimes[ 1 ] - linkEndTimes[ 0 ] );

    // Compare computed against returned transmission time.
    BOOST_CHECK_CLOSE_FRACTION(
                static_cast< double >( transmitterObservationTime ), linkEndTimes[ 0 ],
            std::numeric_limits< double >::epsilon( ) );

    // Recompute observable while fixing transmission time.
    std::vector< double > linkEndTimes2;
    std::vector< Eigen::Vector6d > linkEndStates2;
    Eigen::Vector2d observationFromTransmissionTime = observationModel->computeObservations(
                transmitterObservationTime, transmitter );
    Eigen::Vector2d observationFromTransmissionTime2 = observationModel->computeObservationsWithLinkEndData(
                transmitterObservationTime, transmitter, linkEndTimes2, linkEndStates2 );

    // Compare results against those obtained when keeping reception fixed.
    for( unsigned int i = 0; i < 2; i++ )
    {
        BOOST_CHECK_SMALL( observationFromTransmissionTime( i ) - observationFromTransmissionTime2( i ), 1.0E-15 );
        BOOST_CHECK_SMALL( observationFromTransmissionTime( i ) - observationFromReceptionTime( i ), 1.0E-15 );
        BOOST_CHECK_CLOSE_FRACTION( linkEndTimes2.at( i ), linkEndTimes2.at( i ), 1.0E-15 );
    }


    // Check consistency of link end states and time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates2.at( 0 ), bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( linkEndTimes2.at( 0 ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates2.at( 1 ), bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris( linkEndTimes2.at( 1 ) ),
                std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

