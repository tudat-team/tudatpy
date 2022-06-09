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
#include "tudat/astro/observation_models/relativeAngularPositionObservationModel.h"
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


BOOST_AUTO_TEST_SUITE( test_relative_angular_position_model )


BOOST_AUTO_TEST_CASE( testRelativeAngularPositionModel )
{
//    spice_interface::loadStandardSpiceKernels( );

        spice_interface::loadStandardSpiceKernels( { paths::getSpiceKernelPath( ) + "/de430_mar097_small.bsp" } );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Phobos" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.addSettings( "Phobos" );
    bodySettings.at( "Phobos" )->ephemerisSettings = getDefaultEphemerisSettings( "Phobos" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ receiver ] = std::make_pair( "Earth" , ""  );
    linkEnds[ transmitter ] = std::make_pair( "Mars" , ""  );
    linkEnds[ transmitter2 ] = std::make_pair( "Phobos" , ""  );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                lightTimePerturbingBodies ) );

    // Create observation settings
    std::shared_ptr< ObservationModelSettings > observableSettings = std::make_shared< ObservationModelSettings >
            ( relative_angular_position, linkEnds, lightTimeCorrectionSettings,
              std::make_shared< ConstantObservationBiasSettings >( ( Eigen::Vector2d( ) << 3.2E-9, -1.5E-8 ).finished( ), true ) );

    // Create observation model.
    std::shared_ptr< ObservationModel< 2, double, double > > observationModel =
           ObservationModelCreator< 2, double, double >::createObservationModel( observableSettings, bodies );
    std::shared_ptr< ObservationBias< 2 > > observationBias = observationModel->getObservationBiasCalculator( );


    // Compute observation separately with two functions.
    double receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    Eigen::Vector2d observationFromReceptionTime = observationModel->computeObservations( receiverObservationTime, receiver );
    Eigen::Vector2d observationFromReceptionTime2 = observationModel->computeObservationsWithLinkEndData(
                receiverObservationTime, receiver, linkEndTimes, linkEndStates );
    BOOST_CHECK_EQUAL( linkEndTimes.size( ), 3 );
    BOOST_CHECK_EQUAL( linkEndStates.size( ), 3 );

    // Manually create and compute light time corrections
    std::shared_ptr< LightTimeCorrection > lightTimeCorrectionCalculatorFirstTransmitter =
            createLightTimeCorrections( lightTimeCorrectionSettings.at( 0 ), bodies, linkEnds[ transmitter ], linkEnds[ receiver ] );
    double lightTimeCorrectionFirstTransmitter = lightTimeCorrectionCalculatorFirstTransmitter->calculateLightTimeCorrection(
                linkEndStates.at( 0 ), linkEndStates.at( 2 ), linkEndTimes.at( 0 ), linkEndTimes.at( 2 ) );

    std::shared_ptr< LightTimeCorrection > lightTimeCorrectionCalculatorSecondTransmitter =
            createLightTimeCorrections( lightTimeCorrectionSettings.at( 0 ), bodies, linkEnds[ transmitter2 ], linkEnds[ receiver ] );
    double lightTimeCorrectionSecondTransmitter = lightTimeCorrectionCalculatorSecondTransmitter->calculateLightTimeCorrection(
            linkEndStates.at( 1 ), linkEndStates.at( 2 ), linkEndTimes.at( 1 ), linkEndTimes.at( 2 ) );

    // Check equality of computed observations.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                observationFromReceptionTime, observationFromReceptionTime2, std::numeric_limits< double >::epsilon( ) );

    // Check consistency of link end states and time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 0 ), bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 0 ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 1 ), bodies.at( "Phobos" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 1 ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            linkEndStates.at( 2 ), bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 2 ) ),
            std::numeric_limits< double >::epsilon( ) );

    // Check that reception time is kept fixed.
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( receiverObservationTime ),
                                linkEndTimes[ 2 ], std::numeric_limits< double >::epsilon( ) );

    // Manually compute light time
    Eigen::Vector3d positionDifferenceFirstTransmitter = ( linkEndStates[ 0 ] - linkEndStates[ 2 ] ).segment( 0, 3 );
    BOOST_CHECK_CLOSE_FRACTION(
                positionDifferenceFirstTransmitter.norm( ) / physical_constants::SPEED_OF_LIGHT + lightTimeCorrectionFirstTransmitter,
                linkEndTimes[ 2 ] - linkEndTimes[ 0 ],
                std::numeric_limits< double >::epsilon( ) * 1000.0 );
                        // Poor tolerance due to rounding errors when subtracting times

    Eigen::Vector3d positionDifferenceSecondTransmitter = ( linkEndStates[ 1 ] - linkEndStates[ 2 ] ).segment( 0, 3 );
    BOOST_CHECK_CLOSE_FRACTION(
            positionDifferenceSecondTransmitter.norm( ) / physical_constants::SPEED_OF_LIGHT + lightTimeCorrectionSecondTransmitter,
            linkEndTimes[ 2 ] - linkEndTimes[ 1 ],
            std::numeric_limits< double >::epsilon( ) * 1000.0 );
                        // Poor tolerance due to rounding errors when subtracting times

    // Check computed relative right ascension/declination from link end states
    Eigen::Vector3d sphericalCoordinatesFirstTransmitter = coordinate_conversions::convertCartesianToSpherical( positionDifferenceFirstTransmitter );
    double rightAscensionFirstTransmitter = sphericalCoordinatesFirstTransmitter.z( );
    double declinationFirstTransmitter = mathematical_constants::PI / 2.0 - sphericalCoordinatesFirstTransmitter.y( );

    Eigen::Vector3d sphericalCoordinatesSecondTransmitter = coordinate_conversions::convertCartesianToSpherical( positionDifferenceSecondTransmitter );
    double rightAscensionSecondTransmitter = sphericalCoordinatesSecondTransmitter.z( );
    double declinationSecondTransmitter = mathematical_constants::PI / 2.0 - sphericalCoordinatesSecondTransmitter.y( );

    BOOST_CHECK_CLOSE_FRACTION(
            ( rightAscensionSecondTransmitter -  rightAscensionFirstTransmitter ) + observationBias->getObservationBias(
                    std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ), observationFromReceptionTime( 0 ), std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
            ( declinationSecondTransmitter - declinationFirstTransmitter ) + observationBias->getObservationBias(
                    std::vector< double >( ), std::vector< Eigen::Vector6d>( ) ).y( ), observationFromReceptionTime( 1 ), std::numeric_limits< double >::epsilon( ) );

    // Compute transmission time from light time.
    double firstTransmitterObservationTime = receiverObservationTime - ( linkEndTimes[ 2 ] - linkEndTimes[ 0 ] );
    double secondTransmitterObservationTime = receiverObservationTime - ( linkEndTimes[ 2 ] - linkEndTimes[ 1 ] );

    // Compare computed against returned transmission time.
    BOOST_CHECK_CLOSE_FRACTION(
                static_cast< double >( firstTransmitterObservationTime ), linkEndTimes[ 0 ], std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
            static_cast< double >( secondTransmitterObservationTime ), linkEndTimes[ 1 ], std::numeric_limits< double >::epsilon( ) );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}

