/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <algorithm>
#include <limits>
#include <stdexcept>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/relativity/relativisticLightTimeCorrection.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/observation_models/corrections/firstOrderRelativisticCorrection.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::relativity;
using namespace tudat::ephemerides;

BOOST_AUTO_TEST_SUITE( test_shapiro_delay )

BOOST_AUTO_TEST_CASE( testShapiroDelay )
{
    Eigen::Vector6d groundStationState;
    groundStationState << 0.0, 0.0, 6378.0, 0.0, 0.0, 0.0;
    Eigen::Vector6d satelliteState;
    satelliteState << 0.0, 0.0, 26600.0, 0.0, 0.0, 0.0;
    Eigen::Vector6d centralBodyPosition = Eigen::Vector6d::Zero( );

    std::shared_ptr< ConstantEphemeris > ephemeris = std::make_shared< ConstantEphemeris >( [ & ]( ) { return centralBodyPosition; } );

    double earthGravitationalParameter = 398600.44189E9;

    // Calculate without bending
    double directCalculationNoBending = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                               groundStationState.segment( 0, 3 ),
                                                                                               satelliteState.segment( 0, 3 ),
                                                                                               centralBodyPosition.segment( 0, 3 ),
                                                                                               1.0,
                                                                                               false );

    std::vector< std::function< Eigen::Vector6d( const double ) > > perturbingBodyStateFunctions;
    std::vector< std::function< double( ) > > perturbingBodyGravitationalParameterFunctions;

    perturbingBodyStateFunctions.push_back( std::bind( &Ephemeris::getCartesianState, ephemeris, std::placeholders::_1 ) );
    perturbingBodyGravitationalParameterFunctions.push_back( [ & ]( ) { return earthGravitationalParameter; } );

    // Calculate with bending through calculator (default)
    FirstOrderLightTimeCorrectionCalculator correctionCalculatorNoBending(
            perturbingBodyStateFunctions,
            perturbingBodyGravitationalParameterFunctions,
            std::vector< std::string >{ "Earth" },
            "Satellite",
            "Earth",
            []( ) { return 1.0; },
            false );

    double classInterfaceCalculationNoBending =
            correctionCalculatorNoBending.calculateLightTimeCorrection( groundStationState, satelliteState, 0.0, 0.0 );

    // Living reviews in relativity, GPS
    double expectedResult = 6.3E-3;

    BOOST_CHECK_CLOSE_FRACTION( 0.5 * classInterfaceCalculationNoBending * physical_constants::SPEED_OF_LIGHT, expectedResult, 6.0E-2 );
    BOOST_CHECK_CLOSE_FRACTION( 0.5 * directCalculationNoBending * physical_constants::SPEED_OF_LIGHT, expectedResult, 6.0E-2 );

    // Now calculate with bending enabled
    double directCalculationWithBending = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                                 groundStationState.segment( 0, 3 ),
                                                                                                 satelliteState.segment( 0, 3 ),
                                                                                                 centralBodyPosition.segment( 0, 3 ),
                                                                                                 1.0,
                                                                                                 true );

    FirstOrderLightTimeCorrectionCalculator correctionCalculatorWithBending(
            perturbingBodyStateFunctions,
            perturbingBodyGravitationalParameterFunctions,
            std::vector< std::string >{ "Earth" },
            "Satellite",
            "Earth",
            []( ) { return 1.0; },
            true );

    double classInterfaceCalculationWithBending =
            correctionCalculatorWithBending.calculateLightTimeCorrection( groundStationState, satelliteState, 0.0, 0.0 );

    // Expected result for bending case
    double expectedResultBending = 6.3E-3;

    BOOST_CHECK_CLOSE_FRACTION(
            0.5 * classInterfaceCalculationWithBending * physical_constants::SPEED_OF_LIGHT, expectedResultBending, 6.0E-2 );
    BOOST_CHECK_CLOSE_FRACTION( 0.5 * directCalculationWithBending * physical_constants::SPEED_OF_LIGHT, expectedResultBending, 6.0E-2 );

    // Check that direct and class interface calculations are consistent with each other for both modes
    BOOST_CHECK_CLOSE_FRACTION( directCalculationNoBending, classInterfaceCalculationNoBending, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( directCalculationWithBending, classInterfaceCalculationWithBending, 1.0E-14 );

    // Chenckk bending value
    double bendingEffectMeter =
            ( classInterfaceCalculationWithBending - classInterfaceCalculationNoBending ) * physical_constants::SPEED_OF_LIGHT;

    // Expected value calculated manually on previous runs
    double expectedBendingEffectMeter = -4.69e-09;
    // Verify the value is as expected (using a small tolerance)
    BOOST_CHECK_CLOSE_FRACTION( bendingEffectMeter, expectedBendingEffectMeter, 1.0E-3 );
}
BOOST_AUTO_TEST_CASE( testShapiroDelayGradient )
{
    Eigen::Vector6d groundStationState;
    groundStationState << -513.0E3, 120E3, 6378.0E3, 0.0, 0.0, 0.0;
    Eigen::Vector6d satelliteState;
    satelliteState << 324.0E3, 1434.0E3, 26600.0E3, 0.0, 0.0, 0.0;
    Eigen::Vector6d centralBodyPosition;
    centralBodyPosition << -80.0E3, 120.0E3, 1.0E3, 0.0, 0.0, 0.0;

    double earthGravitationalParameter = 398600.44189E15;
    double positionPerturbation = 1.0E2;

    // Test with bending (default behavior)
    {
        double directCalculation = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                          groundStationState.segment( 0, 3 ),
                                                                                          satelliteState.segment( 0, 3 ),
                                                                                          centralBodyPosition.segment( 0, 3 ),
                                                                                          1.0,
                                                                                          true );

        Eigen::Matrix< double, 1, 3 > numericalLightTimePartialWrtReceiver = Eigen::Matrix< double, 1, 3 >::Zero( );
        Eigen::Matrix< double, 1, 3 > numericalLightTimePartialWrtTransmitter = Eigen::Matrix< double, 1, 3 >::Zero( );

        Eigen::Vector6d perturbedSatelliteState;
        perturbedSatelliteState.setZero( );

        Eigen::Vector6d perturbedStationState;
        perturbedStationState.setZero( );
        for( int i = 0; i < 3; i++ )
        {
            perturbedSatelliteState = satelliteState;
            perturbedSatelliteState( i ) += positionPerturbation;
            double upperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                                 groundStationState.segment( 0, 3 ),
                                                                                                 perturbedSatelliteState.segment( 0, 3 ),
                                                                                                 centralBodyPosition.segment( 0, 3 ),
                                                                                                 1.0,
                                                                                                 true );

            perturbedSatelliteState = satelliteState;
            perturbedSatelliteState( i ) -= positionPerturbation;
            double downperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                                   groundStationState.segment( 0, 3 ),
                                                                                                   perturbedSatelliteState.segment( 0, 3 ),
                                                                                                   centralBodyPosition.segment( 0, 3 ),
                                                                                                   1.0,
                                                                                                   true );

            numericalLightTimePartialWrtReceiver( i ) = ( upperturbedLightTime - downperturbedLightTime ) / ( 2.0 * positionPerturbation );

            perturbedStationState = groundStationState;
            perturbedStationState( i ) += positionPerturbation;
            upperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                          perturbedStationState.segment( 0, 3 ),
                                                                                          satelliteState.segment( 0, 3 ),
                                                                                          centralBodyPosition.segment( 0, 3 ),
                                                                                          1.0,
                                                                                          true );

            perturbedStationState = groundStationState;
            perturbedStationState( i ) -= positionPerturbation;
            downperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                            perturbedStationState.segment( 0, 3 ),
                                                                                            satelliteState.segment( 0, 3 ),
                                                                                            centralBodyPosition.segment( 0, 3 ),
                                                                                            1.0,
                                                                                            true );

            numericalLightTimePartialWrtTransmitter( i ) =
                    ( upperturbedLightTime - downperturbedLightTime ) / ( 2.0 * positionPerturbation );
        }

        Eigen::Matrix< double, 1, 3 > analyticalLightTimePartialWrtReceiver =
                calculateFirstOrderCentralBodyLightTimeCorrectionGradient( earthGravitationalParameter,
                                                                           groundStationState.segment( 0, 3 ),
                                                                           satelliteState.segment( 0, 3 ),
                                                                           centralBodyPosition.segment( 0, 3 ),
                                                                           true,
                                                                           1.0,
                                                                           true );

        Eigen::Matrix< double, 1, 3 > analyticalLightTimePartialWrtTransmitter =
                calculateFirstOrderCentralBodyLightTimeCorrectionGradient( earthGravitationalParameter,
                                                                           groundStationState.segment( 0, 3 ),
                                                                           satelliteState.segment( 0, 3 ),
                                                                           centralBodyPosition.segment( 0, 3 ),
                                                                           false,
                                                                           1.0,
                                                                           true );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalLightTimePartialWrtReceiver, numericalLightTimePartialWrtReceiver, 1.0E-7 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalLightTimePartialWrtTransmitter, numericalLightTimePartialWrtTransmitter, 1.0E-7 );
    }

    // Test without bending
    {
        double directCalculationNoBending = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                                   groundStationState.segment( 0, 3 ),
                                                                                                   satelliteState.segment( 0, 3 ),
                                                                                                   centralBodyPosition.segment( 0, 3 ),
                                                                                                   1.0,
                                                                                                   false );

        Eigen::Matrix< double, 1, 3 > numericalLightTimePartialWrtReceiverNoBending = Eigen::Matrix< double, 1, 3 >::Zero( );
        Eigen::Matrix< double, 1, 3 > numericalLightTimePartialWrtTransmitterNoBending = Eigen::Matrix< double, 1, 3 >::Zero( );

        Eigen::Vector6d perturbedSatelliteState;
        perturbedSatelliteState.setZero( );

        Eigen::Vector6d perturbedStationState;
        perturbedStationState.setZero( );
        for( int i = 0; i < 3; i++ )
        {
            perturbedSatelliteState = satelliteState;
            perturbedSatelliteState( i ) += positionPerturbation;
            double upperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                                 groundStationState.segment( 0, 3 ),
                                                                                                 perturbedSatelliteState.segment( 0, 3 ),
                                                                                                 centralBodyPosition.segment( 0, 3 ),
                                                                                                 1.0,
                                                                                                 false );

            perturbedSatelliteState = satelliteState;
            perturbedSatelliteState( i ) -= positionPerturbation;
            double downperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                                   groundStationState.segment( 0, 3 ),
                                                                                                   perturbedSatelliteState.segment( 0, 3 ),
                                                                                                   centralBodyPosition.segment( 0, 3 ),
                                                                                                   1.0,
                                                                                                   false );

            numericalLightTimePartialWrtReceiverNoBending( i ) =
                    ( upperturbedLightTime - downperturbedLightTime ) / ( 2.0 * positionPerturbation );

            perturbedStationState = groundStationState;
            perturbedStationState( i ) += positionPerturbation;
            upperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                          perturbedStationState.segment( 0, 3 ),
                                                                                          satelliteState.segment( 0, 3 ),
                                                                                          centralBodyPosition.segment( 0, 3 ),
                                                                                          1.0,
                                                                                          false );

            perturbedStationState = groundStationState;
            perturbedStationState( i ) -= positionPerturbation;
            downperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody( earthGravitationalParameter,
                                                                                            perturbedStationState.segment( 0, 3 ),
                                                                                            satelliteState.segment( 0, 3 ),
                                                                                            centralBodyPosition.segment( 0, 3 ),
                                                                                            1.0,
                                                                                            false );

            numericalLightTimePartialWrtTransmitterNoBending( i ) =
                    ( upperturbedLightTime - downperturbedLightTime ) / ( 2.0 * positionPerturbation );
        }

        Eigen::Matrix< double, 1, 3 > analyticalLightTimePartialWrtReceiverNoBending =
                calculateFirstOrderCentralBodyLightTimeCorrectionGradient( earthGravitationalParameter,
                                                                           groundStationState.segment( 0, 3 ),
                                                                           satelliteState.segment( 0, 3 ),
                                                                           centralBodyPosition.segment( 0, 3 ),
                                                                           true,
                                                                           1.0,
                                                                           false );

        Eigen::Matrix< double, 1, 3 > analyticalLightTimePartialWrtTransmitterNoBending =
                calculateFirstOrderCentralBodyLightTimeCorrectionGradient( earthGravitationalParameter,
                                                                           groundStationState.segment( 0, 3 ),
                                                                           satelliteState.segment( 0, 3 ),
                                                                           centralBodyPosition.segment( 0, 3 ),
                                                                           false,
                                                                           1.0,
                                                                           false );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                analyticalLightTimePartialWrtReceiverNoBending, numericalLightTimePartialWrtReceiverNoBending, 1.0E-7 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                analyticalLightTimePartialWrtTransmitterNoBending, numericalLightTimePartialWrtTransmitterNoBending, 1.0E-7 );
    }
}

namespace
{

void checkClosestApproachCase(
        const Eigen::Vector3d& bodyPositionAtTransmission,
        const Eigen::Vector3d& bodyVelocity,
        const Eigen::Vector3d& transmitterPosition,
        const Eigen::Vector3d& receiverPosition,
        const double transmissionTime,
        const double duration,
        const double expectedOffset,
        const std::string& bodyName = "Perturber",
        const std::string& transmitterName = "Transmitter",
        const std::string& receiverName = "Receiver" )
{
    const double gm = 1.266865341960128e17;
    const double expectedEpoch = transmissionTime + expectedOffset;
    Eigen::Vector6d transmitter = Eigen::Vector6d::Zero( );
    Eigen::Vector6d receiver = Eigen::Vector6d::Zero( );
    transmitter.head< 3 >( ) = transmitterPosition;
    receiver.head< 3 >( ) = receiverPosition;
    for( const bool bending : { false, true } )
    {
        double lastBodyEpoch = TUDAT_NAN;
        const auto stateFunction = [ & ]( const double epoch ) -> Eigen::Vector6d {
            lastBodyEpoch = epoch;
            Eigen::Vector6d state;
            state.head< 3 >( ) = bodyPositionAtTransmission + ( epoch - transmissionTime ) * bodyVelocity;
            state.tail< 3 >( ) = bodyVelocity;
            return state;
        };
        FirstOrderLightTimeCorrectionCalculator calculator(
                { stateFunction }, { [ & ]( ) { return gm; } }, { bodyName }, transmitterName, receiverName,
                []( ) { return 1.0; }, bending );
        const double calculated = calculator.calculateLightTimeCorrection(
                transmitter, receiver, transmissionTime, transmissionTime + duration );
        const double epochTolerance = std::max( 1.0e-9, 2.0 * std::numeric_limits< double >::epsilon( ) *
                                               std::max( std::abs( transmissionTime ), std::abs( transmissionTime + duration ) ) );
        BOOST_CHECK_SMALL( lastBodyEpoch - expectedEpoch, epochTolerance );
        const Eigen::Vector3d expectedBodyPosition =
                bodyPositionAtTransmission + ( expectedEpoch - transmissionTime ) * bodyVelocity;
        const double expected = calculateFirstOrderLightTimeCorrectionFromCentralBody(
                gm, transmitterPosition, receiverPosition, expectedBodyPosition, 1.0, bending );
        BOOST_CHECK_CLOSE_FRACTION( calculated, expected, 1.0e-12 );
        if( bodyVelocity.isZero( 0.0 ) )
        {
            // Changing an epoch must not change the stationary-body delay.
            BOOST_CHECK_EQUAL( calculated, expected );
        }
        if( ( receiverPosition - transmitterPosition ).isZero( 0.0 ) )
        {
            continue;  // The pre-existing position gradient is undefined for a zero spatial path.
        }
        // Check the evaluation epoch in the existing fixed-body-epoch gradient,
        // without claiming that it differentiates the selected epoch itself.
        const Eigen::Vector3d gradient = calculator.calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
                transmitter, receiver, transmissionTime, transmissionTime + duration, observation_models::receiver );
        BOOST_CHECK_SMALL( lastBodyEpoch - expectedEpoch, epochTolerance );
        const Eigen::Vector3d expectedGradient = calculateFirstOrderCentralBodyLightTimeCorrectionGradient(
                gm, transmitterPosition, receiverPosition, expectedBodyPosition, true, 1.0, bending );
        BOOST_CHECK_SMALL( ( gradient - expectedGradient ).norm( ), 1.0e-20 );
    }
}

}  // namespace

BOOST_AUTO_TEST_CASE( testClosestApproachStationaryInteriorAndNameIndependence )
{
    const Eigen::Vector3d transmitter = Eigen::Vector3d::Zero( );
    const Eigen::Vector3d receiver( 1.0e9, 0.0, 0.0 );
    const Eigen::Vector3d body( 4.0e8, 2.0e8, 0.0 );
    // Legacy constructor fractions would select transmit, receive or midpoint;
    // all three name combinations must now select the same interior epoch.
    for( const std::string name : { "Earth", "Jupiter", "Other" } )
    {
        checkClosestApproachCase( body, Eigen::Vector3d::Zero( ), transmitter, receiver,
                                  1000.0, 10.0, 4.0, name, "Earth", "Jupiter" );
    }
    checkClosestApproachCase( Eigen::Vector3d( 4.123456789e8, 2.0e8, 0.0 ), Eigen::Vector3d::Zero( ),
                              transmitter, receiver, 5.0e8, 10.0, 4.123456789 );
}

BOOST_AUTO_TEST_CASE( testClosestApproachMovingBodyAnalytic )
{
    const Eigen::Vector3d transmitter = Eigen::Vector3d::Zero( );
    const Eigen::Vector3d receiver( 1.0e9, 0.0, 0.0 );
    const Eigen::Vector3d body( 4.0e8, 2.0e8, 0.0 );
    const Eigen::Vector3d bodyVelocity( 1.0e6, -2.0e6, 3.0e5 );
    const Eigen::Vector3d relativeVelocity = receiver / 10.0 - bodyVelocity;
    const double analyticOffset = body.dot( relativeVelocity ) / relativeVelocity.squaredNorm( );
    checkClosestApproachCase( body, bodyVelocity, transmitter, receiver, 1000.0, 10.0, analyticOffset );
    BOOST_CHECK( std::abs( analyticOffset - 4.0 ) > 0.01 );
}

BOOST_AUTO_TEST_CASE( testClosestApproachClampedLegEndpoints )
{
    const Eigen::Vector3d transmitter = Eigen::Vector3d::Zero( );
    const Eigen::Vector3d receiver( 1.0e9, 0.0, 0.0 );
    checkClosestApproachCase( Eigen::Vector3d( -1.0e8, 2.0e8, 0.0 ), Eigen::Vector3d::Zero( ),
                              transmitter, receiver, 1000.0, 10.0, 0.0 );
    checkClosestApproachCase( Eigen::Vector3d( 1.2e9, 2.0e8, 0.0 ), Eigen::Vector3d::Zero( ),
                              transmitter, receiver, 1000.0, 10.0, 10.0 );
}

BOOST_AUTO_TEST_CASE( testClosestApproachJupiterLikeUplinkAndDownlink )
{
    const Eigen::Vector3d earth = Eigen::Vector3d::Zero( );
    const Eigen::Vector3d spacecraft( 7.0e11, 0.0, 0.0 );
    const Eigen::Vector3d jupiter( 7.004e11, 7.0e7, 0.0 );
    const Eigen::Vector3d velocity( 0.0, 1.3e4, 0.0 );
    // The gravitating body's name deliberately differs from the spacecraft end.
    checkClosestApproachCase( jupiter, velocity, earth, spacecraft, 5.0e8, 2500.0, 2500.0,
                              "Jupiter", "Earth", "JUNO" );
    checkClosestApproachCase( jupiter, velocity, spacecraft, earth, 5.0e8, 2500.0, 0.0,
                              "Jupiter", "JUNO", "Earth" );
}

BOOST_AUTO_TEST_CASE( testClosestApproachZeroDurationAndZeroSpatialPath )
{
    const Eigen::Vector3d transmitter = Eigen::Vector3d::Zero( );
    const Eigen::Vector3d receiver( 1.0e9, 0.0, 0.0 );
    const Eigen::Vector3d body( 4.0e8, 2.0e8, 0.0 );
    checkClosestApproachCase( body, Eigen::Vector3d( 1.0e3, 2.0e3, 0.0 ), transmitter, receiver,
                              1000.0, 0.0, 0.0 );
    checkClosestApproachCase( body, Eigen::Vector3d( 1.0e3, 2.0e3, 0.0 ), transmitter, transmitter,
                              1000.0, 10.0, 5.0 );
}

BOOST_AUTO_TEST_CASE( testClosestApproachInvalidInputAndConvergenceGuards )
{
    Eigen::Vector6d transmitter = Eigen::Vector6d::Zero( );
    Eigen::Vector6d receiver = Eigen::Vector6d::Zero( );
    receiver( 0 ) = 1000.0;
    const auto finiteState = []( const double ) -> Eigen::Vector6d {
        Eigen::Vector6d result = Eigen::Vector6d::Zero( );
        result( 1 ) = 100.0;
        return result;
    };
    FirstOrderLightTimeCorrectionCalculator valid( { finiteState }, { []( ) { return 1.0e14; } },
                                                  { "Body" }, "Transmitter", "Receiver" );
    BOOST_CHECK_THROW( valid.calculateLightTimeCorrection( transmitter, receiver, 1.0, 0.0 ), std::invalid_argument );
    BOOST_CHECK_THROW( valid.calculateLightTimeCorrection( transmitter, receiver, 0.0,
                                                         std::numeric_limits< double >::infinity( ) ), std::invalid_argument );
    const auto invalidState = []( const double ) -> Eigen::Vector6d {
        return Eigen::Vector6d::Constant( std::numeric_limits< double >::quiet_NaN( ) );
    };
    FirstOrderLightTimeCorrectionCalculator invalid( { invalidState }, { []( ) { return 1.0e14; } },
                                                    { "Body" }, "Transmitter", "Receiver" );
    BOOST_CHECK_THROW( invalid.calculateLightTimeCorrection( transmitter, receiver, 0.0, 10.0 ), std::runtime_error );
    const auto nonconvergentState = []( const double epoch ) -> Eigen::Vector6d {
        Eigen::Vector6d result = Eigen::Vector6d::Zero( );
        result( 0 ) = 100.0 * epoch + ( epoch < 5.0 ? 100.0 : -100.0 );
        result( 1 ) = 100.0;
        return result;
    };
    FirstOrderLightTimeCorrectionCalculator nonconvergent( { nonconvergentState }, { []( ) { return 1.0e14; } },
                                                          { "Body" }, "Transmitter", "Receiver" );
    BOOST_CHECK_THROW( nonconvergent.calculateLightTimeCorrection( transmitter, receiver, 0.0, 10.0 ), std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
