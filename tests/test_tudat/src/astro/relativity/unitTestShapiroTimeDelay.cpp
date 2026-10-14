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

BOOST_AUTO_TEST_CASE( testShapiroTimePartialsAndCache )
{
    Eigen::Vector6d transmitterState;
    transmitterState << -513.0E3, 120.0E3, 6378.0E3, -100.0, 300.0, 20.0;
    Eigen::Vector6d receiverState;
    receiverState << 324.0E3, 1434.0E3, 26600.0E3, 500.0, -200.0, 40.0;
    const double transmissionTime = 1000.0;
    const double receptionTime = 1020.0;
    const double timePerturbation = 0.01;
    const std::vector< std::string > bodyNames = { "Transmitter", "Receiver", "ThirdBody" };

    for( const bool useBending : { false, true } )
    {
        for( const bool movingBodies : { false, true } )
        {
            std::vector< std::function< Eigen::Vector6d( const double ) > > stateFunctions;
            std::vector< std::function< double( ) > > gravitationalParameterFunctions;
            for( unsigned int i = 0; i < bodyNames.size( ); i++ )
            {
                stateFunctions.push_back( [ = ]( const double time ) {
                    const double elapsedTime = time - transmissionTime;
                    const double motionScale = movingBodies ? static_cast< double >( i + 1 ) : 0.0;
                    const Eigen::Vector3d velocity = motionScale * Eigen::Vector3d( 30000.0, -15000.0, 8000.0 );
                    const Eigen::Vector3d acceleration = motionScale * Eigen::Vector3d( 0.2, -0.1, 0.3 );
                    Eigen::Vector6d state;
                    state.head< 3 >( ) = Eigen::Vector3d( -80.0E3, 120.0E3, 1.0E3 ) * static_cast< double >( i + 1 ) +
                            elapsedTime * velocity + 0.5 * elapsedTime * elapsedTime * acceleration;
                    state.tail< 3 >( ) = velocity + elapsedTime * acceleration;
                    return state;
                } );
                // Amplify the bending effect so both gradient implementations are exercised meaningfully.
                gravitationalParameterFunctions.push_back( [ = ]( ) { return 398600.44189E15 / ( i + 1 ); } );
            }

            // Test each evaluation epoch separately, then the sum of all three bodies.
            for( unsigned int bodyCase = 0; bodyCase <= bodyNames.size( ); bodyCase++ )
            {
                const bool allBodies = bodyCase == bodyNames.size( );
                const std::vector< std::string > selectedNames =
                        allBodies ? bodyNames : std::vector< std::string >{ bodyNames.at( bodyCase ) };
                const std::vector< std::function< Eigen::Vector6d( const double ) > > selectedStates = allBodies
                        ? stateFunctions
                        : std::vector< std::function< Eigen::Vector6d( const double ) > >{ stateFunctions.at( bodyCase ) };
                const std::vector< std::function< double( ) > > selectedParameters = allBodies
                        ? gravitationalParameterFunctions
                        : std::vector< std::function< double( ) > >{ gravitationalParameterFunctions.at( bodyCase ) };
                FirstOrderLightTimeCorrectionCalculator correction(
                        selectedStates, selectedParameters, selectedNames, "Transmitter", "Receiver", []( ) { return 1.1; }, useBending );

                const double nominalCorrection =
                        correction.calculateLightTimeCorrection( transmitterState, receiverState, transmissionTime, receptionTime );
                std::vector< double > nominalComponents;
                for( unsigned int i = 0; i < selectedNames.size( ); i++ )
                {
                    nominalComponents.push_back( correction.getCurrentLightTimeCorrectionComponent( i ) );
                }

                for( const LinkEndType linkEnd : { transmitter, receiver } )
                {
                    BOOST_TEST_CONTEXT( "bending=" << useBending << ", moving=" << movingBodies << ", body case=" << bodyCase
                                                   << ", link end=" << linkEnd )
                    {
                        const double analyticalPartial = correction.calculateLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
                                transmitterState, receiverState, transmissionTime, receptionTime, linkEnd, nullptr );

                        // Single-body helpers must reproduce the summed partials exactly (same operation order).
                        const Eigen::Vector3d analyticalGradient =
                                correction.calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
                                        transmitterState, receiverState, transmissionTime, receptionTime, linkEnd, nullptr );
                        std::vector< double > singleBodyTimePartials;
                        double summedTimePartial = 0.0;
                        Eigen::Vector3d summedGradient = Eigen::Vector3d::Zero( );
                        for( unsigned int i = 0; i < selectedNames.size( ); i++ )
                        {
                            singleBodyTimePartials.push_back(
                                    correction.calculateSingleBodyLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
                                            transmitterState, receiverState, transmissionTime, receptionTime, linkEnd, i ) );
                            summedTimePartial += singleBodyTimePartials.at( i );
                            summedGradient += correction.calculateSingleBodyLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
                                    transmitterState, receiverState, transmissionTime, receptionTime, linkEnd, i );
                        }
                        BOOST_CHECK_EQUAL( summedTimePartial, analyticalPartial );
                        for( unsigned int j = 0; j < 3; j++ )
                        {
                            BOOST_CHECK_EQUAL( summedGradient( j ), analyticalGradient( j ) );
                        }
                        BOOST_CHECK_EQUAL( correction.getCurrentTotalLightTimeCorrection( ), nominalCorrection );
                        for( unsigned int i = 0; i < selectedNames.size( ); i++ )
                        {
                            BOOST_CHECK_EQUAL( correction.getCurrentLightTimeCorrectionComponent( i ), nominalComponents.at( i ) );
                        }

                        // Only the explicit time dependence is differentiated: both endpoint states remain fixed.
                        const double transmitterStep = ( linkEnd == transmitter ) ? timePerturbation : 0.0;
                        const double receiverStep = ( linkEnd == receiver ) ? timePerturbation : 0.0;
                        const double upCorrection = correction.calculateLightTimeCorrection(
                                transmitterState, receiverState, transmissionTime + transmitterStep, receptionTime + receiverStep );
                        std::vector< double > upComponents;
                        for( unsigned int i = 0; i < selectedNames.size( ); i++ )
                        {
                            upComponents.push_back( correction.getCurrentLightTimeCorrectionComponent( i ) );
                        }
                        const double downCorrection = correction.calculateLightTimeCorrection(
                                transmitterState, receiverState, transmissionTime - transmitterStep, receptionTime - receiverStep );

                        // Compare each single-body time partial with the numerical derivative of its own component.
                        for( unsigned int i = 0; i < selectedNames.size( ); i++ )
                        {
                            const double numericalSingleBodyPartial =
                                    ( upComponents.at( i ) - correction.getCurrentLightTimeCorrectionComponent( i ) ) /
                                    ( 2.0 * timePerturbation );
                            if( numericalSingleBodyPartial == 0.0 )
                            {
                                BOOST_CHECK_EQUAL( singleBodyTimePartials.at( i ), 0.0 );
                            }
                            else
                            {
                                BOOST_CHECK_CLOSE_FRACTION( singleBodyTimePartials.at( i ), numericalSingleBodyPartial, 1.0E-7 );
                            }
                        }
                        const double numericalPartial = ( upCorrection - downCorrection ) / ( 2.0 * timePerturbation );
                        if( !movingBodies || ( bodyCase == 0 && linkEnd == receiver ) || ( bodyCase == 1 && linkEnd == transmitter ) )
                        {
                            BOOST_CHECK_EQUAL( analyticalPartial, 0.0 );
                            BOOST_CHECK_EQUAL( numericalPartial, 0.0 );
                        }
                        else
                        {
                            BOOST_CHECK_CLOSE_FRACTION( analyticalPartial, numericalPartial, 1.0E-7 );
                        }
                        correction.calculateLightTimeCorrection( transmitterState, receiverState, transmissionTime, receptionTime );
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
