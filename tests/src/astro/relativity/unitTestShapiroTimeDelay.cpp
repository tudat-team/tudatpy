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

#include <boost/test/unit_test.hpp>

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
    satelliteState  <<  0.0, 0.0, 26600.0, 0.0, 0.0, 0.0;
    Eigen::Vector6d centralBodyPosition = Eigen::Vector6d::Zero( );

    std::shared_ptr< ConstantEphemeris > ephemeris = std::make_shared< ConstantEphemeris >(
                [ & ]( ){ return centralBodyPosition; } );

    double earthGravitationalParameter = 398600.44189E9;

    double directCalculation = calculateFirstOrderLightTimeCorrectionFromCentralBody(
                earthGravitationalParameter, groundStationState.segment( 0, 3 ),
                satelliteState.segment( 0, 3 ), centralBodyPosition.segment( 0, 3 ) );

    std::vector< std::function< Eigen::Vector6d( const double ) > > perturbingBodyStateFunctions;
    std::vector< std::function< double( ) > > perturbingBodyGravitationalParameterFunctions;

    perturbingBodyStateFunctions.push_back( std::bind( &Ephemeris::getCartesianState, ephemeris, std::placeholders::_1 ) );
    perturbingBodyGravitationalParameterFunctions.push_back( [ & ]( ){ return earthGravitationalParameter; } );

    FirstOrderLightTimeCorrectionCalculator correctionCalculator(
                perturbingBodyStateFunctions, perturbingBodyGravitationalParameterFunctions,
                std::vector< std::string >{ "Earth" }, "Satellite", "Earth" );

    double classInterfaceCalculation = correctionCalculator.calculateLightTimeCorrection(
                groundStationState, satelliteState, 0.0, 0.0 );

    // Living reviews in relativity, GPS.
    double expectedResult = 6.3E-3;

    BOOST_CHECK_CLOSE_FRACTION( 0.5 * classInterfaceCalculation * physical_constants::SPEED_OF_LIGHT, expectedResult, 6.0E-2 );
    BOOST_CHECK_CLOSE_FRACTION( 0.5 * directCalculation * physical_constants::SPEED_OF_LIGHT, expectedResult, 6.0E-2 );
}

BOOST_AUTO_TEST_CASE( testShapiroDelayGradient )
{
    Eigen::Vector6d groundStationState;
    groundStationState << -513.0E3, 120E3, 6378.0E3, 0.0, 0.0, 0.0;
    Eigen::Vector6d satelliteState;
    satelliteState  <<  324.0E3, 1434.0E3, 26600.0E3, 0.0, 0.0, 0.0;
    Eigen::Vector6d centralBodyPosition;
    centralBodyPosition  <<  -80.0E3, 120.0E3, 1.0E3, 0.0, 0.0, 0.0;

    double earthGravitationalParameter = 398600.44189E15;

    double directCalculation = calculateFirstOrderLightTimeCorrectionFromCentralBody(
        earthGravitationalParameter, groundStationState.segment( 0, 3 ),
        satelliteState.segment( 0, 3 ), centralBodyPosition.segment( 0, 3 ) );

    double positionPerturbation = 1.0E2;

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
        double upperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody(
            earthGravitationalParameter, groundStationState.segment( 0, 3 ),
            perturbedSatelliteState.segment( 0, 3 ), centralBodyPosition.segment( 0, 3 ) );
        
        perturbedSatelliteState = satelliteState;
        perturbedSatelliteState( i ) -= positionPerturbation;
        double downperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody(
            earthGravitationalParameter, groundStationState.segment( 0, 3 ),
            perturbedSatelliteState.segment( 0, 3 ), centralBodyPosition.segment( 0, 3 ) );

        numericalLightTimePartialWrtReceiver( i ) = ( upperturbedLightTime - downperturbedLightTime ) / ( 2.0 * positionPerturbation );

        perturbedStationState = groundStationState;
        perturbedStationState( i ) += positionPerturbation;
        upperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody(
            earthGravitationalParameter, perturbedStationState.segment( 0, 3 ),
            satelliteState.segment( 0, 3 ), centralBodyPosition.segment( 0, 3 ) );

        perturbedStationState = groundStationState;
        perturbedStationState( i ) -= positionPerturbation;
        downperturbedLightTime = calculateFirstOrderLightTimeCorrectionFromCentralBody(
            earthGravitationalParameter, perturbedStationState.segment( 0, 3 ),
            satelliteState.segment( 0, 3 ), centralBodyPosition.segment( 0, 3 ) );

        numericalLightTimePartialWrtTransmitter( i ) = ( upperturbedLightTime - downperturbedLightTime ) / ( 2.0 * positionPerturbation );
    }
    Eigen::Matrix< double, 1, 3 > analyticalLightTimePartialWrtReceiver = calculateFirstOrderCentralBodyLightTimeCorrectionGradient(
        earthGravitationalParameter, groundStationState.segment( 0, 3 ),
        satelliteState.segment( 0, 3 ), centralBodyPosition.segment( 0, 3 ), true );

    Eigen::Matrix< double, 1, 3 > analyticalLightTimePartialWrtTransmitter = calculateFirstOrderCentralBodyLightTimeCorrectionGradient(
        earthGravitationalParameter, groundStationState.segment( 0, 3 ),
        satelliteState.segment( 0, 3 ), centralBodyPosition.segment( 0, 3 ), false );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalLightTimePartialWrtReceiver, numericalLightTimePartialWrtReceiver, 1.0E-7 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalLightTimePartialWrtTransmitter, numericalLightTimePartialWrtTransmitter, 1.0E-7 );

}



BOOST_AUTO_TEST_SUITE_END( )

}

}
