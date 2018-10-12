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

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/Relativity/relativisticLightTimeCorrection.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/firstOrderRelativisticLightTimeCorrection.h"

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

BOOST_AUTO_TEST_SUITE_END( )

}

}
