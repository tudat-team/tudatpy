/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"

#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/astro/relativity/metric.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createMetric.h"
#include "tudat/simulation/environment_setup/createRelativisticTimeConverter.h"

#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/propagation_setup/propagationPrintSettings.h"

namespace tudat
{

namespace unit_tests
{

using namespace simulation_setup;
using namespace numerical_integrators;
using namespace propagators;
using namespace spice_interface;
using namespace basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_relativistic_time_propagation )


BOOST_AUTO_TEST_CASE( testCombinedProperTimeAndStateDynamics2 )
{
    loadStandardSpiceKernels( );

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 2.0E7;
    const double maximumTimeStep = 3600.0;
    const double buffer = 5.0 * maximumTimeStep;

    const std::vector< std::string > bodyNames{ "Sun", "Earth", "Moon", "Mars", "Venus" };
    auto bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->bodyDeformationSettings.clear( );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );

    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStations;
    groundStations[ std::make_pair( "Earth", "Graz" ) ] =
            ( Eigen::Vector3d( ) << 4194511.7, 1162789.7, 4647362.5 ).finished( );

    SystemOfBodies bodiesDirect = createSystemOfBodies( bodySettings );
    SystemOfBodies bodiesMulti = createSystemOfBodies( bodySettings );

    createGroundStations( bodiesDirect, groundStations );
    createGroundStations( bodiesMulti, groundStations );

    for( const auto& bodyName : bodyNames )
    {
        if( bodiesDirect.doesBodyExist( bodyName ) )
        {
            bodiesDirect.getBody( bodyName )->setStateFromEphemeris( initialEphemerisTime );
        }
        if( bodiesMulti.doesBodyExist( bodyName ) )
        {
            bodiesMulti.getBody( bodyName )->setStateFromEphemeris( initialEphemerisTime );
        }
    }

    bodiesDirect.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEphemerisTime );
    bodiesMulti.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEphemerisTime );

    auto terminationSettings = std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime );
    auto integratorSettingsDirect = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, 10.0 );
    auto integratorSettingsMulti = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, 10.0 );

    const std::vector< std::string > perturbingBodies{ "Sun", "Moon", "Mars", "Venus" };

    // Direct metric: only Earth at second order.
    auto metricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
            perturbingBodies,
            std::vector< std::string >{ "Earth" },
            std::map< std::string, std::pair< int, int > >( ),
            std::vector< std::string >( ),
            std::make_shared< relativity::PPNParameterSet >( 1.0, 1.0 ),
            false );

    baseMetric = createSpaceTimeMetric( metricSettings, bodiesDirect );
    evaluatedMetricObjects.clear( );
    evaluatedMetricObjects[ std::make_pair( "Earth", "Graz" ) ] = baseMetric->Clone( );

    auto propagationPrintSettings = std::make_shared< PropagationPrintSettings >(
            true,
            true,
            maximumTimeStep,
            0.0,
            true,
            true,
            true,
            true,
            true,
            true );
    auto directOutputSettings = std::make_shared< SingleArcPropagatorProcessingSettings >(
            true,
            true,
            1,
            TUDAT_NAN
            //propagationPrintSettings 
            );

    auto directFromMetricSettings = std::make_shared< DirectRelativisticTimePropagatorSettings< double, double > >(
            std::make_pair( "Earth", "Graz" ),
            initialEphemerisTime,
            integratorSettingsDirect,
            terminationSettings,
            []( const double inputTime ){ return inputTime; },
            1.0,
            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
            directOutputSettings 
        );

    SingleArcDynamicsSimulator< double > directDynamicsSimulator(
            bodiesDirect,
            directFromMetricSettings,
            true );

    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRelativisticState =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 1 );

    auto firstOrderTimeSettings =
            std::make_shared< SecondOrderBodyCenteredRelativisticTimeConverterSettings< double, double > >(
                    "Earth",
                    perturbingBodies,
                    initialEphemerisTime,
                    integratorSettingsMulti,
                    terminationSettings );
    firstOrderTimeSettings->getOutputSettings( )->setIntegratedResult( true );

    auto bodyToTopoSettings =
            std::make_shared< BodycenteredToTopocentricTimePropagatorSettings< double, double > >(
                    std::make_pair( "Earth", "Graz" ),
                    false,
                    0,
                    false,
                    perturbingBodies,
                    initialRelativisticState,
                    initialEphemerisTime,
                    integratorSettingsMulti,
                    terminationSettings );
    bodyToTopoSettings->getOutputSettings( )->setIntegratedResult( true );

    SingleArcDynamicsSimulator< double > barycentricDynamicsSimulator(
            bodiesMulti,
            firstOrderTimeSettings,
            true );

    SingleArcDynamicsSimulator< double > topocentricDynamicsSimulator(
            bodiesMulti,
            bodyToTopoSettings,
            true );

    auto timeEphemerisDirect = bodiesDirect.getBody( "Earth" )->getTimeScaleConverter( );
    auto timeEphemerisCombined = bodiesMulti.getBody( "Earth" )->getTimeScaleConverter( );

    BOOST_REQUIRE_NE( timeEphemerisDirect, nullptr );
    BOOST_REQUIRE_NE( timeEphemerisCombined, nullptr );

    const auto directFunction = timeEphemerisDirect->getTimeDifferenceFunction(
            barycentric_coordinate_time_scale, local_proper_time_scale, "Graz" );
    const auto combinedFunction = timeEphemerisCombined->getTimeDifferenceFunction(
            barycentric_coordinate_time_scale, local_proper_time_scale, "Graz" );

    double currentTime = initialEphemerisTime + 5000.0;
    const double testTimeStep = 1.0E5;
    while( currentTime < finalEphemerisTime - 5000.0 )
    {
        BOOST_CHECK_SMALL( directFunction( currentTime ) - combinedFunction( currentTime ), 1.0E-1 );
        currentTime += testTimeStep;
    }
    BOOST_CHECK_SMALL( directFunction( finalEphemerisTime ) - combinedFunction( finalEphemerisTime ), 1.0E-3 );
}

// Keep one no-op case so the Boost module remains runnable while the
// full propagation test above is disabled.
BOOST_AUTO_TEST_CASE( testRelativisticTimePropagationPlaceholder )
{
    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
