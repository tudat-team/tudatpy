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
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <tuple>

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
    groundStations[ std::make_pair( "Earth", "Graz" ) ] = ( Eigen::Vector3d( ) << 4194511.7, 1162789.7, 4647362.5 ).finished( );

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
    auto metricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >( perturbingBodies,
                                                                                  std::vector< std::string >{ "Earth" },
                                                                                  std::map< std::string, std::pair< int, int > >( ),
                                                                                  std::vector< std::string >( ),
                                                                                  false );

    createBaseMetric( metricSettings, bodiesDirect );

    auto propagationPrintSettings =
            std::make_shared< PropagationPrintSettings >( true, true, maximumTimeStep, 0.0, true, true, true, true, true, true );
    auto directOutputSettings = std::make_shared< SingleArcPropagatorProcessingSettings >( true, true, 1, TUDAT_NAN
                                                                                           // propagationPrintSettings
    );

    auto directFromMetricSettings = std::make_shared< DirectRelativisticTimePropagatorSettings< double, double > >(
            std::make_pair( "Earth", "Graz" ),
            initialEphemerisTime,
            integratorSettingsDirect,
            terminationSettings,
            []( const double inputTime ) { return inputTime; },
            1.0,
            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
            directOutputSettings );

    SingleArcDynamicsSimulator< double > directDynamicsSimulator( bodiesDirect, directFromMetricSettings, true );

    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRelativisticState = Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 1 );

    auto firstOrderTimeSettings = std::make_shared< SecondOrderBodyCenteredRelativisticTimeConverterSettings< double, double > >(
            "Earth", perturbingBodies, initialEphemerisTime, integratorSettingsMulti, terminationSettings );
    firstOrderTimeSettings->getOutputSettings( )->setIntegratedResult( true );

    auto bodyToTopoSettings =
            std::make_shared< BodycenteredToTopocentricTimePropagatorSettings< double, double > >( std::make_pair( "Earth", "Graz" ),
                                                                                                   false,
                                                                                                   0,
                                                                                                   false,
                                                                                                   perturbingBodies,
                                                                                                   initialRelativisticState,
                                                                                                   initialEphemerisTime,
                                                                                                   integratorSettingsMulti,
                                                                                                   terminationSettings );
    bodyToTopoSettings->getOutputSettings( )->setIntegratedResult( true );

    SingleArcDynamicsSimulator< double > barycentricDynamicsSimulator( bodiesMulti, firstOrderTimeSettings, true );

    SingleArcDynamicsSimulator< double > topocentricDynamicsSimulator( bodiesMulti, bodyToTopoSettings, true );

    auto timeEphemerisDirect = bodiesDirect.getBody( "Earth" )->getTimeScaleConverter( );
    auto timeEphemerisCombined = bodiesMulti.getBody( "Earth" )->getTimeScaleConverter( );

    BOOST_REQUIRE_NE( timeEphemerisDirect, nullptr );
    BOOST_REQUIRE_NE( timeEphemerisCombined, nullptr );

    const auto directFunction =
            timeEphemerisDirect->getTimeDifferenceFunction( barycentric_coordinate_time_scale, local_proper_time_scale, "Graz" );
    const auto combinedFunction =
            timeEphemerisCombined->getTimeDifferenceFunction( barycentric_coordinate_time_scale, local_proper_time_scale, "Graz" );

    std::vector< double > sampleTimes, residuals;
    double sumT = 0.0, sumR = 0.0, sumTT = 0.0, sumTR = 0.0;
    double maxAbsRaw = 0.0;
    const double testTimeStep = 1.0E5;
    for( double currentTime = initialEphemerisTime + 5000.0; currentTime < finalEphemerisTime - 5000.0; currentTime += testTimeStep )
    {
        const double residual = directFunction( currentTime ) - combinedFunction( currentTime );
        sampleTimes.push_back( currentTime );
        residuals.push_back( residual );
        sumT += currentTime;
        sumR += residual;
        sumTT += currentTime * currentTime;
        sumTR += currentTime * residual;
        maxAbsRaw = std::max( maxAbsRaw, std::fabs( residual ) );
        BOOST_CHECK_SMALL( residual, 1.0E-1 );
    }
    const double endpointResidual = directFunction( finalEphemerisTime ) - combinedFunction( finalEphemerisTime );
    BOOST_CHECK_SMALL( endpointResidual, 1.0E-3 );

    const double N = static_cast< double >( residuals.size( ) );
    const double den = N * sumTT - sumT * sumT;
    const double slope = std::fabs( den ) > std::numeric_limits< double >::epsilon( ) ? ( N * sumTR - sumT * sumR ) / den : 0.0;
    const double offset = ( sumR - slope * sumT ) / N;
    double detrendedAmp = 0.0;
    for( size_t i = 0; i < residuals.size( ); ++i )
    {
        detrendedAmp = std::max( detrendedAmp, std::fabs( residuals[ i ] - ( offset + slope * sampleTimes[ i ] ) ) );
    }

    std::cout << "[CombinedProperTimeAndStateDynamics2] samples=" << residuals.size( ) << "  slope=" << slope << " s/s  offset=" << offset
              << " s  max|raw|=" << maxAbsRaw << " s  det_amp=" << detrendedAmp << " s  endpoint=" << endpointResidual << " s" << std::endl;
}

// Keep one no-op case so the Boost module remains runnable while the
// full propagation test above is disabled.
BOOST_AUTO_TEST_CASE( testRelativisticTimePropagationPlaceholder )
{
    BOOST_CHECK( true );
}

// Validates the GR/SR-split dependent variables on the post-Newtonian chain
// (FirstOrderBarycentricToBodyCentricTimeStateDerivative): runs
// FirstOrderBodycentricRelativisticTimePropagatorSettings for Earth with both
// proper_time_rate_kinematic_term and proper_time_rate_potential_term saved at
// every step, and reconstructs both terms independently from Earth's BCRS state
// and the Sun's gravitational parameter at each epoch. Each per-step value must
// match the reconstruction to <1e-15 s/s (numerical-noise level) and the sum
// must be of the expected order of magnitude with the expected (non-positive)
// sign.
BOOST_AUTO_TEST_CASE( testProperTimeRateGrSrSplitDependentVariablesPostNewtonian )
{
    loadStandardSpiceKernels( );

    const double initialEphemerisTime = 0.0;
    const double finalEphemerisTime = 5.0 * physical_constants::JULIAN_DAY;
    const double integrationStep = 100.0;
    const double buffer = 5.0 * integrationStep;

    const std::vector< std::string > bodyNames{ "Earth", "Sun" };
    const std::vector< std::string > perturbingBodies{ "Sun" };

    auto bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Sun" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Sun" );
    bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Earth" );
    bodySettings.at( "Earth" )->bodyDeformationSettings.clear( );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );
    bodySettings.at( "Sun" )->bodyDeformationSettings.clear( );
    bodySettings.at( "Sun" )->gravityFieldVariationSettings.clear( );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    for( const auto& bodyName : bodyNames )
    {
        bodies.getBody( bodyName )->setStateFromEphemeris( initialEphemerisTime );
    }
    bodies.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEphemerisTime );

    auto terminationSettings = std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime );
    auto integratorSettings = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, integrationStep );

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables = {
        propagators::properTimeRateKinematicTermDependentVariable( "Earth" ),
        propagators::properTimeRatePotentialTermDependentVariable( "Earth" )
    };

    auto pnSettings = std::make_shared< FirstOrderBodycentricRelativisticTimePropagatorSettings< double, double > >(
            "Earth",
            perturbingBodies,
            initialEphemerisTime,
            integratorSettings,
            terminationSettings,
            std::map< std::string, std::pair< int, int > >( ),
            []( const double t ) { return t; },
            1.0,
            dependentVariables );
    pnSettings->getOutputSettings( )->setIntegratedResult( true );

    SingleArcDynamicsSimulator< double > simulator( bodies, pnSettings, true );
    const std::map< double, Eigen::VectorXd > dependentHistory =
            simulator.getSingleArcPropagationResults( )->getDependentVariableHistory( );

    BOOST_REQUIRE( !dependentHistory.empty( ) );

    const double inv_c2 = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
    const double expectedRateMagnitude =
            ( 0.5 * 30000.0 * 30000.0 + physical_constants::GRAVITATIONAL_CONSTANT * 1.989e30 / 1.496e11 ) * inv_c2;

    auto earthEphemeris = bodies.getBody( "Earth" )->getEphemeris( );
    auto sunEphemeris = bodies.getBody( "Sun" )->getEphemeris( );
    const double sunGm = bodies.getBody( "Sun" )->getGravityFieldModel( )->getGravitationalParameter( );

    double maxAbsRecoveryError = 0.0;
    for( const auto& entry : dependentHistory )
    {
        const double t = entry.first;
        const Eigen::VectorXd& v = entry.second;
        BOOST_REQUIRE_EQUAL( v.size( ), 2 );
        const double kinematic = v( 0 );
        const double potential = v( 1 );

        BOOST_CHECK( std::isfinite( kinematic ) );
        BOOST_CHECK( std::isfinite( potential ) );
        // Both terms have the form -X/c^2 with X >= 0, so non-positive.
        BOOST_CHECK_LE( kinematic, 0.0 );
        BOOST_CHECK_LE( potential, 0.0 );
        // Magnitude check (within an order of magnitude of the expected level).
        BOOST_CHECK_LT( std::fabs( kinematic + potential ), 5.0 * expectedRateMagnitude );
        BOOST_CHECK_GT( std::fabs( kinematic + potential ), 0.1 * expectedRateMagnitude );

        // Independent reconstruction from ephemerides.
        const Eigen::Vector6d earthState = earthEphemeris->getCartesianState( t );
        const double v2 = earthState.segment( 3, 3 ).squaredNorm( );
        const Eigen::Vector6d sunState = sunEphemeris->getCartesianState( t );
        const double r_es = ( earthState.segment( 0, 3 ) - sunState.segment( 0, 3 ) ).norm( );
        const double uExt = sunGm / r_es;

        const double expectedKinematic = -0.5 * v2 * inv_c2;
        const double expectedPotential = -uExt * inv_c2;

        maxAbsRecoveryError =
                std::max( maxAbsRecoveryError, std::fabs( kinematic - expectedKinematic ) + std::fabs( potential - expectedPotential ) );
    }

    std::cout << "[GrSrSplit-PN] samples=" << dependentHistory.size( ) << "  expected_rate~=" << expectedRateMagnitude
              << " s/s  max_abs_recovery_error=" << maxAbsRecoveryError << " s/s" << std::endl;

    BOOST_CHECK_SMALL( maxAbsRecoveryError, 1.0e-15 );
}

// Validates the GR/SR-split dependent variables on the direct-from-metric path
// (DirectProperTimeRateStateDerivative + SolarSystemMetric): runs
// DirectRelativisticTimePropagatorSettings for the ISS with both
// proper_time_rate_kinematic_term and proper_time_rate_potential_term saved.
// Verifies that:
//   - kinematic term tracks -0.5 * v_BCRS^2 / c^2 from the ISS reference state;
//   - potential term tracks -U_total(r_ISS_BCRS) / c^2 from the SolarSystemMetric;
//   - both terms have the expected sign and order of magnitude.
// Also exercises the non-SolarSystemMetric error path: the dependent variable
// creation must throw if the metric backing the direct-from-metric propagator is
// a Schwarzschild metric (no scalar-potential getter exposed).
BOOST_AUTO_TEST_CASE( testProperTimeRateGrSrSplitDependentVariablesDirectFromMetric )
{
    loadStandardSpiceKernels( );

    const double initialEphemerisTime = 0.0;
    const double finalEphemerisTime = 1.0E4;
    const double integrationStep = 10.0;
    const double buffer = 5.0 * integrationStep;

    const std::vector< std::string > bodyNames{ "Earth", "Sun" };

    auto bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Sun" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Sun" );
    bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Earth" );
    bodySettings.at( "Earth" )->bodyDeformationSettings.clear( );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );
    bodySettings.at( "Sun" )->bodyDeformationSettings.clear( );
    bodySettings.at( "Sun" )->gravityFieldVariationSettings.clear( );

    // Earth-only solar-system metric so the metric is well-defined at the ISS BCRS position
    // (Earth's potential there is finite; Sun is excluded to avoid evaluating Sun's potential
    // at the Sun centre via downstream gradient calls).
    bodySettings.setSpaceTimeSettings( std::make_shared< SpaceTimePropertiesSettings >(
            std::make_shared< SolarSystemSpaceTimeMetricSettings >( std::vector< std::string >{ "Earth" } ),
            std::make_shared< relativity::PPNParameterSet >( 1.0, 1.0 ) ) );

    // Add a tabulated ISS-like ephemeris.
    Eigen::Vector6d issState;
    issState << 6.78e6, 0.0, 0.0, 0.0, 7.66e3, 0.0;
    bodySettings.addSettings( "ISS" );
    bodySettings.at( "ISS" )->constantMass = 1.0;
    bodySettings.at( "ISS" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( issState, "Earth", "ECLIPJ2000" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEphemerisTime );

    auto terminationSettings = std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime );
    auto integratorSettings = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, integrationStep );

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables = {
        propagators::properTimeRateKinematicTermDependentVariable( "ISS" ),
        propagators::properTimeRatePotentialTermDependentVariable( "ISS" ),
    };

    auto directSettings = std::make_shared< DirectRelativisticTimePropagatorSettings< double, double > >(
            std::make_pair( "ISS", "" ),
            initialEphemerisTime,
            integratorSettings,
            terminationSettings,
            []( const double t ) { return t; },
            1.0,
            dependentVariables );
    directSettings->getOutputSettings( )->setIntegratedResult( true );

    SingleArcDynamicsSimulator< double > simulator( bodies, directSettings, true );
    const std::map< double, Eigen::VectorXd > dependentHistory =
            simulator.getSingleArcPropagationResults( )->getDependentVariableHistory( );

    BOOST_REQUIRE( !dependentHistory.empty( ) );

    const double inv_c2 = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;

    // Pull the same ISS BCRS state and Earth GM the propagator uses so the comparison is
    // exact (Body::getState returns the BCRS state, which is what the SolarSystemMetric is
    // evaluated against).
    auto issBody = bodies.getBody( "ISS" );
    auto earthBody = bodies.getBody( "Earth" );
    const double earthGm = earthBody->getGravityFieldModel( )->getGravitationalParameter( );

    double maxAbsKinError = 0.0;
    double maxAbsPotError = 0.0;
    for( const auto& entry : dependentHistory )
    {
        const double t = entry.first;
        const Eigen::VectorXd& v = entry.second;
        BOOST_REQUIRE_EQUAL( v.size( ), 2 );
        const double kinematic = v( 0 );
        const double potential = v( 1 );

        BOOST_CHECK( std::isfinite( kinematic ) );
        BOOST_CHECK( std::isfinite( potential ) );
        BOOST_CHECK_LE( kinematic, 0.0 );
        BOOST_CHECK_LE( potential, 0.0 );

        // Tick the body states forward to time t before reading them back.
        issBody->setStateFromEphemeris( t );
        earthBody->setStateFromEphemeris( t );

        const Eigen::Vector6d issBcrsState = issBody->getState( );
        const double v2 = issBcrsState.segment( 3, 3 ).squaredNorm( );
        const double expectedKinematic = -0.5 * v2 * inv_c2;
        maxAbsKinError = std::max( maxAbsKinError, std::fabs( kinematic - expectedKinematic ) );

        // Earth-only metric: U_total at the ISS BCRS position = U_E evaluated relative to Earth's centre.
        const Eigen::Vector6d earthBcrsState = earthBody->getState( );
        const double rIssRelEarth = ( issBcrsState.segment( 0, 3 ) - earthBcrsState.segment( 0, 3 ) ).norm( );
        const double uTotal = earthGm / rIssRelEarth;
        const double expectedPotential = -uTotal * inv_c2;
        maxAbsPotError = std::max( maxAbsPotError, std::fabs( potential - expectedPotential ) );
    }

    std::cout << "[GrSrSplit-Direct] samples=" << dependentHistory.size( ) << "  max_abs_kin_error=" << maxAbsKinError
              << " s/s  max_abs_pot_error=" << maxAbsPotError << " s/s" << std::endl;

    // Tolerances loose enough to absorb ECLIPJ2000<->BCRS rotation noise but tight enough
    // to catch any sign or magnitude regression in the propagator dispatch.
    BOOST_CHECK_SMALL( maxAbsKinError, 1.0e-12 );
    BOOST_CHECK_SMALL( maxAbsPotError, 1.0e-12 );

    // Non-SolarSystemMetric error path: build a separate environment with a Schwarzschild
    // metric and confirm that requesting the potential dependent variable throws when
    // the propagator is created.
    auto bodySettingsSchwarz =
            getDefaultBodySettings( std::vector< std::string >{ "Earth" }, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettingsSchwarz.at( "Earth" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Earth" );
    bodySettingsSchwarz.setSpaceTimeSettings(
            std::make_shared< SpaceTimePropertiesSettings >( std::make_shared< SchwarzschildSpaceTimeMetricSettings >( "Earth" ) ) );
    bodySettingsSchwarz.addSettings( "ISS" );
    bodySettingsSchwarz.at( "ISS" )->constantMass = 1.0;
    bodySettingsSchwarz.at( "ISS" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( issState, "Earth", "ECLIPJ2000" );
    SystemOfBodies bodiesSchwarz = createSystemOfBodies( bodySettingsSchwarz );

    auto schwarzPotentialOnly = std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >{
        propagators::properTimeRatePotentialTermDependentVariable( "ISS" ),
    };
    auto schwarzSettings = std::make_shared< DirectRelativisticTimePropagatorSettings< double, double > >(
            std::make_pair( "ISS", "" ),
            initialEphemerisTime,
            integratorSettings,
            terminationSettings,
            []( const double t ) { return t; },
            1.0,
            schwarzPotentialOnly );
    BOOST_CHECK_THROW( SingleArcDynamicsSimulator< double >( bodiesSchwarz, schwarzSettings, true ), std::runtime_error );
}

// Progressive-complexity diagnostic. Runs several configurations from simplest
// (Sun + Earth point-mass, equator station) up through Earth SH(20) and Moon,
// sampling BOTH the full station residual AND the body-center TCB<->TCG-only
// residual. This separates integrand-level disagreement from station-level
// reconstruction (direct correction -v_E*r/c^2 and topocentric quadrature).
// Output is written to stdout for inspection; the test asserts only finiteness.
BOOST_AUTO_TEST_CASE( testProgressiveComplexityDiagnostic )
{
    loadStandardSpiceKernels( );

    struct Level {
        std::string name;
        std::vector< std::string > bodyNames;
        std::vector< std::string > perturbingBodies;
        std::map< std::string, std::pair< int, int > > metricShOrders;
        std::map< std::string, std::pair< int, int > > tcbToTcgShOrders;
        int topoEarthShDegree;
    };

    const std::vector< Level > levels = {
        { "L1_EarthSun_PointMass_Equator", { "Earth", "Sun" }, { "Sun" }, {}, {}, 0 },
        { "L2_EarthSun_EarthSh20_Equator", { "Earth", "Sun" }, { "Sun" }, { { "Earth", { 20, 20 } } }, {}, 20 },
        { "L3_EarthSunMoon_EarthSh20_MoonSh2_Equator",
          { "Earth", "Sun", "Moon" },
          { "Sun", "Moon" },
          { { "Earth", { 20, 20 } }, { "Moon", { 2, 2 } } },
          { { "Moon", { 2, 2 } } },
          20 }
    };

    const double duration = 90.0 * physical_constants::JULIAN_DAY;
    const double initialEphemerisTime = -0.5 * duration;
    const double finalEphemerisTime = 0.5 * duration;
    const double integrationStep = 100.0;
    const double buffer = 5.0 * integrationStep;
    const std::string stationName = "DiagStation";
    const double equatorialRadius = 6378137.0;
    const Eigen::Vector3d equatorStationPosition =
            ( Eigen::Vector3d( ) << equatorialRadius / std::sqrt( 2.0 ), equatorialRadius / std::sqrt( 2.0 ), 0.0 ).finished( );

    for( const Level& L : levels )
    {
        std::cout << "\n=== " << L.name << " ===" << std::endl;

        auto bodySettings = getDefaultBodySettings( L.bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

        bodySettings.at( "Sun" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Sun" );
        bodySettings.at( "Sun" )->bodyDeformationSettings.clear( );
        bodySettings.at( "Sun" )->gravityFieldVariationSettings.clear( );

        if( L.metricShOrders.count( "Earth" ) )
        {
            bodySettings.at( "Earth" )->gravityFieldSettings =
                    std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( ggm02c, 20 );
        }
        else
        {
            bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Earth" );
        }
        bodySettings.at( "Earth" )->bodyDeformationSettings.clear( );
        bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );

        const bool hasMoon = std::find( L.bodyNames.begin( ), L.bodyNames.end( ), "Moon" ) != L.bodyNames.end( );
        if( hasMoon )
        {
            if( L.metricShOrders.count( "Moon" ) )
            {
                bodySettings.at( "Moon" )->gravityFieldSettings =
                        std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( lpe200, 2 );
            }
            bodySettings.at( "Moon" )->bodyDeformationSettings.clear( );
            bodySettings.at( "Moon" )->gravityFieldVariationSettings.clear( );
        }

        std::map< std::pair< std::string, std::string >, Eigen::Vector3d > stationMap;
        stationMap[ std::make_pair( "Earth", stationName ) ] = equatorStationPosition;

        SystemOfBodies bodiesDirect = createSystemOfBodies( bodySettings );
        SystemOfBodies bodiesPn = createSystemOfBodies( bodySettings );
        createGroundStations( bodiesDirect, stationMap );
        createGroundStations( bodiesPn, stationMap );

        for( const auto& bn : L.bodyNames )
        {
            if( bodiesDirect.doesBodyExist( bn ) )
            {
                bodiesDirect.getBody( bn )->setStateFromEphemeris( initialEphemerisTime );
            }
            if( bodiesPn.doesBodyExist( bn ) )
            {
                bodiesPn.getBody( bn )->setStateFromEphemeris( initialEphemerisTime );
            }
        }
        bodiesDirect.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEphemerisTime );
        bodiesPn.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEphemerisTime );

        auto terminationSettings = std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime );
        auto intSetD = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, integrationStep );
        auto intSetP = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, integrationStep );
        auto intSetT = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, integrationStep );

        auto metricSet = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
                L.bodyNames, std::vector< std::string >( ), L.metricShOrders, std::vector< std::string >( ), false );
        auto baseMetric = createSpaceTimeMetric( metricSet, bodiesDirect );
        evaluatedMetricObjects.clear( );
        evaluatedMetricObjects[ std::make_pair( "Earth", stationName ) ] = baseMetric->Clone( );

        auto directFromMetricSettings = std::make_shared< DirectRelativisticTimePropagatorSettings< double, double > >(
                std::make_pair( "Earth", stationName ), initialEphemerisTime, intSetD, terminationSettings );
        directFromMetricSettings->getOutputSettings( )->setIntegratedResult( true );
        SingleArcDynamicsSimulator< double > simDirect( bodiesDirect, directFromMetricSettings, true );

        auto tcbToTcgSettings = std::make_shared< FirstOrderBodycentricRelativisticTimePropagatorSettings< double, double > >(
                "Earth", L.perturbingBodies, initialEphemerisTime, intSetP, terminationSettings, L.tcbToTcgShOrders );
        tcbToTcgSettings->getOutputSettings( )->setIntegratedResult( true );

        Eigen::Matrix< double, Eigen::Dynamic, 1 > initRel = Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 1 );
        auto topoSettings = std::make_shared< BodycenteredToTopocentricTimePropagatorSettings< double, double > >(
                std::make_pair( "Earth", stationName ),
                false,
                L.topoEarthShDegree,
                false,
                L.perturbingBodies,
                initRel,
                initialEphemerisTime,
                intSetT,
                terminationSettings );
        topoSettings->getOutputSettings( )->setIntegratedResult( true );

        SingleArcDynamicsSimulator< double > simTcb( bodiesPn, tcbToTcgSettings, true );
        SingleArcDynamicsSimulator< double > simTopo( bodiesPn, topoSettings, true );

        auto epD = bodiesDirect.getBody( "Earth" )->getTimeScaleConverter( );
        auto epP = bodiesPn.getBody( "Earth" )->getTimeScaleConverter( );
        BOOST_REQUIRE_NE( epD, nullptr );
        BOOST_REQUIRE_NE( epP, nullptr );

        const auto directStation =
                epD->getTimeDifferenceFunction( barycentric_coordinate_time_scale, local_proper_time_scale, stationName );
        const auto chainStation = epP->getTimeDifferenceFunction( barycentric_coordinate_time_scale, local_proper_time_scale, stationName );
        const auto chainCenter = epP->getTimeDifferenceFunction( barycentric_coordinate_time_scale, body_centered_coordinate_time_scale );

        auto earthEphemeris = bodiesDirect.getBody( "Earth" )->getEphemeris( );
        auto earthRotation = bodiesDirect.getBody( "Earth" )->getRotationalEphemeris( );
        const Eigen::Vector3d stationPcrs = equatorStationPosition;
        const double inv_c2 = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;

        // Missing-correction hypothesis: chain dots v_E_inertial with r_body_fixed
        // (treated as if it were inertial). The "true" direct correction would
        // dot with r_inertial(t) = R_body_to_inertial(t) * r_body_fixed.
        // Missing term = -v_E * ( r_inertial(t) - r_body_fixed ) / c^2.
        auto missingCorrection = [ & ]( double t ) {
            const Eigen::Vector6d stateE = earthEphemeris->getCartesianState( t );
            const Eigen::Vector3d vE = stateE.segment( 3, 3 );
            const Eigen::Matrix3d R = earthRotation->getRotationToBaseFrame( t ).toRotationMatrix( );
            const Eigen::Vector3d rInertial = R * stationPcrs;
            const Eigen::Vector3d delta = rInertial - stationPcrs;
            return -vE.dot( delta ) * inv_c2;
        };

        std::vector< double > times, rTotal, rCorrected;
        times.reserve( static_cast< unsigned int >( duration / integrationStep ) + 1U );
        for( double t = initialEphemerisTime + integrationStep; t <= finalEphemerisTime - integrationStep; t += integrationStep )
        {
            times.push_back( t );
            const double raw = directStation( t ) - chainStation( t );
            rTotal.push_back( raw );
            rCorrected.push_back( raw - missingCorrection( t ) );
        }

        auto detrend = [ & ]( const std::vector< double >& res ) {
            double sumT = 0.0, sumR = 0.0, sumTT = 0.0, sumTR = 0.0;
            for( size_t i = 0; i < times.size( ); ++i )
            {
                sumT += times[ i ];
                sumR += res[ i ];
                sumTT += times[ i ] * times[ i ];
                sumTR += times[ i ] * res[ i ];
            }
            const double N = static_cast< double >( res.size( ) );
            const double den = N * sumTT - sumT * sumT;
            const double slope = std::fabs( den ) > std::numeric_limits< double >::epsilon( ) ? ( N * sumTR - sumT * sumR ) / den : 0.0;
            const double offset = ( sumR - slope * sumT ) / N;
            double amp = 0.0;
            for( size_t i = 0; i < res.size( ); ++i )
            {
                amp = std::max( amp, std::fabs( res[ i ] - ( offset + slope * times[ i ] ) ) );
            }
            return std::make_tuple( slope, offset, amp );
        };

        const auto [ sT, oT, aT ] = detrend( rTotal );
        const auto [ sX, oX, aX ] = detrend( rCorrected );

        std::cout << std::scientific << std::setprecision( 3 );
        std::cout << "  samples=" << rTotal.size( ) << std::endl;
        std::cout << "  raw        (direct - chain):            slope=" << sT << " s/s  offset=" << oT << " s  det_amp=" << aT << " s"
                  << std::endl;
        std::cout << "  corrected  (direct - chain - missing):  slope=" << sX << " s/s  offset=" << oX << " s  det_amp=" << aX << " s"
                  << std::endl;

        BOOST_CHECK( std::isfinite( aT ) );
        BOOST_CHECK( std::isfinite( aX ) );
    }
}

// Original strict one-year test: Earth equator station, 100 s step,
// detrended residual amplitude < 6 ns (matches Dominic's baseline after
// the body-fixed vs inertial station-position fix at
// setNumericallyIntegratedStates.h:1259).
BOOST_AUTO_TEST_CASE( testDirectFromMetricVsChainedPnEquatorOneYearEarthMoonSun )
{
    loadStandardSpiceKernels( );

    const double initialEphemerisTime = -0.5 * 365.0 * physical_constants::JULIAN_DAY;
    const double finalEphemerisTime = 0.5 * 365.0 * physical_constants::JULIAN_DAY;
    const double integrationStep = 100.0;
    const double buffer = 5.0 * integrationStep;

    const std::vector< std::string > bodyNames{ "Earth", "Sun", "Moon" };
    const std::vector< std::string > perturbingBodies{ "Sun", "Moon" };
    const std::string stationName = "DiagStation";
    const double equatorialRadius = 6378137.0;
    const Eigen::Vector3d equatorStationPosition =
            ( Eigen::Vector3d( ) << equatorialRadius / std::sqrt( 2.0 ), equatorialRadius / std::sqrt( 2.0 ), 0.0 ).finished( );

    auto bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Sun" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Sun" );
    bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( ggm02c, 20 );
    bodySettings.at( "Moon" )->gravityFieldSettings = std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( lpe200, 2 );

    bodySettings.at( "Earth" )->bodyDeformationSettings.clear( );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );
    bodySettings.at( "Moon" )->bodyDeformationSettings.clear( );
    bodySettings.at( "Moon" )->gravityFieldVariationSettings.clear( );
    bodySettings.at( "Sun" )->bodyDeformationSettings.clear( );
    bodySettings.at( "Sun" )->gravityFieldVariationSettings.clear( );

    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStations;
    groundStations[ std::make_pair( "Earth", stationName ) ] = equatorStationPosition;

    SystemOfBodies bodiesDirect = createSystemOfBodies( bodySettings );
    SystemOfBodies bodiesPn = createSystemOfBodies( bodySettings );
    createGroundStations( bodiesDirect, groundStations );
    createGroundStations( bodiesPn, groundStations );

    for( const auto& bodyName : bodyNames )
    {
        if( bodiesDirect.doesBodyExist( bodyName ) )
        {
            bodiesDirect.getBody( bodyName )->setStateFromEphemeris( initialEphemerisTime );
        }
        if( bodiesPn.doesBodyExist( bodyName ) )
        {
            bodiesPn.getBody( bodyName )->setStateFromEphemeris( initialEphemerisTime );
        }
    }

    bodiesDirect.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEphemerisTime );
    bodiesPn.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEphemerisTime );

    auto terminationSettings = std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime );
    auto integratorSettingsDirect = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, integrationStep );
    auto integratorSettingsPn = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, integrationStep );
    auto integratorSettingsPnTopo = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, integrationStep );

    const std::map< std::string, std::pair< int, int > > sphericalHarmonicOrders = { { "Earth", std::make_pair( 20, 20 ) },
                                                                                     { "Moon", std::make_pair( 2, 2 ) } };

    auto metricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >( std::vector< std::string >{ "Sun", "Earth", "Moon" },
                                                                                  std::vector< std::string >( ),
                                                                                  sphericalHarmonicOrders,
                                                                                  std::vector< std::string >( ),
                                                                                  false );
    auto baseMetric = createSpaceTimeMetric( metricSettings, bodiesDirect );
    evaluatedMetricObjects.clear( );
    evaluatedMetricObjects[ std::make_pair( "Earth", stationName ) ] = baseMetric->Clone( );

    auto directFromMetricSettings = std::make_shared< DirectRelativisticTimePropagatorSettings< double, double > >(
            std::make_pair( "Earth", stationName ), initialEphemerisTime, integratorSettingsDirect, terminationSettings );
    directFromMetricSettings->getOutputSettings( )->setIntegratedResult( true );

    SingleArcDynamicsSimulator< double > directDynamicsSimulator( bodiesDirect, directFromMetricSettings, true );

    auto firstOrderPnSettings = std::make_shared< FirstOrderBodycentricRelativisticTimePropagatorSettings< double, double > >(
            "Earth",
            perturbingBodies,
            initialEphemerisTime,
            integratorSettingsPn,
            terminationSettings,
            std::map< std::string, std::pair< int, int > >{ { "Moon", std::make_pair( 2, 2 ) } } );
    firstOrderPnSettings->getOutputSettings( )->setIntegratedResult( true );

    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRelativisticState = Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 1 );
    auto bodyToTopoSettings = std::make_shared< BodycenteredToTopocentricTimePropagatorSettings< double, double > >(
            std::make_pair( "Earth", stationName ),
            false,  // useAccelerationTerm OFF (IBP-paired with the -v_E.r/c^2 direct correction)
            20,
            false,
            perturbingBodies,
            initialRelativisticState,
            initialEphemerisTime,
            integratorSettingsPnTopo,
            terminationSettings );
    bodyToTopoSettings->getOutputSettings( )->setIntegratedResult( true );

    SingleArcDynamicsSimulator< double > barycentricPnDynamicsSimulator( bodiesPn, firstOrderPnSettings, true );

    SingleArcDynamicsSimulator< double > topocentricPnDynamicsSimulator( bodiesPn, bodyToTopoSettings, true );

    auto timeEphemerisDirect = bodiesDirect.getBody( "Earth" )->getTimeScaleConverter( );
    auto timeEphemerisPn = bodiesPn.getBody( "Earth" )->getTimeScaleConverter( );
    BOOST_REQUIRE_NE( timeEphemerisDirect, nullptr );
    BOOST_REQUIRE_NE( timeEphemerisPn, nullptr );

    const auto directFunction =
            timeEphemerisDirect->getTimeDifferenceFunction( barycentric_coordinate_time_scale, local_proper_time_scale, stationName );
    const auto pnFunction =
            timeEphemerisPn->getTimeDifferenceFunction( barycentric_coordinate_time_scale, local_proper_time_scale, stationName );

    std::vector< double > sampleTimes, residuals;
    sampleTimes.reserve( static_cast< unsigned int >( ( finalEphemerisTime - initialEphemerisTime ) / integrationStep ) + 1U );
    residuals.reserve( sampleTimes.capacity( ) );

    double sumTime = 0.0, sumResidual = 0.0, sumTimeSquared = 0.0, sumTimeResidual = 0.0;
    const double residualAmplitudeTolerance = 6.0e-9;

    for( double currentTime = initialEphemerisTime + integrationStep; currentTime <= finalEphemerisTime - integrationStep;
         currentTime += integrationStep )
    {
        const double direct = directFunction( currentTime );
        const double pn = pnFunction( currentTime );
        const double residual = direct - pn;

        sampleTimes.push_back( currentTime );
        residuals.push_back( residual );
        sumTime += currentTime;
        sumResidual += residual;
        sumTimeSquared += currentTime * currentTime;
        sumTimeResidual += currentTime * residual;
    }

    BOOST_REQUIRE( !residuals.empty( ) );
    const double N = static_cast< double >( residuals.size( ) );
    const double den = N * sumTimeSquared - sumTime * sumTime;
    const double slope =
            std::fabs( den ) > std::numeric_limits< double >::epsilon( ) ? ( N * sumTimeResidual - sumTime * sumResidual ) / den : 0.0;
    const double offset = ( sumResidual - slope * sumTime ) / N;

    double amp = 0.0, maxAbsRaw = 0.0;
    for( size_t i = 0; i < residuals.size( ); ++i )
    {
        const double detrended = residuals[ i ] - ( offset + slope * sampleTimes[ i ] );
        amp = std::max( amp, std::fabs( detrended ) );
        maxAbsRaw = std::max( maxAbsRaw, std::fabs( residuals[ i ] ) );
    }

    std::cout << "[OneYearEarthMoonSun] samples=" << residuals.size( ) << "  slope=" << slope << " s/s  offset=" << offset
              << " s  max|raw|=" << maxAbsRaw << " s  det_amp=" << amp << " s  tol=" << residualAmplitudeTolerance << " s" << std::endl;

    BOOST_CHECK( std::isfinite( amp ) );
    BOOST_CHECK_SMALL( amp, residualAmplitudeTolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
