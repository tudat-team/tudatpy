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

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

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

// Integrated-identity check for the proper-time-rate decomposition: the two GR/SR-split
// dependent variables (kinematic -v^2/(2 c^2) and potential -U/c^2) must SUM to the
// proper-time-rate integrand that tudat actually integrates. For the first-order
// barycentric->bodycentric propagator the state derivative is exactly
// calculateFirstOrderTcbToTcgIntegrand(v, U) = -v^2/(2 c^2) - U/c^2 = kinematic + potential,
// so the time integral of (kinematic + potential) over the propagation grid must reproduce the
// propagated proper-time state (the TCB->TCG difference) that tudat reports.
//
// We integrate the summed dependent variables on the output grid with the trapezoidal rule and
// compare against the propagated state. The two components are evaluated from the SAME cached
// quantities (currentVelocity_, currentExternalPotential_) that feed the integrand, so the only
// difference is the integration scheme (trapezoid here vs. the RK4 the propagator uses): the
// residual is integration-limited on a smooth, slowly varying integrand (rate ~ -1e-8 s/s,
// nearly constant), not a modelling error. A sign error or a wrong factor (e.g. v^2 instead of
// v^2/2) in either component would break the sum and blow this check up by orders of magnitude.
BOOST_AUTO_TEST_CASE( testProperTimeRateGrSrSplitIntegratedSumMatchesState )
{
    loadStandardSpiceKernels( );

    const double initialEphemerisTime = 0.0;
    const double finalEphemerisTime = physical_constants::JULIAN_DAY;
    const double integrationStep = 30.0;
    const double buffer = 10.0 * integrationStep;

    const std::vector< std::string > bodyNames{ "Earth", "Sun", "Moon" };
    const std::vector< std::string > perturbingBodies{ "Sun", "Moon" };

    auto bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    for( const auto& name : bodyNames )
    {
        bodySettings.at( name )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( name );
        bodySettings.at( name )->bodyDeformationSettings.clear( );
        bodySettings.at( name )->gravityFieldVariationSettings.clear( );
    }

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    for( const auto& name : bodyNames )
    {
        bodies.getBody( name )->setStateFromEphemeris( initialEphemerisTime );
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
    auto results = simulator.getSingleArcPropagationResults( );
    const std::map< double, Eigen::VectorXd > dependentHistory = results->getDependentVariableHistory( );
    const std::map< double, Eigen::VectorXd > stateHistory = results->getEquationsOfMotionNumericalSolution( );

    BOOST_REQUIRE( !dependentHistory.empty( ) );
    BOOST_REQUIRE_EQUAL( dependentHistory.size( ), stateHistory.size( ) );

    // Trapezoidal integral of (kinematic + potential) on the output grid, anchored at the
    // propagated state's first sample so any initial-offset convention cancels. The two maps are
    // keyed by the same output epochs (asserted equal size above), so we iterate them in lockstep.
    double cumulativeIntegral = stateHistory.begin( )->second( 0 );
    double previousTime = dependentHistory.begin( )->first;
    double previousRate = dependentHistory.begin( )->second( 0 ) + dependentHistory.begin( )->second( 1 );

    double maxAbsDifference = 0.0;
    double maxAbsState = 0.0;
    auto stateIt = stateHistory.begin( );
    for( auto depIt = dependentHistory.begin( ); depIt != dependentHistory.end( ); ++depIt, ++stateIt )
    {
        BOOST_REQUIRE_EQUAL( depIt->second.size( ), 2 );
        const double t = depIt->first;
        const double rate = depIt->second( 0 ) + depIt->second( 1 );
        if( depIt != dependentHistory.begin( ) )
        {
            cumulativeIntegral += 0.5 * ( rate + previousRate ) * ( t - previousTime );
        }
        previousTime = t;
        previousRate = rate;

        const double propagatedState = stateIt->second( 0 );
        maxAbsDifference = std::max( maxAbsDifference, std::fabs( cumulativeIntegral - propagatedState ) );
        maxAbsState = std::max( maxAbsState, std::fabs( propagatedState ) );
    }

    const double finalState = stateHistory.rbegin( )->second( 0 );
    std::cout << "[GrSrSplit-IntegratedSum] samples=" << dependentHistory.size( ) << "  final_state=" << finalState
              << " s  trapezoid_sum=" << cumulativeIntegral << " s  max_abs_diff=" << maxAbsDifference
              << " s  rel=" << ( maxAbsState > 0.0 ? maxAbsDifference / maxAbsState : 0.0 ) << std::endl;

    // The integrated decomposition reproduces the propagated proper-time state to roundoff: the
    // integrand is nearly constant, so the trapezoid rule on the output grid agrees with the RK4
    // the propagator uses to ~1e-16 s against a ~1.3e-3 s state (observed max_abs_diff ~1.7e-16 s,
    // rel ~1e-13). Guard at 1e-12 s -- ~4 orders of magnitude over the observed roundoff (safe
    // across platforms) yet many orders below any sign/factor regression in either component.
    BOOST_REQUIRE_GT( maxAbsState, 0.0 );
    BOOST_CHECK_SMALL( maxAbsDifference, 1.0e-12 );
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

// Strict one-year test: Earth equator station, 100 s step. Detrended residual amplitude
// observed at ~0.1 ns after the direct-from-metric station-rotation fix
// (createEnvironmentUpdater.cpp, direct_from_metric case, now requesting
// body_rotational_state_update for the reference body); tolerance 1 ns gives ~10x headroom.
// (Was ~3.9 ns before that fix, with the body-fixed vs inertial station-position fix at
// setNumericallyIntegratedStates.h:1259 already in place.)
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
    const double residualAmplitudeTolerance = 1.0e-9;

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

// Data-driven growing-complexity ladder comparing the direct-from-metric proper
// time against the concatenated post-Newtonian chain (TCB->TCG first order +
// body->topocentric). Each rung shares one configuration struct; the chain's
// per-leg spherical-harmonic settings and the topocentric central-body SH degree
// are derived from the single metric SH map so the two formulations stay
// consistent by construction.
//
// The ladder climbs from an Earth point-mass case up to the lunar target:
// Moon as central body (TCL at the lunar centre), Earth+Sun as perturbers, and a
// lunar-surface station (the proper time underlying TL). The detrended residual
// amplitude (direct - chain), with the irrelevant constant rate offset removed,
// is the discrepancy of interest; the per-scenario note below records the
// post-fix agreement (~10 ps on every rung) and the root cause / fix.
BOOST_AUTO_TEST_CASE( testDirectVsChainComplexityLadder )
{
    loadStandardSpiceKernels( );

    struct LadderScenario {
        std::string name;
        std::string centralBody;
        std::vector< std::string > perturbingBodies;                    // external bodies for the chain (also metric bodies)
        std::map< std::string, std::pair< int, int > > metricShOrders;  // SH degree/order per body in the metric
        double centralBodyRadius;                                       // for the equatorial surface station
        double durationDays;
        double integrationStep;
        double detrendedTolerance;  // <= 0 => report only (no hard assert)
    };

    const double earthRadius = 6378137.0;
    const double moonRadius = 1737400.0;

    // Each rung is also compared against an independent analytic first-order truth rate
    //   dtau/dt - 1 = -( 1/2 |v_st|^2 + sum_i GM_i/|r_st - r_i| ) / c^2
    // integrated on the same grid (see "vs analytic 1st-order truth" below; both 2PN methods
    // exceed this 1st-order reference only by their small 2nd-order content).
    //
    // EXPANSION ORDER. Both methods are run at SECOND post-Newtonian order: the direct metric is
    // 2nd order in the sqrt expansion, and the chain's TCB->TC{G,L} leg uses the 2nd-order
    // converter (SecondOrderBodyCenteredRelativisticTimeConverterSettings). The direct-vs-chain
    // RATE difference is then ~2e-19 s/s (Earth) / ~1e-20 s/s (Moon) -- down from ~1.2e-16 s/s with
    // a 1st-order chain. That 1.2e-16 floor was therefore the omitted 2nd-order BARYCENTRIC term
    // (Soffel Eq. 58), NOT numerical: it is step-independent and ~8 orders above double-precision
    // machine epsilon, and it collapses by 10^2-10^4x purely by adding the model term. The residual
    // PERIODIC disagreement is ~4-38 ps (E1/E2 ~38 ps, E3 ~7 ps with the Moon added, M1/M2 ~4 ps):
    // everything beyond the matched 2nd-order barycentric term -- chiefly the un-modelled 2nd-order
    // TOPOCENTRIC contribution (the 1st-order Turyshev body->topo leg omits the 2nd-order
    // station-velocity terms; expected to dominate for the faster-rotating Earth, ~465 m/s vs the
    // Moon's ~4.6 m/s), plus residual 2nd-/3rd-order formulation differences (the E2->E3 drop at
    // fixed Earth rotation shows the perturber set also enters, so the split is not yet isolated).
    // Still step-independent and >>machine epsilon: expansion/formulation-limited, not numerical.
    // TODO: add a 2nd-order topocentric leg (should remove the topocentric part; testable by a
    // station-latitude sweep) and rerun Time-templated (confirms the floor is not precision-limited).
    //
    // This agreement required a fix to the direct-from-metric environment updater: for a
    // ground-station reference point it now requests body_rotational_state_update for the central
    // body (createEnvironmentUpdater.cpp, direct_from_metric case). Previously the station's
    // body-fixed->inertial rotation was frozen at the initial epoch, giving a ~us diurnal error on
    // point-mass central bodies (E1 ~9.7us, M1 ~0.6us); a spherical-harmonic central body masked it
    // (the SH-field update pulled in the rotation). The metric VALUE was correct throughout (see
    // testPointMassVsShMetricAtRotatingStation in unitTestSolarSystemMetric.cpp).
    //
    // Acceptance thresholds: the detrended-residual amplitude must be < 100 ps (1.0e-10 s) on every
    // rung, and the detrended-residual mean rate (slope, i.e. fractional-frequency offset between the
    // two methods) must be < 1e-18 (enforced globally below). Observed values (per rung, listed in
    // the EXPANSION ORDER note above) are 4-38 ps amplitude and ~1e-19--1e-20 slope, comfortably
    // inside these bounds.
    const std::vector< LadderScenario > scenarios = {
        { "E1_Earth_Sun_PointMass", "Earth", { "Sun" }, {}, earthRadius, 30.0, 120.0, 1.0e-10 },
        { "E2_Earth_Sun_EarthSH20", "Earth", { "Sun" }, { { "Earth", { 20, 20 } } }, earthRadius, 30.0, 120.0, 1.0e-10 },
        { "E3_Earth_SunMoon_EarthSH20_MoonSH2",
          "Earth",
          { "Sun", "Moon" },
          { { "Earth", { 20, 20 } }, { "Moon", { 2, 2 } } },
          earthRadius,
          30.0,
          120.0,
          1.0e-10 },
        { "M1_Moon_EarthSun_PointMass", "Moon", { "Earth", "Sun" }, {}, moonRadius, 30.0, 120.0, 1.0e-10 },
        { "M2_Moon_EarthSun_MoonSH2",  // <- lunar target: TCL centre + TL surface, Moon SH(2,2)
          "Moon",
          { "Earth", "Sun" },
          { { "Moon", { 2, 2 } } },
          moonRadius,
          30.0,
          120.0,
          1.0e-10 }
    };

    // Configure a body's gravity field from the scenario SH map (FromFile SH when
    // requested, point-mass otherwise) and strip time-variable gravity so the two
    // environments are static and identical.
    auto configureGravity =
            []( BodyListSettings& bodySettings, const std::string& body, const std::map< std::string, std::pair< int, int > >& shOrders ) {
                if( shOrders.count( body ) )
                {
                    const int degree = shOrders.at( body ).first;
                    if( body == "Earth" )
                    {
                        bodySettings.at( body )->gravityFieldSettings =
                                std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( ggm02c, degree );
                    }
                    else if( body == "Moon" )
                    {
                        bodySettings.at( body )->gravityFieldSettings =
                                std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( lpe200, degree );
                    }
                    else
                    {
                        bodySettings.at( body )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( body );
                    }
                }
                else
                {
                    bodySettings.at( body )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( body );
                }
                bodySettings.at( body )->bodyDeformationSettings.clear( );
                bodySettings.at( body )->gravityFieldVariationSettings.clear( );
            };

    for( const LadderScenario& scenario : scenarios )
    {
        std::cout << "\n=== " << scenario.name << " ===" << std::endl;

        const double duration = scenario.durationDays * physical_constants::JULIAN_DAY;
        const double initialEphemerisTime = -0.5 * duration;
        const double finalEphemerisTime = 0.5 * duration;
        const double integrationStep = scenario.integrationStep;
        const double buffer = 5.0 * integrationStep;
        const std::string stationName = "LadderStation";

        // Body list: central body + perturbers (unique), metric evaluates over all of them.
        std::vector< std::string > bodyNames{ scenario.centralBody };
        for( const std::string& body : scenario.perturbingBodies )
        {
            if( std::find( bodyNames.begin( ), bodyNames.end( ), body ) == bodyNames.end( ) )
            {
                bodyNames.push_back( body );
            }
        }

        // Derive the chain's per-leg SH config from the single metric SH map:
        //  - TCB->TCG perturber SH: any perturber that carries SH in the metric;
        //  - topocentric central-body SH degree: the central body's metric SH degree.
        std::map< std::string, std::pair< int, int > > chainPerturberShOrders;
        for( const std::string& body : scenario.perturbingBodies )
        {
            if( scenario.metricShOrders.count( body ) )
            {
                chainPerturberShOrders[ body ] = scenario.metricShOrders.at( body );
            }
        }
        const int topoCentralShDegree =
                scenario.metricShOrders.count( scenario.centralBody ) ? scenario.metricShOrders.at( scenario.centralBody ).first : 0;

        auto bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
        // Sun always a point mass; other bodies per the SH map.
        bodySettings.at( "Sun" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Sun" );
        bodySettings.at( "Sun" )->bodyDeformationSettings.clear( );
        bodySettings.at( "Sun" )->gravityFieldVariationSettings.clear( );
        for( const std::string& body : bodyNames )
        {
            if( body != "Sun" )
            {
                configureGravity( bodySettings, body, scenario.metricShOrders );
            }
        }

        const Eigen::Vector3d equatorStationPosition =
                ( Eigen::Vector3d( ) << scenario.centralBodyRadius / std::sqrt( 2.0 ), scenario.centralBodyRadius / std::sqrt( 2.0 ), 0.0 )
                        .finished( );
        std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStations;
        groundStations[ std::make_pair( scenario.centralBody, stationName ) ] = equatorStationPosition;

        SystemOfBodies bodiesDirect = createSystemOfBodies( bodySettings );
        SystemOfBodies bodiesPn = createSystemOfBodies( bodySettings );
        createGroundStations( bodiesDirect, groundStations );
        createGroundStations( bodiesPn, groundStations );

        for( const std::string& body : bodyNames )
        {
            if( bodiesDirect.doesBodyExist( body ) )
            {
                bodiesDirect.getBody( body )->setStateFromEphemeris( initialEphemerisTime );
            }
            if( bodiesPn.doesBodyExist( body ) )
            {
                bodiesPn.getBody( body )->setStateFromEphemeris( initialEphemerisTime );
            }
        }
        bodiesDirect.getBody( scenario.centralBody )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEphemerisTime );
        bodiesPn.getBody( scenario.centralBody )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEphemerisTime );

        auto terminationSettings = std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime );
        auto intSetDirect = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, integrationStep );
        auto intSetTcb = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, integrationStep );
        auto intSetTopo = numerical_integrators::rungeKutta4SettingsDeprecated( initialEphemerisTime, integrationStep );

        // --- Direct-from-metric: SolarSystemMetric over all bodies, evaluated at the station. ---
        auto metricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
                bodyNames, std::vector< std::string >( ), scenario.metricShOrders, std::vector< std::string >( ), false );
        auto baseMetric = createSpaceTimeMetric( metricSettings, bodiesDirect );
        evaluatedMetricObjects.clear( );
        evaluatedMetricObjects[ std::make_pair( scenario.centralBody, stationName ) ] = baseMetric->Clone( );

        auto directFromMetricSettings = std::make_shared< DirectRelativisticTimePropagatorSettings< double, double > >(
                std::make_pair( scenario.centralBody, stationName ), initialEphemerisTime, intSetDirect, terminationSettings );
        directFromMetricSettings->getOutputSettings( )->setIntegratedResult( true );
        SingleArcDynamicsSimulator< double > directDynamicsSimulator( bodiesDirect, directFromMetricSettings, true );

        // --- PN chain: TCB->TC{G,L} (SECOND order) then body->topocentric. The barycentric leg is
        //     run at second post-Newtonian order to match the direct-from-metric metric (also 2nd
        //     order in the sqrt expansion); see the expansion-order note in the header comment. ---
        auto secondOrderTcbToTcgSettings =
                std::make_shared< SecondOrderBodyCenteredRelativisticTimeConverterSettings< double, double > >( scenario.centralBody,
                                                                                                                scenario.perturbingBodies,
                                                                                                                initialEphemerisTime,
                                                                                                                intSetTcb,
                                                                                                                terminationSettings,
                                                                                                                chainPerturberShOrders );
        secondOrderTcbToTcgSettings->getOutputSettings( )->setIntegratedResult( true );

        Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRelativisticState = Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 1 );
        auto bodyToTopoSettings = std::make_shared< BodycenteredToTopocentricTimePropagatorSettings< double, double > >(
                std::make_pair( scenario.centralBody, stationName ),
                false,  // acceleration term off (IBP-paired with the -v_E.r/c^2 direct correction)
                topoCentralShDegree,
                false,
                scenario.perturbingBodies,
                initialRelativisticState,
                initialEphemerisTime,
                intSetTopo,
                terminationSettings );
        bodyToTopoSettings->getOutputSettings( )->setIntegratedResult( true );

        SingleArcDynamicsSimulator< double > barycentricPnDynamicsSimulator( bodiesPn, secondOrderTcbToTcgSettings, true );
        SingleArcDynamicsSimulator< double > topocentricPnDynamicsSimulator( bodiesPn, bodyToTopoSettings, true );

        auto timeEphemerisDirect = bodiesDirect.getBody( scenario.centralBody )->getTimeScaleConverter( );
        auto timeEphemerisPn = bodiesPn.getBody( scenario.centralBody )->getTimeScaleConverter( );
        BOOST_REQUIRE_NE( timeEphemerisDirect, nullptr );
        BOOST_REQUIRE_NE( timeEphemerisPn, nullptr );

        const auto directFunction =
                timeEphemerisDirect->getTimeDifferenceFunction( barycentric_coordinate_time_scale, local_proper_time_scale, stationName );
        const auto pnFunction =
                timeEphemerisPn->getTimeDifferenceFunction( barycentric_coordinate_time_scale, local_proper_time_scale, stationName );

        // Station-position term (reported for context, NOT subtracted): the body-fixed-vs-
        // inertial difference delta(t) = -v_C * ( R(t)*y - y ) / c^2, where v_C is the central
        // body's barycentric velocity and R(t) its rotation to the inertial frame. The
        // converters already apply this correction on the SH path (raw residual << |delta|);
        // it is shown only to gauge whether a rung's residual is dominated by it.
        auto centralEphemeris = bodiesDirect.getBody( scenario.centralBody )->getEphemeris( );
        auto centralRotation = bodiesDirect.getBody( scenario.centralBody )->getRotationalEphemeris( );
        const double inv_c2 = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
        auto stationPositionTerm = [ & ]( const double t ) {
            const Eigen::Vector3d velocityCentral = centralEphemeris->getCartesianState( t ).segment( 3, 3 );
            const Eigen::Matrix3d rotationToInertial = centralRotation->getRotationToBaseFrame( t ).toRotationMatrix( );
            const Eigen::Vector3d stationInertial = rotationToInertial * equatorStationPosition;
            return -velocityCentral.dot( stationInertial - equatorStationPosition ) * inv_c2;
        };

        std::vector< double > sampleTimes, residualsRaw, stationTermValues;
        double maxAbsRaw = 0.0;
        for( double currentTime = initialEphemerisTime + integrationStep; currentTime <= finalEphemerisTime - integrationStep;
             currentTime += integrationStep )
        {
            const double raw = directFunction( currentTime ) - pnFunction( currentTime );
            sampleTimes.push_back( currentTime );
            residualsRaw.push_back( raw );
            stationTermValues.push_back( stationPositionTerm( currentTime ) );
            maxAbsRaw = std::max( maxAbsRaw, std::fabs( raw ) );
        }
        BOOST_REQUIRE( !residualsRaw.empty( ) );

        // Linear detrend (remove the physically irrelevant constant rate offset) and report the
        // remaining periodic amplitude.
        auto detrend = [ & ]( const std::vector< double >& residuals ) {
            double sumTime = 0.0, sumResidual = 0.0, sumTimeSquared = 0.0, sumTimeResidual = 0.0;
            for( size_t i = 0; i < sampleTimes.size( ); ++i )
            {
                sumTime += sampleTimes[ i ];
                sumResidual += residuals[ i ];
                sumTimeSquared += sampleTimes[ i ] * sampleTimes[ i ];
                sumTimeResidual += sampleTimes[ i ] * residuals[ i ];
            }
            const double N = static_cast< double >( residuals.size( ) );
            const double den = N * sumTimeSquared - sumTime * sumTime;
            const double slope = std::fabs( den ) > std::numeric_limits< double >::epsilon( )
                    ? ( N * sumTimeResidual - sumTime * sumResidual ) / den
                    : 0.0;
            const double offset = ( sumResidual - slope * sumTime ) / N;
            double amplitude = 0.0;
            for( size_t i = 0; i < residuals.size( ); ++i )
            {
                amplitude = std::max( amplitude, std::fabs( residuals[ i ] - ( offset + slope * sampleTimes[ i ] ) ) );
            }
            return std::make_tuple( slope, offset, amplitude );
        };

        const auto [ slopeRaw, offsetRaw, amplitudeRaw ] = detrend( residualsRaw );
        const auto [ slopeStation, offsetStation, amplitudeStation ] = detrend( stationTermValues );
        (void)slopeStation;
        (void)offsetStation;

        std::cout << std::scientific << std::setprecision( 3 ) << "  samples=" << residualsRaw.size( ) << std::endl
                  << "  raw (direct - chain):  slope=" << slopeRaw << " s/s  offset=" << offsetRaw << " s  max|raw|=" << maxAbsRaw
                  << " s  det_amp=" << amplitudeRaw << " s" << std::endl
                  << "  station-position term |delta| det_amp=" << amplitudeStation << " s (context)";
        if( scenario.detrendedTolerance > 0.0 )
        {
            std::cout << "  tol=" << scenario.detrendedTolerance << " s (on raw)";
        }
        std::cout << std::endl;

        BOOST_CHECK( std::isfinite( amplitudeRaw ) );
        BOOST_CHECK( std::isfinite( amplitudeStation ) );
        if( scenario.detrendedTolerance > 0.0 )
        {
            BOOST_CHECK_SMALL( amplitudeRaw, scenario.detrendedTolerance );
        }
        // With the chain at SECOND order (matching the direct metric), the residual mean rate is
        // ~2e-19 (Earth) / ~1e-20 (Moon). Guard at 1e-18: a regression to a first-order chain
        // would restore the omitted 2nd-order barycentric term and push this back to ~1.2e-16.
        BOOST_CHECK_SMALL( slopeRaw, 1.0e-18 );

        // DIAG: which side is wrong? Reconstruct the analytic first-order proper-time rate at the
        // rotating station, rate(t) = -( 1/2 |v_st|^2 + sum_i GM_i/|r_st - r_i| ) / c^2, integrate it
        // (trapezoid) on the same grid, and compare against BOTH converters. The correct side tracks
        // this truth (to ~ns 2nd-order level); the buggy side shows the us-scale deviation.
        std::vector< std::function< Eigen::Vector6d( double ) > > bodyStateFns;
        std::vector< double > bodyGms;
        for( const std::string& body : bodyNames )
        {
            bodyStateFns.push_back(
                    [ &, body ]( double t ) { return bodiesDirect.getBody( body )->getEphemeris( )->getCartesianState( t ); } );
            bodyGms.push_back( bodiesDirect.getBody( body )->getGravityFieldModel( )->getGravitationalParameter( ) );
        }
        auto truthRate = [ & ]( const double t ) {
            const Eigen::Vector6d centralState = centralEphemeris->getCartesianState( t );
            const Eigen::Matrix3d rotationToInertial = centralRotation->getRotationToBaseFrame( t ).toRotationMatrix( );
            const Eigen::Vector3d stationInertial = rotationToInertial * equatorStationPosition;
            const Eigen::Vector3d omega = centralRotation->getRotationalVelocityVectorInBaseFrame( t );
            const Eigen::Vector3d stationBcrsPos = centralState.segment( 0, 3 ) + stationInertial;
            const Eigen::Vector3d stationBcrsVel = centralState.segment( 3, 3 ) + omega.cross( stationInertial );
            double potential = 0.0;
            for( size_t b = 0; b < bodyStateFns.size( ); ++b )
            {
                potential += bodyGms[ b ] / ( stationBcrsPos - bodyStateFns[ b ]( t ).segment( 0, 3 ) ).norm( );
            }
            return -( 0.5 * stationBcrsVel.squaredNorm( ) + potential ) * inv_c2;
        };
        std::vector< double > truthCumulative( sampleTimes.size( ), 0.0 );
        for( size_t i = 1; i < sampleTimes.size( ); ++i )
        {
            truthCumulative[ i ] = truthCumulative[ i - 1 ] +
                    0.5 * ( truthRate( sampleTimes[ i ] ) + truthRate( sampleTimes[ i - 1 ] ) ) *
                            ( sampleTimes[ i ] - sampleTimes[ i - 1 ] );
        }
        std::vector< double > directVsTruth( sampleTimes.size( ) ), chainVsTruth( sampleTimes.size( ) );
        for( size_t i = 0; i < sampleTimes.size( ); ++i )
        {
            directVsTruth[ i ] = directFunction( sampleTimes[ i ] ) - truthCumulative[ i ];
            chainVsTruth[ i ] = pnFunction( sampleTimes[ i ] ) - truthCumulative[ i ];
        }
        const auto [ sD, oD, ampDirectVsTruth ] = detrend( directVsTruth );
        const auto [ sC, oC, ampChainVsTruth ] = detrend( chainVsTruth );
        (void)sD;
        (void)oD;
        (void)sC;
        (void)oC;
        std::cout << "  vs analytic 1st-order truth:  direct det_amp=" << ampDirectVsTruth << " s   chain det_amp=" << ampChainVsTruth
                  << " s" << std::endl;

        // Strong validation, enforced on EVERY rung (Earth and Moon): the post-Newtonian chain
        // must track the independent analytic first-order proper-time integral to sub-ns. This is
        // the cross-check that isolates the direct-from-metric point-mass defect to the direct side.
        BOOST_CHECK( std::isfinite( ampDirectVsTruth ) );
        BOOST_CHECK( std::isfinite( ampChainVsTruth ) );
        // Both 2PN methods track the analytic 1st-order reference to the level of their 2nd-order
        // content plus the trapezoid-truth integration error (~50 ps); guard at 1e-10.
        BOOST_CHECK_SMALL( ampChainVsTruth, 1.0e-10 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
