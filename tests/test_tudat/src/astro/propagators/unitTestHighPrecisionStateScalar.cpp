/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <type_traits>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <tudat/config.hpp>

#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/math/integrators/rungeKutta4Integrator.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/math/root_finders/newtonRaphson.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/propagationResults.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"

namespace tudat
{
namespace unit_tests
{

using Scalar = HighPrecisionStateScalar;
using FixedState = Eigen::Matrix< Scalar, 6, 1 >;
using DynamicState = Eigen::Matrix< Scalar, Eigen::Dynamic, 1 >;
using StateHistory = std::map< double, DynamicState >;

#if TUDAT_HIGH_PRECISION_STATE_SCALAR_IS_CPP_BIN_FLOAT_QUAD
static_assert( std::is_same_v< Scalar, boost::multiprecision::cpp_bin_float_quad > );
#else
static_assert( std::is_same_v< Scalar, long double > );
#endif

static_assert( is_state_scalar< Scalar >::value );
static_assert( std::is_same_v< typename Eigen::Vector6ld::Scalar, Scalar > );
static_assert( std::is_same_v< typename Eigen::VectorXld::Scalar, Scalar > );

constexpr bool timeHasExtendedLongDoublePrecision = std::numeric_limits< long double >::digits > std::numeric_limits< double >::digits;

template< typename ValueType >
ValueType scalarFromDecimalString( const char* value )
{
    if constexpr( std::is_same_v< ValueType, double > )
    {
        return std::strtod( value, nullptr );
    }
    else if constexpr( std::is_same_v< ValueType, long double > )
    {
        return std::strtold( value, nullptr );
    }
    else
    {
        return ValueType( value );
    }
}

template< typename ValueType >
ValueType runRepresentativeInterpolationCase( )
{
    using State = Eigen::Matrix< ValueType, 1, 1 >;
    std::map< double, State > history;
    history[ 0.0 ]( 0 ) = scalarFromDecimalString< ValueType >( "1.25" );
    history[ 1.0 ]( 0 ) = scalarFromDecimalString< ValueType >( "2.75" );

    interpolators::LinearInterpolator< double, State > interpolator( history );
    return interpolator.interpolate( 0.25 )( 0 );
}

template< typename ValueType >
ValueType runRepresentativeIntegrationCase( )
{
    using State = Eigen::Matrix< ValueType, Eigen::Dynamic, 1 >;
    State initialState( 2 );
    initialState( 0 ) = scalarFromDecimalString< ValueType >( "1.25" );
    initialState( 1 ) = scalarFromDecimalString< ValueType >( "-0.5" );

    const auto stateDerivative = []( const double, const State& state ) {
        State derivative( 2 );
        derivative( 0 ) = state( 1 );
        derivative( 1 ) = -state( 0 );
        return derivative;
    };

    numerical_integrators::RungeKutta4Integrator< double, State, State, double > integrator( stateDerivative, 0.0, initialState, 0.125 );
    return integrator.performIntegrationStep( 0.125 )( 0 );
}

template< typename ValueType >
ValueType getAbsoluteValue( const ValueType& value )
{
    using std::abs;
    return abs( value );
}

template< typename ValueType >
ValueType getOscillatorRk4Error( const int numberOfSteps )
{
    using State = Eigen::Matrix< ValueType, 2, 1 >;

    State initialState;
    initialState << scalarFromDecimalString< ValueType >( "0" ), scalarFromDecimalString< ValueType >( "1" );

    const auto stateDerivative = []( const double, const State& state ) {
        State derivative;
        derivative << state( 1 ), -state( 0 );
        return derivative;
    };

    const double stepSize = 1.0 / static_cast< double >( numberOfSteps );
    numerical_integrators::RungeKutta4Integrator< double, State, State, double > integrator( stateDerivative, 0.0, initialState, stepSize );
    State finalState = initialState;
    for( int i = 0; i < numberOfSteps; ++i )
    {
        finalState = integrator.performIntegrationStep( stepSize );
    }

    using std::cos;
    using std::sin;
    const ValueType finalTime = scalarFromDecimalString< ValueType >( "1" );
    const ValueType positionError = getAbsoluteValue( finalState( 0 ) - sin( finalTime ) );
    const ValueType velocityError = getAbsoluteValue( finalState( 1 ) - cos( finalTime ) );
    return positionError > velocityError ? positionError : velocityError;
}

template< typename ValueType >
Eigen::Matrix< ValueType, 3, 1 > getPointMassOrbitDynamicsSimulatorError( const int numberOfSteps )
{
    using State = Eigen::Matrix< ValueType, 6, 1 >;

    const double gravitationalParameterAsDouble = 398600441800000.0;
    const ValueType gravitationalParameter = static_cast< ValueType >( gravitationalParameterAsDouble );
    const ValueType semiMajorAxis = scalarFromDecimalString< ValueType >( "8000000" );
    const ValueType eccentricity = scalarFromDecimalString< ValueType >( "0.1" );
    const ValueType pi = mathematical_constants::getPi< ValueType >( );

    State initialKeplerianState;
    initialKeplerianState << semiMajorAxis, eccentricity, pi / static_cast< ValueType >( 6 ), pi / static_cast< ValueType >( 5 ),
            pi / static_cast< ValueType >( 7 ), pi / static_cast< ValueType >( 9 );

    const State initialCartesianState =
            orbital_element_conversions::convertKeplerianToCartesianElements< ValueType >( initialKeplerianState, gravitationalParameter );

    const double propagationTime = 3600.0;
    const Time initialTime( 0.0L );
    const Time stepSize( static_cast< long double >( propagationTime ) / static_cast< long double >( numberOfSteps ) );

    simulation_setup::SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody< ValueType, Time >( "Earth", false );
    bodies.createEmptyBody< ValueType, Time >( "Vehicle", false );
    bodies.at( "Earth" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
    bodies.at( "Earth" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( gravitationalParameterAsDouble ) );
    bodies.processBodyFrameDefinitions< ValueType, Time >( );

    simulation_setup::SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Vehicle" ][ "Earth" ].push_back(
            std::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    const std::vector< std::string > bodiesToPropagate{ "Vehicle" };
    const std::vector< std::string > centralBodies{ "Earth" };
    const basic_astrodynamics::AccelerationMap accelerationModels =
            simulation_setup::createAccelerationModelsMap( bodies, accelerationSettings, bodiesToPropagate, centralBodies );

    const auto integratorSettings = std::make_shared< numerical_integrators::IntegratorSettings< Time > >(
            numerical_integrators::rungeKutta4, initialTime, stepSize );
    const auto outputSettings = std::make_shared< propagators::SingleArcPropagatorProcessingSettings >( false, false, numberOfSteps );
    const auto propagatorSettings = std::make_shared< propagators::TranslationalStatePropagatorSettings< ValueType, Time > >(
            centralBodies,
            accelerationModels,
            bodiesToPropagate,
            initialCartesianState,
            initialTime,
            integratorSettings,
            std::make_shared< propagators::PropagationTimeTerminationSettings >( propagationTime, true ),
            propagators::cowell,
            std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > >( ),
            outputSettings );

    propagators::SingleArcDynamicsSimulator< ValueType, Time > dynamicsSimulator( bodies, propagatorSettings );
    const auto propagationResults = dynamicsSimulator.getSingleArcPropagationResults( );
    const auto& stateHistory = propagationResults->getEquationsOfMotionNumericalSolution( );

    using StateHistoryType = std::decay_t< decltype( stateHistory ) >;
    static_assert( std::is_same_v< typename StateHistoryType::mapped_type::Scalar, ValueType > );
    if( stateHistory.empty( ) )
    {
        throw std::runtime_error( "Point-mass orbit propagation returned no states." );
    }

    const State finalNumericalState = stateHistory.rbegin( )->second;
    const ValueType propagatedDuration =
            static_cast< ValueType >( static_cast< long double >( stateHistory.rbegin( )->first - stateHistory.begin( )->first ) );

    const auto rootFinderTermination = []( const ValueType currentRoot,
                                           const ValueType previousRoot,
                                           const ValueType currentFunctionValue,
                                           const ValueType,
                                           const unsigned int iteration ) {
        using std::abs;
        return currentFunctionValue == static_cast< ValueType >( 0 ) ||
                abs( currentRoot - previousRoot ) < static_cast< ValueType >( 10 ) * std::numeric_limits< ValueType >::epsilon( ) ||
                iteration > 100;
    };
    const std::shared_ptr< root_finders::RootFinder< ValueType > > rootFinder =
            std::make_shared< root_finders::NewtonRaphson< ValueType > >( rootFinderTermination );

    const State analyticalKeplerianState = orbital_element_conversions::propagateKeplerOrbit< ValueType >(
            initialKeplerianState, propagatedDuration, gravitationalParameter, rootFinder );
    const State analyticalCartesianState = orbital_element_conversions::convertKeplerianToCartesianElements< ValueType >(
            analyticalKeplerianState, gravitationalParameter );

    const State stateError = finalNumericalState - analyticalCartesianState;
    const ValueType positionError = stateError.template segment< 3 >( 0 ).norm( );
    const ValueType velocityError = stateError.template segment< 3 >( 3 ).norm( );
    const ValueType characteristicVelocity = sqrt( gravitationalParameter / semiMajorAxis );
    const ValueType normalizedPositionError = positionError / semiMajorAxis;
    const ValueType normalizedVelocityError = velocityError / characteristicVelocity;

    Eigen::Matrix< ValueType, 3, 1 > errorMetrics;
    errorMetrics << ( normalizedPositionError > normalizedVelocityError ? normalizedPositionError : normalizedVelocityError ),
            positionError, velocityError;
    return errorMetrics;
}

struct PropagatedObservationEnvironment {
    simulation_setup::SystemOfBodies bodies;
    std::shared_ptr< propagators::SingleArcSimulationResults< Scalar, Time > > propagationResults;
};

PropagatedObservationEnvironment createPropagatedObservationEnvironment( )
{
    const double gravitationalParameterAsDouble = 1.32712440018e20;
    const Scalar gravitationalParameter = static_cast< Scalar >( gravitationalParameterAsDouble );
    const Scalar pi = mathematical_constants::getPi< Scalar >( );

    FixedState initialKeplerianState;
    initialKeplerianState << scalarFromDecimalString< Scalar >( "778500000000" ), scalarFromDecimalString< Scalar >( "0.0489" ),
            pi * scalarFromDecimalString< Scalar >( "1.3" ) / static_cast< Scalar >( 180 ), pi / static_cast< Scalar >( 5 ),
            pi / static_cast< Scalar >( 7 ), pi / static_cast< Scalar >( 9 );
    const FixedState initialCartesianState =
            orbital_element_conversions::convertKeplerianToCartesianElements< Scalar >( initialKeplerianState, gravitationalParameter );

    simulation_setup::SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody< Scalar, Time >( "Sun", false );
    bodies.createEmptyBody< Scalar, Time >( "Earth", false );
    bodies.createEmptyBody< Scalar, Time >( "Jupiter", false );
    bodies.at( "Sun" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
    bodies.at( "Sun" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( gravitationalParameterAsDouble ) );
    Eigen::Vector6d earthState = Eigen::Vector6d::Zero( );
    earthState( 0 ) = 149597870700.0;
    bodies.at( "Earth" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >( earthState, "SSB", "J2000" ) );
    bodies.at( "Earth" )->setRotationalEphemeris(
            std::make_shared< ephemerides::ConstantRotationalEphemeris >( Eigen::Quaterniond::Identity( ), "J2000", "IAU_Earth" ) );
    const auto jupiterVehicleSystems = std::make_shared< system_models::VehicleSystems >( );
    jupiterVehicleSystems->setDefaultTransponderTurnaroundRatio( );
    bodies.at( "Jupiter" )->setVehicleSystems( jupiterVehicleSystems );
    bodies.processBodyFrameDefinitions< Scalar, Time >( );

    const std::string stationName = "QuadDsnStation";
    simulation_setup::createGroundStation( bodies.at( "Earth" ),
                                           stationName,
                                           ( Eigen::Vector3d( ) << 6378137.0, 0.0, 0.0 ).finished( ),
                                           coordinate_conversions::cartesian_position );
    const std::shared_ptr< ground_stations::StationFrequencyInterpolator > rampedFrequencyInterpolator =
            std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >( std::vector< Time >{ Time( -10000.0L ) },
                                                                                       std::vector< Time >{ Time( 100000.0L ) },
                                                                                       std::vector< double >{ 0.25 },
                                                                                       std::vector< double >{ 7.2e9 } );
    bodies.at( "Earth" )->getGroundStation( stationName )->setTransmittingFrequencyCalculator( rampedFrequencyInterpolator );

    simulation_setup::SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Jupiter" ][ "Sun" ].push_back(
            std::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    const std::vector< std::string > bodiesToPropagate{ "Jupiter" };
    const std::vector< std::string > centralBodies{ "Sun" };
    const basic_astrodynamics::AccelerationMap accelerationModels =
            simulation_setup::createAccelerationModelsMap( bodies, accelerationSettings, bodiesToPropagate, centralBodies );

    const Time initialTime( 0.0L );
    const Time stepSize( 60.0L );
    const auto integratorSettings = std::make_shared< numerical_integrators::IntegratorSettings< Time > >(
            numerical_integrators::rungeKutta4, initialTime, stepSize );
    const auto outputSettings = std::make_shared< propagators::SingleArcPropagatorProcessingSettings >(
            false, true, 1, TUDAT_NAN, std::make_shared< propagators::PropagationPrintSettings >( ) );
    const auto propagatorSettings = std::make_shared< propagators::TranslationalStatePropagatorSettings< Scalar, Time > >(
            centralBodies,
            accelerationModels,
            bodiesToPropagate,
            initialCartesianState,
            initialTime,
            integratorSettings,
            std::make_shared< propagators::PropagationTimeTerminationSettings >( 86400.0, true ),
            propagators::cowell,
            std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > >( ),
            outputSettings );

    propagators::SingleArcDynamicsSimulator< Scalar, Time > dynamicsSimulator( bodies, propagatorSettings );
    return { bodies, dynamicsSimulator.getSingleArcPropagationResults( ) };
}

Scalar computeEarthJupiterRange( const std::shared_ptr< ephemerides::Ephemeris >& earthEphemeris,
                                 const std::shared_ptr< ephemerides::Ephemeris >& jupiterEphemeris,
                                 const Time& time )
{
    return ( jupiterEphemeris->getTemplatedStateFromEphemeris< Scalar, Time >( time ) -
             earthEphemeris->getTemplatedStateFromEphemeris< Scalar, Time >( time ) )
            .template segment< 3 >( 0 )
            .norm( );
}

Scalar computeTwoWayRangeReference( const std::shared_ptr< ephemerides::Ephemeris >& earthEphemeris,
                                    const std::shared_ptr< ephemerides::Ephemeris >& jupiterEphemeris,
                                    const Time& receptionTime )
{
    const Scalar speedOfLight = physical_constants::getSpeedOfLight< Scalar >( );
    Time reflectionTime = receptionTime;
    for( int i = 0; i < 50; ++i )
    {
        const Scalar range = computeEarthJupiterRange( earthEphemeris, jupiterEphemeris, reflectionTime );
        const Time updatedReflectionTime = receptionTime - static_cast< long double >( range / speedOfLight );
        if( updatedReflectionTime == reflectionTime )
        {
            break;
        }
        reflectionTime = updatedReflectionTime;
    }

    return static_cast< Scalar >( 2 ) * computeEarthJupiterRange( earthEphemeris, jupiterEphemeris, reflectionTime );
}

BOOST_AUTO_TEST_SUITE( test_high_precision_state_scalar )

BOOST_AUTO_TEST_CASE( testEigenStateVectorsAndArithmetic )
{
    const Scalar one = scalarFromDecimalString< Scalar >( "1" );
    const Scalar increment = scalarFromDecimalString< Scalar >( "0.125" );

    FixedState fixedState = FixedState::Constant( one );
    DynamicState dynamicState = DynamicState::Constant( 6, increment );
    fixedState += dynamicState;

    BOOST_CHECK( fixedState( 0 ) == one + increment );
    BOOST_CHECK_EQUAL( fixedState.rows( ), 6 );
    BOOST_CHECK_EQUAL( dynamicState.rows( ), 6 );
}

BOOST_AUTO_TEST_CASE( testStateInterpolationStorageAndRetrieval )
{
    const Scalar preciseValue = scalarFromDecimalString< Scalar >( "1.000000000000000000000000000000001" );

    std::map< double, FixedState > fixedStateHistory;
    for( int i = 0; i < 8; ++i )
    {
        fixedStateHistory[ static_cast< double >( i ) ] = FixedState::Constant( preciseValue );
    }

    const auto stateInterpolator = propagators::createStateInterpolator( fixedStateHistory );
    const FixedState interpolatedAtNode = stateInterpolator->interpolate( 3.0 );
    BOOST_CHECK( interpolatedAtNode( 0 ) == fixedStateHistory.at( 3.0 )( 0 ) );

    const auto tabulatedEphemeris =
            std::make_shared< ephemerides::TabulatedCartesianEphemeris< Scalar, double > >( stateInterpolator, "SSB", "J2000" );
    BOOST_CHECK( tabulatedEphemeris->getCartesianLongState( 3.0 )( 0 ) == fixedStateHistory.at( 3.0 )( 0 ) );

    simulation_setup::Body body;
    body.setLongState( interpolatedAtNode );
    BOOST_CHECK( body.getLongState( )( 0 ) == fixedStateHistory.at( 3.0 )( 0 ) );

#if TUDAT_HIGH_PRECISION_STATE_SCALAR_IS_CPP_BIN_FLOAT_QUAD
    const Scalar residual = body.getLongState( )( 0 ) - static_cast< Scalar >( 1 );
    BOOST_CHECK( residual > scalarFromDecimalString< Scalar >( "1e-34" ) );
    BOOST_CHECK( residual < scalarFromDecimalString< Scalar >( "2e-33" ) );
#endif
}

BOOST_AUTO_TEST_CASE( testResetTabulatedEphemerisResolvesPicosecondOffset )
{
    const Time firstEpoch( 100000000.0L );
    const Scalar basePosition = scalarFromDecimalString< Scalar >( "1000000000000" );
    const Scalar velocity = scalarFromDecimalString< Scalar >( "1000000" );

    std::map< Time, FixedState > stateHistory;
    for( int i = 0; i < 8; ++i )
    {
        FixedState state = FixedState::Zero( );
        state( 0 ) = basePosition + velocity * static_cast< Scalar >( i );
        state( 3 ) = velocity;
        stateHistory[ firstEpoch + static_cast< long double >( i ) ] = state;
    }

    simulation_setup::SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody( "Vehicle", false );
    propagators::addEmptyTabulatedEphemeris< Scalar, Time >( bodies, "Vehicle", "SSB" );

    // Exercise the same reset function used when propagated results are selected
    // as the body's post-propagation ephemeris.
    propagators::resetIntegratedEphemerisOfBody( bodies, stateHistory, "Vehicle" );

    const Time nodeTime = firstEpoch + 3.0L;
    const Time queryTime = nodeTime + 1.0e-12L;
    const Scalar representedTimeIncrement = static_cast< Scalar >( static_cast< long double >( queryTime - nodeTime ) );
    const Scalar expectedPositionChange = velocity * representedTimeIncrement;

    const auto resetEphemeris = std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< Scalar, Time > >(
            bodies.at( "Vehicle" )->getEphemeris( ) );
    BOOST_REQUIRE( resetEphemeris != nullptr );

    const Scalar nodePosition = resetEphemeris->getCartesianLongStateFromExtendedTime( nodeTime )( 0 );
    const Scalar offsetPosition = resetEphemeris->getCartesianLongStateFromExtendedTime( queryTime )( 0 );
    const Scalar interpolatedPositionChange = offsetPosition - nodePosition;

    BOOST_CHECK( representedTimeIncrement > static_cast< Scalar >( 0 ) );
    BOOST_CHECK( interpolatedPositionChange > static_cast< Scalar >( 0 ) );
    BOOST_CHECK( getAbsoluteValue( interpolatedPositionChange - expectedPositionChange ) <
                 expectedPositionChange / scalarFromDecimalString< Scalar >( "4" ) );

#if TUDAT_HIGH_PRECISION_STATE_SCALAR_IS_CPP_BIN_FLOAT_QUAD
    // Quad state arithmetic should retain substantially more than merely the
    // fact that the position changed.
    BOOST_CHECK( getAbsoluteValue( interpolatedPositionChange - expectedPositionChange ) < scalarFromDecimalString< Scalar >( "1e-18" ) );

    // At this state magnitude, double cannot resolve the expected change
    // accurately. Depending on the platform and interpolation rounding, it
    // may remain unchanged or jump by one double-precision spacing.
    std::map< Time, Eigen::Vector6d > doubleStateHistory;
    for( const auto& stateEntry : stateHistory )
    {
        doubleStateHistory[ stateEntry.first ] = stateEntry.second.template cast< double >( );
    }
    const auto doubleInterpolator = propagators::createStateInterpolator( doubleStateHistory );
    ephemerides::TabulatedCartesianEphemeris< double, Time > doubleEphemeris( doubleInterpolator, "SSB", "J2000" );
    const double doublePositionChange = doubleEphemeris.getCartesianStateFromExtendedTime( queryTime )( 0 ) -
            doubleEphemeris.getCartesianStateFromExtendedTime( nodeTime )( 0 );
    const Scalar doublePositionError = getAbsoluteValue( static_cast< Scalar >( doublePositionChange ) - expectedPositionChange );
    BOOST_CHECK( doublePositionError > expectedPositionChange / scalarFromDecimalString< Scalar >( "2" ) );
    BOOST_CHECK( getAbsoluteValue( interpolatedPositionChange - expectedPositionChange ) * scalarFromDecimalString< Scalar >( "1e9" ) <
                 doublePositionError );
    std::cout << "Expected 1 ps position change [m]: " << expectedPositionChange
              << ", quad position change [m]: " << interpolatedPositionChange << ", double position change [m]: " << doublePositionChange
              << ", double position-change error [m]: " << doublePositionError << std::endl;
#endif
}

BOOST_AUTO_TEST_CASE( testPropagationInterfacesAndStateHistory )
{
    const Scalar initialValue = scalarFromDecimalString< Scalar >( "1.000000000000000000000000000000001" );
    DynamicState initialState( 1 );
    initialState( 0 ) = initialValue;

    const auto derivative = []( const double, const DynamicState& state ) { return DynamicState::Zero( state.rows( ) ); };
    const auto terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( 1.0 );
    propagators::CustomStatePropagatorSettings< Scalar, double > propagatorSettings( derivative, initialState, terminationSettings );
    BOOST_CHECK( propagatorSettings.getInitialStates( )( 0 ) == initialValue );

    StateHistory stateHistory;
    stateHistory[ 0.0 ] = initialState;
    stateHistory[ 1.0 ] = initialState;

    const auto outputSettings = std::make_shared< propagators::SingleArcPropagatorProcessingSettings >( );
    const auto rawSolutionConversion = []( StateHistory& processedHistory, const StateHistory& rawHistory ) {
        processedHistory = rawHistory;
    };
    propagators::SingleArcSimulationResults< Scalar, double > simulationResults( {}, outputSettings, rawSolutionConversion, nullptr );
    simulationResults.reset( stateHistory,
                             {},
                             {},
                             {},
                             std::make_shared< propagators::PropagationTerminationDetails >( propagators::termination_condition_reached ) );
    simulationResults.finalizePropagation( {} );

    const auto& retrievedHistory = simulationResults.getEquationsOfMotionNumericalSolution( );
    BOOST_CHECK( retrievedHistory.at( 1.0 )( 0 ) == initialValue );

    const Scalar propagatedValue = runRepresentativeIntegrationCase< Scalar >( );
    BOOST_CHECK( propagatedValue != initialValue );
}

BOOST_AUTO_TEST_CASE( testRk4FourthOrderConvergence )
{
    const std::vector< int > coarseStepCounts = { 16, 32, 64, 128 };
    std::vector< Scalar > configuredErrors;
    for( const int stepCount : coarseStepCounts )
    {
        configuredErrors.push_back( getOscillatorRk4Error< Scalar >( stepCount ) );
    }

    for( std::size_t i = 1; i < configuredErrors.size( ); ++i )
    {
        const Scalar errorReduction = configuredErrors.at( i - 1 ) / configuredErrors.at( i );
        BOOST_CHECK( errorReduction > scalarFromDecimalString< Scalar >( "14" ) );
        BOOST_CHECK( errorReduction < scalarFromDecimalString< Scalar >( "18" ) );
    }

#if TUDAT_HIGH_PRECISION_STATE_SCALAR_IS_CPP_BIN_FLOAT_QUAD
    const std::vector< int > fineStepCounts = { 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576 };
    std::vector< Scalar > quadErrors;
    std::vector< double > doubleErrors;
    for( const int stepCount : fineStepCounts )
    {
        quadErrors.push_back( getOscillatorRk4Error< Scalar >( stepCount ) );
        doubleErrors.push_back( getOscillatorRk4Error< double >( stepCount ) );
    }

    for( std::size_t i = 1; i < quadErrors.size( ); ++i )
    {
        const Scalar errorReduction = quadErrors.at( i - 1 ) / quadErrors.at( i );
        BOOST_CHECK( errorReduction > scalarFromDecimalString< Scalar >( "14" ) );
        BOOST_CHECK( errorReduction < scalarFromDecimalString< Scalar >( "18" ) );
    }

    // At these smaller steps the double result has stopped showing
    // fourth-order improvement, while quad continues to converge.
    const double finalDoubleErrorReduction = doubleErrors.at( doubleErrors.size( ) - 2 ) / doubleErrors.back( );
    BOOST_CHECK_LT( finalDoubleErrorReduction, 4.0 );

    BOOST_CHECK( quadErrors.back( ) * scalarFromDecimalString< Scalar >( "1e10" ) < static_cast< Scalar >( doubleErrors.back( ) ) );

    BOOST_TEST_MESSAGE( "RK4 finest step count: " << fineStepCounts.back( ) << ", quad error: " << quadErrors.back( )
                                                  << ", double error: " << doubleErrors.back( ) );
#endif
}

BOOST_AUTO_TEST_CASE( testPointMassOrbitDynamicsSimulatorNoiseFloor )
{
#if TUDAT_HIGH_PRECISION_STATE_SCALAR_IS_CPP_BIN_FLOAT_QUAD
    // On platforms with a 53-bit long double, this end-to-end comparison
    // enters its platform-dependent time/reference cancellation regime one
    // halving earlier. Stop before that regime and retain the same
    // fourth-order checks for every step that is run.
    const std::vector< int > stepCounts = timeHasExtendedLongDoublePrecision
            ? std::vector< int >{ 256, 512, 1024, 2048, 4096, 8192, 16384, 32768 }
            : std::vector< int >{ 256, 512, 1024, 2048, 4096, 8192, 16384 };

    std::vector< Eigen::Matrix< Scalar, 3, 1 > > quadErrors;
    std::vector< Eigen::Matrix< double, 3, 1 > > doubleErrors;
    for( const int stepCount : stepCounts )
    {
        quadErrors.push_back( getPointMassOrbitDynamicsSimulatorError< Scalar >( stepCount ) );
        doubleErrors.push_back( getPointMassOrbitDynamicsSimulatorError< double >( stepCount ) );
    }

    for( std::size_t i = 1; i < quadErrors.size( ); ++i )
    {
        const Scalar quadErrorReduction = quadErrors.at( i - 1 )( 0 ) / quadErrors.at( i )( 0 );
        const double doubleErrorReduction = doubleErrors.at( i - 1 )( 0 ) / doubleErrors.at( i )( 0 );
        BOOST_CHECK( quadErrorReduction > scalarFromDecimalString< Scalar >( "12" ) );
        BOOST_CHECK( quadErrorReduction < scalarFromDecimalString< Scalar >( "20" ) );

        std::cout << "Point-mass orbit step count: " << stepCounts.at( i ) << ", quad normalized error: " << quadErrors.at( i )( 0 )
                  << ", double normalized error: " << doubleErrors.at( i )( 0 ) << ", quad RK4 reduction: " << quadErrorReduction
                  << ", double RK4 reduction: " << doubleErrorReduction << std::endl;
    }

    const double finalDoubleErrorReduction = doubleErrors.at( doubleErrors.size( ) - 2 )( 0 ) / doubleErrors.back( )( 0 );
    BOOST_CHECK_LT( finalDoubleErrorReduction, 8.0 );

    BOOST_CHECK( quadErrors.back( )( 0 ) * scalarFromDecimalString< Scalar >( "10" ) < static_cast< Scalar >( doubleErrors.back( )( 0 ) ) );

    BOOST_TEST_MESSAGE( "Point-mass orbit finest step count: " << stepCounts.back( )
                                                               << ", quad position error [m]: " << quadErrors.back( )( 1 )
                                                               << ", quad velocity error [m/s]: " << quadErrors.back( )( 2 )
                                                               << ", double position error [m]: " << doubleErrors.back( )( 1 )
                                                               << ", double velocity error [m/s]: " << doubleErrors.back( )( 2 ) );
#endif
}

BOOST_AUTO_TEST_CASE( testResetPropagatedEphemerisInQuadRangeObservations )
{
#if TUDAT_HIGH_PRECISION_STATE_SCALAR_IS_CPP_BIN_FLOAT_QUAD
    PropagatedObservationEnvironment environment = createPropagatedObservationEnvironment( );
    const auto& stateHistory = environment.propagationResults->getEquationsOfMotionNumericalSolution( );
    BOOST_REQUIRE( !stateHistory.empty( ) );

    const auto jupiterEphemeris = std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< Scalar, Time > >(
            environment.bodies.at( "Jupiter" )->getEphemeris( ) );
    BOOST_REQUIRE( jupiterEphemeris != nullptr );
    const auto earthEphemeris = environment.bodies.at( "Earth" )->getEphemeris( );

    // Verify that the simulator, rather than this test, reset Jupiter's
    // ephemeris from the quad propagation results.
    const FixedState finalStoredState = jupiterEphemeris->getCartesianLongStateFromExtendedTime( stateHistory.rbegin( )->first );
    BOOST_CHECK( ( finalStoredState - stateHistory.rbegin( )->second ).norm( ) == static_cast< Scalar >( 0 ) );

    const auto convergenceCriteria = std::make_shared< observation_models::LightTimeConvergenceCriteria >( );
    BOOST_CHECK( observation_models::getDefaultLightTimeTolerance< Scalar >( ) < scalarFromDecimalString< Scalar >( "1e-25" ) );

    observation_models::LinkEnds oneWayLinkEnds;
    oneWayLinkEnds[ observation_models::transmitter ] = observation_models::linkEndId( "Earth" );
    oneWayLinkEnds[ observation_models::receiver ] = observation_models::linkEndId( "Jupiter" );
    const auto oneWaySettings = std::make_shared< observation_models::ObservationModelSettings >(
            observation_models::one_way_range,
            observation_models::LinkDefinition( oneWayLinkEnds ),
            std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > >( ),
            nullptr,
            convergenceCriteria );
    const auto oneWayRangeModel =
            observation_models::ObservationModelCreator< 1, Scalar, Time >::createObservationModel( oneWaySettings, environment.bodies );
    BOOST_REQUIRE( ( std::dynamic_pointer_cast< observation_models::OneWayRangeObservationModel< Scalar, Time > >( oneWayRangeModel ) !=
                     nullptr ) );

    observation_models::LinkEnds nWayLinkEnds;
    nWayLinkEnds[ observation_models::transmitter ] = observation_models::linkEndId( "Earth" );
    nWayLinkEnds[ observation_models::reflector1 ] = observation_models::linkEndId( "Jupiter" );
    nWayLinkEnds[ observation_models::receiver ] = observation_models::linkEndId( "Earth" );
    const auto nWaySettings = std::make_shared< observation_models::NWayRangeObservationModelSettings >(
            observation_models::LinkDefinition( nWayLinkEnds ),
            std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > >( ),
            nullptr,
            convergenceCriteria );
    const auto nWayRangeModel =
            observation_models::ObservationModelCreator< 1, Scalar, Time >::createObservationModel( nWaySettings, environment.bodies );
    BOOST_REQUIRE(
            ( std::dynamic_pointer_cast< observation_models::NWayRangeObservationModel< Scalar, Time > >( nWayRangeModel ) != nullptr ) );

    const Time observationTime( 43200.125L );
    const Time offsetObservationTime = observationTime + 1.0e-12L;
    const Scalar representedTimeIncrement = static_cast< Scalar >( static_cast< long double >( offsetObservationTime - observationTime ) );
    BOOST_REQUIRE( representedTimeIncrement > static_cast< Scalar >( 0 ) );

    const auto computeRange = []( const auto& observationModel, const Time& time, const observation_models::LinkEndType referenceLinkEnd ) {
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        return observationModel->computeIdealObservationsWithLinkEndData(
                time, referenceLinkEnd, linkEndTimes, linkEndStates, nullptr )( 0 );
    };

    const Scalar oneWayRange = computeRange( oneWayRangeModel, observationTime, observation_models::receiver );
    const Scalar offsetOneWayRange = computeRange( oneWayRangeModel, offsetObservationTime, observation_models::receiver );
    const Scalar oneWayReference = computeEarthJupiterRange( earthEphemeris, jupiterEphemeris, observationTime );
    const Scalar offsetOneWayReference = computeEarthJupiterRange( earthEphemeris, jupiterEphemeris, offsetObservationTime );

    const Scalar nWayRange = computeRange( nWayRangeModel, observationTime, observation_models::receiver );
    const Scalar offsetNWayRange = computeRange( nWayRangeModel, offsetObservationTime, observation_models::receiver );
    const Scalar nWayReference = computeTwoWayRangeReference( earthEphemeris, jupiterEphemeris, observationTime );
    const Scalar offsetNWayReference = computeTwoWayRangeReference( earthEphemeris, jupiterEphemeris, offsetObservationTime );

    const Scalar oneWayRangeChange = offsetOneWayRange - oneWayRange;
    const Scalar expectedOneWayRangeChange = offsetOneWayReference - oneWayReference;
    const Scalar nWayRangeChange = offsetNWayRange - nWayRange;
    const Scalar expectedNWayRangeChange = offsetNWayReference - nWayReference;

    const Scalar nWayRangeTolerance =
            timeHasExtendedLongDoublePrecision ? scalarFromDecimalString< Scalar >( "1e-11" ) : scalarFromDecimalString< Scalar >( "1e-9" );
    std::cout << "long double digits: " << std::numeric_limits< long double >::digits << ", sizeof(long double): " << sizeof( long double )
              << ", one-way range reference error [m]: " << getAbsoluteValue( oneWayRange - oneWayReference )
              << ", n-way range reference error [m]: " << getAbsoluteValue( nWayRange - nWayReference )
              << ", one-way range-change error [m]: " << getAbsoluteValue( oneWayRangeChange - expectedOneWayRangeChange )
              << ", n-way range-change error [m]: " << getAbsoluteValue( nWayRangeChange - expectedNWayRangeChange )
              << ", n-way range tolerance [m]: " << nWayRangeTolerance << std::endl;
    BOOST_CHECK( getAbsoluteValue( oneWayRange - oneWayReference ) < scalarFromDecimalString< Scalar >( "1e-18" ) );
    BOOST_CHECK( getAbsoluteValue( nWayRange - nWayReference ) < nWayRangeTolerance );
    BOOST_CHECK( oneWayRangeChange != static_cast< Scalar >( 0 ) );
    BOOST_CHECK( nWayRangeChange != static_cast< Scalar >( 0 ) );
    BOOST_CHECK( getAbsoluteValue( oneWayRangeChange - expectedOneWayRangeChange ) < scalarFromDecimalString< Scalar >( "1e-18" ) );
    BOOST_CHECK( getAbsoluteValue( nWayRangeChange - expectedNWayRangeChange ) < nWayRangeTolerance );

    // The test remains entirely quad precision. Double spacing is used only as
    // a benchmark for the resolution unavailable at Jupiter range.
    const double oneWayDoubleSpacing = std::nextafter( static_cast< double >( oneWayRange ), std::numeric_limits< double >::infinity( ) ) -
            static_cast< double >( oneWayRange );
    const double nWayDoubleSpacing = std::nextafter( static_cast< double >( nWayRange ), std::numeric_limits< double >::infinity( ) ) -
            static_cast< double >( nWayRange );
    BOOST_CHECK( getAbsoluteValue( oneWayRangeChange ) * static_cast< Scalar >( 1000 ) < static_cast< Scalar >( oneWayDoubleSpacing ) );
    BOOST_CHECK( getAbsoluteValue( nWayRangeChange ) * static_cast< Scalar >( 1000 ) < static_cast< Scalar >( nWayDoubleSpacing ) );

    BOOST_TEST_MESSAGE( "Represented epoch increment [s]: " << representedTimeIncrement << ", one-way range change [m]: "
                                                            << oneWayRangeChange << ", one-way double spacing [m]: " << oneWayDoubleSpacing
                                                            << ", n-way range change [m]: " << nWayRangeChange
                                                            << ", n-way double spacing [m]: " << nWayDoubleSpacing );

    const std::string stationName = "QuadDsnStation";
    const Eigen::Vector3d stationPosition =
            environment.bodies.at( "Earth" )->getGroundStation( stationName )->getNominalStationState( )->getNominalCartesianPosition( );
    const auto timeScaleConverter = earth_orientation::createDefaultTimeConverter( );
    const Time utcTime = timeScaleConverter->getCurrentTime< Time >(
            basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, observationTime, stationPosition );
    const Time offsetUtcTime = timeScaleConverter->getCurrentTime< Time >(
            basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, offsetObservationTime, stationPosition );
    const Scalar utcTimeIncrement = convertIndependentVariableToScalar< Scalar >( offsetUtcTime - utcTime );
    const Time roundTripTdbTime = timeScaleConverter->getCurrentTime< Time >(
            basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, utcTime, stationPosition );
    const Time offsetRoundTripTdbTime = timeScaleConverter->getCurrentTime< Time >(
            basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, offsetUtcTime, stationPosition );
    const Scalar roundTripTdbTimeIncrement = convertIndependentVariableToScalar< Scalar >( offsetRoundTripTdbTime - roundTripTdbTime );
    const Scalar tdbToUtcCorrectionChange = timeScaleConverter->getTimeScaleConversionCorrectionDifference< Scalar, Time >(
            basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, observationTime, offsetObservationTime, stationPosition );
    const Scalar utcToTdbCorrectionChange = timeScaleConverter->getTimeScaleConversionCorrectionDifference< Scalar, Time >(
            basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, utcTime, offsetUtcTime, stationPosition );

    // The SOFA TDB-TT model is evaluated in double, but it is added to a Time.
    // At epochs away from a leap-second boundary, the common correction must
    // therefore not erase a picosecond separation in either direction.
    const Scalar timeIncrementTolerance = timeHasExtendedLongDoublePrecision ? scalarFromDecimalString< Scalar >( "2e-16" )
                                                                             : scalarFromDecimalString< Scalar >( "7.5e-13" );
    std::cout << "TDB->UTC increment error [s]: " << getAbsoluteValue( utcTimeIncrement - representedTimeIncrement )
              << ", UTC->TDB round-trip increment error [s]: " << getAbsoluteValue( roundTripTdbTimeIncrement - representedTimeIncrement )
              << ", time increment tolerance [s]: " << timeIncrementTolerance
              << ", TDB->UTC correction change [s]: " << tdbToUtcCorrectionChange
              << ", UTC->TDB correction change [s]: " << utcToTdbCorrectionChange << std::endl;
    BOOST_CHECK( utcTimeIncrement > static_cast< Scalar >( 0 ) );
    BOOST_CHECK( roundTripTdbTimeIncrement > static_cast< Scalar >( 0 ) );
    BOOST_CHECK( getAbsoluteValue( utcTimeIncrement - representedTimeIncrement ) < timeIncrementTolerance );
    BOOST_CHECK( getAbsoluteValue( roundTripTdbTimeIncrement - representedTimeIncrement ) < timeIncrementTolerance );
    BOOST_CHECK( getAbsoluteValue( tdbToUtcCorrectionChange ) < scalarFromDecimalString< Scalar >( "1e-20" ) );
    BOOST_CHECK( getAbsoluteValue( utcToTdbCorrectionChange ) < scalarFromDecimalString< Scalar >( "1e-20" ) );

    observation_models::LinkEnds dsnLinkEnds;
    dsnLinkEnds[ observation_models::transmitter ] = observation_models::linkEndId( "Earth", stationName );
    dsnLinkEnds[ observation_models::retransmitter ] = observation_models::linkEndId( "Jupiter" );
    dsnLinkEnds[ observation_models::receiver ] = observation_models::linkEndId( "Earth", stationName );
    const auto dsnSettings = std::make_shared< observation_models::DsnNWayAveragedDopplerObservationModelSettings >(
            observation_models::LinkDefinition( dsnLinkEnds ),
            std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > >( ),
            nullptr,
            convergenceCriteria );
    const auto quadDsnModel =
            observation_models::ObservationModelCreator< 1, Scalar, Time >::createObservationModel( dsnSettings, environment.bodies );
    BOOST_REQUIRE( ( std::dynamic_pointer_cast< observation_models::DsnNWayAveragedDopplerObservationModel< Scalar, Time > >(
                             quadDsnModel ) != nullptr ) );

    const double integrationTime = 0.1;
    const double referenceFrequency = 7.2e9;
    const auto ancillarySettings = observation_models::getDsnNWayAveragedDopplerAncillarySettings(
            { observation_models::x_band, observation_models::x_band }, observation_models::x_band, referenceFrequency, integrationTime );
    const auto computeDsnDoppler = [ & ]( const auto& model, const Time& epoch ) {
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        return model->computeIdealObservationsWithLinkEndData(
                epoch, observation_models::receiver, linkEndTimes, linkEndStates, ancillarySettings )( 0 );
    };

    const Scalar quadDoppler = computeDsnDoppler( quadDsnModel, observationTime );
    const Scalar offsetQuadDoppler = computeDsnDoppler( quadDsnModel, offsetObservationTime );
    const Scalar quadDopplerChange = offsetQuadDoppler - quadDoppler;
    BOOST_CHECK( quadDopplerChange != static_cast< Scalar >( 0 ) );
    BOOST_CHECK( getAbsoluteValue( quadDopplerChange ) < scalarFromDecimalString< Scalar >( "1e-6" ) );

    const Time localSlopeEndTime = observationTime + 1.0e-3L;
    const Scalar localSlopeTimeInterval = convertIndependentVariableToScalar< Scalar >( localSlopeEndTime - observationTime );
    const Scalar localDopplerSlope = ( computeDsnDoppler( quadDsnModel, localSlopeEndTime ) - quadDoppler ) / localSlopeTimeInterval;

    const auto computeDopplerSequence = [ & ]( const Time& startTime, const long double stepSize, const unsigned int numberOfSteps ) {
        std::vector< Time > epochs;
        std::vector< Scalar > observations;
        epochs.reserve( numberOfSteps + 1 );
        observations.reserve( numberOfSteps + 1 );
        for( unsigned int i = 0; i <= numberOfSteps; ++i )
        {
            const Time epoch = startTime + static_cast< long double >( i ) * stepSize;
            epochs.push_back( epoch );
            observations.push_back( computeDsnDoppler( quadDsnModel, epoch ) );
        }
        return std::make_pair( epochs, observations );
    };

    const auto picosecondSequence = computeDopplerSequence( observationTime, 1.0e-12L, 100 );
    const Time nanosecondSequenceStart = observationTime + 100.0e-12L;
    const auto nanosecondSequence = computeDopplerSequence( nanosecondSequenceStart, 1.0e-9L, 100 );
    BOOST_REQUIRE( picosecondSequence.first.back( ) == nanosecondSequence.first.front( ) );
    BOOST_CHECK( picosecondSequence.second.back( ) == nanosecondSequence.second.front( ) );

    struct DopplerSequenceMetrics {
        Scalar minimumIncrement;
        Scalar maximumIncrement;
        Scalar maximumAbsoluteIncrement;
        Scalar maximumIncrementResidual;
        Scalar maximumTrendResidual;
        unsigned int directionReversals;
    };
    const auto getDopplerSequenceMetrics = [ & ]( const std::pair< std::vector< Time >, std::vector< Scalar > >& sequence ) {
        const Scalar firstIncrement = sequence.second.at( 1 ) - sequence.second.at( 0 );
        DopplerSequenceMetrics metrics{
            firstIncrement, firstIncrement, getAbsoluteValue( firstIncrement ), static_cast< Scalar >( 0 ), static_cast< Scalar >( 0 ), 0
        };
        for( std::size_t i = 1; i < sequence.first.size( ); ++i )
        {
            const Scalar representedIncrement =
                    convertIndependentVariableToScalar< Scalar >( sequence.first.at( i ) - sequence.first.at( i - 1 ) );
            const Scalar observationIncrement = sequence.second.at( i ) - sequence.second.at( i - 1 );
            const Scalar expectedIncrement = localDopplerSlope * representedIncrement;
            const Scalar absoluteIncrement = getAbsoluteValue( observationIncrement );
            const Scalar incrementResidual = getAbsoluteValue( observationIncrement - expectedIncrement );
            const Scalar representedOffset = convertIndependentVariableToScalar< Scalar >( sequence.first.at( i ) - observationTime );
            const Scalar trendResidual = getAbsoluteValue( sequence.second.at( i ) - quadDoppler - localDopplerSlope * representedOffset );

            metrics.minimumIncrement = observationIncrement < metrics.minimumIncrement ? observationIncrement : metrics.minimumIncrement;
            metrics.maximumIncrement = observationIncrement > metrics.maximumIncrement ? observationIncrement : metrics.maximumIncrement;
            metrics.maximumAbsoluteIncrement =
                    absoluteIncrement > metrics.maximumAbsoluteIncrement ? absoluteIncrement : metrics.maximumAbsoluteIncrement;
            metrics.maximumIncrementResidual =
                    incrementResidual > metrics.maximumIncrementResidual ? incrementResidual : metrics.maximumIncrementResidual;
            metrics.maximumTrendResidual = trendResidual > metrics.maximumTrendResidual ? trendResidual : metrics.maximumTrendResidual;
            if( observationIncrement * localDopplerSlope <= static_cast< Scalar >( 0 ) )
            {
                ++metrics.directionReversals;
            }
        }
        return metrics;
    };

    const DopplerSequenceMetrics picosecondMetrics = getDopplerSequenceMetrics( picosecondSequence );
    const DopplerSequenceMetrics nanosecondMetrics = getDopplerSequenceMetrics( nanosecondSequence );
    const Scalar maximumDopplerIncrement =
            timeHasExtendedLongDoublePrecision ? scalarFromDecimalString< Scalar >( "1e-9" ) : scalarFromDecimalString< Scalar >( "2e-6" );
    const Scalar maximumDopplerTrendResidual =
            timeHasExtendedLongDoublePrecision ? scalarFromDecimalString< Scalar >( "1e-9" ) : scalarFromDecimalString< Scalar >( "1e-6" );
    std::cout << "1 ps maximum absolute increment [Hz]: " << picosecondMetrics.maximumAbsoluteIncrement
              << ", 1 ns maximum absolute increment [Hz]: " << nanosecondMetrics.maximumAbsoluteIncrement
              << ", 1 ps maximum trend residual [Hz]: " << picosecondMetrics.maximumTrendResidual
              << ", 1 ns maximum trend residual [Hz]: " << nanosecondMetrics.maximumTrendResidual
              << ", 1 ps maximum increment residual [Hz]: " << picosecondMetrics.maximumIncrementResidual
              << ", 1 ns maximum increment residual [Hz]: " << nanosecondMetrics.maximumIncrementResidual
              << ", increment tolerance [Hz]: " << maximumDopplerIncrement << ", trend tolerance [Hz]: " << maximumDopplerTrendResidual
              << std::endl;
    BOOST_CHECK( picosecondMetrics.maximumAbsoluteIncrement < maximumDopplerIncrement );
    BOOST_CHECK( nanosecondMetrics.maximumAbsoluteIncrement < maximumDopplerIncrement );
    BOOST_CHECK( picosecondMetrics.maximumTrendResidual < maximumDopplerTrendResidual );
    BOOST_CHECK( nanosecondMetrics.maximumTrendResidual < maximumDopplerTrendResidual );
    BOOST_CHECK( picosecondMetrics.maximumIncrementResidual < maximumDopplerIncrement );
    BOOST_CHECK( nanosecondMetrics.maximumIncrementResidual < maximumDopplerIncrement );

    const Scalar nanosecondSweepDuration =
            convertIndependentVariableToScalar< Scalar >( nanosecondSequence.first.back( ) - nanosecondSequence.first.front( ) );
    const Scalar nanosecondSweepChange = nanosecondSequence.second.back( ) - nanosecondSequence.second.front( );
    const Scalar expectedNanosecondSweepChange = localDopplerSlope * nanosecondSweepDuration;
    std::cout << "100 ns sweep change [Hz]: " << nanosecondSweepChange
              << ", expected 100 ns sweep change [Hz]: " << expectedNanosecondSweepChange
              << ", sweep residual [Hz]: " << getAbsoluteValue( nanosecondSweepChange - expectedNanosecondSweepChange ) << std::endl;
    if constexpr( timeHasExtendedLongDoublePrecision )
    {
        BOOST_CHECK( nanosecondSweepChange * expectedNanosecondSweepChange > static_cast< Scalar >( 0 ) );
        BOOST_CHECK( getAbsoluteValue( nanosecondSweepChange - expectedNanosecondSweepChange ) <
                     getAbsoluteValue( expectedNanosecondSweepChange ) * scalarFromDecimalString< Scalar >( "0.02" ) );
    }
    else
    {
        BOOST_CHECK( getAbsoluteValue( nanosecondSweepChange - expectedNanosecondSweepChange ) < maximumDopplerTrendResidual );
    }

    const double dopplerDoubleSpacing = std::nextafter( static_cast< double >( quadDoppler ), std::numeric_limits< double >::infinity( ) ) -
            static_cast< double >( quadDoppler );
    const auto doubleDsnModel =
            observation_models::ObservationModelCreator< 1, double, Time >::createObservationModel( dsnSettings, environment.bodies );
    const double doubleDoppler = computeDsnDoppler( doubleDsnModel, observationTime );
    const double offsetDoubleDoppler = computeDsnDoppler( doubleDsnModel, offsetObservationTime );
    const double doubleDopplerChange = offsetDoubleDoppler - doubleDoppler;
    const Scalar doubleObservableError = getAbsoluteValue( static_cast< Scalar >( doubleDoppler ) - quadDoppler );
    BOOST_CHECK_EQUAL( doubleDopplerChange, 0.0 );
    if constexpr( timeHasExtendedLongDoublePrecision )
    {
        BOOST_CHECK( doubleObservableError > getAbsoluteValue( quadDopplerChange ) * scalarFromDecimalString< Scalar >( "1e5" ) );
    }

    BOOST_TEST_MESSAGE( "TDB->UTC picosecond increment [s]: "
                        << utcTimeIncrement << ", UTC->TDB round-trip increment [s]: " << roundTripTdbTimeIncrement
                        << ", 0.1 s quad DSN Doppler [Hz]: " << quadDoppler << ", picosecond Doppler change [Hz]: " << quadDopplerChange
                        << ", double Doppler change [Hz]: " << doubleDopplerChange << ", double-vs-quad observable error [Hz]: "
                        << doubleObservableError << ", double spacing at observable [Hz]: " << dopplerDoubleSpacing
                        << ", local quad Doppler slope [Hz/s]: " << localDopplerSlope << ", 1 ps increment range [Hz]: ["
                        << picosecondMetrics.minimumIncrement << ", " << picosecondMetrics.maximumIncrement
                        << "], 1 ps maximum trend residual [Hz]: " << picosecondMetrics.maximumTrendResidual
                        << ", 1 ps direction reversals: " << picosecondMetrics.directionReversals << ", 1 ns increment range [Hz]: ["
                        << nanosecondMetrics.minimumIncrement << ", " << nanosecondMetrics.maximumIncrement
                        << "], 1 ns maximum trend residual [Hz]: " << nanosecondMetrics.maximumTrendResidual
                        << ", 1 ns direction reversals: " << nanosecondMetrics.directionReversals << ", 100 ns sweep change [Hz]: "
                        << nanosecondSweepChange << ", expected 100 ns sweep change [Hz]: " << expectedNanosecondSweepChange );
#endif
}

BOOST_AUTO_TEST_CASE( testCrossConfigurationConsistency )
{
    const long double longDoubleInterpolation = runRepresentativeInterpolationCase< long double >( );
    const Scalar configuredInterpolation = runRepresentativeInterpolationCase< Scalar >( );
    const long double configuredInterpolationAsLongDouble = static_cast< long double >( configuredInterpolation );

    const long double interpolationTolerance = 32.0L * std::numeric_limits< long double >::epsilon( );
    BOOST_CHECK_SMALL( configuredInterpolationAsLongDouble - longDoubleInterpolation, interpolationTolerance );

    const long double longDoublePropagation = runRepresentativeIntegrationCase< long double >( );
    const Scalar configuredPropagation = runRepresentativeIntegrationCase< Scalar >( );
    const long double configuredPropagationAsLongDouble = static_cast< long double >( configuredPropagation );
    const long double propagationTolerance = 128.0L * std::numeric_limits< long double >::epsilon( );
    BOOST_CHECK_SMALL( configuredPropagationAsLongDouble - longDoublePropagation, propagationTolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
