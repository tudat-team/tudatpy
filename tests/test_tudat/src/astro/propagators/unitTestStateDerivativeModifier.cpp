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

#include <cmath>
#include <iterator>

#include <Eigen/Core>

#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::numerical_integrators;
using namespace tudat::propagators;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_state_derivative_modifier )

namespace
{

constexpr double testGravitationalParameter = 3.986004418E14;
const Eigen::Vector3d testAdditionalAcceleration = ( Eigen::Vector3d( ) << 0.1, 0.2, 0.3 ).finished( );

SystemOfBodies createTestBodies( )
{
    BodyListSettings bodySettings( "SSB", "J2000" );
    bodySettings.addSettings( "Earth" );
    bodySettings.addSettings( "Vehicle" );

    bodySettings.at( "Earth" )->gravityFieldSettings = centralGravitySettings( testGravitationalParameter );
    bodySettings.at( "Earth" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings.at( "Vehicle" )->constantMass = 1000.0;

    return createSystemOfBodies< double, double >( bodySettings );
}

std::shared_ptr< SingleArcDynamicsSimulator< double, double > > createTestSimulator( const SystemOfBodies& bodies )
{
    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Vehicle" ][ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    const std::vector< std::string > bodiesToPropagate = { "Vehicle" };
    const std::vector< std::string > centralBodies = { "Earth" };
    const AccelerationMap accelerationModels =
            createAccelerationModelsMap( bodies, accelerationSettings, bodiesToPropagate, centralBodies );

    Eigen::VectorXd initialState( 6 );
    initialState << 7.0E6, 0.0, 0.0, 0.0, 7.5E3, 0.0;

    const double initialTime = 0.0;
    const double timeStep = 1.0;

    const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables = { totalAccelerationDependentVariable(
            "Vehicle" ) };

    const std::shared_ptr< IntegratorSettings< double > > integratorSettings = eulerSettings( timeStep );
    const std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            translationalStatePropagatorSettings( centralBodies,
                                                  accelerationModels,
                                                  bodiesToPropagate,
                                                  initialState,
                                                  initialTime,
                                                  integratorSettings,
                                                  propagationTimeTerminationSettings( initialTime + 3.0 * timeStep, true ),
                                                  cowell,
                                                  dependentVariables );

    return std::dynamic_pointer_cast< SingleArcDynamicsSimulator< double, double > >(
            createDynamicsSimulator< double, double >( bodies, propagatorSettings, false ) );
}

Eigen::Vector3d computePointMassGravityAcceleration( const Eigen::Vector6d& state )
{
    const Eigen::Vector3d position = state.segment( 0, 3 );
    return -testGravitationalParameter * position / std::pow( position.norm( ), 3 );
}

Eigen::Vector6d computeExpectedModifiedDerivative( const Eigen::Vector6d& state )
{
    Eigen::Vector6d expectedDerivative;
    expectedDerivative.segment( 0, 3 ) = state.segment( 3, 3 );
    expectedDerivative.segment( 3, 3 ) = computePointMassGravityAcceleration( state ) + testAdditionalAcceleration;
    return expectedDerivative;
}

void checkVectorEntryClose( const double computedValue, const double expectedValue )
{
    constexpr double tolerance = 1.0E-10;
    if( std::fabs( expectedValue ) < 1.0E-14 )
    {
        BOOST_CHECK_SMALL( computedValue - expectedValue, tolerance );
    }
    else
    {
        BOOST_CHECK_CLOSE_FRACTION( computedValue, expectedValue, tolerance );
    }
}

}  // namespace

BOOST_AUTO_TEST_CASE( testStateDerivativeModifierPropagation )
{
    const SystemOfBodies bodies = createTestBodies( );
    const std::shared_ptr< SingleArcDynamicsSimulator< double, double > > dynamicsSimulator = createTestSimulator( bodies );

    dynamicsSimulator->getDynamicsStateDerivative( )->setStateDerivativeModifierFunction(
            []( const Eigen::MatrixXd&, const Eigen::MatrixXd& nominalStateDerivative ) {
                Eigen::MatrixXd modifiedStateDerivative = nominalStateDerivative;
                modifiedStateDerivative.block( 3, 0, 3, modifiedStateDerivative.cols( ) ).colwise( ) += testAdditionalAcceleration;
                return modifiedStateDerivative;
            } );

    Eigen::Vector6d initialState;
    initialState << 7.0E6, 0.0, 0.0, 0.0, 7.5E3, 0.0;
    dynamicsSimulator->integrateEquationsOfMotion( initialState );

    const std::map< double, Eigen::VectorXd > stateHistory =
            dynamicsSimulator->getSingleArcPropagationResults( )->getEquationsOfMotionNumericalSolution( );
    BOOST_REQUIRE_EQUAL( stateHistory.size( ), 4 );

    for( auto stateIterator = stateHistory.begin( ); std::next( stateIterator ) != stateHistory.end( ); stateIterator++ )
    {
        const auto nextStateIterator = std::next( stateIterator );
        const double stepSize = nextStateIterator->first - stateIterator->first;
        const Eigen::Vector6d numericalDerivative = ( nextStateIterator->second - stateIterator->second ) / stepSize;
        const Eigen::Vector6d expectedDerivative = computeExpectedModifiedDerivative( stateIterator->second );

        for( int i = 0; i < 6; i++ )
        {
            checkVectorEntryClose( numericalDerivative( i ), expectedDerivative( i ) );
        }
    }

    const std::map< double, Eigen::VectorXd >& dependentVariableHistory = dynamicsSimulator->getDependentVariableHistory( );
    BOOST_REQUIRE_EQUAL( dependentVariableHistory.size( ), 4 );

    for( const auto& dependentVariableEntry : dependentVariableHistory )
    {
        const Eigen::Vector6d currentState = stateHistory.at( dependentVariableEntry.first );
        const Eigen::Vector3d expectedNominalAcceleration = computePointMassGravityAcceleration( currentState );
        for( int i = 0; i < 3; i++ )
        {
            checkVectorEntryClose( dependentVariableEntry.second( i ), expectedNominalAcceleration( i ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testStateDerivativeModifierWrongSize )
{
    const SystemOfBodies bodies = createTestBodies( );
    const std::shared_ptr< SingleArcDynamicsSimulator< double, double > > dynamicsSimulator = createTestSimulator( bodies );

    dynamicsSimulator->getDynamicsStateDerivative( )->setStateDerivativeModifierFunction(
            []( const Eigen::MatrixXd&, const Eigen::MatrixXd& ) { return Eigen::MatrixXd::Zero( 5, 1 ); } );

    Eigen::Vector6d initialState;
    initialState << 7.0E6, 0.0, 0.0, 0.0, 7.5E3, 0.0;

    BOOST_CHECK_THROW( dynamicsSimulator->integrateEquationsOfMotion( initialState ), std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
