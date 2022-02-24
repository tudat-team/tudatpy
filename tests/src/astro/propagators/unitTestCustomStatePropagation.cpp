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

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

#include <boost/make_shared.hpp>
#include <memory>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/math/integrators/rungeKutta4Integrator.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::propagators;
using namespace tudat::simulation_setup;
using namespace tudat::numerical_integrators;

BOOST_AUTO_TEST_SUITE( test_body_custom_state_propagation )

double getDummyCustomState1(
        const double currentTime, const double currentCustomState )
{
    return -0.02;
}

double getDummyCustomState2(
        const double currentTime, const double currentCustomState )
{
    return -0.02 * currentTime;
}

double getDummyCustomState3(
        const double currentTime, const double currentCustomState )
{
    return -0.002 * currentCustomState;
}

double getDummyCustomState4(
        const double currentTime, const double currentCustomState )
{
    return -0.00002 * currentCustomState * currentTime;
}


// Test custom state propagation, linearly decreasing with time
BOOST_AUTO_TEST_CASE( testSingleCustomStatePropagation )
{
    // Crate bodies
    SystemOfBodies bodies;

    // Create settings for propagation
    double initialCustomState = 500.0;
    std::shared_ptr< CustomStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< CustomStatePropagatorSettings< double > >(
                std::bind( &getDummyCustomState1, std::placeholders::_1, std::placeholders::_2 ), initialCustomState,
                std::make_shared< PropagationTimeTerminationSettings >( 1000.0 ) );

    // Define numerical integrator settings.
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, 1.0 );

    // Create dynamics simulation object.
    SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, false );

    // Test propagated solution.
    std::map< double, Eigen::VectorXd > integratedState = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integratedState.begin( );
         stateIterator != integratedState.end( ); stateIterator++ )
    {
        BOOST_CHECK_EQUAL( stateIterator->second.rows( ), 1 );
        BOOST_CHECK_SMALL( std::fabs( stateIterator->second( 0 ) - ( 500.0 - 0.02 * stateIterator->first ) ), 1.0E-9 );
    }
}

// Test custom state propagation, quadratically decreasing with time
BOOST_AUTO_TEST_CASE( testSingleCustomStatePropagation2 )
{
    // Crate bodies
    SystemOfBodies bodies;

    // Create settings for propagation
    double initialCustomState = 500.0;
    std::shared_ptr< CustomStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< CustomStatePropagatorSettings< double > >(
                std::bind( &getDummyCustomState2, std::placeholders::_1, std::placeholders::_2 ), initialCustomState,
                std::make_shared< PropagationTimeTerminationSettings >( 1000.0 ) );

    // Define numerical integrator settings.
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, 1.0 );

    // Create dynamics simulation object.
    SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, false );

    // Test propagated solution.
    std::map< double, Eigen::VectorXd > integratedState = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integratedState.begin( );
         stateIterator != integratedState.end( ); stateIterator++ )
    {
        BOOST_CHECK_EQUAL( stateIterator->second.rows( ), 1 );
        BOOST_CHECK_SMALL( std::fabs( stateIterator->second( 0 ) -
                                      ( 500.0 - 0.5 * 0.02 * stateIterator->first  * stateIterator->first ) ), 1.0E-9 );
    }
}

// Test custom state propagation, exponentially decreasing with time
BOOST_AUTO_TEST_CASE( testSingleCustomStatePropagation3 )
{
    // Crate bodies
    SystemOfBodies bodies;

    // Create settings for propagation
    double initialCustomState = 500.0;
    std::shared_ptr< CustomStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< CustomStatePropagatorSettings< double > >(
                std::bind( &getDummyCustomState3, std::placeholders::_1, std::placeholders::_2 ), initialCustomState,
                std::make_shared< PropagationTimeTerminationSettings >( 1000.0 ) );

    // Define numerical integrator settings.
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, 1.0 );

    // Create dynamics simulation object.
    SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, false );

    // Test propagated solution.
    std::map< double, Eigen::VectorXd > integratedState = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integratedState.begin( );
         stateIterator != integratedState.end( ); stateIterator++ )
    {
        BOOST_CHECK_EQUAL( stateIterator->second.rows( ), 1 );
        BOOST_CHECK_CLOSE_FRACTION( stateIterator->second( 0 ), 500.0 * std::exp( -stateIterator->first * 0.002 ), 1.0E-9 );
    }
}

// Test custom state propagation, exponentially decreasing with square of time
BOOST_AUTO_TEST_CASE( testSingleCustomStatePropagation4 )
{
    // Crate bodies
    SystemOfBodies bodies;

    // Create settings for propagation
    double initialCustomState = 500.0;
    std::shared_ptr< CustomStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< CustomStatePropagatorSettings< double > >(
                std::bind( &getDummyCustomState4, std::placeholders::_1, std::placeholders::_2 ), initialCustomState,
                std::make_shared< PropagationTimeTerminationSettings >( 100.0 ) );

    // Define numerical integrator settings.
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, 0.01 );

    // Create dynamics simulation object.
    SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, false );

    // Test propagated solution.
    std::map< double, Eigen::VectorXd > integratedState = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integratedState.begin( );
         stateIterator != integratedState.end( ); stateIterator++ )
    {

        BOOST_CHECK_EQUAL( stateIterator->second.rows( ), 1 );
        BOOST_CHECK_CLOSE_FRACTION( stateIterator->second( 0 ),
                                    500.0 * std::exp( -0.5 * stateIterator->first * stateIterator->first * 0.00002 ), 1.0E-9 );
    }
}

// Test if custom mass rate is properly computed in multi-type propagation
// (majority of test copied form unitTestMultiTypeStatePropagation).
BOOST_AUTO_TEST_CASE( testMultiTypeCustomStatePropagation )
{
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace gravitation;
    using namespace numerical_integrators;
    using namespace unit_conversions;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 10.0;

    // Define body settings for simulation.
    BodyListSettings bodySettings = BodyListSettings( "SSB", "J2000" );
    bodySettings.addSettings( "Earth" );
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );

    // Create Earth object
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "Asterix" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::point_mass_gravity ) );
    accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );


    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Convert Asterix state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements,
                earthGravitationalParameter );


    std::shared_ptr< SingleArcPropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              std::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );

    // Create mass rate model and mass propagation settings
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Asterix" ] = std::make_shared< basic_astrodynamics::CustomMassRateModel >(
                [ ]( const double ){ return -0.01; } );
    Eigen::VectorXd initialMass = Eigen::VectorXd( 1 );
    initialMass( 0 ) = 500.0;
    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                std::vector< std::string >{ "Asterix" }, massRateModels, initialMass,
                std::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );

    // Create custom state derivative model settings
    double initialCustomState = 500.0;
    std::shared_ptr< SingleArcPropagatorSettings< double > > customPropagatorSettings =
            std::make_shared< CustomStatePropagatorSettings< double > >(
                std::bind( &getDummyCustomState1, std::placeholders::_1, std::placeholders::_2 ), initialCustomState,
                std::make_shared< PropagationTimeTerminationSettings >( 1000.0 ) );

    // Create total propagator settings, depending on current case.
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings;

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );
    propagatorSettingsList.push_back( massPropagatorSettings );
    propagatorSettingsList.push_back( customPropagatorSettings );

    propagatorSettings = std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsList,
                std::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );


    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, true );

    std::map< double, Eigen::VectorXd > integratedState = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    // Test propagated solution.
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integratedState.begin( );
         stateIterator != integratedState.end( ); stateIterator++ )
    {
        BOOST_CHECK_EQUAL( stateIterator->second.rows( ), 8 );
        BOOST_CHECK_SMALL( std::fabs( stateIterator->second( 6 ) - ( 500.0 - 0.01 * stateIterator->first ) ), 1.0E-9 );
        BOOST_CHECK_SMALL( std::fabs( stateIterator->second( 7 ) - ( 500.0 - 0.02 * stateIterator->first ) ), 1.0E-9 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
