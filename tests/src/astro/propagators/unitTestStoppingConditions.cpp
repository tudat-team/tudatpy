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
#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

#include <boost/make_shared.hpp>
#include <memory>

#include "tudat/astro/aerodynamics/testApolloCapsuleCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/io/basicInputOutput.h"
#include <limits>
#include <string>

#include <Eigen/Core>

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_propagation_stopping_conditions )

std::shared_ptr< propagators::PropagationTerminationSettings > getTerminationSettings( const int testType )
{

    // Define stopping conditions, depending on test case.
    std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings;
    switch( testType )
    {
    // Stop at given time.
    case 0:
        terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( 3200.0 );
        break;
        // Stop at given Mach number
    case 1:
        terminationSettings = std::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
                    std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                        propagators::mach_number_dependent_variable, "Apollo" ), 3.0, 1 );
        break;
        // Stop at given altitude
    case 2:
        terminationSettings = std::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
                    std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                        propagators::altitude_dependent_variable, "Apollo" ), 10.0E3, 1 );
        break;
        // Stop at given density
    case 3:
        terminationSettings = std::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
                    std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                        propagators::local_density_dependent_variable, "Apollo" ), 1.1, 0);
        break;
        // Stop when a single of the conditions 0-3 is fulfilled.
    case 4:
    {
        std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > constituentSettings;
        constituentSettings.push_back( getTerminationSettings( 0 ) );
        constituentSettings.push_back( getTerminationSettings( 1 ) );
        constituentSettings.push_back( getTerminationSettings( 2 ) );
        constituentSettings.push_back( getTerminationSettings( 3 ) );

        terminationSettings = std::make_shared< propagators::PropagationHybridTerminationSettings >(
                    constituentSettings, 1 );
        break;
    }
        // Stop when all of the conditions 0-3 is fulfilled.
    case 5:
    {
        std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > constituentSettings;
        constituentSettings.push_back( getTerminationSettings( 0 ) );
        constituentSettings.push_back( getTerminationSettings( 1 ) );
        constituentSettings.push_back( getTerminationSettings( 2 ) );
        constituentSettings.push_back( getTerminationSettings( 3 ) );

        terminationSettings = std::make_shared< propagators::PropagationHybridTerminationSettings >(
                    constituentSettings, 0 );
        break;
    }
    }
    return terminationSettings;
}

void performSimulation( const int testType )
{
    using namespace ephemerides;
    using namespace interpolators;
    using namespace numerical_integrators;
    using namespace spice_interface;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace orbital_element_conversions;
    using namespace propagators;
    using namespace aerodynamics;
    using namespace basic_mathematics;
    using namespace input_output;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 3300.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 1.0;


    // Set Keplerian elements for Capsule.
    Vector6d apolloInitialStateInKeplerianElements;
    apolloInitialStateInKeplerianElements( semiMajorAxisIndex ) = spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    apolloInitialStateInKeplerianElements( eccentricityIndex ) = 0.005;
    apolloInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    apolloInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    apolloInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    apolloInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    // Convert apollo state from Keplerian elements to Cartesian elements.
    const Eigen::Vector6d apolloInitialState = convertKeplerianToCartesianElements(
                apolloInitialStateInKeplerianElements,
                getBodyGravitationalParameter( "Earth" ) );


    // Define simulation body settings.
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize, "SSB", "J2000"  );
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings.at( "Earth" )->gravityFieldSettings =
            std::make_shared< simulation_setup::GravityFieldSettings >(
                central_spice );
    bodySettings.at( "Earth" )->rotationModelSettings->resetOriginalFrame( "J2000" );

    // Create Earth object
    simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Create vehicle objects.
    bodies.createEmptyBody( "Apollo" );

    // Create vehicle aerodynamic coefficients
    bodies.at( "Apollo" )->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );
    bodies.at( "Apollo" )->setConstantBodyMass( 5.0E3 );


    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
    accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[ "Apollo" ] = accelerationsOfApollo;

    bodiesToPropagate.push_back( "Apollo" );
    centralBodies.push_back( "Earth" );

    // Set initial state
    Eigen::Vector6d systemInitialState = apolloInitialState;

    // Create acceleration models and propagation settings, using current test case to retrieve stop settings..
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              getTerminationSettings( testType ) );
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, false );

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    // Check whether propagation has stopped at given conditions
    switch( testType )
    {
    // Check whether propagation stopped when passing t=3200.0
    case 0:
    {

        BOOST_CHECK_EQUAL( ( (--( numericalSolution.end( ) ) )->first >= 3200.0 ) &&
                           ( (--(--( numericalSolution.end( ) ) ) )->first < 3200.0 ), true );
        break;
    }
    case 1:
    {
        // Compute Mach number for last two time steps, and check whether last step was first to pass below Mach=3
        std::shared_ptr< DynamicsStateDerivativeModel< double, double > > stateDerivativeModel =
                dynamicsSimulator.getDynamicsStateDerivative( );
        stateDerivativeModel->computeStateDerivative(
                    (--(--( numericalSolution.end( ) ) ) )->first, (--(--( numericalSolution.end( ) ) ) )->second );
        std::shared_ptr< AtmosphericFlightConditions > flightConditions =
                std::dynamic_pointer_cast< AtmosphericFlightConditions >( bodies.at( "Apollo" )->getFlightConditions( ) );
        double secondToLastMachNumber = flightConditions->getCurrentAirspeed( ) / flightConditions->getCurrentSpeedOfSound( );
        stateDerivativeModel->computeStateDerivative(
                    (--( numericalSolution.end( ) ) )->first, (--( numericalSolution.end( ) ) )->second );
        double lastMachNumber = flightConditions->getCurrentAirspeed( ) / flightConditions->getCurrentSpeedOfSound( );

        BOOST_CHECK_EQUAL( ( lastMachNumber <= 3.0 ) && ( secondToLastMachNumber > 3.0 ), true );

        break;
    }
    case 2:
    {
        // Compute altitude for last two time steps, and check whether last step was first to pass below altitude = 10 km
        std::shared_ptr< DynamicsStateDerivativeModel< double, double > > stateDerivativeModel =
                dynamicsSimulator.getDynamicsStateDerivative( );
        stateDerivativeModel->computeStateDerivative(
                    (--(--( numericalSolution.end( ) ) ) )->first, (--(--( numericalSolution.end( ) ) ) )->second );
        std::shared_ptr< FlightConditions > flightConditions = bodies.at( "Apollo" )->getFlightConditions( );
        double secondToLastAltitude = flightConditions->getCurrentAltitude( );
        stateDerivativeModel->computeStateDerivative(
                    (--( numericalSolution.end( ) ) )->first, (--( numericalSolution.end( ) ) )->second );
        double lastAltitude = flightConditions->getCurrentAltitude( );

        BOOST_CHECK_EQUAL( ( lastAltitude <= 10.0E3 ) && ( secondToLastAltitude > 10.0E3 ), true );

        break;
    }
    case 3:
    {
        // Compute density for last two time steps, and check whether last step was first to pass below density = 1.1 kg/m3
        std::shared_ptr< DynamicsStateDerivativeModel< double, double > > stateDerivativeModel =
                dynamicsSimulator.getDynamicsStateDerivative( );
        stateDerivativeModel->computeStateDerivative(
                    (--(--( numericalSolution.end( ) ) ) )->first, (--(--( numericalSolution.end( ) ) ) )->second );
        std::shared_ptr< AtmosphericFlightConditions > flightConditions =
                std::dynamic_pointer_cast< AtmosphericFlightConditions >( bodies.at( "Apollo" )->getFlightConditions( ) );
        double secondToLastDensity = flightConditions->getCurrentDensity( );
        stateDerivativeModel->computeStateDerivative(
                    (--( numericalSolution.end( ) ) )->first, (--( numericalSolution.end( ) ) )->second );
        double lastDensity = flightConditions->getCurrentDensity( );

        BOOST_CHECK_EQUAL( ( lastDensity >= 1.1 ) && ( secondToLastDensity < 1.1 ), true );

        break;
    }
        // Check whether at least a single of the conditions for case 0-3 was first reached at last time step.
    case 4:
    {
        std::shared_ptr< DynamicsStateDerivativeModel< double, double > > stateDerivativeModel =
                dynamicsSimulator.getDynamicsStateDerivative( );
        stateDerivativeModel->computeStateDerivative(
                    (--(--( numericalSolution.end( ) ) ) )->first, (--(--( numericalSolution.end( ) ) ) )->second );
        std::shared_ptr< AtmosphericFlightConditions > flightConditions =
                std::dynamic_pointer_cast< AtmosphericFlightConditions >( bodies.at( "Apollo" )->getFlightConditions( ) );

        double secondToLastMachNumber = flightConditions->getCurrentAirspeed( ) / flightConditions->getCurrentSpeedOfSound( );
        double secondToLastAltitude = flightConditions->getCurrentAltitude( );
        double secondToLastDensity = flightConditions->getCurrentDensity( );

        stateDerivativeModel->computeStateDerivative(
                    (--( numericalSolution.end( ) ) )->first, (--( numericalSolution.end( ) ) )->second );

        double lastMachNumber = flightConditions->getCurrentAirspeed( ) / flightConditions->getCurrentSpeedOfSound( );
        double lastAltitude = flightConditions->getCurrentAltitude( );
        double lastDensity = flightConditions->getCurrentDensity( );

        BOOST_CHECK_EQUAL( ( ( (--( numericalSolution.end( ) ) )->first >= 3200.0 ) && ( (--(--( numericalSolution.end( ) ) ) )->first < 3200.0 ) ) ||
                           ( ( lastMachNumber <= 3.0 ) && ( secondToLastMachNumber > 3.0 ) ) ||
                           ( ( lastAltitude <= 10.0E3 ) && ( secondToLastAltitude > 10.0E3 ) ) ||
                           ( ( lastDensity >= 1.1 ) && ( secondToLastDensity < 1.1 ) ), 1 );

        break;
    }
        // Check whether all conditions for case 0-3 was first reached at last time step.
    case 5:
    {
        std::shared_ptr< DynamicsStateDerivativeModel< double, double > > stateDerivativeModel =
                dynamicsSimulator.getDynamicsStateDerivative( );
        stateDerivativeModel->computeStateDerivative(
                    (--(--( numericalSolution.end( ) ) ) )->first, (--(--( numericalSolution.end( ) ) ) )->second );
        std::shared_ptr< AtmosphericFlightConditions > flightConditions =
                std::dynamic_pointer_cast< AtmosphericFlightConditions >( bodies.at( "Apollo" )->getFlightConditions( ) );

        double secondToLastMachNumber = flightConditions->getCurrentAirspeed( ) / flightConditions->getCurrentSpeedOfSound( );
        double secondToLastAltitude = flightConditions->getCurrentAltitude( );
        double secondToLastDensity = flightConditions->getCurrentDensity( );

        stateDerivativeModel->computeStateDerivative(
                    (--( numericalSolution.end( ) ) )->first, (--( numericalSolution.end( ) ) )->second );

        double lastMachNumber = flightConditions->getCurrentAirspeed( ) / flightConditions->getCurrentSpeedOfSound( );
        double lastAltitude = flightConditions->getCurrentAltitude( );
        double lastDensity = flightConditions->getCurrentDensity( );

        BOOST_CHECK_EQUAL( ( ( (--( numericalSolution.end( ) ) )->first >= 3200.0 ) ) &&
                           ( ( lastMachNumber <= 3.0 ) ) &&
                           ( ( lastAltitude <= 10.0E3 ) ) &&
                           ( ( lastDensity >= 1.1 ) ), 1 );

        BOOST_CHECK_EQUAL( ( (--(--( numericalSolution.end( ) ) ) )->first < 3200.0 )||
                           ( ( secondToLastMachNumber > 3.0 ) ) ||
                           ( ( secondToLastAltitude > 10.0E3 ) ) ||
                           ( ( secondToLastDensity < 1.1 ) ), 1 );

        break;
    }

    }
}

//! Test to perform propagation of Apollo capsule for various stopping conditions, and checking whether the final state
//! corresponds to condition that was given.
BOOST_AUTO_TEST_CASE( testPropagationStoppingConditions )
{
    for( unsigned int i = 0; i < 6; i++ )
    {
        performSimulation( i );
    }
}

//! Test to see if environment is updated to evaluate stopping conditions if required
//! (test uses altitude above Earth for Kepler orbit dynamics)
BOOST_AUTO_TEST_CASE( testPropagationStoppingConditionsWithDependentVariableUpdate )
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::unit_conversions;

    spice_interface::loadStandardSpiceKernels( );

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodiesToCreate );
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ) );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Asterix" );
    

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::point_mass_gravity ) );
    accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to Cartesian
    // elements.

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements =Eigen::Vector6d::Zero( );
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7000.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.5;
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 180.0 );

    // Convert Asterix state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements,
                earthGravitationalParameter );

    // Define termination conditions
    std::shared_ptr< PropagationTerminationSettings > stoppingCondition =
            std::make_shared< PropagationDependentVariableTerminationSettings >(
                std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                propagators::altitude_dependent_variable, "Asterix", "Earth" ), 25.0E3, 1 );

    // Create simulation object and propagate dynamics.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, stoppingCondition );
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, 10.0 );
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );


    // Check whether termination done correctly
    auto stateIterator = integrationResult.rbegin( );
    BOOST_CHECK_EQUAL(
               ( stateIterator->second.segment( 0, 3 ).norm( ) -
                 tudat::spice_interface::getAverageRadius( "Earth" ) ) < 25.0E3, true );
    stateIterator++;
    BOOST_CHECK_EQUAL(
               ( stateIterator->second.segment( 0, 3 ).norm( ) -
                 tudat::spice_interface::getAverageRadius( "Earth" ) ) < 25.0E3, false );
}



BOOST_AUTO_TEST_SUITE_END( )


}

}

