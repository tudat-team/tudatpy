/*    Copyright (c) 2010-2016, Delft University of Technology
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
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>

#include <Tudat/Astrodynamics/Propagators/dynamicsSimulator.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/External/SpiceInterface/spiceInterface.h>
#include <Tudat/SimulationSetup/body.h>
#include <Tudat/SimulationSetup/createAccelerationModels.h>
#include <Tudat/SimulationSetup/defaultBodies.h>
#include <Tudat/InputOutput/basicInputOutput.h>

#include <iostream>
#include <limits>
#include <string>

#include <Eigen/Core>

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_propagation_stopping_conditions )

boost::shared_ptr< propagators::PropagationTerminationSettings > getTerminationSettings( const int testType )
{
    boost::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings;
    switch( testType )
    {
    case 0:
        terminationSettings = boost::make_shared< propagators::PropagationTimeTerminationSettings >( 3200.0 );
        break;
    case 1:
        terminationSettings = boost::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
                    boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                        propagators::mach_number_dependent_variable, "Apollo" ), 3.0, 1 );
        break;
    case 2:
        terminationSettings = boost::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
                    boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                        propagators::altitude_dependent_variable, "Apollo" ), 10.0E3, 1 );
        break;
    case 3:
        terminationSettings = boost::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
                    boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                        propagators::local_density_dependent_variable, "Apollo" ), 1.1, 0);
        break;

    case 4:
    {
        std::vector< boost::shared_ptr< propagators::PropagationTerminationSettings > > constituentSettings;
        constituentSettings.push_back( getTerminationSettings( 0 ) );
        constituentSettings.push_back( getTerminationSettings( 1 ) );
        constituentSettings.push_back( getTerminationSettings( 2 ) );
        constituentSettings.push_back( getTerminationSettings( 3 ) );

        terminationSettings = boost::make_shared< propagators::PropagationHybridTerminationSettings >(
                    constituentSettings, 1 );
        break;
    }
    case 5:
    {
        std::vector< boost::shared_ptr< propagators::PropagationTerminationSettings > > constituentSettings;
        constituentSettings.push_back( getTerminationSettings( 0 ) );
        constituentSettings.push_back( getTerminationSettings( 1 ) );
        constituentSettings.push_back( getTerminationSettings( 2 ) );
        constituentSettings.push_back( getTerminationSettings( 3 ) );

        terminationSettings = boost::make_shared< propagators::PropagationHybridTerminationSettings >(
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
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );


    ///////////////////////////////////////////////////////////////////////////

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

    // --------------------------------------------------------------------------------------------
    // CONVERT INITIAL STATE FROM KEPLERIAN TO CARTESIAN ELEMENTS
    // --------------------------------------------------------------------------------------------

    // Convert apollo state from Keplerian elements to Cartesian elements.
    const Vector6d apolloInitialState = convertKeplerianToCartesianElements(
                apolloInitialStateInKeplerianElements,
                getBodyGravitationalParameter( "Earth" ) );



    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< simulation_setup::ConstantEphemerisSettings >(
                basic_mathematics::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->gravityFieldSettings =
            boost::make_shared< simulation_setup::GravityFieldSettings >(
                central_spice );

    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    Eigen::VectorXd systemInitialState = Eigen::VectorXd( 6 );

    // Create vehicle objects.
    bodyMap[ "Apollo" ] = boost::make_shared< simulation_setup::Body >( );

    // Create vehicle aerodynamic coefficients
    bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );
    bodyMap[ "Apollo" ]->setConstantBodyMass( 5.0E3 );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
    accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[  "Apollo" ] = accelerationsOfApollo;

    bodiesToPropagate.push_back( "Apollo" );
    centralBodies.push_back( "Earth" );

    // Set initial state
    systemInitialState.segment( 0, 6 ) = apolloInitialState;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              getTerminationSettings( testType ) );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    switch( testType )
    {
    case 0:
    {

        BOOST_CHECK_EQUAL( ( (--( numericalSolution.end( ) ) )->first >= 3200.0 ) && ( (--(--( numericalSolution.end( ) ) ) )->first < 3200.0 ), true );
        break;
    }
    case 1:
    {
        boost::shared_ptr< DynamicsStateDerivativeModel< double, double > > stateDerivativeModel =
                dynamicsSimulator.getDynamicsStateDerivative( );
        stateDerivativeModel->computeStateDerivative(
                    (--(--( numericalSolution.end( ) ) ) )->first, (--(--( numericalSolution.end( ) ) ) )->second );
        boost::shared_ptr< FlightConditions > flightConditions = bodyMap.at( "Apollo" )->getFlightConditions( );
        double secondToLastMachNumber = flightConditions->getCurrentAirspeed( ) / flightConditions->getCurrentSpeedOfSound( );
        stateDerivativeModel->computeStateDerivative(
                    (--( numericalSolution.end( ) ) )->first, (--( numericalSolution.end( ) ) )->second );
        double lastMachNumber = flightConditions->getCurrentAirspeed( ) / flightConditions->getCurrentSpeedOfSound( );

        BOOST_CHECK_EQUAL( ( lastMachNumber <= 3.0 ) && ( secondToLastMachNumber > 3.0 ), true );

        break;
    }
    case 2:
    {
        boost::shared_ptr< DynamicsStateDerivativeModel< double, double > > stateDerivativeModel =
                dynamicsSimulator.getDynamicsStateDerivative( );
        stateDerivativeModel->computeStateDerivative(
                    (--(--( numericalSolution.end( ) ) ) )->first, (--(--( numericalSolution.end( ) ) ) )->second );
        boost::shared_ptr< FlightConditions > flightConditions = bodyMap.at( "Apollo" )->getFlightConditions( );
        double secondToLastAltitude = flightConditions->getCurrentAltitude( );
        stateDerivativeModel->computeStateDerivative(
                    (--( numericalSolution.end( ) ) )->first, (--( numericalSolution.end( ) ) )->second );
        double lastAltitude = flightConditions->getCurrentAltitude( );

        BOOST_CHECK_EQUAL( ( lastAltitude <= 10.0E3 ) && ( secondToLastAltitude > 10.0E3 ), true );

        break;
    }
    case 3:
    {
        boost::shared_ptr< DynamicsStateDerivativeModel< double, double > > stateDerivativeModel =
                dynamicsSimulator.getDynamicsStateDerivative( );
        stateDerivativeModel->computeStateDerivative(
                    (--(--( numericalSolution.end( ) ) ) )->first, (--(--( numericalSolution.end( ) ) ) )->second );
        boost::shared_ptr< FlightConditions > flightConditions = bodyMap.at( "Apollo" )->getFlightConditions( );
        double secondToLastDensity = flightConditions->getCurrentDensity( );
        stateDerivativeModel->computeStateDerivative(
                    (--( numericalSolution.end( ) ) )->first, (--( numericalSolution.end( ) ) )->second );
        double lastDensity = flightConditions->getCurrentDensity( );

        BOOST_CHECK_EQUAL( ( lastDensity >= 1.1 ) && ( secondToLastDensity < 1.1 ), true );

        break;
    }
    case 4:
    {
        boost::shared_ptr< DynamicsStateDerivativeModel< double, double > > stateDerivativeModel =
                dynamicsSimulator.getDynamicsStateDerivative( );
        stateDerivativeModel->computeStateDerivative(
                    (--(--( numericalSolution.end( ) ) ) )->first, (--(--( numericalSolution.end( ) ) ) )->second );
        boost::shared_ptr< FlightConditions > flightConditions = bodyMap.at( "Apollo" )->getFlightConditions( );

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
    case 5:
    {
        boost::shared_ptr< DynamicsStateDerivativeModel< double, double > > stateDerivativeModel =
                dynamicsSimulator.getDynamicsStateDerivative( );
        stateDerivativeModel->computeStateDerivative(
                    (--(--( numericalSolution.end( ) ) ) )->first, (--(--( numericalSolution.end( ) ) ) )->second );
        boost::shared_ptr< FlightConditions > flightConditions = bodyMap.at( "Apollo" )->getFlightConditions( );

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


BOOST_AUTO_TEST_CASE( testPropagationStoppingConditions )
{
    for( unsigned int i = 0; i < 6; i++ )
    {
        performSimulation( i );
    }
}


BOOST_AUTO_TEST_SUITE_END( )


}

}

