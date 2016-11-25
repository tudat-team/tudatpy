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
#include <Tudat/Basics/testMacros.h>
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/External/SpiceInterface/spiceInterface.h>
#include <Tudat/SimulationSetup/EnvironmentSetup/body.h>
#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
#include <Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h>
#include <Tudat/InputOutput/basicInputOutput.h>

#include <iostream>
#include <limits>
#include <string>

#include <Eigen/Core>

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_dependent_variable_output )

//! Propagate entry of Apollo capsule, and save a list of dependent variables during entry. The saved dependent variables
//! are compared against theoretical/manual values in this test.
BOOST_AUTO_TEST_CASE( testDependentVariableOutput )
{
    using namespace tudat;
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
    const Vector6d apolloInitialState = convertKeplerianToCartesianElements(
                apolloInitialStateInKeplerianElements,
                getBodyGravitationalParameter( "Earth" ) );

    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings[ "Earth" ]->gravityFieldSettings =
            boost::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects.
    bodyMap[ "Apollo" ] = boost::make_shared< simulation_setup::Body >( );

    // Create vehicle aerodynamic coefficients
    bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );
    bodyMap[ "Apollo" ]->setConstantBodyMass( 5.0E3 );
    bodyMap[ "Apollo" ]->setEphemeris(
                boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, basic_mathematics::Vector6d  > >( ),
                    "Earth" ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
    accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationsOfApollo[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Apollo" ] = accelerationsOfApollo;

    bodiesToPropagate.push_back( "Apollo" );
    centralBodies.push_back( "Earth" );

    // Set initial state
    basic_mathematics::Vector6d systemInitialState = apolloInitialState;

    // Define list of dependent variables to save.
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "Apollo" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable,
                                                                           "Apollo", "Earth" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >( relative_distance_dependent_variable,
                                                                           "Apollo", "Earth" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >( relative_speed_dependent_variable,
                                                                           "Apollo", "Earth" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    central_gravity, "Apollo", "Earth", 1 ) );


    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >( relative_position_dependent_variable,
                                                                           "Apollo", "Earth" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >( relative_velocity_dependent_variable,
                                                                           "Apollo", "Earth" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    central_gravity, "Apollo", "Earth", 0 ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    total_acceleration_dependent_variable, "Apollo" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_moment_coefficients_dependent_variable, "Apollo" ) );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    setTrimmedConditions( bodyMap.at( "Apollo" ) );

    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              boost::make_shared< propagators::PropagationTimeTerminationSettings >( 3200.0 ), cowell,
              boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );

    // Retrieve numerical solutions for state and dependent variables
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableSolution =
            dynamicsSimulator.getDependentVariableHistory( );

    // Iterate over results for dependent variables, and check against manually retrieved values.
    basic_mathematics::Vector6d currentStateDerivative;
    Eigen::Vector3d manualCentralGravity;
    for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
         variableIterator != dependentVariableSolution.end( ); variableIterator++ )
    {
        currentStateDerivative = dynamicsSimulator.getDynamicsStateDerivative( )->computeStateDerivative(
                    variableIterator->first, numericalSolution.at( variableIterator->first ) );

        // Manually compute central gravity.
        manualCentralGravity =
                -bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ) *
                variableIterator->second.segment( 5, 3 ) /
                std::pow( variableIterator->second.segment( 5, 3 ).norm( ), 3 );

        // Check output time consistency
        BOOST_CHECK_EQUAL( numericalSolution.count( variableIterator->first ), 1 );

        // Check relative position and velocity against state
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( numericalSolution.at( variableIterator->first )( i ) -
                                   variableIterator->second( 5 + i ) ), 2.0E-5 );
            BOOST_CHECK_SMALL(
                        std::fabs( numericalSolution.at( variableIterator->first )( 3 + i ) -
                                   variableIterator->second( 8 + i ) ), 5.0E-11 );
        }

        // Check central gravity acceleration
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    manualCentralGravity.segment( 0, 3 ),
                    variableIterator->second.segment( 11, 3 ), ( 5.0 * std::numeric_limits< double >::epsilon( ) ) );

        // Check total acceleration (tolerance is not epsilon due to numerical root finding for trim)
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( currentStateDerivative( 3 + i ) -
                                   variableIterator->second( 14 + i ) ), 1.0E-13 );
        }

        // Check relative position and velocity norm.
        BOOST_CHECK_SMALL(
                    std::fabs( ( numericalSolution.at( variableIterator->first ).segment( 0, 3 ) ).norm( ) -
                               variableIterator->second( 2 ) ), 2.0E-5 );
        BOOST_CHECK_SMALL(
                    std::fabs( ( numericalSolution.at( variableIterator->first ).segment( 3, 3 ) ).norm( ) -
                               variableIterator->second( 3 ) ), 2.0E-11 );

        // Check central gravity acceleration norm
        BOOST_CHECK_CLOSE_FRACTION(
                    manualCentralGravity.norm( ),
                    variableIterator->second( 4 ), 5.0 * std::numeric_limits< double >::epsilon( ) );

        // Check Mach number
        BOOST_CHECK_CLOSE_FRACTION(
                    bodyMap.at( "Apollo" )->getFlightConditions( )->getCurrentAirspeed( ) /
                    bodyMap.at( "Apollo" )->getFlightConditions( )->getCurrentSpeedOfSound( ),
                    variableIterator->second( 0 ), std::numeric_limits< double >::epsilon( ) );

        // Check altitude.
        BOOST_CHECK_CLOSE_FRACTION(
                    bodyMap.at( "Apollo" )->getFlightConditions( )->getCurrentAltitude( ),
                    variableIterator->second( 1 ), std::numeric_limits< double >::epsilon( ) );

        // Check trimmed condition (y-term)/symmetric vehicle shape (x- and z-term).
        BOOST_CHECK_SMALL(
                    std::fabs( variableIterator->second( 17 ) ), 1.0E-14 );
        BOOST_CHECK_SMALL(
                    std::fabs( variableIterator->second( 18 ) ), 1.0E-10 );
        BOOST_CHECK_SMALL(
                    std::fabs( variableIterator->second( 19 ) ), 1.0E-14 );


    }

}

BOOST_AUTO_TEST_SUITE_END( )


}

}


