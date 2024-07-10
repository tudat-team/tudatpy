/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN





#include <memory>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/simulation.h"
//
//namespace tudat
//{
//
//namespace unit_tests
//{

using namespace tudat;
using namespace tudat::propagators;
using namespace tudat::simulation_setup;
using namespace tudat::numerical_integrators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;

//BOOST_AUTO_TEST_SUITE( test_propagation_radiation_pressure )
//
//// Test custom state propagation, linearly decreasing with time
//BOOST_AUTO_TEST_CASE( testMultiTypeRadiationPressurePropagation )
//{
int main( )
{
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames = { "Earth", "Moon", "Sun" };

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 3.0 * 3600.0;

    // Get initial state vector as input to integration.
    // Set Keplerian elements for Capsule.
    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = spice_interface::getAverageRadius( "Earth" ) + 800.0E3;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.005;
    initialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex )
        = unit_conversions::convertDegreesToRadians( 235.7 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
        = unit_conversions::convertDegreesToRadians( 23.4 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );


    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings(
        bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.addSettings( "Vehicle" );
    bodySettings.get( "Vehicle" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Earth", "ECLIPJ2000", "Vehicle_Fixed" );
    bodySettings.get( "Vehicle" )->bodyExteriorPanelSettings_ = bodyWingPanelledGeometry( 2.0, 1.0, 4.0, 20.0, 0.4, 0.3, 0.5, 0.1, true, true );
    bodySettings.get( "Vehicle" )->ephemerisSettings = keplerEphemerisSettings(
        initialStateInKeplerianElements, initialEphemerisTime, spice_interface::getBodyGravitationalParameter( "Earth" ), "Earth", "ECLIPJ2000" );
    bodySettings.get( "Vehicle" )->rigidBodyPropertiesSettings = constantRigidBodyPropertiesSettings( 300.0 );
    bodySettings.get("Earth")->radiationSourceModelSettings =
        extendedRadiationSourceModelSettings( {
            albedoPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::albedo_knocke, "Sun"),
            delayedThermalPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke, "Sun" ) }, {6, 12} );
    // Crate bodies
    std::map<std::string, std::vector<std::string>> sourceToTargetOccultingBodies;
    sourceToTargetOccultingBodies[ "Sun" ] = { "Earth" };
    sourceToTargetOccultingBodies[ "Moon" ] = { "Earth" };

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    addRadiationPressureTargetModel( bodies, "Vehicle", paneledRadiationPressureTargetModelSettingsWithOccultationMap( sourceToTargetOccultingBodies ) );
    addRadiationPressureTargetModel( bodies, "Vehicle", cannonballRadiationPressureTargetModelSettingsWithOccultationMap( 3.0, 1.5, sourceToTargetOccultingBodies ) );

    // Set accelerations between bodies that are to be taken into account (mutual point mass gravity between all bodies).
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Sun" ].push_back( radiationPressureAcceleration( paneled_target ) );
    accelerationsOfVehicle[ "Earth" ].push_back( radiationPressureAcceleration( paneled_target ) );
    accelerationsOfVehicle[ "Earth" ].push_back( pointMassGravityAcceleration(  ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Define list of bodies to propagate
    std::vector< std::string > bodiesToIntegrate = { "Vehicle" };
    std::vector< std::string > centralBodies = { "Earth" };

    // Create acceleration models and propagation settings.
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
        bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Define numerical integrator settings.
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
        std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialEphemerisTime, 60.0 );


    // Convert apollo state from Keplerian elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = convertKeplerianToCartesianElements(
        initialStateInKeplerianElements, bodies.at( "Earth" )->getGravitationalParameter( ) );

    // Create dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( singleAccelerationDependentVariable( radiation_pressure, "Vehicle", "Sun" ) );
    dependentVariablesList.push_back( singleAccelerationDependentVariable( radiation_pressure, "Vehicle", "Earth" ) );
    dependentVariablesList.push_back( vehiclePanelInertialSurfaceNormals( "Vehicle" ) );
    dependentVariablesList.push_back( vehiclePanelBodyFixedSurfaceNormals( "Vehicle" ) );
    dependentVariablesList.push_back( perVehiclePanelRadiationPressureForce( "Vehicle", "Sun" ) );
    dependentVariablesList.push_back( perVehiclePanelRadiationPressureForce( "Vehicle", "Earth" ) );
    dependentVariablesList.push_back( paneledRadiationSourcePerPanelIrradiance( "Vehicle", "Earth" ) );
    dependentVariablesList.push_back( paneledRadiationSourceGeometry( "Vehicle", "Earth" ) );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
        std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, initialEphemerisTime, integratorSettings,
              std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ), cowell, dependentVariablesList );

    propagatorSettings->getPrintSettings( )->setPrintDependentVariableData( true );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
        bodies, propagatorSettings );
//
//    std::cout<<dynamicsSimulator.getDependentVariableHistory( ).size( )<<std::endl;
//
//    input_output::writeDataMapToTextFile(
//        dynamicsSimulator.getDependentVariableHistory( ),
//        "radiationPressureDependentVariables.dat",
//        "/home/dominic/Downloads/",
//        "",
//        std::numeric_limits< double >::digits10,
//        std::numeric_limits< double >::digits10,
//        "," );

}
//
//BOOST_AUTO_TEST_SUITE_END( )
//
//} // namespace unit_tests
//
//} // namespace tudat
https://github.com/tudat-team/tudatpy/tags