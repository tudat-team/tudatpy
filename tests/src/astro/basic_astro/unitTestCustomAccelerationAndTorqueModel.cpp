/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *
 *    Notes:
 *      Test tolerance was set at 5.0e-15 (or 5.0e-7 for floats) instead of epsilon due to
 *      rounding errors in Eigen types with entries over a number of orders of magnitude,
 *      presumably causing the observed larger than epsilon relative differences between
 *      expected and computed values.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/simulation.h"

namespace tudat
{
namespace unit_tests
{


BOOST_AUTO_TEST_SUITE( test_customAccelerationAndTorqueModels )

Eigen::Vector3d customAcceleration( const double time )
{
    return ( Eigen::Vector3d( ) << 1.0E-6 * std::sin ( 2.0 * mathematical_constants::PI * time / 1.0E4 + 0.1 ),
             -2.0E-6 * std::sin ( 2.0 * mathematical_constants::PI * time / 3.0E4 + 0.5 ),
             -8.0E-7 * std::sin ( 2.0 * mathematical_constants::PI * time / 5.0E4 + 1.5 ) ).finished( );
}

BOOST_AUTO_TEST_CASE( test_customAccelerationModelCreation )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::unit_conversions;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Create body objects.
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun" };
    BodyListSettings bodySettings = getDefaultBodySettings(
                bodiesToCreate, "Sun", "ECLIPJ2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back(
                std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Sun" ].push_back(
                std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );


    std::map< double, Eigen::Vector3d > customAccelerationMap;
    double timeStep = 30.0;
    double currentTime = -120.0;
    while( currentTime < simulationEndEpoch + 120.0 )
    {
        customAccelerationMap[ currentTime ] = customAcceleration( currentTime );
        currentTime += timeStep;
    }
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > customAccelerationInterpolator =
            interpolators::createOneDimensionalInterpolator(
                customAccelerationMap, std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) );

    std::function< Eigen::Vector3d( const double ) > customAccelerationFunction =
            std::bind( static_cast< Eigen::Vector3d( interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d >::* )
                       ( const double ) >( &interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d >::interpolate ),
                       customAccelerationInterpolator, std::placeholders::_1 );

    std::function< double( const double) > customAccelerationScalingFunction =
            tudat::simulation_setup::getOccultationFunction(
                bodies, "Sun", "Earth", "Vehicle" );

    accelerationsOfVehicle[ "Earth" ].push_back(
                std::make_shared< CustomAccelerationSettings >(
                    customAccelerationFunction, customAccelerationScalingFunction ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set initial conditions for the Vehicle satellite that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to Cartesian
    // elements.

    // Set Keplerian elements for Vehicle.
    Eigen::Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 0.0 );
    vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) =
            convertDegreesToRadians( 235.7 );
    vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) =
            convertDegreesToRadians( 23.4 );
    vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Convert Vehicle state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                vehicleInitialStateInKeplerianElements,
                earthGravitationalParameter );

    // Create propagator settings.
    // Set variables to save
    std::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings;
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::custom_acceleration, "Vehicle", "Earth", 0 ) );
    dependentVariableSaveSettings = std::make_shared< DependentVariableSaveSettings >( dependentVariables );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch,
              cowell, dependentVariableSaveSettings );

    // Create numerical integrator settings.
    double simulationStartEpoch = 0.0;
    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, integratorSettings, propagatorSettings );
    std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );
    std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    double sourceBodyRadius = bodies.at( "Sun" )->getShapeModel( )->getAverageRadius( );
    double occultingBodyRadius = bodies.at( "Earth" )->getShapeModel( )->getAverageRadius( );

    for( auto it : dependentVariableHistory )
    {
        double shadowFunction =
                mission_geometry::computeShadowFunction(
                    bodies.at( "Sun" )->getStateInBaseFrameFromEphemeris( it.first ).segment( 0, 3 ), sourceBodyRadius,
                    bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( it.first ).segment( 0, 3 ), occultingBodyRadius,
                    stateHistory.at( it.first ).segment( 0, 3 ) +
                    bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( it.first ).segment( 0, 3 ) );

        Eigen::Vector3d expectedAcceleration = shadowFunction * customAcceleration( it.first );
        Eigen::Vector3d savedAcceleration = it.second;

//        std::cout<<savedAcceleration.transpose( )<<std::endl;
//        std::cout<<expectedAcceleration.transpose( )<<std::endl<<std::endl;

        BOOST_CHECK_SMALL( savedAcceleration( 0 ) - expectedAcceleration( 0 ), 1.0E-19 );
        BOOST_CHECK_SMALL( savedAcceleration( 1 ) - expectedAcceleration( 1 ), 1.0E-19 );
        BOOST_CHECK_SMALL( savedAcceleration( 2 ) - expectedAcceleration( 2 ), 1.0E-19 );
    }
}

Eigen::Vector3d customTorque( const double time )
{
    return ( Eigen::Vector3d( ) << 1.0E-6 * std::sin ( 2.0 * mathematical_constants::PI * time / 1.0E4 + 0.1 ),
             -2.0E-6 * std::sin ( 2.0 * mathematical_constants::PI * time / 3.0E4 + 0.5 ),
             -8.0E-7 * std::sin ( 2.0 * mathematical_constants::PI * time / 5.0E4 + 1.5 ) ).finished( );
}

BOOST_AUTO_TEST_CASE( test_customTorqueModelCreation )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::unit_conversions;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Create body objects.
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun" };
    BodyListSettings bodySettings = getDefaultBodySettings(
                bodiesToCreate, "Sun", "ECLIPJ2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    Eigen::Matrix3d bodyInertiaTensor = 100.0 * Eigen::Matrix3d::Identity( );
    bodyInertiaTensor( 2, 2 ) *= 2.0;
    bodyInertiaTensor( 1, 1 ) *= 1.5;

    bodies.at( "Vehicle" )->setBodyInertiaTensor( bodyInertiaTensor );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back(
                std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Sun" ].push_back(
                std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set initial conditions for the Vehicle satellite that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to Cartesian
    // elements.

    // Set Keplerian elements for Vehicle.
    Eigen::Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 0.0 );
    vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) =
            convertDegreesToRadians( 235.7 );
    vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) =
            convertDegreesToRadians( 23.4 );
    vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Convert Vehicle state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialTranslationalState = convertKeplerianToCartesianElements(
                vehicleInitialStateInKeplerianElements,
                earthGravitationalParameter );

    // Create propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialTranslationalState, simulationEndEpoch );


    Eigen::VectorXd systemInitialRotationalState = Eigen::VectorXd::Zero( 7 );
    systemInitialRotationalState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
    systemInitialRotationalState.segment( 4, 3 ) = Eigen::Vector3d::Constant( 1.0E-5 );

    // Create torque models
    SelectedTorqueMap torqueMap;
    torqueMap[ "Vehicle" ][ "Earth" ].push_back( std::make_shared< TorqueSettings >(
                                                     basic_astrodynamics::second_order_gravitational_torque ) );

    std::map< double, Eigen::Vector3d > customTorqueMap;
    double timeStep = 30.0;
    double currentTime = -120.0;
    while( currentTime < simulationEndEpoch + 120.0 )
    {
        customTorqueMap[ currentTime ] = customTorque( currentTime );
        currentTime += timeStep;
    }
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > customTorqueInterpolator =
            interpolators::createOneDimensionalInterpolator(
                customTorqueMap, std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) );

    std::function< Eigen::Vector3d( const double ) > customTorqueFunction =
            std::bind( static_cast< Eigen::Vector3d( interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d >::* )
                       ( const double ) >( &interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d >::interpolate ),
                       customTorqueInterpolator, std::placeholders::_1 );

    std::function< double( const double) > customTorqueScalingFunction =
            tudat::simulation_setup::getOccultationFunction(
                bodies, "Sun", "Earth", "Vehicle" );

    torqueMap[ "Vehicle" ][ "Earth" ].push_back( std::make_shared< CustomTorqueSettings >(
                                                     customTorqueFunction, customTorqueScalingFunction ) );

    basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                bodies, torqueMap, bodiesToPropagate );


    // Define propagator settings.
    std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
            std::make_shared< RotationalStatePropagatorSettings< double > >
            ( torqueModelMap, bodiesToPropagate, systemInitialRotationalState, std::make_shared< PropagationTimeTerminationSettings >(
                  simulationEndEpoch ) );

    // Set variables to save
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                std::make_shared< SingleTorqueDependentVariableSaveSettings >(
                    custom_torque, "Vehicle", "Earth" ) );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariables );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );
    propagatorSettingsList.push_back( rotationalPropagatorSettings );

    std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsList,
                std::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ),
                dependentVariablesToSave );

    // Create numerical integrator settings.
    double simulationStartEpoch = 0.0;
    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, integratorSettings, propagatorSettings );
    std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );
    std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    double sourceBodyRadius = bodies.at( "Sun" )->getShapeModel( )->getAverageRadius( );
    double occultingBodyRadius = bodies.at( "Earth" )->getShapeModel( )->getAverageRadius( );

    for( auto it : dependentVariableHistory )
    {
        double shadowFunction =
                mission_geometry::computeShadowFunction(
                    bodies.at( "Sun" )->getStateInBaseFrameFromEphemeris( it.first ).segment( 0, 3 ), sourceBodyRadius,
                    bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( it.first ).segment( 0, 3 ), occultingBodyRadius,
                    stateHistory.at( it.first ).segment( 0, 3 ) +
                    bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( it.first ).segment( 0, 3 ) );

        Eigen::Vector3d expectedTorque = shadowFunction * customTorque( it.first );
        Eigen::Vector3d savedTorque = it.second;
        BOOST_CHECK_SMALL( savedTorque( 0 ) - expectedTorque( 0 ), 1.0E-19 );
        BOOST_CHECK_SMALL( savedTorque( 1 ) - expectedTorque( 1 ), 1.0E-19 );
        BOOST_CHECK_SMALL( savedTorque( 2 ) - expectedTorque( 2 ), 1.0E-19 );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
