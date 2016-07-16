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

#include <Tudat/Basics/testMacros.h>
#include <Tudat/Astrodynamics/Propagators/dynamicsSimulator.h>
#include <Tudat/External/SpiceInterface/spiceEphemeris.h>
#include <Tudat/External/SpiceInterface/spiceRotationalEphemeris.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/SimulationSetup/body.h>
#include <Tudat/SimulationSetup/createAccelerationModels.h>
#include <Tudat/SimulationSetup/createMassRateModels.h>
#include <Tudat/SimulationSetup/defaultBodies.h>

#include <iostream>
#include <limits>
#include <string>

#include <Eigen/Core>

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_thrust_acceleration )

BOOST_AUTO_TEST_CASE( testConstantThrustAcceleration )
{
    using namespace tudat;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap;

    // Create vehicle objects.
    double vehicleMass = 5.0E3;
    bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );
    bodyMap[ "Vehicle" ]->setEphemeris(
                boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, basic_mathematics::Vector6d  > >( ),
                    "SSB" ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    Eigen::Vector3d thrustDirection;
    thrustDirection << -1.4, 2.4, 5,6;

    double thrustMagnitude = 1.0E3;
    double specificImpulse = 250.0;

    double massRate = thrustMagnitude / ( specificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );
    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                       boost::make_shared< CustomThrustDirectionSettings >(
                                                           boost::lambda::constant( thrustDirection ) ),
                                                       boost::make_shared< ConstantThrustMagnitudeSettings >(
                                                           thrustMagnitude, specificImpulse ) ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "SSB" );

    // Set initial state
    basic_mathematics::Vector6d systemInitialState = basic_mathematics::Vector6d::Zero( );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, 0.1 );

    for( unsigned int i = 0; i < 2; i++ )
    {
        if( i == 0 )
        {
            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, translationalPropagatorSettings, true, false, false );

            // Retrieve numerical solutions for state and dependent variables
            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
                    dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            Eigen::Vector3d constantAcceleration = thrustDirection.normalized( ) * thrustMagnitude / vehicleMass;
            for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
                 numericalSolution.begin( ); outputIterator != numericalSolution.end( ); outputIterator++ )
            {
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( outputIterator->second.segment( 0, 3 ) ),
                            ( 0.5 * constantAcceleration * std::pow( outputIterator->first, 2.0 ) ), 1.0E-12 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( outputIterator->second.segment( 3, 3 ) ),
                            ( constantAcceleration * outputIterator->first ), 1.0E-12 );
            }
        }
        else if( i == 1 )
        {
            std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
            massRateModels[ "Vehicle" ] = (
                        createMassRateModel( "Vehicle", boost::make_shared< FromThrustMassModelSettings >( 1 ),
                                             bodyMap, accelerationModelMap ) );

            boost::shared_ptr< PropagatorSettings< double > > massPropagatorSettings =
                    boost::make_shared< MassPropagatorSettings< double > >(
                        boost::assign::list_of( "Vehicle" ), massRateModels,
                        ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

            std::vector< boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
            propagatorSettingsVector.push_back( translationalPropagatorSettings );
            propagatorSettingsVector.push_back( massPropagatorSettings );

            boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, false, false );

            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
                    dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
                 numericalSolution.begin( ); outputIterator != numericalSolution.end( ); outputIterator++ )
            {
                double currentMass = vehicleMass - outputIterator->first * massRate;

                Eigen::Vector3d currentVelocity = thrustDirection.normalized( ) *
                        specificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION *
                        std::log( vehicleMass / currentMass );
                BOOST_CHECK_CLOSE_FRACTION( outputIterator->second( 6 ), currentMass, 1.0E-12 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( outputIterator->second.segment( 3, 3 ) ), currentVelocity, 1.0E-11 );

            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testFromEngineThrustAcceleration )
{
    using namespace tudat;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap;

    // Create vehicle objects.
    double vehicleMass = 5.0E3;
    double dryVehicleMass = 2.0E3;

    bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );
    bodyMap[ "Vehicle" ]->setEphemeris(
                boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, basic_mathematics::Vector6d  > >( ),
                    "SSB" ) );


    double thrustMagnitude1 = 1.0E3;
    double specificImpulse1 = 250.0;
    double massFlow1 = propulsion::computePropellantMassRateFromSpecificImpulse(
                thrustMagnitude1, specificImpulse1 );


    double thrustMagnitude2 = 2.0E3;
    double specificImpulse2 = 300.0;
    double massFlow2 = propulsion::computePropellantMassRateFromSpecificImpulse(
                thrustMagnitude2, specificImpulse2 );

    boost::shared_ptr< system_models::VehicleSystems > vehicleSystems = boost::make_shared<
            system_models::VehicleSystems >( dryVehicleMass );
    boost::shared_ptr< system_models::EngineModel > vehicleEngineModel1 =
            boost::make_shared< system_models::DirectEngineModel >( specificImpulse1, boost::lambda::constant( massFlow1 ) );
    boost::shared_ptr< system_models::EngineModel > vehicleEngineModel2 =
            boost::make_shared< system_models::DirectEngineModel >( specificImpulse2, boost::lambda::constant( massFlow2 ) );
    vehicleSystems->setEngineModel( vehicleEngineModel1, "Engine1" );
    vehicleSystems->setEngineModel( vehicleEngineModel2, "Engine2" );
    bodyMap.at( "Vehicle" )->setVehicleSystems( vehicleSystems );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    for( unsigned int i = 0; i < 4; i++ )
    {
        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        Eigen::Vector3d thrustDirection;
        thrustDirection << -1.4, 2.4, 5,6;

        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        // Define acceleration model settings.
        switch( i )
        {
        case 0:
        {
            accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                               boost::make_shared< CustomThrustDirectionSettings >(
                                                                   boost::lambda::constant( thrustDirection ) ),
                                                               boost::make_shared< FromEngineThrustMagnitudeSettings >(
                                                                   1, "" ) ) );
            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
            break;
        }
        case 1:
        {
            accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                               boost::make_shared< CustomThrustDirectionSettings >(
                                                                   boost::lambda::constant( thrustDirection ) ),
                                                               boost::make_shared< FromEngineThrustMagnitudeSettings >(
                                                                   0, "Engine1" ) ) );
            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
            break;
        }
        case 2:
        {
            accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                               boost::make_shared< CustomThrustDirectionSettings >(
                                                                   boost::lambda::constant( thrustDirection ) ),
                                                               boost::make_shared< FromEngineThrustMagnitudeSettings >(
                                                                   0, "Engine2" ) ) );
            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
            break;
        }
        case 3:
        {
            accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                               boost::make_shared< CustomThrustDirectionSettings >(
                                                                   boost::lambda::constant( thrustDirection ) ),
                                                               boost::make_shared< FromEngineThrustMagnitudeSettings >(
                                                                   0, "Engine1" ) ) );
            accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                               boost::make_shared< CustomThrustDirectionSettings >(
                                                                   boost::lambda::constant( thrustDirection ) ),
                                                               boost::make_shared< FromEngineThrustMagnitudeSettings >(
                                                                   0, "Engine2" ) ) );
            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
            break;
        }
        }

        bodiesToPropagate.push_back( "Vehicle" );
        centralBodies.push_back( "SSB" );

        // Set initial state
        basic_mathematics::Vector6d systemInitialState = basic_mathematics::Vector6d::Zero( );

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                boost::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings );
        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, 0.1 );

        std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;

        double totalMassRate, totalThrust;
        switch( i )
        {
        case 0:
        {
            massRateModels[ "Vehicle" ] = (
                        createMassRateModel( "Vehicle", boost::make_shared< FromThrustMassModelSettings >( 1 ),
                                             bodyMap, accelerationModelMap ) );
            totalMassRate = massFlow1 + massFlow2;
            totalThrust = thrustMagnitude1 + thrustMagnitude2;
            break;
        }
        case 1:
        {
            massRateModels[ "Vehicle" ] = (
                        createMassRateModel( "Vehicle", boost::make_shared< FromThrustMassModelSettings >( 0, "Engine1" ),
                                             bodyMap, accelerationModelMap ) );
            totalMassRate = massFlow1;
            totalThrust = thrustMagnitude1;
            break;
        }
        case 2:
        {
            massRateModels[ "Vehicle" ] = (
                        createMassRateModel( "Vehicle", boost::make_shared< FromThrustMassModelSettings >( 0, "Engine2" ),
                                             bodyMap, accelerationModelMap ) );
            totalMassRate = massFlow2;
            totalThrust = thrustMagnitude2;
            break;
        }
        case 3:
        {
            massRateModels[ "Vehicle" ] = (
                        createMassRateModel( "Vehicle", boost::make_shared< FromThrustMassModelSettings >( 0, "Engine1" ),
                                             bodyMap, accelerationModelMap ) );
            totalMassRate = massFlow1;
            totalThrust = thrustMagnitude1 + thrustMagnitude2;
            break;
        }
        }

        double totalSpecificImpulse = totalThrust / ( physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * totalMassRate );

        boost::shared_ptr< PropagatorSettings< double > > massPropagatorSettings =
                boost::make_shared< MassPropagatorSettings< double > >(
                    boost::assign::list_of( "Vehicle" ), massRateModels,
                    ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

        std::vector< boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( translationalPropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                boost::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );

        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
                dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
             numericalSolution.begin( ); outputIterator != numericalSolution.end( ); outputIterator++ )
        {
            double currentMass = vehicleMass - outputIterator->first * totalMassRate;

            Eigen::Vector3d currentVelocity = thrustDirection.normalized( ) *
                    totalSpecificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION *
                    std::log( vehicleMass / currentMass );
            BOOST_CHECK_CLOSE_FRACTION( outputIterator->second( 6 ), currentMass, 1.0E-12 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        ( outputIterator->second.segment( 3, 3 ) ), currentVelocity, 1.0E-11 );

        }
    }
}

BOOST_AUTO_TEST_CASE( testRadialAndVelocityThrustAcceleration )
{
    using namespace tudat;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;

    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");

    double thrustMagnitude = 1.0E3;
    double specificImpulse = 250.0;

    for( unsigned int i = 0; i < 2; i++ )
    {
        // Create Earth object
        simulation_setup::NamedBodyMap bodyMap;

        // Create vehicle objects.
        double vehicleMass = 5.0E3;
        bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
        bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );
        bodyMap[ "Vehicle" ]->setEphemeris(
                    boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, basic_mathematics::Vector6d  > >( ),
                        "Earth" ) );
        bodyMap[ "Earth" ] = boost::make_shared< Body >( );

        bodyMap[ "Earth" ]->setEphemeris(
                    boost::make_shared< ephemerides::SpiceEphemeris >( "Sun", "SSB", false, false ) );
        bodyMap[ "Earth" ]->setGravityFieldModel( boost::make_shared< gravitation::GravityFieldModel >(
                                                      spice_interface::getBodyGravitationalParameter( "Earth" ) ) );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        bool isThurstInVelocityDirection;

        if( i == 0 )
        {
            isThurstInVelocityDirection = 0;
        }
        else
        {
            isThurstInVelocityDirection = 1;
        }

        // Define acceleration model settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                           boost::make_shared< ThrustDirectionFromStateGuidanceSettings >(
                                                               "Earth", isThurstInVelocityDirection, 1  ),
                                                           boost::make_shared< ConstantThrustMagnitudeSettings >(
                                                               thrustMagnitude, specificImpulse ) ) );
        if( i == 1 )
        {
            accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        }

        accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

        bodiesToPropagate.push_back( "Vehicle" );
        centralBodies.push_back( "Earth" );

        // Set initial state
        double radius = 1.0E3;
        double circularVelocity = std::sqrt( radius * thrustMagnitude / vehicleMass );
        basic_mathematics::Vector6d systemInitialState = basic_mathematics::Vector6d::Zero( );

        if( i == 0 )
        {
            systemInitialState( 0 ) = radius;
            systemInitialState( 4 ) = circularVelocity;
        }
        else
        {
            systemInitialState( 0 ) = 8.0E6;
            systemInitialState( 4 ) = 7.5E3;

        }

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        boost::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings;
        if( i == 1 )
        {
            std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
            dependentVariables.push_back(
                        boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
            dependentVariableSaveSettings = boost::make_shared< DependentVariableSaveSettings >( dependentVariables );

        }
        boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                boost::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings,
                  cowell, dependentVariableSaveSettings );
        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, 0.1 );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, translationalPropagatorSettings, true, false, false );

        // Retrieve numerical solutions for state and dependent variables
        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
                dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableSolution =
                dynamicsSimulator.getDependentVariableHistory( );

        if( i == 0 )
        {
            double angularVelocity = circularVelocity / radius;

            for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
                 numericalSolution.begin( ); outputIterator != numericalSolution.end( ); outputIterator++ )
            {
                double currentAngle = angularVelocity * outputIterator->first;

                BOOST_CHECK_CLOSE_FRACTION(
                            ( outputIterator->second.segment( 0, 3 ).norm( ) ), radius, 1.0E-10 * radius );
                BOOST_CHECK_CLOSE_FRACTION(
                            ( outputIterator->second.segment( 3, 3 ).norm( ) ), circularVelocity, 1.0E-10 * circularVelocity );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 0 ) - radius * std::cos( currentAngle ) ), 1.0E-10 * radius );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 1 ) - radius * std::sin( currentAngle ) ), 1.0E-10 * radius  );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 2 ) ), 1.0E-15 );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 3 ) + circularVelocity * std::sin( currentAngle ) ), 1.0E-10 * circularVelocity  );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 4 ) - circularVelocity * std::cos( currentAngle ) ), 1.0E-10 * circularVelocity  );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 5 ) ), 1.0E-15 );
            }
        }
        else if( i == 1 )
        {
            for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
                 numericalSolution.begin( ); outputIterator != numericalSolution.end( ); outputIterator++ )
            {
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( -1.0 * thrustMagnitude / vehicleMass * outputIterator->second.segment( 3, 3 ).normalized( ) ),
                                                   ( dependentVariableSolution.at( outputIterator->first ) ), 1.0E-14 );

            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testThrustAccelerationFromExistingRotation )
{
    using namespace tudat;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;

    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");

    double thrustMagnitude = 1.0E3;
    double specificImpulse = 250.0;

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap;

    // Create vehicle objects.
    double vehicleMass = 5.0E3;
    bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );
    bodyMap[ "Vehicle" ]->setEphemeris(
                boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, basic_mathematics::Vector6d  > >( ),
                    "Earth" ) );
    bodyMap[ "Vehicle" ]->setRotationalEphemeris(
                boost::make_shared< ephemerides::SpiceRotationalEphemeris >( "ECLIPJ2000", "IAU_MOON" ) );
    bodyMap[ "Earth" ] = boost::make_shared< Body >( );

    bodyMap[ "Earth" ]->setEphemeris(
                boost::make_shared< ephemerides::SpiceEphemeris >( "Sun", "SSB", false, false ) );
    bodyMap[ "Earth" ]->setGravityFieldModel( boost::make_shared< gravitation::GravityFieldModel >(
                                                  spice_interface::getBodyGravitationalParameter( "Earth" ) ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                       boost::make_shared< ThrustDirectionGuidanceSettings >(
                                                           thrust_direction_from_existing_body_orientation, "Earth" ),
                                                       boost::make_shared< ConstantThrustMagnitudeSettings >(
                                                           thrustMagnitude, specificImpulse ) ) );
    accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Set initial state
    basic_mathematics::Vector6d systemInitialState = basic_mathematics::Vector6d::Zero( );


    systemInitialState( 0 ) = 8.0E6;
    systemInitialState( 4 ) = 7.5E3;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    boost::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings;


    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
    dependentVariableSaveSettings = boost::make_shared< DependentVariableSaveSettings >( dependentVariables );


    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings,
              cowell, dependentVariableSaveSettings );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, 2.5 );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, translationalPropagatorSettings, true, false, false );

    // Retrieve numerical solutions for state and dependent variables
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput =
            dynamicsSimulator.getDependentVariableHistory( );

    Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( );
    double thrustAcceleration = thrustMagnitude / vehicleMass;
    Eigen::Quaterniond rotationToInertialFrame;
    for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
         dependentVariableOutput.begin( ); outputIterator != dependentVariableOutput.end( ); outputIterator++ )
    {
        rotationToInertialFrame = bodyMap.at( "Vehicle" )->getRotationalEphemeris( )->getRotationToBaseFrame(
                    outputIterator->first );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( thrustAcceleration * ( rotationToInertialFrame * bodyFixedThrustDirection ) ),
                    outputIterator->second, 1.0E-15 );

    }
}



BOOST_AUTO_TEST_SUITE_END( )

}

}


