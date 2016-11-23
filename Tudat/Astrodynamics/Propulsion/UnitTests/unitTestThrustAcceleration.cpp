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

#include "Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include <Tudat/Basics/testMacros.h>
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include <Tudat/External/SpiceInterface/spiceEphemeris.h>
#include <Tudat/External/SpiceInterface/spiceRotationalEphemeris.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/SimulationSetup/EnvironmentSetup/body.h>
#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
#include <Tudat/SimulationSetup/PropagationSetup/createMassRateModels.h>
#include <Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h>

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
                                                       boost::make_shared< ConstantThrustEngineSettings >(
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

BOOST_AUTO_TEST_CASE( testFromEngineThrustAcceleration )
{
    using namespace tudat;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;


    for( unsigned int i = 0; i < 4; i++ )
    {
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
                boost::make_shared< system_models::DirectEngineModel >(
                    boost::lambda::constant( specificImpulse1 ), boost::lambda::constant( massFlow1 ) );
        boost::shared_ptr< system_models::EngineModel > vehicleEngineModel2 =
                boost::make_shared< system_models::DirectEngineModel >(
                    boost::lambda::constant( specificImpulse2 ), boost::lambda::constant( massFlow2 ) );
        vehicleSystems->setEngineModel( vehicleEngineModel1, "Engine1" );
        vehicleSystems->setEngineModel( vehicleEngineModel2, "Engine2" );
        bodyMap.at( "Vehicle" )->setVehicleSystems( vehicleSystems );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

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
                                                               boost::make_shared< FromBodyThrustEngineSettings >(
                                                                   1, "" ) ) );
            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
            break;
        }
        case 1:
        {
            accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                               boost::make_shared< CustomThrustDirectionSettings >(
                                                                   boost::lambda::constant( thrustDirection ) ),
                                                               boost::make_shared< FromBodyThrustEngineSettings >(
                                                                   0, "Engine1" ) ) );
            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
            break;
        }
        case 2:
        {
            accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                               boost::make_shared< CustomThrustDirectionSettings >(
                                                                   boost::lambda::constant( thrustDirection ) ),
                                                               boost::make_shared< FromBodyThrustEngineSettings >(
                                                                   0, "Engine2" ) ) );
            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
            break;
        }
        case 3:
        {
            accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                               boost::make_shared< CustomThrustDirectionSettings >(
                                                                   boost::lambda::constant( thrustDirection ) ),
                                                               boost::make_shared< FromBodyThrustEngineSettings >(
                                                                   0, "Engine1" ) ) );
            accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                               boost::make_shared< CustomThrustDirectionSettings >(
                                                                   boost::lambda::constant( thrustDirection ) ),
                                                               boost::make_shared< FromBodyThrustEngineSettings >(
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
                                                           boost::make_shared< ConstantThrustEngineSettings >(
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
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( -1.0 * thrustMagnitude / vehicleMass * outputIterator->second.segment( 3, 3 ).normalized( ) ),
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
    Eigen::Vector3d bodyFixedThrustDirection = ( Eigen::Vector3d( ) << 1.4, 3.1, -0.5 ).finished( ).normalized( );

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                       boost::make_shared< ThrustDirectionGuidanceSettings >(
                                                           thrust_direction_from_existing_body_orientation, "Earth" ),
                                                       boost::make_shared< ConstantThrustEngineSettings >(
                                                           thrustMagnitude, specificImpulse, bodyFixedThrustDirection ) ) );
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

    double thrustAcceleration = thrustMagnitude / vehicleMass;
    Eigen::Quaterniond rotationToInertialFrame;
    for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
         dependentVariableOutput.begin( ); outputIterator != dependentVariableOutput.end( ); outputIterator++ )
    {
        rotationToInertialFrame = bodyMap.at( "Vehicle" )->getRotationalEphemeris( )->getRotationToBaseFrame(
                    outputIterator->first );
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION(
                        ( thrustAcceleration * ( rotationToInertialFrame * bodyFixedThrustDirection )( i ) ),
                        outputIterator->second( i ), 2.0E-15 );
        }
    }
}

BOOST_AUTO_TEST_CASE( testConcurrentThrustAndAerodynamicAcceleration )
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
    double vehicleMass = 5.0E3;
    bodyMap[ "Apollo" ]->setConstantBodyMass( vehicleMass );

    // Create vehicle aerodynamic coefficients
    bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );
    bodyMap[ "Apollo" ]->setEphemeris(
                boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    boost::shared_ptr< interpolators::OneDimensionalInterpolator<
                    double, basic_mathematics::Vector6d  > >( ), "Earth" ) );

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

    double thrustMagnitude = 1.0E-3;
    double specificImpulse = 250.0;
    //double massRate = thrustMagnitude / ( specificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );
    accelerationsOfApollo[ "Apollo" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                      boost::make_shared< ThrustDirectionGuidanceSettings >(
                                                          thrust_direction_from_existing_body_orientation, "Earth" ),
                                                      boost::make_shared< ConstantThrustEngineSettings >(
                                                          thrustMagnitude, specificImpulse ) ) );

    accelerationMap[ "Apollo" ] = accelerationsOfApollo;

    bodiesToPropagate.push_back( "Apollo" );
    centralBodies.push_back( "Earth" );

    // Set initial state
    basic_mathematics::Vector6d systemInitialState = apolloInitialState;


    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    setTrimmedConditions( bodyMap.at( "Apollo" ) );

    // Define list of dependent variables to save.
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "Apollo" ) );
    dependentVariables.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "Apollo", reference_frames::angle_of_attack ) );
    dependentVariables.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "Apollo", reference_frames::angle_of_sideslip ) );
    dependentVariables.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "Apollo", reference_frames::bank_angle ) );
    dependentVariables.push_back(
                boost::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                    "Apollo", reference_frames::inertial_frame, reference_frames::body_frame ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    rotation_matrix_to_body_fixed_frame_variable, "Apollo" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic, "Apollo", "Earth", 0 ) );
    dependentVariables.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    thrust_acceleration, "Apollo", "Apollo", 0 ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_force_coefficients_dependent_variable, "Apollo" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_moment_coefficients_dependent_variable, "Apollo" ) );

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

    // Iterate over results for dependent variables, and check against computed values.
    Eigen::Matrix3d rotationToBodyFixedFrame1, rotationToBodyFixedFrame2;
    Eigen::Vector3d expectedThrustDirection, computedThrustDirection;
    Eigen::Vector3d aerodynamicCoefficients;

    Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( );

    boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > vehicelCoefficientInterface =
            bodyMap.at( "Apollo" )->getAerodynamicCoefficientInterface( );

    for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
         variableIterator != dependentVariableSolution.end( ); variableIterator++ )
    {
        vehicelCoefficientInterface->updateCurrentCoefficients(
                    boost::assign::list_of( variableIterator->second( 0 ) )(
                        variableIterator->second( 1 ) )( variableIterator->second( 2 ) ) );
        aerodynamicCoefficients = vehicelCoefficientInterface->getCurrentForceCoefficients( );

        rotationToBodyFixedFrame1 = getMatrixFromVectorRotationRepresentation( variableIterator->second.segment( 4, 9 ) );
        rotationToBodyFixedFrame2 = getMatrixFromVectorRotationRepresentation( variableIterator->second.segment( 13, 9 ) );

        expectedThrustDirection = rotationToBodyFixedFrame1.transpose( ) * bodyFixedThrustDirection;
        computedThrustDirection = variableIterator->second.segment( 25, 3 ).normalized( );

        // Check thrust magnitude
        BOOST_CHECK_CLOSE_FRACTION(
                    ( variableIterator->second.segment( 25, 3 ) ).norm( ), thrustMagnitude / vehicleMass,
                    2.0 * std::numeric_limits< double >::epsilon( ) );
        for( unsigned int i = 0; i < 3; i ++ )
        {
            // Check rotation matrices
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( rotationToBodyFixedFrame1( i, j ) - rotationToBodyFixedFrame2( i, j ) ),
                                    8.0 * std::numeric_limits< double >::epsilon( ) );
            }
            // Check thrust direction
            BOOST_CHECK_SMALL( std::fabs( expectedThrustDirection( i ) - computedThrustDirection( i ) ),
                                15.0 * std::numeric_limits< double >::epsilon( ) );

            // Check aerodynamic coefficients
            BOOST_CHECK_SMALL(
                        std::fabs( variableIterator->second( 28 ) - aerodynamicCoefficients( 0 ) ), 1.0E-10 );
            BOOST_CHECK_SMALL(
                        std::fabs( variableIterator->second( 29 ) - aerodynamicCoefficients( 1 ) ), 1.0E-10 );
            BOOST_CHECK_SMALL(
                        std::fabs( variableIterator->second( 30 ) - aerodynamicCoefficients( 2 ) ), 1.0E-10 );

            // Check trimmed condition (y-term)/symmetric vehicle shape (x- and z-term).
            BOOST_CHECK_SMALL(
                        std::fabs( variableIterator->second( 31 ) ), 1.0E-14 );
            BOOST_CHECK_SMALL(
                        std::fabs( variableIterator->second( 32 ) ), 1.0E-10 );
            BOOST_CHECK_SMALL(
                        std::fabs( variableIterator->second( 33 ) ), 1.0E-14 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
