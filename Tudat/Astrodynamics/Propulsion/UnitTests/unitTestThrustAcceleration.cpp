/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/External/SpiceInterface/spiceRotationalEphemeris.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/multiDimensionalArrayReader.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
#include "Tudat/SimulationSetup/PropagationSetup/createMassRateModels.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
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
                    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d  > >( ),
                    "SSB" ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    Eigen::Vector3d thrustDirection;
    thrustDirection << -1.4, 2.4, 5.6;

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
    Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );

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

        boost::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings =
                boost::make_shared< MassPropagatorSettings< double > >(
                    boost::assign::list_of( "Vehicle" ), massRateModels,
                    ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

        std::vector< boost::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( translationalPropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        boost::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
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


//! In this unit test, the thrust acceleration is tested for the case where the thrust force is taken from a (set of)
//! engine objects stored in the VehicleSystems object. This is tested for a single engine (out of two), two engines, as
//! well as the (unrealistic) case of thrust from 2 engines and mass flow from only of these engines (for the purposes of
//! mass propagation).
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
                        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d  > >( ),
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
        thrustDirection << -1.4, 2.4, 5.6;

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
        Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );

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

        boost::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings =
                boost::make_shared< MassPropagatorSettings< double > >(
                    boost::assign::list_of( "Vehicle" ), massRateModels,
                    ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

        std::vector< boost::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( translationalPropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        boost::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
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

//! In this unit test, the thrust force set to be colinear with the position and velocity vectors is checked.
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
    spice_interface::loadStandardSpiceKernels( );

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
                        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d  > >( ),
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
        Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );

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

            // Check if the spacecraft is in a circular orbit, with the orbit being maintained exactly by the radial thrust.
            for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
                 numericalSolution.begin( ); outputIterator != numericalSolution.end( ); outputIterator++ )
            {
                double currentAngle = angularVelocity * outputIterator->first;

                // Check constancy of position and velocoty scalars.
                BOOST_CHECK_CLOSE_FRACTION(
                            ( outputIterator->second.segment( 0, 3 ).norm( ) ), radius, 1.0E-10 * radius );
                BOOST_CHECK_CLOSE_FRACTION(
                            ( outputIterator->second.segment( 3, 3 ).norm( ) ),
                            circularVelocity, 1.0E-10 * circularVelocity );

                // Check whether orbit is planar (in xy-plane) with a constant angular motion at the prescribed mean motion.
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 0 ) - radius * std::cos( currentAngle ) ), 1.0E-10 * radius );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 1 ) - radius * std::sin( currentAngle ) ), 1.0E-10 * radius );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 2 ) ), 1.0E-15 );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 3 ) + circularVelocity * std::sin( currentAngle ) ),
                            1.0E-10 * circularVelocity  );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 4 ) - circularVelocity * std::cos( currentAngle ) ),
                            1.0E-10 * circularVelocity  );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 5 ) ), 1.0E-15 );
            }
        }
        else if( i == 1 )
        {
            for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
                 numericalSolution.begin( ); outputIterator != numericalSolution.end( ); outputIterator++ )
            {
                // Check if the thrust acceleration is of the correct magnitude, and in the same direction as the velocity.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( -1.0 * thrustMagnitude / vehicleMass * outputIterator->second.segment( 3, 3 ).normalized( ) ),
                            ( dependentVariableSolution.at( outputIterator->first ) ), 1.0E-14 );

            }
        }
    }

}

//! Test the thrust force direction when it is taken from a predetermined vehicle rotation (RotationalEphemeris)
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
    spice_interface::loadStandardSpiceKernels( );



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
                    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d  > >( ),
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
    Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );


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
            // Test thrust direction during propagation against rotational ephemeris of vehicle.
            BOOST_CHECK_CLOSE_FRACTION(
                        ( thrustAcceleration * ( rotationToInertialFrame * bodyFixedThrustDirection )( i ) ),
                        outputIterator->second( i ), 2.0E-15 );
        }
    }
}

//! Test whether the concurrent use of aerodynamic and thrust forces is handled correctly. In particular, this test checks
//! whether the aerodynamic orientation is correctly used to determine the inertial thrust direction.
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
    const Vector6d apolloInitialState = convertKeplerianToCartesianElements(
                apolloInitialStateInKeplerianElements,
                getBodyGravitationalParameter( "Earth" ) );

    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings[ "Earth" ]->gravityFieldSettings =
            boost::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

    for( unsigned int testCase = 0; testCase < 2; testCase++ )
    {
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
                        double, Eigen::Vector6d  > >( ), "Earth" ) );

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
        accelerationsOfApollo[ "Apollo" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                         boost::make_shared< ThrustDirectionGuidanceSettings >(
                                                             thrust_direction_from_existing_body_orientation, "Earth" ),
                                                         boost::make_shared< ConstantThrustEngineSettings >(
                                                             thrustMagnitude, specificImpulse ) ) );

        accelerationMap[ "Apollo" ] = accelerationsOfApollo;

        bodiesToPropagate.push_back( "Apollo" );
        centralBodies.push_back( "Earth" );

        // Set initial state
        Eigen::Vector6d systemInitialState = apolloInitialState;


        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        // Set trim conditions for the vehicle.
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
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        airspeed_dependent_variable, "Apollo", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        local_density_dependent_variable, "Apollo", "Earth" ) );
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
        dependentVariables.push_back(
                    boost::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                        "Apollo", reference_frames::inertial_frame, reference_frames::aerodynamic_frame ) );
        dependentVariables.push_back(
                    boost::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                        "Apollo", reference_frames::aerodynamic_frame, reference_frames::body_frame ) );

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
        Eigen::Vector3d expectedAerodynamicForceDirection;

        Eigen::Vector3d aerodynamicCoefficients;

        Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( );

        boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > vehicelCoefficientInterface =
                bodyMap.at( "Apollo" )->getAerodynamicCoefficientInterface( );

        for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
             variableIterator != dependentVariableSolution.end( ); variableIterator++ )
        {
            double currentMachNumber = variableIterator->second( 0 );
            double currentAngleOfAttack = variableIterator->second( 1 );
            double currentAngleOfSideSlip = variableIterator->second( 2 );

            double currentBankAngle = variableIterator->second( 3 );
            double currentAirspeed = variableIterator->second( 4 );
            double currentDensity = variableIterator->second( 5 );

            BOOST_CHECK_SMALL(
                        std::fabs( currentAngleOfSideSlip ), std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( currentBankAngle ), std::numeric_limits< double >::epsilon( ) );

            rotationToBodyFixedFrame1 =
                    getMatrixFromVectorRotationRepresentation( variableIterator->second.segment( 6, 9 ) );
            rotationToBodyFixedFrame2 =
                    getMatrixFromVectorRotationRepresentation( variableIterator->second.segment( 15, 9 ) );

            Eigen::Vector3d currentAerodynamicAcceleration = variableIterator->second.segment( 24, 3 );
            Eigen::Vector3d currentThrustAcceleration = variableIterator->second.segment( 27, 3 );

            Eigen::Vector3d currentAerodynamicForceCoefficients = variableIterator->second.segment( 30, 3 );
            Eigen::Vector3d currentAerodynamicMomentCoefficients = variableIterator->second.segment( 33, 3 );

            Eigen::Matrix3d currentRotationFromInertialToAerodynamicFrame =
                    getMatrixFromVectorRotationRepresentation( variableIterator->second.segment( 36, 9 ) );

            Eigen::Matrix3d currentRotationFromAerodynamicToBodyFixedFrame =
                    getMatrixFromVectorRotationRepresentation( variableIterator->second.segment( 45, 9 ) );

            vehicelCoefficientInterface->updateCurrentCoefficients(
                        boost::assign::list_of( currentMachNumber )(
                            currentAngleOfAttack )( currentAngleOfSideSlip ) );
            aerodynamicCoefficients = vehicelCoefficientInterface->getCurrentForceCoefficients( );



            expectedThrustDirection = rotationToBodyFixedFrame1.transpose( ) * bodyFixedThrustDirection;
            computedThrustDirection = currentThrustAcceleration.normalized( );

            // Check thrust magnitude
            BOOST_CHECK_CLOSE_FRACTION(
                        ( currentThrustAcceleration ).norm( ), thrustMagnitude / vehicleMass,
                        2.0 * std::numeric_limits< double >::epsilon( ) );




            Eigen::Vector3d expectedAerodynamicAcceleration = 0.5 * currentDensity * currentAirspeed * currentAirspeed *
                    bodyMap.at( "Apollo" )->getAerodynamicCoefficientInterface( )->getReferenceArea( ) *
                    ( currentRotationFromInertialToAerodynamicFrame.inverse( ) ) *
                    currentAerodynamicForceCoefficients / vehicleMass;


            // Check aerodynamic coefficients
            BOOST_CHECK_SMALL(
                        std::fabs( currentAerodynamicForceCoefficients( 0 ) - aerodynamicCoefficients( 0 ) ), 1.0E-10 );
            BOOST_CHECK_SMALL(
                        std::fabs( currentAerodynamicForceCoefficients( 1 ) - aerodynamicCoefficients( 1 ) ), 1.0E-10 );
            BOOST_CHECK_SMALL(
                        std::fabs( currentAerodynamicForceCoefficients( 2 ) - aerodynamicCoefficients( 2 ) ), 1.0E-10 );

            // Check trimmed condition (y-term)/symmetric vehicle shape (x- and z-term).
            BOOST_CHECK_SMALL(
                        std::fabs( currentAerodynamicMomentCoefficients( 0 ) ), 1.0E-14 );
            BOOST_CHECK_SMALL(
                        std::fabs( currentAerodynamicMomentCoefficients( 1 ) ), 1.0E-10 );
            BOOST_CHECK_SMALL(
                        std::fabs( currentAerodynamicMomentCoefficients( 2 ) ), 1.0E-14 );
            Eigen::Matrix3d currentRotationFromAerodynamicToBodyFixedFrame2 =
                    rotationToBodyFixedFrame1 * ( currentRotationFromInertialToAerodynamicFrame.inverse( ) );
            Eigen::Matrix3d computedAerodynamicToBodyFixedFrame =
                    reference_frames::getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
                        currentAngleOfAttack, 0.0 );


            for( unsigned int i = 0; i < 3; i ++ )
            {
                // Check rotation matrices
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( rotationToBodyFixedFrame1( i, j ) - rotationToBodyFixedFrame2( i, j ) ),
                                10.0 * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL(
                                std::fabs( currentRotationFromAerodynamicToBodyFixedFrame( i, j ) -
                                           computedAerodynamicToBodyFixedFrame( i, j ) ),
                                10.0 * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL(
                                std::fabs( currentRotationFromAerodynamicToBodyFixedFrame2( i, j ) -
                                           computedAerodynamicToBodyFixedFrame( i, j ) ),
                                10.0 * std::numeric_limits< double >::epsilon( ) );
                }
                // Check thrust direction
                BOOST_CHECK_SMALL( std::fabs( expectedThrustDirection( i ) - computedThrustDirection( i ) ),
                                   15.0 * std::numeric_limits< double >::epsilon( ) );

                BOOST_CHECK_SMALL( std::fabs( expectedAerodynamicAcceleration( i ) - currentAerodynamicAcceleration( i ) ),
                                   std::max( 10.0 * currentAerodynamicAcceleration.norm( ), 1.0 ) *
                                   std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}

//! Test whether the thrust is properly computed if the full vector (*magnitude and direction) is imposed in either the
//! inertial or LVLH frame.
BOOST_AUTO_TEST_CASE( testInterpolatedThrustVector )
{

    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace gravitation;
    using namespace numerical_integrators;
    using namespace unit_conversions;
    using namespace interpolators;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 60.0;

    // Define body settings for simulation.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = boost::make_shared< BodySettings >( );
    bodySettings[ "Earth" ]->ephemerisSettings = getDefaultEphemerisSettings( "Earth" );
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );

    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    double bodyMass = 1.0;
    bodyMap[ "Asterix" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Asterix" ]->setConstantBodyMass( bodyMass );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

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

    std::map< double, Eigen::Vector3d > randomThrustMap;
    randomThrustMap[ 0 ] = Eigen::MatrixXd::Random( 3, 1 );
    randomThrustMap[ 1.0E4 ] = 20.0 * Eigen::MatrixXd::Random( 3, 1 );
    randomThrustMap[ 2.0E4 ] = 20.0 * Eigen::MatrixXd::Random( 3, 1 );
    randomThrustMap[ 3.0E4 ] = 20.0 * Eigen::MatrixXd::Random( 3, 1 );
    randomThrustMap[ 4.0E4 ] = 20.0 * Eigen::MatrixXd::Random( 3, 1 );
    randomThrustMap[ 5.0E4 ] = 20.0 * Eigen::MatrixXd::Random( 3, 1 );
    randomThrustMap[ 6.0E4 ] = 20.0 * Eigen::MatrixXd::Random( 3, 1 );
    randomThrustMap[ 7.0E4 ] = 20.0 * Eigen::MatrixXd::Random( 3, 1 );
    randomThrustMap[ 8.0E4 ] = 20.0 * Eigen::MatrixXd::Random( 3, 1 );
    randomThrustMap[ 9.0E4 ] = 20.0 * Eigen::MatrixXd::Random( 3, 1 );

    boost::shared_ptr< DataInterpolationSettings< double, Eigen::Vector3d > > thrustDataInterpolation =
            boost::make_shared< DataInterpolationSettings< double, Eigen::Vector3d > >(
                boost::make_shared< DataMapSettings< double, Eigen::Vector3d > >( randomThrustMap ),
                boost::make_shared< InterpolatorSettings >( linear_interpolator ) );

    for( unsigned int testCase = 0; testCase < 2; testCase++ )
    {

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define propagation settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
        accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );

        boost::shared_ptr< ThrustAccelerationSettings > thrustSettings =
                boost::make_shared< ThrustAccelerationSettings >(
                    thrustDataInterpolation, boost::lambda::constant( 300.0 ),
                    testCase == 0 ? inertial_thurst_frame : lvlh_thrust_frame, "Earth" );

        accelerationsOfAsterix[ "Asterix" ].push_back( thrustSettings );
        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );



        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
        // The initial conditions are given in Keplerian elements and later on converted to Cartesian
        // elements.


        // Convert Asterix state from Keplerian elements to Cartesian elements.
        double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                    asterixInitialStateInKeplerianElements,
                    earthGravitationalParameter );

        // Define list of dependent variables to save.
        std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back(
                    boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Asterix", "Asterix", 0 ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        relative_position_dependent_variable, "Asterix", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        relative_velocity_dependent_variable, "Asterix", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        lvlh_to_inertial_frame_rotation_dependent_variable, "Asterix", "Earth" ) );


        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch,
                  cowell, boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );
        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, fixedStepSize );


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            /////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );


        boost::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector3d > > thrustInterpolator =
                thrustSettings->interpolatorInterface_->getThrustInterpolator( );

        Eigen::Vector3d thrustDifference;

        if( testCase == 0 )
        {
            // Test whether thrust direction is indeed the imposed direction (in inertial frame).
            for( std::map< double, Eigen::VectorXd >::iterator outputIterator = dependentVariableResult.begin( );
                 outputIterator != dependentVariableResult.end( ); outputIterator++ )
            {
                thrustDifference =  outputIterator->second.segment( 0, 3 ) -
                        thrustInterpolator->interpolate( outputIterator->first );
                for( unsigned int i = 0; i < 3; i++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( thrustDifference( i ) ), 1.0E-14 );
                }
            }
        }
        else if( testCase == 1 )
        {
            // Test whether thrust direction is indeed the imposed direction (in lvlh frame).
            for( std::map< double, Eigen::VectorXd >::iterator outputIterator = dependentVariableResult.begin( );
                 outputIterator != dependentVariableResult.end( ); outputIterator++ )
            {
                Eigen::Matrix3d manualRotationMatrix =
                        reference_frames::getVelocityBasedLvlhToInertialRotation(
                            outputIterator->second.segment( 3, 6 ), Eigen::Vector6d::Zero( ) );
                Eigen::Matrix3d currentRotationMatrix =
                        getMatrixFromVectorRotationRepresentation( outputIterator->second.segment( 9, 9 ) );
                thrustDifference = manualRotationMatrix
                        * thrustInterpolator->interpolate( outputIterator->first ) - outputIterator->second.segment( 0, 3 );
                double currentThrustMagnitude = thrustInterpolator->interpolate( outputIterator->first ).norm( ) ;

                for( unsigned int i = 0; i < 3; i++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( thrustDifference( i ) ), 1.0E-12 * currentThrustMagnitude );
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( manualRotationMatrix( i, j ) -
                                                      currentRotationMatrix( i, j ) ), 1.0E-12
                                           );
                    }
                }
            }
        }
    }
}

//! Class that is used in the following unit test to check whether the parameterized thrust guidance is working correctly.
//! This class sets (using an algorithm with no physical basis). It computes a 'regular' guidance input, as well as a thrust
//! multiplier (throttle) if required. The regular guidance input is used as the Mach number when interpolating the thrust
//! tables.
class ThrustMultiplierComputation: public simulation_setup::ThrustInputParameterGuidance
{
public:
    ThrustMultiplierComputation( const double startTime, const double endTime,
                                 const bool useDummyMachNumber, const bool useThrustMultiplier ):
        ThrustInputParameterGuidance(
            static_cast< int >( useDummyMachNumber ) + static_cast< int >( useThrustMultiplier ), 0,
            useThrustMultiplier, ( useThrustMultiplier == true ) ? static_cast< int >( useDummyMachNumber ) : -1 ),
        startTime_( startTime ), endTime_( endTime ),
        useDummyMachNumber_( useDummyMachNumber ), useThrustMultiplier_( useThrustMultiplier ){ }

    void updateGuidanceParameters(  )
    {
        if( useDummyMachNumber_ )
        {
            currentThrustGuidanceParameters_[ 0 ] = ( currentTime_ - startTime_ ) / ( endTime_ - startTime_ ) * 15.0;
            if( useThrustMultiplier_ )
            {
                currentThrustGuidanceParameters_[ 1 ] = 1.0 - ( currentTime_ - startTime_ ) / ( endTime_ - startTime_ );
            }
        }
        else if( useThrustMultiplier_ )
        {
            currentThrustGuidanceParameters_[ 0 ] = 1.0 - ( currentTime_ - startTime_ ) / ( endTime_ - startTime_ );
        }
    }

private:

    double startTime_;

    double endTime_;

    bool useDummyMachNumber_;

    bool useThrustMultiplier_;
};

//! Test to check whether the parameterized thrust settings (i.e. thrust as a function of environment and/or guidance input)
//! is working correctly
BOOST_AUTO_TEST_CASE( testConcurrentThrustAndAerodynamicAccelerationWithEnvironmentDependentThrust )
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
    spice_interface::loadStandardSpiceKernels( );


    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 200.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 1.0;


    // Set spherical elements for Apollo.
    Vector6d apolloSphericalEntryState;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 50.0E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 6.0E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            1.0 * mathematical_constants::PI / 180.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;

    // Convert apollo state from spherical elements to Cartesian elements.
    Vector6d apolloInitialState = orbital_element_conversions::convertSphericalOrbitalToCartesianState(
                apolloSphericalEntryState );

    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 1.0E4,
                                    simulationEndEpoch + 1.0E4 );
    bodySettings[ "Earth" ]->gravityFieldSettings =
            boost::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects.
    bodyMap[ "Apollo" ] = boost::make_shared< simulation_setup::Body >( );
    double vehicleMass = 5.0E5;
    bodyMap[ "Apollo" ]->setConstantBodyMass( vehicleMass );

    // Create vehicle aerodynamic coefficients
    bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface( unit_tests::getApolloCoefficientInterface( ) );
    bodyMap[ "Apollo" ]->setEphemeris(
                boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    boost::shared_ptr< interpolators::OneDimensionalInterpolator<
                    double, Eigen::Vector6d  > >( ), "Earth" ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    int numberOfCasesPerSet = 4;
    for( int i = 0; i < numberOfCasesPerSet * 2; i++ )
    {
        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define acceleration model settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
        accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
        accelerationsOfApollo[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );

        // Define specific impulse dependencies.
        std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependencies;
        specificImpulseDependencies.push_back( propulsion::mach_number_dependent_thrust );
        specificImpulseDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );

        // Define variables (thrust dependencies and guidance object) that are different per case.
        std::vector< propulsion::ThrustIndependentVariables > thrustDependencies;
        boost::shared_ptr< ThrustMultiplierComputation > thrustInputParameterGuidance;

        // Use no guidance input
        if( ( i % numberOfCasesPerSet == 0 ) )
        {
            thrustDependencies.push_back( propulsion::mach_number_dependent_thrust );
            thrustDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );
        }
        // Use guidance for throttling the thrust
        if( ( i % numberOfCasesPerSet == 1 ) )
        {
            thrustDependencies.push_back( propulsion::mach_number_dependent_thrust );
            thrustDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );
            thrustDependencies.push_back( propulsion::throttle_dependent_thrust );

            thrustInputParameterGuidance =
                    boost::make_shared< ThrustMultiplierComputation >( simulationStartEpoch, simulationEndEpoch, 0, 1 );
        }
        // Use guidance input to generate a 'fake' Mach number
        else if( ( i % numberOfCasesPerSet == 2 ) )
        {
            thrustDependencies.push_back( propulsion::guidance_input_dependent_thrust );
            thrustDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );

            thrustInputParameterGuidance =
                    boost::make_shared< ThrustMultiplierComputation >( simulationStartEpoch, simulationEndEpoch, 1, 0 );
        }
        // Use guidance input to generate a 'fake' Mach number and a throttle
        else if( ( i % numberOfCasesPerSet == 3 ) )
        {
            thrustDependencies.push_back( propulsion::guidance_input_dependent_thrust );
            thrustDependencies.push_back( propulsion::throttle_dependent_thrust );
            thrustDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );

            thrustInputParameterGuidance =
                    boost::make_shared< ThrustMultiplierComputation >( simulationStartEpoch, simulationEndEpoch, 1, 1 );

        }

        // Load thrust tables from files (as a function of Mach number and specific impulse) and create interpolators.
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > thrustValues =
                MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Propulsion/UnitTests/Tmax_test.txt" );
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > specificImpulseValues =
                MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Propulsion/UnitTests/Isp_test.txt" );

        boost::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator =
                boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    thrustValues.second, thrustValues.first );
        boost::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator =
                boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    specificImpulseValues.second, specificImpulseValues.first );

        double constantSpecificImpulse = 1000.0;

        // Test with interpolated specific impulse
        if( i < numberOfCasesPerSet )
        {
            // Thrust requires no guidance.
            if( !( i % numberOfCasesPerSet == 0 ) )
            {
                accelerationsOfApollo[ "Apollo" ].push_back(
                            boost::make_shared< ThrustAccelerationSettings >(
                                boost::make_shared< ThrustDirectionGuidanceSettings >(
                                    thrust_direction_from_existing_body_orientation, "Earth" ),
                                createParameterizedThrustMagnitudeSettings(
                                    thrustInputParameterGuidance, thrustMagnitudeInterpolator, thrustDependencies,
                                    specificImpulseInterpolator, specificImpulseDependencies ) ) );
            }
            // Thrust requires guidance.
            else
            {
                accelerationsOfApollo[ "Apollo" ].push_back(
                            boost::make_shared< ThrustAccelerationSettings >(
                                boost::make_shared< ThrustDirectionGuidanceSettings >(
                                    thrust_direction_from_existing_body_orientation, "Earth" ),
                                boost::make_shared< ParameterizedThrustMagnitudeSettings >(
                                    thrustMagnitudeInterpolator, thrustDependencies,
                                    specificImpulseInterpolator, specificImpulseDependencies ) ) );
            }



        }
        // Test with constant specific impulse
        else
        {
            // Thrust requires no guidance.
            if( !( i % numberOfCasesPerSet == 0 ) )
            {
                accelerationsOfApollo[ "Apollo" ].push_back(
                            boost::make_shared< ThrustAccelerationSettings >(
                                boost::make_shared< ThrustDirectionGuidanceSettings >(
                                    thrust_direction_from_existing_body_orientation, "Earth" ),
                                createParameterizedThrustMagnitudeSettings(
                                    thrustInputParameterGuidance, thrustMagnitudeInterpolator, thrustDependencies,
                                    constantSpecificImpulse ) ) );
            }
            // Thrust requires guidance.
            else
            {
                accelerationsOfApollo[ "Apollo" ].push_back(
                            boost::make_shared< ThrustAccelerationSettings >(
                                boost::make_shared< ThrustDirectionGuidanceSettings >(
                                    thrust_direction_from_existing_body_orientation, "Earth" ),
                                boost::make_shared< ParameterizedThrustMagnitudeSettings >(
                                    thrustMagnitudeInterpolator, thrustDependencies,
                                    constantSpecificImpulse ) ) );
            }
        }

        // Finalize acceleration settings.
        accelerationMap[ "Apollo" ] = accelerationsOfApollo;
        bodiesToPropagate.push_back( "Apollo" );
        centralBodies.push_back( "Earth" );

        // Set initial state
        Eigen::Vector6d systemInitialState = apolloInitialState;


        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        // Set trimmed conditions for body orientation.
        setTrimmedConditions( bodyMap.at( "Apollo" ) );

        // Define list of dependent variables to save.
        std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        mach_number_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        airspeed_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        local_density_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        thrust_acceleration, "Apollo", "Apollo", 1 ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        total_mass_rate_dependent_variables, "Apollo" ) );

        // Define propagation settings.
        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                  boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ), cowell,
                  boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

        std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
        massRateModels[ "Apollo" ] = createMassRateModel( "Apollo", boost::make_shared< FromThrustMassModelSettings >( 1 ),
                                                          bodyMap, accelerationModelMap );

        boost::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
                boost::make_shared< MassPropagatorSettings< double > >(
                    boost::assign::list_of( "Apollo" ), massRateModels,
                    ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ),
                    boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ) );

        std::vector< boost::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( translationalPropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        boost::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                boost::make_shared< MultiTypePropagatorSettings< double > >(
                    propagatorSettingsVector, boost::make_shared< propagators::PropagationTimeTerminationSettings >(
                        simulationEndEpoch ),
                    boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

        // Define integration settings.
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

        // Declare test variables
        double currentDynamicPressure, currentMachNumber, currentMass;
        double currentThrustForce, currentMassRate, expectedThrust, expectedMassRate;
        double currentSpecificImpulse;
        std::vector< double > currentThrustInput;
        std::vector< double > specificImpulseInput;

        for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
             variableIterator != dependentVariableSolution.end( ); variableIterator++ )
        {
            // Retrieve data from dependent variables/propagated dynamics
            currentMass = numericalSolution.at( variableIterator->first )( 6 );

            currentDynamicPressure =
                    0.5 * variableIterator->second( 2 ) * variableIterator->second( 1 ) * variableIterator->second( 1 );
            currentMachNumber = variableIterator->second( 0 );
            currentThrustForce = variableIterator->second( 3 ) * currentMass;

            currentMassRate = -variableIterator->second( 4 );

            if( !( i % numberOfCasesPerSet == 0 ) )
            {
                thrustInputParameterGuidance->update( variableIterator->first );
            }
            currentThrustInput.clear( );

            if( ( i % numberOfCasesPerSet == 0 )  || ( i % numberOfCasesPerSet == 1 ) )
            {
                currentThrustInput.push_back( currentMachNumber );
            }
            else
            {
                currentThrustInput.push_back( thrustInputParameterGuidance->getThrustInputGuidanceParameter( 0 ) );
            }
            currentThrustInput.push_back( currentDynamicPressure );

            expectedThrust = thrustMagnitudeInterpolator->interpolate( currentThrustInput );
            if( ( i % numberOfCasesPerSet == 1 ) || ( i % numberOfCasesPerSet == 3 ) || ( i % numberOfCasesPerSet == 4 ) )
            {
                expectedThrust *= thrustInputParameterGuidance->getThrustInputGuidanceParameter(
                            thrustInputParameterGuidance->getThrottleSettingIndex( ) );
            }

            specificImpulseInput.clear( );
            specificImpulseInput.push_back( currentMachNumber );
            specificImpulseInput.push_back( currentDynamicPressure );

            if( !( i < numberOfCasesPerSet ) )
            {
                currentSpecificImpulse = constantSpecificImpulse;
            }
            else
            {
                currentSpecificImpulse = specificImpulseInterpolator->interpolate( specificImpulseInput );
            }

            expectedMassRate = expectedThrust /
                    ( currentSpecificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );

            BOOST_CHECK_CLOSE_FRACTION( expectedThrust, currentThrustForce,
                                        5.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( currentMassRate, expectedMassRate,
                                        5.0 * std::numeric_limits< double >::epsilon( ) );
        }
    }
}


//! Test to check whether the acceleration-limited thrust guidance (i.e throttle from a maximum impoised acceleration)
//! works correctly
BOOST_AUTO_TEST_CASE( testAccelerationLimitedGuidedThrust )
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
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 200.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 1.0;


    // Set spherical elements for Apollo.
    Vector6d apolloSphericalEntryState;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 50.0E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 6.0E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            1.0 * mathematical_constants::PI / 180.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;

    // Convert apollo state from spherical elements to Cartesian elements.
    Vector6d apolloInitialState = orbital_element_conversions::convertSphericalOrbitalToCartesianState(
                apolloSphericalEntryState );

    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 1.0E4,
                                    simulationEndEpoch + 1.0E4 );
    bodySettings[ "Earth" ]->gravityFieldSettings =
            boost::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects.
    bodyMap[ "Apollo" ] = boost::make_shared< simulation_setup::Body >( );
    double vehicleMass = 5.0E5;
    bodyMap[ "Apollo" ]->setConstantBodyMass( vehicleMass );

    // Create vehicle aerodynamic coefficients
    bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface( unit_tests::getApolloCoefficientInterface( ) );
    bodyMap[ "Apollo" ]->setEphemeris(
                boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    boost::shared_ptr< interpolators::OneDimensionalInterpolator<
                    double, Eigen::Vector6d  > >( ), "Earth" ) );

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

    std::vector< propulsion::ThrustIndependentVariables > thrustDependencies;
    thrustDependencies.push_back( propulsion::mach_number_dependent_thrust );
    thrustDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );
    thrustDependencies.push_back( propulsion::throttle_dependent_thrust );

    std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependencies;
    specificImpulseDependencies.push_back( propulsion::mach_number_dependent_thrust );
    specificImpulseDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );

    std::string thrustFile =
            tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Propulsion/UnitTests/Tmax_test.txt";
    std::string specificImpulseFile =
            tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Propulsion/UnitTests/Isp_test.txt";

    accelerationsOfApollo[ "Apollo" ].push_back(
                boost::make_shared< ThrustAccelerationSettings >(
                    boost::make_shared< ThrustDirectionGuidanceSettings >(
                        thrust_direction_from_existing_body_orientation, "Earth" ),
                    createAccelerationLimitedParameterizedThrustMagnitudeSettings(
                        bodyMap, "Apollo", physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION,
                        thrustFile, thrustDependencies,
                        specificImpulseFile, specificImpulseDependencies, "Earth" ) ) );

    accelerationMap[ "Apollo" ] = accelerationsOfApollo;

    bodiesToPropagate.push_back( "Apollo" );
    centralBodies.push_back( "Earth" );

    // Set initial state
    Eigen::Vector6d systemInitialState = apolloInitialState;


    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    setTrimmedConditions( bodyMap.at( "Apollo" ) );

    // Define list of dependent variables to save.
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    thrust_acceleration, "Apollo", "Apollo", 1 ) );


    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ), cowell,
              boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

    std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Apollo" ] = createMassRateModel( "Apollo", boost::make_shared< FromThrustMassModelSettings >( 1 ),
                                                      bodyMap, accelerationModelMap );

    boost::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
            boost::make_shared< MassPropagatorSettings< double > >(
                boost::assign::list_of( "Apollo" ), massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ),
                boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ) );

    std::vector< boost::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( translationalPropagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    boost::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
            boost::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector, boost::make_shared< propagators::PropagationTimeTerminationSettings >(
                    simulationEndEpoch ),
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

    for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
         variableIterator != dependentVariableSolution.end( ); variableIterator++ )
    {
        BOOST_CHECK_EQUAL(
                    ( variableIterator->second( 0 ) ) <=
                    physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * (
                        1.0 + 4.0 * std::numeric_limits< double >::epsilon( ) ), true );
    }
}

//! Test to check whether the mee-costate based thrust guidance is working correctly
BOOST_AUTO_TEST_CASE( testMeeCostateBasedThrust )
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
    using namespace unit_conversions;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 24.0 * 3600.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 15.0;


    // Set spherical elements for Apollo.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = 0.5;
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = 0.5;
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.5;
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Convert asterix state from spherical elements to Cartesian elements.
    Vector6d asterixInitialState = orbital_element_conversions::convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, getBodyGravitationalParameter( "Earth" ) );

    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 1.0E4,
                                    simulationEndEpoch + 1.0E4 );
    bodySettings[ "Earth" ]->gravityFieldSettings =
            boost::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects.
    bodyMap[ "Asterix" ] = boost::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "Earth", "ECLIPJ2000" );

    double vehicleMass = 5.0E5;

    // Run simulations for a single MEE costate not equal to zero, for each of the first 5 elements
    for( unsigned int i = 0; i < 5; i++ )
    {
        bodyMap[ "Asterix" ]->setConstantBodyMass( vehicleMass );

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define acceleration model settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
        accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );

        Eigen::VectorXd costates = Eigen::VectorXd::Zero( 5 );
        costates( i ) = 100.0;

        boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
                boost::make_shared< MeeCostateBasedThrustDirectionSettings >(
                    "Asterix", "Earth", boost::lambda::constant( costates ) );
        boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings =
                boost::make_shared< ConstantThrustEngineSettings >( 1.0E4, 30000.0 );

        accelerationsOfAsterix[ "Asterix" ].push_back(
                    boost::make_shared< ThrustAccelerationSettings >(
                        thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );

        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );

        // Set initial state
        Eigen::Vector6d systemInitialState = asterixInitialState;


        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        // Define list of dependent variables to save.
        std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        modified_equinocial_state_dependent_variable, "Asterix", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        thrust_acceleration, "Asterix", "Asterix" ) );

        // Create propagator/integrator settings
        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                  boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ), cowell,
                  boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );
        std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
        massRateModels[ "Asterix" ] = createMassRateModel(
                    "Asterix", boost::make_shared< FromThrustMassModelSettings >( 1 ),
                    bodyMap, accelerationModelMap );
        boost::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
                boost::make_shared< MassPropagatorSettings< double > >(
                    boost::assign::list_of( "Asterix" ), massRateModels,
                    ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ),
                    boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ) );
        std::vector< boost::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( translationalPropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );
        boost::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                boost::make_shared< MultiTypePropagatorSettings< double > >(
                    propagatorSettingsVector, boost::make_shared< propagators::PropagationTimeTerminationSettings >(
                        simulationEndEpoch ),
                    boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4, simulationStartEpoch, fixedStepSize );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );

        // Retrieve change in Modified Equinoctial Elements
        std::map< double, Eigen::VectorXd > dependentVariableSolution =
                dynamicsSimulator.getDependentVariableHistory( );
        Eigen::Vector6d finalModifiedEquinoctialElementsError =
                dependentVariableSolution.rbegin( )->second.segment( 0, 6 ) -
                dependentVariableSolution.begin( )->second.segment( 0, 6 );

        // Test whether MEE rates are within reasonable bounds (values are determined empirically)
        BOOST_CHECK_EQUAL( ( finalModifiedEquinoctialElementsError( i ) < 0 ), true );
        if( i == 0 )
        {
            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 0 ) ) > 2.5E6 ), true );
        }
        else if( i < 3 )
        {
            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 0 ) ) < 2E5 ), true );
        }
        else
        {
            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 0 ) ) < 0.1 ), true );
        }

        if( i == 1 )
        {
            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 1 ) ) > 0.1 ), true );
        }
        else
        {
            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 1 ) ) < 0.025 ), true );
        }

        if( i == 2 )
        {
            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 2 ) ) > 0.1 ), true );
        }
        else
        {
            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 2 ) ) < 0.025 ), true );
        }

        if( i == 3 )
        {
            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 3 ) ) > 0.075 ), true );
        }
        else
        {
            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 3 ) ) < 0.005 ), true );
        }

        if( i == 4 )
        {
            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 4 ) ) > 0.075 ), true );
        }
        else
        {
            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 4 ) ) < 0.005 ), true );
        }

//        std::cout<<finalModifiedEquinoctialElementsError.transpose( )<<std::endl;
//        input_output::writeDataMapToTextFile( dependentVariableSolution, "meeThrustDep_" + std::to_string( i ) + ".dat" );
//        input_output::writeDataMapToTextFile( numericalSolution, "meeThrustState_" + std::to_string( i ) + ".dat" );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
