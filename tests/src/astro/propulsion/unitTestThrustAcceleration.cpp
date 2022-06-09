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
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/basics/testMacros.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/multiDimensionalArrayReader.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/propagation_setup/createMassRateModels.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
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
    using namespace ephemerides;

    // Create Earth object
    simulation_setup::SystemOfBodies bodies;
    // Create vehicle objects.
    double vehicleMass = 5.0E3;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );
    Eigen::Vector3d thrustDirection;
    thrustDirection << -1.4, 2.4, 5.6;

    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
            [=](const double){ return thrustDirection; };
    bodies.at( "Vehicle" )->setRotationalEphemeris(
                createRotationModel(
                    std::make_shared< BodyFixedDirectionBasedRotationSettings >(
                        thrustDirectionFunction, "ECLIPJ2000", "VehicleFixed" ),
                    "Vehicle", bodies ) );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    double thrustMagnitude = 1.0E3;
    double specificImpulse = 250.0;
    double massRate = thrustMagnitude / ( specificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );

    addEngineModel(
                "Vehicle", "MainEngine",
                std::make_shared< ConstantThrustMagnitudeSettings >(
                    thrustMagnitude, specificImpulse ), bodies );



    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                       "MainEngine" ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "SSB" );

    // Set initial state
    Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            std::make_shared< propagators::PropagationTimeTerminationSettings >( 10.0 );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings );
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, 0.1 );
    {
        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, translationalPropagatorSettings, true, false, false );

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
        std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
        massRateModels[ "Vehicle" ] = (
                    createMassRateModel( "Vehicle", std::make_shared< FromThrustMassRateSettings >(1 ),
                                         bodies, accelerationModelMap ) );

        std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings =
                std::make_shared< MassPropagatorSettings< double > >(
                    std::vector< std::string >{ "Vehicle" }, massRateModels,
                    ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( translationalPropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings, true, false, false );

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

Eigen::Vector3d getEarthJupiterVector( const double time )
{
    return spice_interface::getBodyCartesianPositionAtEpoch(
                "Jupiter", "Earth", "J2000", "None", time );
}

Eigen::Vector3d getEarthJupiterVectorFromEnvironment(
        const double time, const simulation_setup::SystemOfBodies& bodies )
{
    return bodies.at( "Jupiter" )->getPosition( ) - bodies.at( "Earth" )->getPosition( );
}

double getFreeRotationAngle( const double time )
{
    return 2.0 * std::sin( 2.0 * mathematical_constants::PI * time / 3600.0 );
}


BOOST_AUTO_TEST_CASE( testDirectionBasedRotationWithThrustAcceleration )
{
    using namespace tudat;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace ephemerides;
    using namespace orbital_element_conversions;

    spice_interface::loadStandardSpiceKernels( );

    // Create Earth object
    simulation_setup::SystemOfBodies bodies;

    // Set simulation time settings.
    const double simulationEndEpoch = 4.0 * 3600.0;

    double thrustMagnitude1 = 1.0E3;
    double bodyMass = 400.0;
    double specificImpulse1 = 250.0;

    std::vector< std::map< double, Eigen::Matrix3d > > rotationMatrixHistoriesWithoutFreeRotationAngle;
    for( int test = 0; test < 8; test++ )
    {
        // Define body settings for simulation.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Sun" );
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Jupiter" );

        // Create body objects.
        BodyListSettings bodySettings =
                getDefaultBodySettings( bodiesToCreate, "Earth", "J2000" );
        SystemOfBodies bodies;

        std::function< double( const double ) > freeRotationAngleFunction = nullptr;
        if( test > 3 )
        {
            freeRotationAngleFunction = &getFreeRotationAngle;
        }

        if( ( test % 4 ) == 3 )
        {
            bodySettings.addSettings( "Vehicle" );
            bodySettings.at( "Vehicle" )->rotationModelSettings =
                    std::make_shared< BodyFixedDirectionBasedRotationSettings >(
                        std::function< Eigen::Vector3d( const double ) >( ), "J2000", "VehicleFixed" );
            bodies = createSystemOfBodies( bodySettings );
            std::dynamic_pointer_cast< tudat::ephemerides::DirectionBasedRotationalEphemeris >(
                        bodies.at( "Vehicle")->getRotationalEphemeris( ) )->setInertialBodyAxisDirectionCalculator(
                        std::make_shared< CustomBodyFixedDirectionCalculator >(
                        std::bind( &getEarthJupiterVectorFromEnvironment, std::placeholders::_1, bodies ) ) );
            std::dynamic_pointer_cast< tudat::ephemerides::DirectionBasedRotationalEphemeris >(
                        bodies.at( "Vehicle")->getRotationalEphemeris( ) )->setFreeRotationAngleFunction(
                        freeRotationAngleFunction );
        }
        else if( ( test % 4 ) == 2 )
        {
            bodies = createSystemOfBodies( bodySettings );

            // Create spacecraft object.
            bodies.createEmptyBody( "Vehicle" );
            bodies.at( "Vehicle" )->setRotationalEphemeris(
                        createRotationModel(
                            std::make_shared< BodyFixedDirectionBasedRotationSettings >(
                                &getEarthJupiterVector, "J2000", "VehicleFixed", freeRotationAngleFunction ),
                            "Vehicle", bodies ) );
        }
        else
        {
            bodySettings.addSettings( "Vehicle" );
            bodySettings.at( "Vehicle" )->rotationModelSettings =
                    std::make_shared< BodyFixedDirectionBasedRotationSettings >(
                        &getEarthJupiterVector, "J2000", "VehicleFixed", freeRotationAngleFunction );
            bodies = createSystemOfBodies( bodySettings );

        }
        bodies.at( "Vehicle" )->setConstantBodyMass( bodyMass );


        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define propagation settings.
        Eigen::Vector3d bodyFixedThrustDirection;
        if( ( test % 4 ) == 1 )
        {
            bodyFixedThrustDirection = Eigen::Vector3d::UnitY( );
        }
        else
        {
            bodyFixedThrustDirection = Eigen::Vector3d::UnitX( );
        }

        addEngineModel(
                    "Vehicle", "MainEngine",
                    std::make_shared< ConstantThrustMagnitudeSettings >(
                        thrustMagnitude1, specificImpulse1 ), bodies, bodyFixedThrustDirection );


        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfVehicle[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >( "MainEngine" ) );

        accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
        bodiesToPropagate.push_back( "Vehicle" );
        centralBodies.push_back( "Earth" );
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );


        // Set Keplerian elements for Vehicle.
        Eigen::Vector6d vehicleInitialStateInKeplerianElements;
        vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 8000.0E3;
        vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
        vehicleInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
        vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = unit_conversions::convertDegreesToRadians( 235.7 );
        vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = unit_conversions::convertDegreesToRadians( 23.4 );
        vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

        double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        const Eigen::Vector6d vehicleInitialState = convertKeplerianToCartesianElements(
                    vehicleInitialStateInKeplerianElements, earthGravitationalParameter );


        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;

        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        inertial_to_body_fixed_rotation_matrix_variable, "Vehicle" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
//        dependentVariables.push_back(
//                    std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
//                        "Vehicle", reference_frames::inertial_frame, reference_frames::body_frame ) );

        // Define propagator settings (Cowell)
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, vehicleInitialState, simulationEndEpoch,
                  cowell, std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

        // Define integrator settings.
        const double fixedStepSize = 5.0;
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, fixedStepSize );

        // Propagate orbit with Cowell method
        SingleArcDynamicsSimulator< double > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings );

        std::map< double, Eigen::Matrix3d > currentRotationMatrixHistory;
        for( auto it : dynamicsSimulator.getDependentVariableHistory( ) )
        {
            double currentTime = it.first;
            Eigen::Matrix3d currentRotationMatrixToBodyFixedFrame =
                    propagators::getMatrixFromVectorRotationRepresentation( it.second.segment( 0, 9 ) );
            if( test < 4 )
            {
                currentRotationMatrixHistory[ currentTime ] = currentRotationMatrixToBodyFixedFrame;
            }
            else
            {
                Eigen::Matrix3d remainingRotationMatrix =
                        rotationMatrixHistoriesWithoutFreeRotationAngle[ test - 4 ][ currentTime ] *
                        currentRotationMatrixToBodyFixedFrame.inverse( );
                double freeRotationAngle = getFreeRotationAngle( currentTime );
                Eigen::Matrix3d testRemainingRotationMatrix = Eigen::AngleAxisd(
                            freeRotationAngle, Eigen::Vector3d::UnitX( ) ).toRotationMatrix( );
                for( int i = 0; i < 3; i++ )
                {
                    for( int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( remainingRotationMatrix( i, j ) - testRemainingRotationMatrix( i, j ) ),
                                           8.0 * std::numeric_limits< double >::epsilon( ) );
                    }
                }
//                std::cout<<test<<" "<<freeRotationAngle<<std::endl<<remainingRotationMatrix - testRemainingRotationMatrix<<std::endl<<std::endl;
//                std::cout<<testRemainingRotationMatrix<<std::endl<<std::endl<<std::endl;

            }

            Eigen::Vector3d inertialThrustVector = it.second.segment( 9, 3 );
            if( test % 4 != 1 )
            {
                Eigen::Vector3d inertialThrustDirection = inertialThrustVector.normalized( );
                Eigen::Vector3d testInertialThrustDirection = getEarthJupiterVector( currentTime ).normalized( );

                for( int i = 0; i < 3; i++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( inertialThrustDirection( i ) - testInertialThrustDirection( i ) ),
                                       5.0 * std::numeric_limits< double >::epsilon( ) );
                }
            }

            Eigen::Vector3d bodyFixedThrustVector = currentRotationMatrixToBodyFixedFrame * inertialThrustVector;

            int zeroIndex = 1;
            int nonZeroIndex = 0;
            if( ( test % 4 ) == 1 )
            {
                zeroIndex = 0;
                nonZeroIndex = 1;
            }
            if( it.first < 50.0 )
            {
                std::cout<<"Testt "<<test<<" Time "<<it.first<<std::endl;
            }
            BOOST_CHECK_CLOSE_FRACTION( bodyFixedThrustVector( nonZeroIndex ), thrustMagnitude1 / bodyMass,
                                        5.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( bodyFixedThrustVector( zeroIndex ), thrustMagnitude1 / bodyMass *
                               5.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( bodyFixedThrustVector( 2 ), thrustMagnitude1 / bodyMass *
                               10.0 * std::numeric_limits< double >::epsilon( ) );
        }
        if( test < 4 )
        {
            rotationMatrixHistoriesWithoutFreeRotationAngle.push_back( currentRotationMatrixHistory );
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
        simulation_setup::SystemOfBodies bodies;

        // Create vehicle objects.
        double vehicleMass = 5.0E3;
        double dryVehicleMass = 2.0E3;

        bodies.createEmptyBody( "Vehicle" );
        bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );

        double thrustMagnitude1 = 1.0E3;
        double specificImpulse1 = 250.0;
        double massFlow1 = propulsion::computePropellantMassRateFromSpecificImpulse(
                    thrustMagnitude1, specificImpulse1 );

        double thrustMagnitude2 = 2.0E3;
        double specificImpulse2 = 300.0;
        double massFlow2 = propulsion::computePropellantMassRateFromSpecificImpulse(
                    thrustMagnitude2, specificImpulse2 );

        addEngineModel(
                    "Vehicle", "Engine1",
                    std::make_shared< ConstantThrustMagnitudeSettings >(
                        thrustMagnitude1, specificImpulse1 ), bodies );

        addEngineModel(
                    "Vehicle", "Engine2",
                    std::make_shared< ConstantThrustMagnitudeSettings >(
                        thrustMagnitude2, specificImpulse2 ), bodies );


        Eigen::Vector3d thrustDirection;
        thrustDirection << -1.4, 2.4, 5.6;

        std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
                [=](const double){ return thrustDirection; };
        bodies.at( "Vehicle" )->setRotationalEphemeris(
                    createRotationModel(
                        std::make_shared< BodyFixedDirectionBasedRotationSettings >(
                            thrustDirectionFunction, "ECLIPJ2000", "VehicleFixed" ),
                        "Vehicle", bodies ) );

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        // Define acceleration model settings.
        switch( i )
        {
        case 0:
        {
            accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
            std::vector< std::string >( { "Engine1", "Engine2" } ) ) );
            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
            break;
        }
        case 1:
        {
            accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                                   "Engine1" ) );
            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
            break;
        }
        case 2:
        {
            accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                               "Engine2" ) );
            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
            break;
        }
        case 3:
        {
            accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                                   "Engine1" ) );
            accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                                    "Engine2" ) );
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
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );

        std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                std::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings );
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, 0.1 );

        std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;

        double totalMassRate, totalThrust;
        switch( i )
        {
        case 0:
        {
            massRateModels[ "Vehicle" ] = (
                        createMassRateModel( "Vehicle", std::make_shared< FromThrustMassRateSettings >(1 ),
                                             bodies, accelerationModelMap ) );
            totalMassRate = massFlow1 + massFlow2;
            totalThrust = thrustMagnitude1 + thrustMagnitude2;
            break;
        }
        case 1:
        {
            massRateModels[ "Vehicle" ] = (
                        createMassRateModel( "Vehicle", std::make_shared< FromThrustMassRateSettings >(0, "Engine1" ),
                                             bodies, accelerationModelMap ) );
            totalMassRate = massFlow1;
            totalThrust = thrustMagnitude1;
            break;
        }
        case 2:
        {
            massRateModels[ "Vehicle" ] = (
                        createMassRateModel( "Vehicle", std::make_shared< FromThrustMassRateSettings >(0, "Engine2" ),
                                             bodies, accelerationModelMap ) );
            totalMassRate = massFlow2;
            totalThrust = thrustMagnitude2;
            break;
        }
        case 3:
        {
            massRateModels[ "Vehicle" ] = (
                        createMassRateModel( "Vehicle", std::make_shared< FromThrustMassRateSettings >(0, "Engine1" ),
                                             bodies, accelerationModelMap ) );
            totalMassRate = massFlow1;
            totalThrust = thrustMagnitude1 + thrustMagnitude2;
            break;
        }
        }

        double totalSpecificImpulse = totalThrust / ( physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * totalMassRate );

        std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings =
                std::make_shared< MassPropagatorSettings< double > >(
                    std::vector< std::string >{ "Vehicle" }, massRateModels,
                    ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( translationalPropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings, true, false, false );

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

    double thrustMagnitude1 = 1.0E3;
    double specificImpulse1 = 250.0;

    for( unsigned int i = 0; i < 2; i++ )
    {
        // Create Earth object
        simulation_setup::SystemOfBodies bodies;

        // Create vehicle objects.
        double vehicleMass = 5.0E3;
        bodies.createEmptyBody( "Earth" );

        bodies.at( "Earth" )->setEphemeris(
                    std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "ECLIPJ2000" ) );

//                    std::make_shared< ephemerides::SpiceEphemeris >( "Sun", "SSB", false, false ) );
        bodies.at( "Earth" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >(
                                                      spice_interface::getBodyGravitationalParameter( "Earth" ) ) );

        bodies.createEmptyBody( "Vehicle" );
        bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );

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

        bodies.at( "Vehicle" )->setRotationalEphemeris(
                    createRotationModel(
                        std::make_shared< OrbitalStateBasedRotationSettings >(
                            "Earth", isThurstInVelocityDirection, true, "ECLIPJ2000", "BodyFixed" ),
                        "Vehicle", bodies ) );

        addEngineModel(
                    "Vehicle", "MainEngine",
                    std::make_shared< ConstantThrustMagnitudeSettings >(
                        thrustMagnitude1, specificImpulse1 ), bodies );

        // Define acceleration model settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                           "MainEngine" ) );
        if( i == 1 )
        {
            accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        }

        accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

        bodiesToPropagate.push_back( "Vehicle" );
        centralBodies.push_back( "Earth" );

        // Set initial state
        double radius = 1.0E3;
        double circularVelocity = std::sqrt( radius * thrustMagnitude1 / vehicleMass );
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
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );

        std::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings;
        if( i == 1 )
        {
            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
            dependentVariables.push_back(
                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
            dependentVariableSaveSettings = std::make_shared< DependentVariableSaveSettings >( dependentVariables );

        }
        std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                std::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings,
                  cowell, dependentVariableSaveSettings );
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, 0.1 );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, translationalPropagatorSettings, true, false, false );

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
                            std::fabs( outputIterator->second( 0 ) - radius * std::cos( currentAngle ) ), 1.0E-9 * radius );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 1 ) - radius * std::sin( currentAngle ) ), 1.0E-9 * radius );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 2 ) ), 1.0E-15 );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 3 ) + circularVelocity * std::sin( currentAngle ) ),
                            1.0E-9 * circularVelocity  );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 4 ) - circularVelocity * std::cos( currentAngle ) ),
                            1.0E-9 * circularVelocity  );
                BOOST_CHECK_SMALL(
                            std::fabs( outputIterator->second( 5 ) ), 1.0E-15 );
            }
        }
        else if( i == 1 )
        {
            for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
                 numericalSolution.begin( ); outputIterator != numericalSolution.end( ); outputIterator++ )
            {
                Eigen::Vector3d vectorDifference =
                        ( -1.0 * thrustMagnitude1 / vehicleMass * outputIterator->second.segment( 3, 3 ).normalized( ) ) -
                        ( dependentVariableSolution.at( outputIterator->first ).normalized( ) );

//                for( int j = 0; j < 3; j++ )
//                {
//                    BOOST_CHECK_SMALL(
//                                std::fabs( vectorDifference( j ) ) / dependentVariableSolution.at( outputIterator->first ).norm( ), 1.0E-14 );
//                }

                // Check if the thrust acceleration is of the correct magnitude, and in the same direction as the velocity.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( -1.0 * thrustMagnitude1 / vehicleMass * outputIterator->second.segment( 3, 3 ).normalized( ) ),
                            ( dependentVariableSolution.at( outputIterator->first ) ), 1.0E-14 );
//                }

            }
        }
    }

}

//! Test the thrust force direction when it is taken from a predetermined vehicle rotation (RotationalEphemeris)
BOOST_AUTO_TEST_CASE( testThrustAccelerationFromExistingRotation )
{
    std::cout<<"testThrustAccelerationFromExistingRotation"<<std::endl;
    using namespace tudat;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );



    double thrustMagnitude1 = 1.0E3;
    double specificImpulse1 = 250.0;

    // Create Earth object
    simulation_setup::SystemOfBodies bodies;

    // Create vehicle objects.
    bodies.createEmptyBody( "Earth" );
    bodies.at( "Earth" )->setEphemeris(
                std::make_shared< ephemerides::SpiceEphemeris >( "Sun", "SSB", false, false ) );
    bodies.at( "Earth" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >(
                                                    spice_interface::getBodyGravitationalParameter( "Earth" ) ) );

    double vehicleMass = 5.0E3;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );
    bodies.at( "Vehicle" )->setRotationalEphemeris(
                std::make_shared< ephemerides::SpiceRotationalEphemeris >( "ECLIPJ2000", "IAU_MOON" ) );
    

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    Eigen::Vector3d bodyFixedThrustDirection = ( Eigen::Vector3d( ) << 1.4, 3.1, -0.5 ).finished( ).normalized( );

    addEngineModel(
                "Vehicle", "MainEngine",
                std::make_shared< ConstantThrustMagnitudeSettings >(
                    thrustMagnitude1, specificImpulse1 ), bodies, bodyFixedThrustDirection );

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >( "MainEngine" ) );
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );


    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Set initial state
    Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );


    systemInitialState( 0 ) = 8.0E6;
    systemInitialState( 4 ) = 7.5E3;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    std::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings;


    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
    dependentVariableSaveSettings = std::make_shared< DependentVariableSaveSettings >( dependentVariables );


    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            std::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings,
              cowell, dependentVariableSaveSettings );
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, 2.5 );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, translationalPropagatorSettings, true, false, false );

    // Retrieve numerical solutions for state and dependent variables
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput =
            dynamicsSimulator.getDependentVariableHistory( );

    double thrustAcceleration = thrustMagnitude1 / vehicleMass;
    Eigen::Quaterniond rotationToInertialFrame;
    for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
         dependentVariableOutput.begin( ); outputIterator != dependentVariableOutput.end( ); outputIterator++ )
    {
        rotationToInertialFrame = bodies.at( "Vehicle" )->getRotationalEphemeris( )->getRotationToBaseFrame(
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
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings.at( "Earth" )->gravityFieldSettings =
            std::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

    for( unsigned int testCase = 0; testCase < 2; testCase++ )
    {
        // Create Earth object
        simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );

        // Create vehicle objects.
        bodies.createEmptyBody( "Apollo" );
        double vehicleMass = 5.0E3;
        bodies.at( "Apollo" )->setConstantBodyMass( vehicleMass );

        // Create vehicle aerodynamic coefficients
        bodies.at( "Apollo" )->setAerodynamicCoefficientInterface(
                    unit_tests::getApolloCoefficientInterface( ) );
        bodies.at( "Apollo" )->setRotationalEphemeris(
                    createRotationModel(
                        std::make_shared< PitchTrimRotationSettings >( "Earth", "ECLIPJ2000", "VehicleFixed" ),
                        "Apollo", bodies ) );



        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define acceleration model settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
        accelerationsOfApollo[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

        double thrustMagnitude1 = 1.0E-3;
        double specificImpulse1 = 250.0;

        addEngineModel(
                    "Apollo", "MainEngine",
                    std::make_shared< ConstantThrustMagnitudeSettings >(
                        thrustMagnitude1, specificImpulse1 ), bodies );


        accelerationsOfApollo[ "Apollo" ].push_back( std::make_shared< ThrustAccelerationSettings >( "MainEngine" ) );

        accelerationMap[ "Apollo" ] = accelerationsOfApollo;

        bodiesToPropagate.push_back( "Apollo" );
        centralBodies.push_back( "Earth" );

        // Set initial state
        Eigen::Vector6d systemInitialState = apolloInitialState;


        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );


        // Define list of dependent variables to save.
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "Apollo", reference_frames::angle_of_attack ) );
        dependentVariables.push_back(
                    std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "Apollo", reference_frames::angle_of_sideslip ) );
        dependentVariables.push_back(
                    std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "Apollo", reference_frames::bank_angle ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        airspeed_dependent_variable, "Apollo", "Earth" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        local_density_dependent_variable, "Apollo", "Earth" ) );
        dependentVariables.push_back(
                    std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                        "Apollo", reference_frames::inertial_frame, reference_frames::body_frame ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                            inertial_to_body_fixed_rotation_matrix_variable, "Apollo" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        aerodynamic, "Apollo", "Earth", 0 ) );
        dependentVariables.push_back(
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        thrust_acceleration, "Apollo", "Apollo", 0 ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        aerodynamic_force_coefficients_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        aerodynamic_moment_coefficients_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                        "Apollo", reference_frames::inertial_frame, reference_frames::aerodynamic_frame ) );
        dependentVariables.push_back(
                    std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                        "Apollo", reference_frames::aerodynamic_frame, reference_frames::body_frame ) );

        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                  std::make_shared< propagators::PropagationTimeTerminationSettings >( 3200.0 ), cowell,
                  std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, simulationStartEpoch, fixedStepSize );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings, true, false, false );

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

        std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > vehicelCoefficientInterface =
                bodies.at( "Apollo" )->getAerodynamicCoefficientInterface( );

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

            vehicelCoefficientInterface->updateCurrentCoefficients( { currentMachNumber, currentAngleOfAttack, currentAngleOfSideSlip } );
            aerodynamicCoefficients = vehicelCoefficientInterface->getCurrentForceCoefficients( );



            expectedThrustDirection = rotationToBodyFixedFrame1.transpose( ) * bodyFixedThrustDirection;
            computedThrustDirection = currentThrustAcceleration.normalized( );

            // Check thrust magnitude
            BOOST_CHECK_CLOSE_FRACTION(
                        ( currentThrustAcceleration ).norm( ), thrustMagnitude1 / vehicleMass,
                        20.0 * std::numeric_limits< double >::epsilon( ) );




            Eigen::Vector3d expectedAerodynamicAcceleration = 0.5 * currentDensity * currentAirspeed * currentAirspeed *
                    bodies.at( "Apollo" )->getAerodynamicCoefficientInterface( )->getReferenceArea( ) *
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
                                20.0 * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL(
                                std::fabs( currentRotationFromAerodynamicToBodyFixedFrame( i, j ) -
                                           computedAerodynamicToBodyFixedFrame( i, j ) ),
                                10.0 * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL(
                                std::fabs( currentRotationFromAerodynamicToBodyFixedFrame2( i, j ) -
                                           computedAerodynamicToBodyFixedFrame( i, j ) ),
                                20.0 * std::numeric_limits< double >::epsilon( ) );
                }
                // Check thrust direction
                BOOST_CHECK_SMALL( std::fabs( expectedThrustDirection( i ) - computedThrustDirection( i ) ),
                                   40.0 * std::numeric_limits< double >::epsilon( ) );

                BOOST_CHECK_SMALL( std::fabs( expectedAerodynamicAcceleration( i ) - currentAerodynamicAcceleration( i ) ),
                                   std::max( 10.0 * currentAerodynamicAcceleration.norm( ), 1.0 ) *
                                   std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}

void createRotationAndEngineModelForFullThrust(
        const std::string& bodyName,
        const std::function< Eigen::Vector3d( const double ) > thrustFunction,
        const std::string& engineName,
        const std::string& bodyFixedDirectionName,
        const simulation_setup::SystemOfBodies& bodies,
        const ephemerides::SatelliteBasedFrames directionFrame = ephemerides::SatelliteBasedFrames::inertial_satellite_based_frame )
{
    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
            [=](const double time) -> Eigen::Vector3d
    {
        if( time == time )
        {
            return thrustFunction( time ).normalized( );
        }
        else
        {
            return Eigen::Vector3d::Constant( TUDAT_NAN );
        }
    };
    std::function< double( const double ) > thrustMagnitudeFunction =
            [=](const double time)
    {
        if( time == time )
        {
            return thrustFunction( time ).norm( );
        }
        else
        {
            return TUDAT_NAN;
        }
    };

    bodies.at( bodyName )->setRotationalEphemeris(
                createRotationModel(
                    std::make_shared< simulation_setup::BodyFixedDirectionBasedRotationSettings >(
                        thrustDirectionFunction, bodies.getFrameOrientation( ), bodyFixedDirectionName, nullptr,
                        std::make_pair( directionFrame, "Earth" ) ),
                    bodyName, bodies ) );
    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustSettings =
            std::make_shared< simulation_setup::CustomThrustMagnitudeSettings >( thrustMagnitudeFunction , [=](const double){return TUDAT_NAN; } );

    addEngineModel( bodyName,  engineName, thrustSettings,
                    bodies );
}

//! Test whether the thrust is properly computed if the full vector (*magnitude and direction) is imposed in either the
//! inertial or TNW frame.
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
    using namespace ephemerides;


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
    BodyListSettings bodySettings;
    bodySettings.addSettings( "Earth" );
    bodySettings.at( "Earth" )->ephemerisSettings = getDefaultEphemerisSettings( "Earth" );
    bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );

    for( unsigned int testCase = 0; testCase < 2; testCase++ )
    {
        // Create Earth object
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Create spacecraft object.
        double bodyMass = 1.0;
        bodies.createEmptyBody( "Asterix" );
        bodies.at( "Asterix" )->setConstantBodyMass( bodyMass );





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

        std::shared_ptr< DataInterpolationSettings< double, Eigen::Vector3d > > thrustDataInterpolation =
                std::make_shared< DataInterpolationSettings< double, Eigen::Vector3d > >(
                    std::make_shared< DataMapSettings< double, Eigen::Vector3d > >( randomThrustMap ),
                    std::make_shared< InterpolatorSettings >( linear_interpolator ) );

        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector3d > > thrustInterpolator =
                interpolators::createOneDimensionalInterpolator( thrustDataInterpolation );
        typedef interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > LocalInterpolator;
        std::function< Eigen::Vector3d( const double ) > thrustFunction =
                std::bind(
                    static_cast< Eigen::Vector3d( LocalInterpolator::* )( const double ) >
                    ( &LocalInterpolator::interpolate ), thrustInterpolator, std::placeholders::_1 );
//        std::function< Eigen::Vector3d( const double ) > thrustFunction =
//                [=](const double time){
//            std::cout<<time<<" "<<thrustFunction2( time ).transpose( )<<std::endl<<std::endl;
//            return thrustFunction2( time );

//        };

        createRotationAndEngineModelForFullThrust(
                    "Asterix", thrustFunction, "MainEngine", "Vehicle_Fixed", bodies,
                    testCase == 0 ? inertial_satellite_based_frame : tnw_satellite_based_frame );

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define propagation settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
        accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::point_mass_gravity ) );


        std::shared_ptr< ThrustAccelerationSettings > thrustSettings =
                std::make_shared< ThrustAccelerationSettings >( "MainEngine" );

        accelerationsOfAsterix[ "Asterix" ].push_back( thrustSettings );
        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );



        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
        // The initial conditions are given in Keplerian elements and later on converted to Cartesian
        // elements.


        // Convert Asterix state from Keplerian elements to Cartesian elements.
        double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                    asterixInitialStateInKeplerianElements,
                    earthGravitationalParameter );

        // Define list of dependent variables to save.
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back(
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Asterix", "Asterix", 0 ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        relative_position_dependent_variable, "Asterix", "Earth" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        relative_velocity_dependent_variable, "Asterix", "Earth" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        tnw_to_inertial_frame_rotation_dependent_variable, "Asterix", "Earth" ) );


        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch,
                  cowell, std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, fixedStepSize );


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            /////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings, true, false, false );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );


        Eigen::Vector3d thrustDifference;

        if( testCase == 0 )
        {
            // Test whether thrust direction is indeed the imposed direction (in inertial frame).
            for( std::map< double, Eigen::VectorXd >::iterator outputIterator = dependentVariableResult.begin( );
                 outputIterator != dependentVariableResult.end( ); outputIterator++ )
            {
                thrustDifference =  outputIterator->second.segment( 0, 3 ) -
                        thrustFunction( outputIterator->first );

                for( unsigned int i = 0; i < 3; i++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( thrustDifference( i ) ), 1.0E-14 );
                }
            }
        }
        else if( testCase == 1 )
        {
            // Test whether thrust direction is indeed the imposed direction (in tnw frame).
            for( std::map< double, Eigen::VectorXd >::iterator outputIterator = dependentVariableResult.begin( );
                 outputIterator != dependentVariableResult.end( ); outputIterator++ )
            {
                Eigen::Matrix3d manualRotationMatrix =
                        reference_frames::getTnwToInertialRotation(
                            outputIterator->second.segment( 3, 6 ), Eigen::Vector6d::Zero( ) );
                Eigen::Matrix3d currentRotationMatrix =
                        getMatrixFromVectorRotationRepresentation( outputIterator->second.segment( 9, 9 ) );
                thrustDifference = manualRotationMatrix
                        * thrustFunction( outputIterator->first ) - outputIterator->second.segment( 0, 3 );
                double currentThrustMagnitude = thrustFunction( outputIterator->first ).norm( ) ;

                for( unsigned int i = 0; i < 3; i++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( thrustDifference( i ) ), 1.0E-12 * currentThrustMagnitude );
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( manualRotationMatrix( i, j ) -
                                                      currentRotationMatrix( i, j ) ), 1.0E-12 );
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
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 1.0E4,
                                    simulationEndEpoch + 1.0E4 );
    bodySettings.at( "Earth" )->gravityFieldSettings =
            std::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

    // Create Earth object
    simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );

    // Create vehicle objects.
    bodies.createEmptyBody( "Apollo" );
    double vehicleMass = 5.0E5;
    bodies.at( "Apollo" )->setConstantBodyMass( vehicleMass );

    // Create vehicle aerodynamic coefficients
    bodies.at( "Apollo" )->setAerodynamicCoefficientInterface( unit_tests::getApolloCoefficientInterface( ) );
    bodies.at( "Apollo" )->setRotationalEphemeris(
                createRotationModel(
                    std::make_shared< PitchTrimRotationSettings >( "Earth", "ECLIPJ2000", "VehicleFixed" ),
                    "Apollo", bodies ) );

    int numberOfCasesPerSet = 4;
    for( int i = 0; i < numberOfCasesPerSet * 2; i++ )
    {
        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define acceleration model settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
        accelerationsOfApollo[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

        // Define specific impulse dependencies.
        std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependencies;
        specificImpulseDependencies.push_back( propulsion::mach_number_dependent_thrust );
        specificImpulseDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );

        // Define variables (thrust dependencies and guidance object) that are different per case.
        std::vector< propulsion::ThrustIndependentVariables > thrustDependencies;
        std::shared_ptr< ThrustMultiplierComputation > thrustInputParameterGuidance;

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
                    std::make_shared< ThrustMultiplierComputation >( simulationStartEpoch, simulationEndEpoch, 0, 1 );
        }
        // Use guidance input to generate a 'fake' Mach number
        else if( ( i % numberOfCasesPerSet == 2 ) )
        {
            thrustDependencies.push_back( propulsion::guidance_input_dependent_thrust );
            thrustDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );

            thrustInputParameterGuidance =
                    std::make_shared< ThrustMultiplierComputation >( simulationStartEpoch, simulationEndEpoch, 1, 0 );
        }
        // Use guidance input to generate a 'fake' Mach number and a throttle
        else if( ( i % numberOfCasesPerSet == 3 ) )
        {
            thrustDependencies.push_back( propulsion::guidance_input_dependent_thrust );
            thrustDependencies.push_back( propulsion::throttle_dependent_thrust );
            thrustDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );

            thrustInputParameterGuidance =
                    std::make_shared< ThrustMultiplierComputation >( simulationStartEpoch, simulationEndEpoch, 1, 1 );

        }

        // Load thrust tables from files (as a function of Mach number and specific impulse) and create interpolators.
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > thrustValues =
                MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    tudat::paths::getTudatTestDataPath( ) + "/Tmax_test.txt" );
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > specificImpulseValues =
                MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    tudat::paths::getTudatTestDataPath( ) + "/Isp_test.txt" );

        std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator =
                std::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    thrustValues.second, thrustValues.first );
        std::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator =
                std::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    specificImpulseValues.second, specificImpulseValues.first );

        double constantSpecificImpulse = 1000.0;

        // Test with interpolated specific impulse
        if( i < numberOfCasesPerSet )
        {
            // Thrust requires no guidance.
            if( !( i % numberOfCasesPerSet == 0 ) )
            {
                addEngineModel(
                            "Apollo", "MainEngine",
                            createParameterizedThrustMagnitudeSettings(
                                thrustInputParameterGuidance, thrustMagnitudeInterpolator, thrustDependencies,
                                specificImpulseInterpolator, specificImpulseDependencies ), bodies );
                accelerationsOfApollo[ "Apollo" ].push_back(
                            std::make_shared< ThrustAccelerationSettings >(
                                 "MainEngine" ) );
            }
            // Thrust requires guidance.
            else
            {
                addEngineModel(
                            "Apollo", "MainEngine",
                            std::make_shared< ParameterizedThrustMagnitudeSettings >(
                                thrustMagnitudeInterpolator, thrustDependencies,
                                specificImpulseInterpolator, specificImpulseDependencies ), bodies );
                accelerationsOfApollo[ "Apollo" ].push_back(
                            std::make_shared< ThrustAccelerationSettings >(
                                 "MainEngine" ) );
            }



        }
        // Test with constant specific impulse
        else
        {
            // Thrust requires no guidance.
            if( !( i % numberOfCasesPerSet == 0 ) )
            {
                addEngineModel(
                            "Apollo", "MainEngine",
                            createParameterizedThrustMagnitudeSettings(
                                thrustInputParameterGuidance, thrustMagnitudeInterpolator, thrustDependencies,
                                constantSpecificImpulse ), bodies );
                accelerationsOfApollo[ "Apollo" ].push_back(
                            std::make_shared< ThrustAccelerationSettings >(
                                  "MainEngine") );
            }
            // Thrust requires guidance.
            else
            {
                addEngineModel(
                            "Apollo", "MainEngine",
                            std::make_shared< ParameterizedThrustMagnitudeSettings >(
                                thrustMagnitudeInterpolator, thrustDependencies,
                                constantSpecificImpulse ), bodies );
                accelerationsOfApollo[ "Apollo" ].push_back(
                            std::make_shared< ThrustAccelerationSettings >(
                                 "MainEngine") );
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
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );

//        // Set trimmed conditions for body orientation.
//        setTrimmedConditions( bodies.at( "Apollo" ) );

        // Define list of dependent variables to save.
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        mach_number_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        airspeed_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        local_density_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        thrust_acceleration, "Apollo", "Apollo", 1 ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        total_mass_rate_dependent_variables, "Apollo" ) );

        // Define propagation settings.
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                  std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ), cowell,
                  std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

        std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
        massRateModels[ "Apollo" ] = createMassRateModel( "Apollo", std::make_shared< FromThrustMassRateSettings >(1 ),
                                                          bodies, accelerationModelMap );

        std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
                std::make_shared< MassPropagatorSettings< double > >(
                    std::vector< std::string >{ "Apollo" }, massRateModels,
                    ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ),
                    std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ) );

        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( translationalPropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                std::make_shared< MultiTypePropagatorSettings< double > >(
                    propagatorSettingsVector, std::make_shared< propagators::PropagationTimeTerminationSettings >(
                        simulationEndEpoch ),
                    std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

        // Define integration settings.
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, simulationStartEpoch, fixedStepSize );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings, true, false, false );

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
                                        24.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( currentMassRate, expectedMassRate,
                                        12.0 * std::numeric_limits< double >::epsilon( ) );
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
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 1.0E4,
                                    simulationEndEpoch + 1.0E4 );
    bodySettings.at( "Earth" )->gravityFieldSettings =
            std::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

    // Create Earth object
    simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );

    // Create vehicle objects.
    bodies.createEmptyBody( "Apollo" );
    double vehicleMass = 5.0E5;
    bodies.at( "Apollo" )->setConstantBodyMass( vehicleMass );

    // Create vehicle aerodynamic coefficients
    bodies.at( "Apollo" )->setAerodynamicCoefficientInterface( unit_tests::getApolloCoefficientInterface( ) );
    bodies.at( "Apollo" )->setRotationalEphemeris(
                createRotationModel(
                    std::make_shared< PitchTrimRotationSettings >( "Earth", "ECLIPJ2000", "VehicleFixed" ),
                    "Apollo", bodies ) );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
    accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationsOfApollo[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    std::vector< propulsion::ThrustIndependentVariables > thrustDependencies;
    thrustDependencies.push_back( propulsion::mach_number_dependent_thrust );
    thrustDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );
    thrustDependencies.push_back( propulsion::throttle_dependent_thrust );

    std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependencies;
    specificImpulseDependencies.push_back( propulsion::mach_number_dependent_thrust );
    specificImpulseDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );

    std::string thrustFile =
            tudat::paths::getTudatTestDataPath( ) + "/Tmax_test.txt";
    std::string specificImpulseFile =
            tudat::paths::getTudatTestDataPath( ) + "/Isp_test.txt";

    addEngineModel(
                "Apollo", "MainEngine",
                createAccelerationLimitedParameterizedThrustMagnitudeSettings(
                    bodies, "Apollo", physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION,
                    thrustFile, thrustDependencies,
                    specificImpulseFile, specificImpulseDependencies, "Earth" ), bodies );

    accelerationsOfApollo[ "Apollo" ].push_back(
                std::make_shared< ThrustAccelerationSettings >(
                     "MainEngine" ) );

    accelerationMap[ "Apollo" ] = accelerationsOfApollo;

    bodiesToPropagate.push_back( "Apollo" );
    centralBodies.push_back( "Earth" );

    // Set initial state
    Eigen::Vector6d systemInitialState = apolloInitialState;


    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    thrust_acceleration, "Apollo", "Apollo", 1 ) );


    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ), cowell,
              std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Apollo" ] = createMassRateModel( "Apollo", std::make_shared< FromThrustMassRateSettings >(1 ),
                                                      bodies, accelerationModelMap );

    std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                std::vector< std::string >{ "Apollo" }, massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ),
                std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ) );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( translationalPropagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector, std::make_shared< propagators::PropagationTimeTerminationSettings >(
                    simulationEndEpoch ),
                std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );


    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, false );

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

//////! Test to check whether the mee-costate based thrust guidance is working correctly
////BOOST_AUTO_TEST_CASE( testMeeCostateBasedThrust )
////{
////    using namespace tudat;
////    using namespace ephemerides;
////    using namespace interpolators;
////    using namespace numerical_integrators;
////    using namespace spice_interface;
////    using namespace simulation_setup;
////    using namespace basic_astrodynamics;
////    using namespace orbital_element_conversions;
////    using namespace propagators;
////    using namespace aerodynamics;
////    using namespace basic_mathematics;
////    using namespace input_output;
////    using namespace unit_conversions;

////    // Load Spice kernels.
////    spice_interface::loadStandardSpiceKernels( );

////    // Set simulation start epoch.
////    const double simulationStartEpoch = 0.0;

////    // Set simulation end epoch.
////    const double simulationEndEpoch = 24.0 * 3600.0;

////    // Set numerical integration fixed step size.
////    const double fixedStepSize = 15.0;


////    // Set spherical elements for Apollo.
////    Eigen::Vector6d asterixInitialStateInKeplerianElements;
////    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
////    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
////    asterixInitialStateInKeplerianElements( inclinationIndex ) = 0.5;
////    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = 0.5;
////    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.5;
////    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

////    // Convert asterix state from spherical elements to Cartesian elements.
////    Vector6d asterixInitialState = orbital_element_conversions::convertKeplerianToCartesianElements(
////                asterixInitialStateInKeplerianElements, getBodyGravitationalParameter( "Earth" ) );

////    // Define simulation body settings.
////    BodyListSettings bodySettings =
////            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 1.0E4,
////                                    simulationEndEpoch + 1.0E4, "Earth", "ECLIPJ2000" );
////    bodySettings.at( "Earth" )->gravityFieldSettings =
////            std::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

////    // Create Earth object
////    simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );

////    // Create vehicle objects.
////    bodies.createEmptyBody( "Asterix" );


////    double vehicleMass = 5.0E5;

////    // Run simulations for a single MEE costate not equal to zero, for each of the first 5 elements
////    for( unsigned int i = 0; i < 5; i++ )
////    {
////        bodies.at( "Asterix" )->setConstantBodyMass( vehicleMass );

////        // Define propagator settings variables.
////        SelectedAccelerationMap accelerationMap;
////        std::vector< std::string > bodiesToPropagate;
////        std::vector< std::string > centralBodies;

////        // Define acceleration model settings.
////        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
////        accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

////        Eigen::VectorXd costates = Eigen::VectorXd::Zero( 5 );
////        costates( i ) = 100.0;

////        std::shared_ptr< ThrustDirectionSettings > thrustDirectionGuidanceSettings =
////                std::make_shared< MeeCostateBasedThrustDirectionSettings >(
////                    "Asterix", "Earth", [ & ]( const double ){ return costates; } );
////        std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings =
////                std::make_shared< ConstantThrustMagnitudeSettings >( 1.0E4, 30000.0 );

////        accelerationsOfAsterix[ "Asterix" ].push_back(
////                    std::make_shared< ThrustAccelerationSettings >(
////                        thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );

////        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

////        bodiesToPropagate.push_back( "Asterix" );
////        centralBodies.push_back( "Earth" );

////        // Set initial state
////        Eigen::Vector6d systemInitialState = asterixInitialState;


////        // Create acceleration models and propagation settings.
////        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
////                    bodies, accelerationMap, bodiesToPropagate, centralBodies );

////        // Define list of dependent variables to save.
////        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
////        dependentVariables.push_back(
////                    std::make_shared< SingleDependentVariableSaveSettings >(
////                        modified_equinocial_state_dependent_variable, "Asterix", "Earth" ) );
////        dependentVariables.push_back(
////                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
////                        thrust_acceleration, "Asterix", "Asterix" ) );

////        // Create propagator/integrator settings
////        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
////                std::make_shared< TranslationalStatePropagatorSettings< double > >
////                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
////                  std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ), cowell,
////                  std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );
////        std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
////        massRateModels[ "Asterix" ] = createMassRateModel(
////                    "Asterix", std::make_shared< FromThrustMassRateSettings >(1 ),
////                    bodies, accelerationModelMap );
////        std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
////                std::make_shared< MassPropagatorSettings< double > >(
////                    std::vector< std::string >{ "Asterix" }, massRateModels,
////                    ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ),
////                    std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ) );
////        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
////        propagatorSettingsVector.push_back( translationalPropagatorSettings );
////        propagatorSettingsVector.push_back( massPropagatorSettings );
////        std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
////                std::make_shared< MultiTypePropagatorSettings< double > >(
////                    propagatorSettingsVector, std::make_shared< propagators::PropagationTimeTerminationSettings >(
////                        simulationEndEpoch ),
////                    std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

////        std::shared_ptr< IntegratorSettings< > > integratorSettings =
////                std::make_shared< IntegratorSettings< > >
////                ( rungeKutta4, simulationStartEpoch, fixedStepSize );

////        // Create simulation object and propagate dynamics.
////        SingleArcDynamicsSimulator< > dynamicsSimulator(
////                    bodies, integratorSettings, propagatorSettings, true, false, false );

////        // Retrieve change in Modified Equinoctial Elements
////        std::map< double, Eigen::VectorXd > dependentVariableSolution =
////                dynamicsSimulator.getDependentVariableHistory( );
////        Eigen::Vector6d finalModifiedEquinoctialElementsError =
////                dependentVariableSolution.rbegin( )->second.segment( 0, 6 ) -
////                dependentVariableSolution.begin( )->second.segment( 0, 6 );

////        // Test whether MEE rates are within reasonable bounds (values are determined empirically)
////        BOOST_CHECK_EQUAL( ( finalModifiedEquinoctialElementsError( i ) < 0 ), true );
////        if( i == 0 )
////        {
////            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 0 ) ) > 2.5E6 ), true );
////        }
////        else if( i < 3 )
////        {
////            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 0 ) ) < 2E5 ), true );
////        }
////        else
////        {
////            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 0 ) ) < 0.1 ), true );
////        }

////        if( i == 1 )
////        {
////            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 1 ) ) > 0.1 ), true );
////        }
////        else
////        {
////            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 1 ) ) < 0.025 ), true );
////        }

////        if( i == 2 )
////        {
////            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 2 ) ) > 0.1 ), true );
////        }
////        else
////        {
////            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 2 ) ) < 0.025 ), true );
////        }

////        if( i == 3 )
////        {
////            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 3 ) ) > 0.075 ), true );
////        }
////        else
////        {
////            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 3 ) ) < 0.005 ), true );
////        }

////        if( i == 4 )
////        {
////            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 4 ) ) > 0.075 ), true );
////        }
////        else
////        {
////            BOOST_CHECK_EQUAL( ( std::fabs( finalModifiedEquinoctialElementsError( 4 ) ) < 0.005 ), true );
////        }
////    }
////}

//! Test to check whether the mee-costate based thrust guidance is working correctly
BOOST_AUTO_TEST_CASE( testMomentumWheelDesaturationThrust )
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
    using namespace estimatable_parameters;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 4.0 * 3600.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 2.0;


    // Create vehicle objects.
    simulation_setup::SystemOfBodies bodies;
    bodies.createEmptyBody( "Asterix" );




    double vehicleMass = 5.0E5;

    bodies.at( "Asterix" )->setConstantBodyMass( vehicleMass );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define times and deltaV magnitudes for momentum wheel desaturation maneuvers.
    std::vector< double > thrustMidTimes = { 1.0 * 3600.0, 2.0 * 3600.0, 3.0 * 3600.0 };
    std::vector< Eigen::Vector3d > deltaVValues =
    { 1.0E-3 * ( Eigen::Vector3d( ) << 0.3, -2.5, 3.4 ).finished( ),
      1.0E-3 * ( Eigen::Vector3d( ) << 2.0, 5.9, -0.5 ).finished( ),
      1.0E-3 * ( Eigen::Vector3d( ) << -1.6, 4.4, -5.8 ).finished( ) };
    double totalManeuverTime = 90.0;
    double maneuverRiseTime = 15.0;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Asterix" ].push_back(
                std::make_shared< MomentumWheelDesaturationAccelerationSettings >(
                    thrustMidTimes, deltaVValues, totalManeuverTime, maneuverRiseTime ) );
    accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "SSB" );

    // Set initial state
    Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    momentum_wheel_desaturation_acceleration, "Asterix", "Asterix" ) );

    // Create propagator/integrator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ), cowell,
              std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

    std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
            translationalPropagatorSettings;

    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch + fixedStepSize / 9.0, fixedStepSize );


    // Define list of parameters to estimate.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Asterix", desaturation_delta_v_values ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, propagatorSettings );

    // Create simulation object and propagate dynamics.
    SingleArcVariationalEquationsSolver< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, parametersToEstimate,
                true, std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
                false, true, false );

    auto stateHistory = dynamicsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
    auto dependentVariableResult = dynamicsSimulator.getDynamicsSimulator( )->getDependentVariableHistory( );

    auto stateTransitionHistory = dynamicsSimulator.getNumericalVariationalEquationsSolution( )[ 0 ];
    auto sensitivityHistory = dynamicsSimulator.getNumericalVariationalEquationsSolution( )[ 1 ];

    // Compute thrust start times from maneuvers mid-times.
    std::vector< double > thrustStartTimes;
    for( unsigned int i = 0; i < thrustMidTimes.size( ); i++ )
    {
        thrustStartTimes.push_back( thrustMidTimes.at( i ) - totalManeuverTime / 2.0 );
    }
    thrustStartTimes.push_back( std::numeric_limits< double >::max( ) );

    // Create interpolator to look up maneuvers start times.
    std::shared_ptr< tudat::interpolators::LookUpScheme< double > > timeLookup =
            std::make_shared< tudat::interpolators::HuntingAlgorithmLookupScheme< double > >(
                thrustStartTimes );

    for( auto variableIterator : dependentVariableResult )
    {
        // Identify maneuver start time closest to current time.
        double currentTime = variableIterator.first;
        int currentNearestNeighbour = timeLookup->findNearestLowerNeighbour( currentTime );

        double currentStartTime = thrustStartTimes.at( currentNearestNeighbour );

        Eigen::Vector3d expectedAcceleration = Eigen::Vector3d::Zero( );
        double scalingNorm = 0.0;

        // If maneuver still ongoing at current time.
        if( ( std::fabs( currentTime - currentStartTime ) < totalManeuverTime ) && ( currentTime > currentStartTime )  )
        {
            // Compute peak desaturation acceleration.
            Eigen::Vector3d peakAcceleration = deltaVValues.at( currentNearestNeighbour ) /
                    ( totalManeuverTime - maneuverRiseTime );
            scalingNorm = peakAcceleration.norm( );

            // Compute time elapsed since maneuver start.
            double timeSinceStart = currentTime - currentStartTime;

            // Compute expected acceleration from peak acceleration and time elapsed since maneuver initiation.
            if( timeSinceStart < maneuverRiseTime )
            {
                double timeRatio = timeSinceStart / maneuverRiseTime;
                expectedAcceleration = peakAcceleration * timeRatio * timeRatio * (
                            3.0 - 2.0 * timeRatio );
            }
            else if( timeSinceStart < totalManeuverTime - maneuverRiseTime )
            {
                expectedAcceleration = peakAcceleration;
            }
            else
            {
                double timeRatio = ( totalManeuverTime - timeSinceStart ) / maneuverRiseTime;
                expectedAcceleration = peakAcceleration * timeRatio * timeRatio * (
                            3.0 - 2.0 * timeRatio );
            }
        }

        // If maneuver already completed at current time.
        else if( currentTime > currentStartTime )
        {
            Eigen::Vector3d expectedDeltaV = Eigen::Vector3d::Zero( );
            for( int i = 0; i <= currentNearestNeighbour; i++ )
            {
                // Compute expected deltaV.
                expectedDeltaV += deltaVValues.at( i );

                // Check that the sensivity matrix blocks which describe the velocity partials w.r.t. the deltaV values
                // of all the maneuvers encountered until current time are almost identity blocks.
                if( currentTime - currentStartTime > totalManeuverTime )
                {
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                sensitivityHistory.at( currentTime ).block( 3, i * 3, 3, 3 ), Eigen::Matrix3d::Identity( ), 1.0E-4 );
                }
            }
            for( int i = currentNearestNeighbour + 1; i <= 2; i++ )
            {
                // Check that the sensitivity matrix blocks which describe the velocity partials w.r.t. the deltaV values
                // of the upcoming maneuvers are filled with zeros.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            sensitivityHistory.at( currentTime ).block( 3, i * 3, 3, 3 ), Eigen::Matrix3d::Zero( ),
                            std::numeric_limits< double >::epsilon( ) );
            }

            Eigen::Vector3d currentVelocity = stateHistory.at( variableIterator.first ).segment( 3, 3 );

            // Check deltaV values consistency.
            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( expectedDeltaV( i ) - currentVelocity( i ) ), 1.0E-5 * currentVelocity.norm( ) );
            }

        }

        // Check accelerations consistency.
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( expectedAcceleration( i ) - variableIterator.second( i ) ),
                               5.0 * std::numeric_limits< double >::epsilon( ) * scalingNorm );
        }


        // Check state transition matrix consistency.
        // The state transition matrix is expected to be equal to the identity matrix, expect for the current
        // position partials w.r.t. the initial velocity, expected to show a linear time-dependence.
        Eigen::Matrix6d stateTransitionMatrix = stateTransitionHistory.at( currentTime );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    stateTransitionMatrix.block( 0, 3, 3, 3 ),
                    ( ( currentTime - integratorSettings->initialTime_ ) * Eigen::Matrix3d::Identity( ) ),
                    1.0E-8 );
        stateTransitionMatrix.block( 0, 3, 3, 3 ) = Eigen::Matrix3d::Zero( );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    stateTransitionMatrix, Eigen::Matrix6d::Identity( ),
                    std::numeric_limits< double >::epsilon( ) );
    }
}
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
