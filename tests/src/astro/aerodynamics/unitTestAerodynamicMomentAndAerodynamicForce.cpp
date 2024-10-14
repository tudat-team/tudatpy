/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      The unit tests here are based off of expected values that are internally computed. Ideally,
 *      these should be based off of published values that can be found in literature. This will
 *      have to be updated in future for the code to be considered completely tested. In addition,
 *      more test values are required, as more the unit tests are benchmarked off of one set of
 *      data.
 *
 *      The class objects aerodynamicCoefficientInterface, aerodynamicForce, and aerodynamicMoment
 *      are declared multiple times in local scope currently since their member variables aren't
 *      set at construction, but rather through set-functions. Once these are adapted to be set
 *      through the constructor, const class objects can be declared that can be shared between the
 *      unit tests.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>


#include <memory>
#include <boost/lambda/lambda.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
//#include <Eigen/Geometry>

#include "tudat/basics/testMacros.h"
//#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
#include "tudat/astro/aerodynamics/testApolloCapsuleCoefficients.h"

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_aerodynamic_acceleration_force_moment_models )

using namespace aerodynamics;
using namespace numerical_integrators;
using namespace simulation_setup;
using namespace propagators;
using namespace basic_mathematics;
using namespace basic_astrodynamics;
using namespace orbital_element_conversions;
using namespace spice_interface;

//! Test implementation of aerodynamic force and acceleration models.
BOOST_AUTO_TEST_CASE( testAerodynamicForceAndAcceleration )
{
    // Set force coefficients.
    const Eigen::Vector3d forceCoefficients( 1.1, 1.2, 1.3 );

    // Set dynamical model parameters.
    const double density = 3.5e-5;
    const double airSpeed = 3.491e3;
    const double dynamicPressure = 0.5 * density * airSpeed * airSpeed;
    const double referenceArea = 2.2;
    const double referenceLength = 3.2;
    const double mass = 1.93;

    // Compute expected force.
    const Eigen::Vector3d expectedForce = forceCoefficients * dynamicPressure * referenceArea;

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    // Test 1: test the force model implemented as free function with primitive arguments.
    {
        // Compute force.
        Eigen::Vector3d force = computeAerodynamicForce( dynamicPressure,
                                                         referenceArea, forceCoefficients );

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 2: test the force model implemented as free function, with coefficient interface
    //         argument.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        AerodynamicCoefficientInterfacePointer aerodynamicCoefficientInterface =
                createConstantCoefficientAerodynamicCoefficientInterface(
                    forceCoefficients, Eigen::Vector3d::Zero( ),
                    referenceLength, referenceArea, Eigen::Vector3d::Zero( ) );

        // Compute aerodynamic force using free function with coefficient interface argument.
        Eigen::Vector3d force = computeAerodynamicForce( dynamicPressure,
                                                         aerodynamicCoefficientInterface );

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );

        // Test aerodynamic coefficient interface properties
        BOOST_CHECK_EQUAL(
                    aerodynamicCoefficientInterface->getIndependentVariableNames( ).size( ), 0 );

        bool isVariableIndexTooHigh = 0;
        try
        {
            aerodynamicCoefficientInterface->getIndependentVariableName( 0 );
        }
        catch( std::runtime_error const& )

        {
            isVariableIndexTooHigh = 1;
        }
        BOOST_CHECK_EQUAL( isVariableIndexTooHigh, 1 );
    }

    // Test 3: test the acceleration model implemented as free function with primitive arguments,
    //         based on the force that can be derived from the computed acceleration.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        AerodynamicCoefficientInterfacePointer aerodynamicCoefficientInterface =
                createConstantCoefficientAerodynamicCoefficientInterface(
                    forceCoefficients, Eigen::Vector3d::Zero( ),
                    referenceLength, referenceArea, Eigen::Vector3d::Zero( ) );

        // Compute aerodynamic force from aerodynamic acceleration free function with primitive
        // arguments.
        Eigen::Vector3d force = computeAerodynamicAcceleration(
                    dynamicPressure, aerodynamicCoefficientInterface, mass ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 4: test the acceleration model implemented as free function with coefficient interface
    //         argument, based on the force that can be derived from the computed acceleration.
    {
        // Compute aerodynamic force from aerodynamic acceleration free function with
        // coefficient interface argument.
        Eigen::Vector3d force = computeAerodynamicAcceleration(
                    dynamicPressure, referenceArea, forceCoefficients, mass ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 5: Test the acceleration model class without inverted coefficients.
    {
        // Create aaerodynamic acceleration model class, no inverted coefficients, direct mass
        // and reference area.
        AerodynamicAccelerationPointer accelerationClass
                = std::make_shared< AerodynamicAcceleration >(
                    [ & ]( Eigen::Vector3d& input ){ input = forceCoefficients; },
                [ & ]( ){ return density; },
        [ & ]( ){ return airSpeed; },
        mass, referenceArea, false );
        accelerationClass->updateMembers( );
        Eigen::Vector3d force = accelerationClass->getAcceleration( ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );

        // Create aerodynamic acceleration model class, no inverted coefficients, mass and
        // reference area set through std::functions.
        AerodynamicAccelerationPointer accelerationClass2 =
                std::make_shared< AerodynamicAcceleration >(
                    [ & ]( Eigen::Vector3d& input ){ input = forceCoefficients; },
                [ & ]( ){ return density; },
        [ & ]( ){ return airSpeed; },
        [ & ]( ){ return mass; },
        [ & ]( ){ return referenceArea; },
        false );
        accelerationClass2->updateMembers( );
        force = accelerationClass2->getAcceleration( ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 6: Test the acceleration model class with inverted coefficients
    {
        // Create aaerodynamic acceleration model class, inverted coefficients, direct mass
        // and reference area.
        AerodynamicAccelerationPointer accelerationClass =
                std::make_shared< AerodynamicAcceleration >(
                    [ & ]( Eigen::Vector3d& input ){ input = -forceCoefficients; },
                [ & ]( ){ return density; },
        [ & ]( ){ return airSpeed; },
        mass, referenceArea, true );
        accelerationClass->updateMembers( );
        Eigen::Vector3d force = accelerationClass->getAcceleration( ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );

        // Create aerodynamic acceleration model class, inverted coefficients, mass and
        // reference area set through std::functions.
        AerodynamicAccelerationPointer accelerationClass2 =
                std::make_shared< AerodynamicAcceleration >(
                    [ & ]( Eigen::Vector3d& input ){ input = -forceCoefficients; },
                [ & ]( ){ return density; },
        [ & ]( ){ return airSpeed; },
        [ & ]( ){ return mass; },
        [ & ]( ){ return referenceArea; },
        true );
        accelerationClass2->updateMembers( );
        force = accelerationClass2->getAcceleration( ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }
}

//! Test implementation of aerodynamic moment and rotational acceleration models.
BOOST_AUTO_TEST_CASE( testAerodynamicMomentAndRotationalAcceleration )
{
    // Set moment coefficients.
    const Eigen::Vector3d momentCoefficients( -3.2, 1.0, 8.4 );

    // Set dynamical model parameters.
    const double dynamicPressure = 123.6;
    const double referenceArea = 1.7;
    const double referenceLength = 2.6;

    // Calculate expected moment.
    const Eigen::Vector3d expectedMoment = dynamicPressure * referenceArea *
            referenceLength * momentCoefficients;

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    // Test 1: test the moment model implemented as free function with primitive arguments.
    {
        // Compute aerodynamic moment using free function with primitive arguments.
        Eigen::Vector3d moment = computeAerodynamicMoment( dynamicPressure, referenceArea,
                                                           referenceLength, momentCoefficients );

        // Check if computed moment matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMoment, moment, tolerance );
    }

    // Test 2: test the moment moment implemented as free function with coefficient interface
    //         argument.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        AerodynamicCoefficientInterfacePointer aerodynamicCoefficientInterface =
                createConstantCoefficientAerodynamicCoefficientInterface(
                    Eigen::Vector3d::Zero( ), momentCoefficients,
                    referenceLength, referenceArea, Eigen::Vector3d::Zero( ) );

        // Compute aerodynamic moment using free function with coefficient interface argument.
        Eigen::Vector3d moment = computeAerodynamicMoment( dynamicPressure,
                                                           aerodynamicCoefficientInterface );

        // Check if computed moment matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMoment, moment, tolerance );
    }
}

class DummyAngleCalculator
{
public:

    Eigen::Vector3d getAerodynamicAngles( const double time )
    {
        time_ = time;

        double currentAngleOfAttack = getDummyAngleOfAttack( );
        double currentAngleOfSideslip = getDummyAngleOfSideslip( );
        double currentBankAngle = getDummyBankAngle( );

        return ( Eigen::Vector3d( ) << currentAngleOfAttack, currentAngleOfSideslip, currentBankAngle ).finished( );

    }

    double getDummyAngleOfAttack( )
    {
        return 0.2 - time_ / 1000.0;
    }

    double getDummyAngleOfSideslip( )
    {
        return 0.6 + 0.5 * time_ / 1000.0;
    }

    double getDummyBankAngle( )
    {
        return 1.3 + 0.24 * time_ / 1000.0;
    }

    double time_;
};

void testAerodynamicForceDirection( const bool includeThrustForce,
                                    const bool useSpiceRotationModel,
                                    const bool swapCreationOrder )
{
    std::cout<<includeThrustForce<<" "<<useSpiceRotationModel<<" "<<swapCreationOrder<<std::endl;
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

    for( unsigned int i = 0; i < 4; i++ )
    {
        // Create Earth object
        BodyListSettings defaultBodySettings =
                getDefaultBodySettings( std::vector< std::string >{ "Earth" }, -1.0E6, 1.0E6 );
        defaultBodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                    Eigen::Vector6d::Zero( ) );
        SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

        // Create vehicle objects.
        double vehicleMass = 5.0E3;
        bodies.createEmptyBody( "Vehicle" );

        bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );

        std::shared_ptr< DummyAngleCalculator > testAngles =
                std::make_shared< DummyAngleCalculator >( );

        if( useSpiceRotationModel )
        {
            bodies.at( "Vehicle" )->setRotationalEphemeris(
                        std::make_shared< ephemerides::SpiceRotationalEphemeris >( "ECLIPJ2000", "IAU_MOON" ) );
        }
        else
        {
            std::shared_ptr< ephemerides::RotationalEphemeris > vehicleRotationModel =
                    createRotationModel(
                        std::make_shared< AerodynamicAngleRotationSettings >(
                            "Earth", "ECLIPJ2000", "VehicleFixed",
                            std::bind( &DummyAngleCalculator::getAerodynamicAngles, testAngles, std::placeholders::_1 ) ),
                        "Vehicle", bodies );

            bodies.at( "Vehicle" )->setRotationalEphemeris( vehicleRotationModel );
        }


        bool areCoefficientsInAerodynamicFrame;
        if( ( i % 2 ) == 0 )
        {
            areCoefficientsInAerodynamicFrame = 1;
        }
        else
        {
            areCoefficientsInAerodynamicFrame = 0;
        }

        Eigen::Vector3d aerodynamicCoefficients;
        if( ( i / 2 ) % 2 == 0 )
        {
            aerodynamicCoefficients = Eigen::Vector3d::UnitX( );
        }
        else
        {
            aerodynamicCoefficients = ( Eigen::Vector3d( ) << 1.0, -0.1, 0.5 ).finished( );

        }

        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                std::make_shared< ConstantAerodynamicCoefficientSettings >(
                    2.0, 4.0, Eigen::Vector3d::Zero( ), aerodynamicCoefficients, Eigen::Vector3d::Zero( ),
                    aerodynamics::getAerodynamicCoefficientFrame( areCoefficientsInAerodynamicFrame, 1 ),
                    aerodynamics::getAerodynamicCoefficientFrame( areCoefficientsInAerodynamicFrame, 1 ) );
        bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ) );
        Eigen::Vector3d aerodynamicCoefficientsDirection = aerodynamicCoefficients.normalized( );



        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define acceleration model settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

        if( swapCreationOrder )
        {
            accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
        }

        Eigen::Vector3d bodyFixedThrustDirection = ( Eigen::Vector3d( ) << 1.4, 3.1, -0.5 ).finished( ).normalized( );

        if( includeThrustForce )
        {
            addEngineModel(
                        "Vehicle", "MainEngine",
                        std::make_shared< ConstantThrustMagnitudeSettings >(
                            thrustMagnitude, specificImpulse ),
                        bodies, bodyFixedThrustDirection );

            accelerationsOfVehicle[ "Vehicle" ].push_back(
                        std::make_shared< ThrustAccelerationSettings >(
                             "MainEngine" ) );
        }

        if( !swapCreationOrder )
        {
            accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
        }



        accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

        bodiesToPropagate.push_back( "Vehicle" );
        centralBodies.push_back( "Earth" );

        // Set initial state
        Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );

        systemInitialState( 0 ) = 6.8E6;
        systemInitialState( 4 ) = 7.5E3;

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );

        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;

        dependentVariables.push_back(
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        aerodynamic, "Vehicle", "Earth", 0 ) );
        dependentVariables.push_back(
                    std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                        "Vehicle", reference_frames::inertial_frame, reference_frames::aerodynamic_frame ) );
        dependentVariables.push_back(
                    std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                        "Vehicle", reference_frames::inertial_frame, reference_frames::body_frame ) );
        dependentVariables.push_back(
                    std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                        "Vehicle", reference_frames::inertial_frame, reference_frames::corotating_frame ) );
        dependentVariables.push_back(
                    std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                        "Vehicle", reference_frames::inertial_frame, reference_frames::trajectory_frame ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        body_fixed_airspeed_based_velocity_variable, "Vehicle" ) );
        if( includeThrustForce )
        {
            dependentVariables.push_back(
                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
        }


        std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                std::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, 5.0 );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, 0.0, integratorSettings, terminationSettings,
                  cowell, dependentVariables );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, translationalPropagatorSettings );

        // Retrieve numerical solutions for state and dependent variables
        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput =
                dynamicsSimulator.getDependentVariableHistory( );

        std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris =
                bodies.at( "Vehicle" )->getRotationalEphemeris( );

        double thrustAccelerationMagnitude = thrustMagnitude / vehicleMass;
        Eigen::Matrix3d matrixDifference;
        for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
             dependentVariableOutput.begin( ); outputIterator != dependentVariableOutput.end( ); outputIterator++ )
        {


            // Retrieve dependent variables from output;
            Eigen::Matrix3d rotationToAerodynamicFrame =
                    getMatrixFromVectorRotationRepresentation(
                        outputIterator->second.segment( 3, 9 ) );
            Eigen::Matrix3d rotationToBodyFrame =
                    getMatrixFromVectorRotationRepresentation(
                        outputIterator->second.segment( 12, 9 ) );
            Eigen::Matrix3d rotationToCorotatingFrame = getMatrixFromVectorRotationRepresentation(
                        outputIterator->second.segment( 21, 9 ) );
            Eigen::Matrix3d rotationToTrajectoryFrame = getMatrixFromVectorRotationRepresentation(
                        outputIterator->second.segment( 30, 9 ) );
            Eigen::Vector3d bodyFixedAirspeed = outputIterator->second.segment( 39, 3 );

            // Velocity vector in aerodynamic and trajectory frames should have component in positivie x-direction only.
            Eigen::Vector3d bodyFixedAirspeedInAerodynamicFrame = rotationToAerodynamicFrame *
                    rotationToCorotatingFrame.transpose( ) * bodyFixedAirspeed;
            Eigen::Vector3d bodyFixedAirspeedInTrajectoryFrame = rotationToTrajectoryFrame *
                    rotationToCorotatingFrame.transpose( ) * bodyFixedAirspeed;

            // Check velocity in aerodynamic frame
            BOOST_CHECK_CLOSE_FRACTION( bodyFixedAirspeedInAerodynamicFrame.x( ),
                                        bodyFixedAirspeedInAerodynamicFrame.norm( ),
                                        std::numeric_limits< double >::epsilon( ) );

            // Check velocity in trajectory frame
            BOOST_CHECK_CLOSE_FRACTION( bodyFixedAirspeedInTrajectoryFrame.x( ),
                                        bodyFixedAirspeedInTrajectoryFrame.norm( ),
                                        std::numeric_limits< double >::epsilon( ) );

            // For drag-only aerodynamics
            if( ( i % 4 ) == 0 )
            {
                Eigen::Vector3d aerodynamicForceInAerodynamicFrame = rotationToAerodynamicFrame * outputIterator->second.segment( 0, 3 );
                BOOST_CHECK_CLOSE_FRACTION( -aerodynamicForceInAerodynamicFrame.x( ),
                                            aerodynamicForceInAerodynamicFrame.norm( ),
                                            std::numeric_limits< double >::epsilon( ) );
            }

            // For C_{x}-only aerodynamics
            else if( ( i % 4 ) == 1 )
            {
                Eigen::Vector3d aerodynamicForceInBodyFrame = rotationToBodyFrame * outputIterator->second.segment( 0, 3 );
                BOOST_CHECK_CLOSE_FRACTION( -aerodynamicForceInBodyFrame.x( ),
                                            aerodynamicForceInBodyFrame.norm( ),
                                            std::numeric_limits< double >::epsilon( ) );
            }

            // Check that aerodynamic force is in correct direction (in aerodynamic frame).
            else if( ( i % 4 ) == 2 )
            {
                Eigen::Vector3d aerodynamicForceDirectionInAerodynamicFrame =
                        ( rotationToAerodynamicFrame * outputIterator->second.segment( 0, 3 ) ).normalized( );
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_CLOSE_FRACTION( -aerodynamicForceDirectionInAerodynamicFrame( j ),
                                                aerodynamicCoefficientsDirection( j ), 5.0E-14 );
                }
            }

            // Check that aerodynamic force is in correct direction (in body frame).
            else if( ( i % 4 ) == 3 )
            {
                Eigen::Vector3d aerodynamicForceDirectionInBodyFrame =
                        ( rotationToBodyFrame * outputIterator->second.segment( 0, 3 ) ).normalized( );
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_CLOSE_FRACTION( -aerodynamicForceDirectionInBodyFrame( j ),
                                                aerodynamicCoefficientsDirection( j ), 5.0E-14 );
                }
            }

            // Check if thrust force is in correct direction.
            if( includeThrustForce && !useSpiceRotationModel )
            {
                Eigen::Vector3d thrustForceInBodyFrame = rotationToBodyFrame * outputIterator->second.segment( 42, 3 );

                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( thrustForceInBodyFrame( j ) - bodyFixedThrustDirection( j ) *
                                           thrustAccelerationMagnitude ), 1.0E-14 );
                }
            }
            else if( includeThrustForce && useSpiceRotationModel )
            {
                Eigen::Vector3d thrustForceInPropagationFrame = ( outputIterator->second.segment( 42, 3 ) );

                Eigen::Matrix3d rotationToBodyFrameFromEphemeris = spice_interface::computeRotationQuaternionBetweenFrames(
                            "ECLIPJ2000", "IAU_MOON", outputIterator->first ).toRotationMatrix( );
                Eigen::Vector3d imposedThrustForceInPropagationFrame =
                        thrustAccelerationMagnitude * ( rotationToBodyFrameFromEphemeris.transpose( ) * bodyFixedThrustDirection );



                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( thrustForceInPropagationFrame( j ) - imposedThrustForceInPropagationFrame( j ) ),
                                1.0E-15 );
                }

                matrixDifference = rotationToBodyFrameFromEphemeris - rotationToBodyFrame;

                for( unsigned int j = 0; j < 3; j++ )
                {
                    for( unsigned int k = 0; k < 3; k++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( j, k ) ), 1.0E-13 );
                    }
                }

            }

            if( useSpiceRotationModel )
            {
                // Check if imposed and indirectly obtained rotation matrices are equal.
                Eigen::Matrix3d rotationToBodyFrameFromEphemeris = rotationalEphemeris->getRotationToTargetFrame(
                            outputIterator->first ).toRotationMatrix( );
                matrixDifference = rotationToBodyFrameFromEphemeris - rotationToBodyFrame;
                for( unsigned int j = 0; j < 3; j++ )
                {
                    for( unsigned int k = 0; k < 3; k++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( j, k ) ), 1.0E-13 );
                    }
                }
            }
            else
            {
                testAngles->getAerodynamicAngles( outputIterator->first );

                Eigen::Matrix3d aerodynamicToBodyFrame = rotationToBodyFrame * rotationToAerodynamicFrame.inverse( );
                Eigen::Matrix3d aerodynamicToTrajectoryFrame = rotationToTrajectoryFrame * rotationToAerodynamicFrame.inverse( );

                Eigen::Matrix3d manualAerodynamicToBodyFrame =
                        reference_frames::getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
                            testAngles->getDummyAngleOfAttack( ), testAngles->getDummyAngleOfSideslip( ) );
                Eigen::Matrix3d manualAerodynamicToTrajectoryFrame =
                        reference_frames::getAerodynamicToTrajectoryFrameTransformationMatrix(
                            testAngles->getDummyBankAngle( ) );

                matrixDifference = aerodynamicToBodyFrame - manualAerodynamicToBodyFrame;

                for( unsigned int j = 0; j < 3; j++ )
                {
                    for( unsigned int k = 0; k < 3; k++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( j, k ) ), 1.0E-14 );
                    }
                }

                matrixDifference = aerodynamicToTrajectoryFrame - manualAerodynamicToTrajectoryFrame;
                for( unsigned int j = 0; j < 3; j++ )
                {
                    for( unsigned int k = 0; k < 3; k++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( j, k ) ), 1.0E-14 );
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testAerodynamicForceDirectionInPropagation )
{
    testAerodynamicForceDirection( 0, 0, 0 );
    testAerodynamicForceDirection( 1, 0, 0 );
    testAerodynamicForceDirection( 1, 1, 0 );
    testAerodynamicForceDirection( 1, 0, 1 );
    testAerodynamicForceDirection( 1, 1, 1 );
}

Eigen::Vector2d sideslipBankAngleFunction( const double time )
{
    double bankAngle = 2.0 * tudat::mathematical_constants::PI * time / 1000.0;
    double sideslipAngle = 0.01 * std::sin( 2.0 * tudat::mathematical_constants::PI * time / 60.0 );
    return ( Eigen::Vector2d( ) << sideslipAngle, bankAngle ).finished( );
}

BOOST_AUTO_TEST_CASE( testAerodynamicTrimWithFreeAngles )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 3300.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 1.0;


    // Set Keplerian elements for Capsule.
    Eigen::Vector6d apolloInitialStateInKeplerianElements;
    apolloInitialStateInKeplerianElements( semiMajorAxisIndex ) = spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    apolloInitialStateInKeplerianElements( eccentricityIndex ) = 0.005;
    apolloInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    apolloInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    apolloInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    apolloInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    // Convert apollo state from Keplerian elements to Cartesian elements.
    const double earthGravitationalParameter = getBodyGravitationalParameter( "Earth" );
    const Eigen::Vector6d apolloInitialState = convertKeplerianToCartesianElements(
                apolloInitialStateInKeplerianElements,
                getBodyGravitationalParameter( "Earth" ) );


    // Define simulation body settings.
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->gravityFieldSettings =
            std::make_shared< simulation_setup::GravityFieldSettings >( central_spice );


    // Create Earth object
    simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );

    // Create vehicle objects.
    bodies.createEmptyBody( "Apollo" );

    // Create vehicle aerodynamic coefficients
    bodies.at( "Apollo" )->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );
    bodies.at( "Apollo" )->setConstantBodyMass( 5.0E3 );
    bodies.at( "Apollo" )->setRotationalEphemeris(
                createRotationModel(
                    std::make_shared< PitchTrimRotationSettings >(
                        "Earth", "ECLIPJ2000", "VehicleFixed", &sideslipBankAngleFunction ),
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
    accelerationMap[ "Apollo" ] = accelerationsOfApollo;

    bodiesToPropagate.push_back( "Apollo" );
    centralBodies.push_back( "Earth" );

    // Set initial state
    Eigen::Vector6d systemInitialState = apolloInitialState;

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
                    aerodynamic_moment_coefficients_dependent_variable, "Apollo" ) );
    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_force_coefficients_dependent_variable, "Apollo" ) );


    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    std::shared_ptr< TranslationalStatePropagatorSettings < double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationStartEpoch, integratorSettings,
              std::make_shared< propagators::PropagationTimeTerminationSettings >( 3200.0 ), cowell,
              dependentVariables );


    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, propagatorSettings );

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput =
            dynamicsSimulator.getDependentVariableHistory( );

    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > aerodynamicCoefficientInterface =
            bodies.at( "Apollo" )->getAerodynamicCoefficientInterface( );

    for( auto it : dependentVariableOutput )
    {
        double machNumber = it.second( 0 );
        double angleOfAttack = it.second( 1 );
        double sideslipAngle = it.second( 2 );
        double bankAngle = it.second( 3 );
        aerodynamicCoefficientInterface->updateFullCurrentCoefficients(
        { machNumber, angleOfAttack,sideslipAngle } );

        Eigen::Vector3d testMomentCoefficients = aerodynamicCoefficientInterface->getCurrentMomentCoefficients( );
        Eigen::Vector3d testForceCoefficients = aerodynamicCoefficientInterface->getCurrentForceCoefficients( );
        Eigen::Vector2d testAngles = sideslipBankAngleFunction( it.first );
        Eigen::Vector3d momentCoefficients = it.second.segment( 4, 3 );
        Eigen::Vector3d forceCoefficients = it.second.segment( 7, 3 );

        BOOST_CHECK_EQUAL( testAngles( 0 ), sideslipAngle );
        BOOST_CHECK_EQUAL( testAngles( 1 ), bankAngle );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( momentCoefficients, testMomentCoefficients, std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( forceCoefficients, testForceCoefficients, std::numeric_limits< double >::epsilon( ) );
    }
}

Eigen::Matrix3d massDependentInertiaTensor( const double mass )
{
    return 5.0E3 * Eigen::Matrix3d::Identity( ) - ( mass - 5.0E3 ) * Eigen::Vector3d::Ones( ) * Eigen::Vector3d::Ones( ).transpose( );
}

Eigen::Vector3d massDependentCenterOfMass( const double mass )
{
    Eigen::Vector3d com = ( Eigen::Vector3d( ) <<  -0.6624, 0.0, 0.1369 ).finished( ) - Eigen::Vector3d::Ones( ) * ( mass - 5.0E3 ) / 1000.0;
    return com;
}

double customMassDerivativeFunction( const double time )
{
    return -1.0;
}

BOOST_AUTO_TEST_CASE( testCombinedAerodynamicForceAndMoment )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 10.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 1000.0;


    // Set Keplerian elements for Capsule.
    Eigen::Vector6d apolloInitialStateInKeplerianElements;
    apolloInitialStateInKeplerianElements( semiMajorAxisIndex ) = spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    apolloInitialStateInKeplerianElements( eccentricityIndex ) = 0.005;
    apolloInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    apolloInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    apolloInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    apolloInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    // Convert apollo state from Keplerian elements to Cartesian elements.
    const double earthGravitationalParameter = getBodyGravitationalParameter( "Earth" );
    const Eigen::Vector6d apolloInitialState = convertKeplerianToCartesianElements(
                apolloInitialStateInKeplerianElements,
                getBodyGravitationalParameter( "Earth" ) );

    for( unsigned propagationType = 0; propagationType < 2 ; propagationType++ )
    {
        for( unsigned int i = 0; i < 3; i++ )
        {

            // Define simulation body settings.
            BodyListSettings bodySettings =
                    getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                            simulationEndEpoch + 10.0 * fixedStepSize, "Earth", "ECLIPJ2000" );
            bodySettings.at( "Earth" )->gravityFieldSettings =
                    std::make_shared< simulation_setup::GravityFieldSettings >( central_spice );


            // Create Earth object
            simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );

            // Create vehicle objects.
            bodies.createEmptyBody( "Apollo" );

            // Create vehicle aerodynamic coefficients
            auto aerodynamicCoefficients = unit_tests::getApolloCoefficientInterface( ) ;
            if( i == 1 )
            {
                aerodynamicCoefficients->resetForceCoefficientsFrame( body_fixed_frame_coefficients );
            }
            else if( i == 2 )
            {
                aerodynamicCoefficients->resetForceCoefficientsFrame( aerodynamics::negative_aerodynamic_frame_coefficients );
                aerodynamicCoefficients->resetMomentCoefficientsFrame( body_fixed_frame_coefficients );

            }

            bodies.at( "Apollo" )->setAerodynamicCoefficientInterface( aerodynamicCoefficients );

            double initialVehicleMass = 5.0E3;

            if( propagationType == 0 )
            {
                bodies.at( "Apollo" )->setConstantBodyMass( initialVehicleMass );
                bodies.at( "Apollo" )->setBodyInertiaTensor( 5.0E3 * Eigen::Matrix3d::Identity( ) );
            }
            else
            {
                auto rigidBodyProperties = std::make_shared< MassDependentMassDistributionSettings >(
                    initialVehicleMass,
                    &massDependentCenterOfMass,
                    &massDependentInertiaTensor );

                addRigidBodyProperties(
                    bodies, "Apollo", rigidBodyProperties );

                addFlightConditions(
                    bodies, "Apollo", "Earth" );
                std::shared_ptr< aerodynamics::AerodynamicMomentContributionInterface > momentCoefficientInterface =
                    createMomentContributionInterface(
                        aerodynamicCoefficients->getForceCoefficientsFrame( ),
                        aerodynamicCoefficients->getMomentCoefficientsFrame( ),
                        bodies.at( "Apollo" ) );
                bodies.at( "Apollo" )->getAerodynamicCoefficientInterface( )->setMomentContributionInterface(
                    momentCoefficientInterface );
            }


            bodies.at( "Apollo" )->setRotationalEphemeris(
                        createRotationModel(
                            constantRotationModelSettings( "ECLIPJ2000", "VehicleFixed", Eigen::Matrix3d::Identity( ) ),
//                            std::make_shared< PitchTrimRotationSettings >(
//                                "Earth", "ECLIPJ2000", "VehicleFixed", &sideslipBankAngleFunction ),
                            "Apollo", bodies ) );

            // Define propagator settings variables.
            SelectedAccelerationMap accelerationMap;
            SelectedTorqueMap torqueMap;
            std::vector< std::string > bodiesToPropagate;
            std::vector< std::string > centralBodies;

            // Define acceleration model settings.
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
            accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
            accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
            accelerationsOfApollo[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
            accelerationMap[ "Apollo" ] = accelerationsOfApollo;

            std::map< std::string, std::vector< std::shared_ptr< TorqueSettings > > > torquesOfApollo;
            torquesOfApollo[ "Earth" ].push_back( std::make_shared< TorqueSettings >( aerodynamic_torque ) );
            torqueMap[ "Apollo" ] = torquesOfApollo;


            bodiesToPropagate.push_back( "Apollo" );
            centralBodies.push_back( "Earth" );

            // Set initial state
            Eigen::Vector6d systemInitialState = apolloInitialState;
            Eigen::VectorXd systemInitialRotationalState = Eigen::VectorXd::Zero( 7 );
            systemInitialRotationalState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
            systemInitialRotationalState.segment( 4, 3 ) = Eigen::Vector3d::Constant( 1.0E-5 );

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
                std::make_shared< SingleDependentVariableSaveSettings >( airspeed_dependent_variable, "Apollo" ) );
            dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >( local_density_dependent_variable, "Apollo" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            aerodynamic_moment_coefficients_dependent_variable, "Apollo" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            aerodynamic_force_coefficients_dependent_variable, "Apollo" ) );
            dependentVariables.push_back(
                std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                    "Apollo", reference_frames::body_frame, reference_frames::inertial_frame, "Earth" ) );
            dependentVariables.push_back(
                std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                    "Apollo", reference_frames::aerodynamic_frame, reference_frames::inertial_frame, "Earth" ) );
            dependentVariables.push_back(
                std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                    "Apollo", reference_frames::aerodynamic_frame, reference_frames::body_frame, "Earth" ) );
            dependentVariables.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic, "Apollo", "Earth", 0 ) );
            dependentVariables.push_back(
                std::make_shared< SingleTorqueDependentVariableSaveSettings >(
                    aerodynamic_torque, "Apollo", "Earth", 0 ) );
            if( propagationType == 1 )
            {
                dependentVariables.push_back(
                    std::make_shared<SingleDependentVariableSaveSettings>(
                        body_center_of_mass, "Apollo" ));
                dependentVariables.push_back(
                    std::make_shared<SingleDependentVariableSaveSettings>(
                        body_inertia_tensor, "Apollo" ));
            }

            // Create acceleration models and propagation settings.
            basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodies, accelerationMap, bodiesToPropagate, centralBodies );
            basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                bodies, torqueMap, bodiesToPropagate );

            std::shared_ptr< IntegratorSettings< > > integratorSettings =
                    std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, simulationStartEpoch, fixedStepSize );

            auto terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch );
            std::shared_ptr< TranslationalStatePropagatorSettings < double > > translationalPropagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationStartEpoch, integratorSettings,
                      terminationSettings, cowell,
                      dependentVariables );

            std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
                std::make_shared< RotationalStatePropagatorSettings< double > >(
                    torqueModelMap, bodiesToPropagate, systemInitialRotationalState, simulationStartEpoch, integratorSettings,
                    terminationSettings, quaternions,
                    dependentVariables );

            std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > propagatorSettingsVector;
            propagatorSettingsVector.push_back( translationalPropagatorSettings );
            propagatorSettingsVector.push_back( rotationalPropagatorSettings );

            if( propagationType == 1 )
            {
                std::map<std::string, std::vector<std::shared_ptr<basic_astrodynamics::MassRateModel> > > massRateModels;
                massRateModels[ "Apollo" ].push_back(
                    createMassRateModel( "Apollo", std::make_shared< CustomMassRateSettings >( &customMassDerivativeFunction ),
                                         bodies, accelerationModelMap ));

                std::shared_ptr<SingleArcPropagatorSettings<double> > massPropagatorSettings =
                    std::make_shared<MassPropagatorSettings<double> >(
                        bodiesToPropagate, massRateModels,
                        ( Eigen::Matrix<double, 1, 1>( ) << initialVehicleMass ).finished( ), simulationStartEpoch, integratorSettings,
                        terminationSettings,
                        dependentVariables );

                propagatorSettingsVector.push_back( massPropagatorSettings );
            }

            std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings =
                std::make_shared< MultiTypePropagatorSettings< double > >(
                    propagatorSettingsVector, integratorSettings, simulationStartEpoch,
                    terminationSettings,
                    dependentVariables );

//            if( propagationType == 1  && i == 0 )
            {
                propagatorSettings->getPrintSettings()->setPrintDependentVariableData( true );
            }

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodies, propagatorSettings );

            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > stateOutput =
                dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > rawStateOutput =
                dynamicsSimulator.getEquationsOfMotionNumericalSolutionRaw( );
            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput =
                    dynamicsSimulator.getDependentVariableHistory( );

            std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > aerodynamicCoefficientInterface =
                    bodies.at( "Apollo" )->getAerodynamicCoefficientInterface( );

            for( auto it : dependentVariableOutput )
            {
                // Extract data from dependent variables
                double mass = propagationType == 0 ? 5000.0 : stateOutput.at( it.first )( 13 );
                double machNumber = it.second( 0 );
                double angleOfAttack = it.second( 1 );
                double sideslipAngle = it.second( 2 );
                double bankAngle = it.second( 3 );
                double airspeed = it.second( 4 );
                double density = it.second( 5 );
                Eigen::Vector3d momentCoefficients = it.second.segment( 6, 3 );
                Eigen::Vector3d forceCoefficients = it.second.segment( 9, 3 );
                Eigen::Matrix3d rotationFromBodyFrame = getMatrixFromVectorRotationRepresentation(
                    it.second.segment( 12, 9 ) );
                Eigen::Matrix3d rotationFromAerodynamicFrame = getMatrixFromVectorRotationRepresentation(
                    it.second.segment( 21, 9 ) );
                Eigen::Matrix3d rotationFromAerodynamicToBodyFrame = getMatrixFromVectorRotationRepresentation(
                    it.second.segment( 30, 9 ) );
                Eigen::Vector3d aerodynamicForce = it.second.segment( 39, 3 ) * mass;
                Eigen::Vector3d aerodynamicMoment = it.second.segment( 42, 3 );

                // Update environment to current state for comparison
                dynamicsSimulator.getDynamicsStateDerivative( )->computeStateDerivative(
                    it.first, rawStateOutput.at( it.first ) );
                aerodynamicCoefficientInterface->updateFullCurrentCoefficients(
                { machNumber, angleOfAttack, sideslipAngle } );

                // Extract aerodynamic coefficients
                Eigen::Vector3d testMomentCoefficients = aerodynamicCoefficientInterface->getCurrentMomentCoefficients( );
                Eigen::Vector3d testForceCoefficients = aerodynamicCoefficientInterface->getCurrentForceCoefficients( );

                // Check aerodynamic coefficients
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( momentCoefficients, testMomentCoefficients, std::numeric_limits< double >::epsilon( ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( forceCoefficients, testForceCoefficients, std::numeric_limits< double >::epsilon( ) );

                // Set relevant rotation matrices
                Eigen::Matrix3d bodyToMomentFrame = rotationFromAerodynamicToBodyFrame.transpose( );
                if( i == 2 )
                {
                    bodyToMomentFrame = Eigen::Matrix3d::Identity( );
                }

                Eigen::Matrix3d inertialToForceFrame = rotationFromAerodynamicFrame;
                if( i == 1 )
                {
                    inertialToForceFrame = rotationFromBodyFrame;
                }
                else if( i == 2 )
                {
                    inertialToForceFrame = -rotationFromAerodynamicFrame;
                }

                Eigen::Matrix3d forceToMomentFrame = Eigen::Matrix3d::Identity( );
                if( i == 1 )
                {
                    forceToMomentFrame = rotationFromAerodynamicToBodyFrame.transpose( );
                }
                else if( i == 2 )
                {
                    forceToMomentFrame = -rotationFromAerodynamicToBodyFrame;
                }


                // Convert aerodynamic coefficients to frame used in propagation
                Eigen::Vector3d bodyFixedMomentCoefficients =
                    bodyToMomentFrame.transpose( )  * momentCoefficients;
                Eigen::Vector3d inertialForceCoefficients =
                    inertialToForceFrame * forceCoefficients;

                // Manually compute acceleration and torque
                Eigen::Vector3d testAerodynamicForce =
                    0.5 * density * airspeed * airspeed * aerodynamicCoefficients->getReferenceArea( ) * inertialForceCoefficients;
                Eigen::Vector3d testAerodynamicMoment =
                    0.5 * density * airspeed * airspeed *
                    aerodynamicCoefficients->getReferenceArea( ) * aerodynamicCoefficients->getReferenceLength( ) *
                    bodyFixedMomentCoefficients;

                // Compare manual and saved force and torque
                for( unsigned int index = 0; index < 3; index++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( aerodynamicForce( index ) - testAerodynamicForce( index ) ),
                                       10.0 * std::numeric_limits< double >::epsilon( ) * aerodynamicForce.norm( ) );
                    BOOST_CHECK_SMALL( std::fabs( aerodynamicMoment( index ) - testAerodynamicMoment( index ) ),
                                       10.0 * std::numeric_limits< double >::epsilon( ) * aerodynamicForce.norm( ) );
                }

                if( propagationType == 1  )
                {
                    Eigen::Vector3d centerOfMass = it.second.segment( 45, 3 );

                    // Get force contribution to moment calculation
                    aerodynamicCoefficientInterface->updateFullCurrentCoefficients(
                        { machNumber, angleOfAttack,sideslipAngle },
                        std::map< std::string, std::vector< double > > ( ), TUDAT_NAN, false );
                    Eigen::Vector3d testMomentCoefficientsWithoutForceContribution =
                        aerodynamicCoefficientInterface->getCurrentMomentCoefficients( );
                    Eigen::Vector3d directMomentContribution = testMomentCoefficients - testMomentCoefficientsWithoutForceContribution;

                    // Manually compute force contribution to moment calculation
                    Eigen::Vector3d testForceContribution =
                        ( ( bodyToMomentFrame * ( aerodynamicCoefficients->getMomentReferencePoint( ) - centerOfMass ) )
                        ).cross( forceToMomentFrame * forceCoefficients ) / aerodynamicCoefficients->getReferenceLength( );

                    // Check force contribution to moment calculation
                    for( unsigned int index = 0; index < 3; index++ )
                    {
                        BOOST_CHECK_SMALL(
                                std::fabs( directMomentContribution( index ) - testForceContribution( index )),
                                20.0 * std::numeric_limits<double>::epsilon( ) * testForceContribution.norm( ));
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
