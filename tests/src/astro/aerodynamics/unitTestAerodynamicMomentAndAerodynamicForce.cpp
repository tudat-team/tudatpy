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

#include <boost/make_shared.hpp>
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
//#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/body.h"
//#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_aerodynamic_acceleration_force_moment_models )

using namespace aerodynamics;
using namespace simulation_setup;

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
                    referenceLength, referenceArea, referenceLength, Eigen::Vector3d::Zero( ) );

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
                    referenceLength, referenceArea, referenceLength, Eigen::Vector3d::Zero( ) );

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
            referenceLength, referenceArea, referenceLength, Eigen::Vector3d::Zero( ) );

        // Compute aerodynamic moment using free function with coefficient interface argument.
        Eigen::Vector3d moment = computeAerodynamicMoment( dynamicPressure,
                                                           aerodynamicCoefficientInterface );

        // Check if computed moment matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMoment, moment, tolerance );
    }
}

class DummyAngleCalculator: public AerodynamicGuidance
{
public:

    void updateGuidance( const double time )
    {
        time_ = time;

        currentAngleOfAttack_ = getDummyAngleOfAttack( );
        currentAngleOfSideslip_ = getDummyAngleOfSideslip( );
        currentBankAngle_ = getDummyBankAngle( );

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
                                    const bool imposeThrustDirection,
                                    const bool swapCreationOrder,
                                    const bool setAngleGuidanceManually )
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

    unsigned int maximumIndex = 8;
    if( imposeThrustDirection )
    {
        maximumIndex = 4;
    }
    for( unsigned int i = 0; i < maximumIndex; i++ )
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
        if( i < 4 && !imposeThrustDirection )
        {
            bodies.at( "Vehicle" )->setRotationalEphemeris(
                        std::make_shared< ephemerides::SpiceRotationalEphemeris >( "ECLIPJ2000", "IAU_MOON" ) );
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
                    2.0, 4.0, 1.5, Eigen::Vector3d::Zero( ), aerodynamicCoefficients, Eigen::Vector3d::Zero( ),
                    areCoefficientsInAerodynamicFrame, 1 );
        bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );
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
            if( !imposeThrustDirection )
            {
                accelerationsOfVehicle[ "Vehicle" ].push_back(
                            std::make_shared< ThrustAccelerationSettings >(
                                std::make_shared< ThrustDirectionSettings >(
                                    thrust_direction_from_existing_body_orientation, "Earth" ),
                                std::make_shared< ConstantThrustMagnitudeSettings >(
                                    thrustMagnitude, specificImpulse, bodyFixedThrustDirection ) ) );
            }
            else
            {
                accelerationsOfVehicle[ "Vehicle" ].push_back(
                            std::make_shared< ThrustAccelerationSettings >(
                                std::make_shared< CustomThrustOrientationSettings >(
                                    std::bind( spice_interface::computeRotationQuaternionBetweenFrames,
                                                     "IAU_Mars", "IAU_Earth", std::placeholders::_1 ) ),
                                std::make_shared< ConstantThrustMagnitudeSettings >(
                                    thrustMagnitude, specificImpulse, bodyFixedThrustDirection ) ) );
            }
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

        std::shared_ptr< DummyAngleCalculator > testAngles =
                std::make_shared< DummyAngleCalculator >( );

        if( !( i < 4 ) )
        {
            if( !setAngleGuidanceManually )
            {
                setGuidanceAnglesFunctions( testAngles, bodies.at( "Vehicle" ) );
            }
            else
            {
                std::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator =
                        bodies.at( "Vehicle" )->getFlightConditions( )->getAerodynamicAngleCalculator( );
                angleCalculator->setOrientationAngleFunctions(
                            std::bind( &DummyAngleCalculator::getDummyAngleOfAttack, testAngles ),
                            std::bind( &DummyAngleCalculator::getDummyAngleOfSideslip, testAngles ),
                            std::bind( &DummyAngleCalculator::getDummyBankAngle, testAngles ),
                            std::bind( &DummyAngleCalculator::updateGuidance, testAngles, std::placeholders::_1 ) );
            }
        }


        std::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings;

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

        dependentVariableSaveSettings = std::make_shared< DependentVariableSaveSettings >( dependentVariables );


        std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                std::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings,
                  cowell, dependentVariableSaveSettings );
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, 5.0 );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, translationalPropagatorSettings, true, false, false );

        // Retrieve numerical solutions for state and dependent variables
        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput =
                dynamicsSimulator.getDependentVariableHistory( );

        std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris =
                bodies.at( "Vehicle" )->getRotationalEphemeris( );

        double thrustAcceleration = thrustMagnitude / vehicleMass;
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
            if( includeThrustForce && !imposeThrustDirection )
            {
                Eigen::Vector3d thrustForceInBodyFrame = rotationToBodyFrame * outputIterator->second.segment( 42, 3 );

                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( thrustForceInBodyFrame( j ) - bodyFixedThrustDirection( j ) *
                                           thrustAcceleration ), 1.0E-14 );
                }
            }
            else if( includeThrustForce && imposeThrustDirection )
            {
                Eigen::Vector3d thrustForceInPropagationFrame = ( outputIterator->second.segment( 42, 3 ) );

                Eigen::Matrix3d rotationToBodyFrameFromEphemeris = spice_interface::computeRotationQuaternionBetweenFrames(
                            "IAU_Earth", "IAU_Mars", outputIterator->first ).toRotationMatrix( );
                Eigen::Vector3d imposedThrustForceInPropagationFrame =
                        thrustAcceleration * ( rotationToBodyFrameFromEphemeris.transpose( ) * bodyFixedThrustDirection );



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
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( j, k ) ), 5.0E-14 );
                    }
                }

            }

            if( i < 4 && !imposeThrustDirection )
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
            else if( !( i < 4 ) )
            {
                testAngles->updateGuidance( outputIterator->first );

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
    testAerodynamicForceDirection( 0, 0, 0, 0 );
    testAerodynamicForceDirection( 1, 0, 0, 0 );
    testAerodynamicForceDirection( 1, 1, 0, 0 );
    testAerodynamicForceDirection( 1, 0, 1, 0 );
    testAerodynamicForceDirection( 1, 1, 1, 0 );

    testAerodynamicForceDirection( 0, 0, 0, 1 );
    testAerodynamicForceDirection( 1, 0, 0, 1 );
    testAerodynamicForceDirection( 1, 1, 0, 1 );
    testAerodynamicForceDirection( 1, 0, 1, 1 );
    testAerodynamicForceDirection( 1, 1, 1, 1 );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
