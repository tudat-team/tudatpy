/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
//
#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include <memory>

#include "tudat/astro/aerodynamics/testApolloCapsuleCoefficients.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/basics/testMacros.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/multiDimensionalArrayReader.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/propagation_setup/createMassRateModels.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/gravitation/polynomialGravityFieldVariations.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include <limits>
#include <string>

#include <Eigen/Core>

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_rigid_body_properties )

double dummyMassFunction( const double time )
{
    return 5.0E3 - 0.5 * ( time - 1.0E7 ) - 0.5 * std::sin( ( time - 1.0E7 ) / 300.0 );
}

Eigen::Vector3d dummyCenterOfMassFunction( const double time )
{
    return ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( ) -
            ( Eigen::Vector3d( ) << 0.1, 0.3, 0.5 ).finished( ) / 1000.0 * ( time - 1.0E7 );
}

Eigen::Matrix3d dummyInertiaTensorFunction( const double time )
{
    return 3.0 * Eigen::Matrix3d::Identity( ) -
            ( Eigen::Vector3d::Ones( ) * Eigen::Vector3d::Ones( ).transpose( ) ) / 1000.0 * ( time - 1.0E7 );
}

Eigen::Vector3d dummyCenterOfMassFunction2( const double mass )
{
    return ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( ) -
            ( Eigen::Vector3d( ) << 0.1, 0.3, 0.5 ).finished( ) * ( 5.0E3 - mass ) / 1000.0;
}

Eigen::Matrix3d dummyInertiaTensorFunction2( const double mass )
{
    return 3.0 * Eigen::Matrix3d::Identity( ) -
            ( Eigen::Vector3d::Ones( ) * Eigen::Vector3d::Ones( ).transpose( ) ) * ( 5.0E3 - mass ) / 1000.0;
}

//! Test inertia tensor derived from gravity field (if available), and that explicitly set rigid body properties are not replaced by
//! properties derived from the gravity field.
BOOST_AUTO_TEST_CASE( testGravityLinkedInertiaAvailabilityAndOwnership )
{
    using namespace gravitation;
    using namespace simulation_setup;

    std::shared_ptr< Body > body = std::make_shared< Body >( );
    body->setBodyName( "TestBody" );

    // Check that a point-mass gravity field defines no inertia tensor
    body->setGravityFieldModel( std::make_shared< GravityFieldModel >( 4.0E5 ) );
    BOOST_CHECK( !body->getMassProperties( )->isInertiaTensorAvailable( ) );

    // Check inertia tensor (and zero derivative) from static degree two field and scaled mean moment of inertia, against principal
    // moments computed from the unnormalized coefficients
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients( 0, 0 ) = 1.0;
    cosineCoefficients( 2, 0 ) = -1.0E-3;
    cosineCoefficients( 2, 2 ) = 2.0E-4;
    std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicField = std::make_shared< SphericalHarmonicsGravityField >(
            4.0E5, 2.0E3, cosineCoefficients, Eigen::MatrixXd::Zero( 3, 3 ), "BodyFixed" );
    body->setGravityFieldModel( sphericalHarmonicField );
    body->setMassProperties( std::make_shared< FromGravityFieldRigidBodyProperties >( sphericalHarmonicField, 0.4 ) );

    const double c20 = -1.0E-3 * std::sqrt( 5.0 );
    const double c22 = 2.0E-4 * std::sqrt( 5.0 / 12.0 );
    const Eigen::Vector3d scaledPrincipalMoments( 0.4 + c20 / 3.0 - 2.0 * c22, 0.4 + c20 / 3.0 + 2.0 * c22, 0.4 - 2.0 * c20 / 3.0 );
    const Eigen::Matrix3d expectedInertiaTensor =
            ( 4.0E5 / physical_constants::GRAVITATIONAL_CONSTANT ) * ( 2.0E3 * 2.0E3 ) * scaledPrincipalMoments.asDiagonal( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( body->getBodyInertiaTensor( ), expectedInertiaTensor, 5.0E-15 );
    BOOST_CHECK_SMALL( body->getBodyInertiaTensorDerivative( ).norm( ), 1.0E-30 );

    // Check that explicitly set rigid body properties are retained when resetting the gravity field
    const Eigen::Matrix3d explicitInertiaTensor = 7.0 * Eigen::Matrix3d::Identity( );
    std::shared_ptr< TimeDependentRigidBodyProperties > explicitProperties =
            std::make_shared< TimeDependentRigidBodyProperties >( 12.0, Eigen::Vector3d::Zero( ), explicitInertiaTensor );
    body->setMassProperties( explicitProperties );
    body->setGravityFieldModel( std::make_shared< GravityFieldModel >( 9.0E5 ) );
    BOOST_CHECK( body->getMassProperties( ) == explicitProperties );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( body->getBodyInertiaTensor( ), explicitInertiaTensor, 5.0E-15 );
}

//! Test creation of rigid body properties from deprecated (gravity field) settings, gravity-derived settings and explicit settings,
//! and the update of gravity-derived properties when the gravity field changes.
BOOST_AUTO_TEST_CASE( testGravityDerivedSettingsCompatibilityAndSynchronization )
{
    using namespace gravitation;
    using namespace simulation_setup;

    const double gravitationalParameter = 4.0e5;
    const double referenceRadius = 2.0e3;
    const double scaledMeanMomentOfInertia = 0.4;
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients( 0, 0 ) = 1.0;
    cosineCoefficients( 1, 0 ) = 1.0e-5;
    cosineCoefficients( 1, 1 ) = -2.0e-5;
    sineCoefficients( 1, 1 ) = 3.0e-5;
    cosineCoefficients( 2, 0 ) = -1.0e-3;
    cosineCoefficients( 2, 2 ) = 2.0e-4;
    sineCoefficients( 2, 1 ) = -3.0e-4;

    const auto makeGravitySettings = [ & ]( const double legacyScaledMeanMoment = TUDAT_NAN ) {
        return std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                gravitationalParameter, referenceRadius, cosineCoefficients, sineCoefficients, "BodyFixed", legacyScaledMeanMoment );
    };

    BodyListSettings settings;
    for( const std::string bodyName : { "Legacy", "FromGravityField", "NoInertia", "Explicit" } )
    {
        settings.addSettings( bodyName );
    }
    settings.at( "Legacy" )->gravityFieldSettings = makeGravitySettings( scaledMeanMomentOfInertia );
    settings.at( "FromGravityField" )->gravityFieldSettings = makeGravitySettings( );
    settings.at( "FromGravityField" )->rigidBodyPropertiesSettings =
            fromGravityFieldRigidBodyPropertiesSettings( scaledMeanMomentOfInertia );
    settings.at( "NoInertia" )->gravityFieldSettings = makeGravitySettings( );

    const double explicitMass = 123.0;
    const Eigen::Vector3d explicitCenterOfMass = ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( );
    const Eigen::Matrix3d explicitInertia = 7.0 * Eigen::Matrix3d::Identity( );
    settings.at( "Explicit" )->gravityFieldSettings = makeGravitySettings( scaledMeanMomentOfInertia );
    settings.at( "Explicit" )->rigidBodyPropertiesSettings =
            constantRigidBodyPropertiesSettings( explicitMass, explicitCenterOfMass, explicitInertia );

    SystemOfBodies bodies = createSystemOfBodies( settings );
    const std::shared_ptr< Body > legacyBody = bodies.at( "Legacy" );
    const std::shared_ptr< Body > fromGravityFieldBody = bodies.at( "FromGravityField" );
    const std::shared_ptr< Body > noInertiaBody = bodies.at( "NoInertia" );
    const std::shared_ptr< Body > explicitBody = bodies.at( "Explicit" );

    // Check that deprecated and gravity-derived settings result in identical rigid body properties
    const std::shared_ptr< FromGravityFieldRigidBodyProperties > legacyProperties =
            std::dynamic_pointer_cast< FromGravityFieldRigidBodyProperties >( legacyBody->getMassProperties( ) );
    const std::shared_ptr< FromGravityFieldRigidBodyProperties > fromGravityFieldProperties =
            std::dynamic_pointer_cast< FromGravityFieldRigidBodyProperties >( fromGravityFieldBody->getMassProperties( ) );
    BOOST_REQUIRE( legacyProperties != nullptr );
    BOOST_REQUIRE( fromGravityFieldProperties != nullptr );
    BOOST_CHECK_CLOSE_FRACTION( legacyProperties->getScaledMeanMomentOfInertia( ), scaledMeanMomentOfInertia, 5.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( fromGravityFieldProperties->getScaledMeanMomentOfInertia( ), scaledMeanMomentOfInertia, 5.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( legacyBody->getBodyMass( ), fromGravityFieldBody->getBodyMass( ), 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            legacyBody->getBodyFixedCenterOfMass( ), fromGravityFieldBody->getBodyFixedCenterOfMass( ), 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( legacyBody->getBodyInertiaTensor( ), fromGravityFieldBody->getBodyInertiaTensor( ), 5.0e-15 );

    // Check center of mass from degree one coefficients
    const Eigen::Vector3d expectedCenterOfMass =
            referenceRadius * std::sqrt( 3.0 ) * ( Eigen::Vector3d( ) << -2.0e-5, 3.0e-5, 1.0e-5 ).finished( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fromGravityFieldBody->getBodyFixedCenterOfMass( ), expectedCenterOfMass, 5.0e-15 );

    // Check that no inertia tensor is defined without scaled mean moment of inertia
    BOOST_REQUIRE( noInertiaBody->getMassProperties( ) != nullptr );
    BOOST_CHECK( !noInertiaBody->getMassProperties( )->isInertiaTensorAvailable( ) );

    // Check that explicit rigid body settings take precedence over deprecated gravity field settings
    BOOST_CHECK( std::dynamic_pointer_cast< FromGravityFieldRigidBodyProperties >( explicitBody->getMassProperties( ) ) == nullptr );
    BOOST_CHECK_CLOSE_FRACTION( explicitBody->getBodyMass( ), explicitMass, 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( explicitBody->getBodyFixedCenterOfMass( ), explicitCenterOfMass, 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( explicitBody->getBodyInertiaTensor( ), explicitInertia, 5.0e-15 );

    // Check update of mass and inertia tensor after changing gravitational parameter
    const double initialMass = fromGravityFieldBody->getBodyMass( );
    const Eigen::Matrix3d initialInertia = fromGravityFieldBody->getBodyInertiaTensor( );
    fromGravityFieldBody->getGravityFieldModel( )->resetGravitationalParameter( 2.0 * gravitationalParameter );
    fromGravityFieldBody->getMassProperties( )->update( 0.0 );
    BOOST_CHECK_CLOSE_FRACTION( fromGravityFieldBody->getBodyMass( ), 2.0 * initialMass, 5.0e-15 );
    const Eigen::Matrix3d expectedInertiaAfterMassUpdate = 2.0 * initialInertia;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fromGravityFieldBody->getBodyInertiaTensor( ), expectedInertiaAfterMassUpdate, 5.0e-15 );

    const std::shared_ptr< SphericalHarmonicsGravityField > fromGravityFieldGravityField =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( fromGravityFieldBody->getGravityFieldModel( ) );
    BOOST_REQUIRE( fromGravityFieldGravityField != nullptr );
    Eigen::MatrixXd modifiedCosineCoefficients = cosineCoefficients;
    // Check update of center of mass and inertia tensor after changing degree one and two coefficients
    modifiedCosineCoefficients( 1, 1 ) *= 2.0;
    modifiedCosineCoefficients( 2, 0 ) *= 1.5;
    fromGravityFieldGravityField->setCosineCoefficients( modifiedCosineCoefficients );
    fromGravityFieldBody->getMassProperties( )->update( 0.0 );
    Eigen::Vector3d expectedUpdatedCenterOfMass = expectedCenterOfMass;
    expectedUpdatedCenterOfMass.x( ) *= 2.0;
    const double inertiaChangeScale = fromGravityFieldBody->getBodyMass( ) * referenceRadius * referenceRadius * std::sqrt( 5.0 ) *
            ( modifiedCosineCoefficients( 2, 0 ) - cosineCoefficients( 2, 0 ) ) / 3.0;
    const Eigen::Matrix3d expectedCoefficientInertiaChange = inertiaChangeScale * Eigen::Vector3d( 1.0, 1.0, -2.0 ).asDiagonal( );
    const Eigen::Matrix3d expectedUpdatedInertia = expectedInertiaAfterMassUpdate + expectedCoefficientInertiaChange;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fromGravityFieldBody->getBodyFixedCenterOfMass( ), expectedUpdatedCenterOfMass, 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fromGravityFieldBody->getBodyInertiaTensor( ), expectedUpdatedInertia, 5.0e-15 );

    // Check that explicit rigid body properties are not changed by changing the gravity field
    explicitBody->getGravityFieldModel( )->resetGravitationalParameter( 3.0 * gravitationalParameter );
    explicitBody->getMassProperties( )->update( 0.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( explicitBody->getBodyInertiaTensor( ), explicitInertia, 5.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( explicitBody->getBodyMass( ), explicitMass, 5.0e-15 );
}

//! Test inertia tensor derivative from linearly varying degree two coefficients against numerical derivatives, and its contribution to
//! the inertial torque.
BOOST_AUTO_TEST_CASE( testPrescribedGravityVariationUpdatesInertiaAndDerivative )
{
    using namespace gravitation;
    using namespace simulation_setup;

    Eigen::MatrixXd nominalCosine = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd nominalSine = Eigen::MatrixXd::Zero( 3, 3 );
    nominalCosine( 0, 0 ) = 1.0;
    nominalCosine( 2, 0 ) = -1.0e-3;
    nominalCosine( 2, 2 ) = 2.0e-4;

    Eigen::MatrixXd cosineRate = Eigen::MatrixXd::Zero( 1, 3 );
    Eigen::MatrixXd sineRate = Eigen::MatrixXd::Zero( 1, 3 );
    cosineRate << 2.0e-8, -3.0e-8, 5.0e-8;
    sineRate << 0.0, 7.0e-8, -11.0e-8;
    const std::shared_ptr< PolynomialGravityFieldVariations > prescribedVariation = std::make_shared< PolynomialGravityFieldVariations >(
            std::map< int, Eigen::MatrixXd >( { { 1, cosineRate } } ), std::map< int, Eigen::MatrixXd >( { { 1, sineRate } } ), 0.0, 2, 0 );
    const std::shared_ptr< GravityFieldVariationsSet > variationSet = std::make_shared< GravityFieldVariationsSet >(
            std::vector< std::shared_ptr< GravityFieldVariations > >( { prescribedVariation } ),
            std::vector< BodyDeformationTypes >( { polynomial_variation } ),
            std::vector< std::string >( { "" } ) );
    const std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > gravityField =
            std::make_shared< TimeDependentSphericalHarmonicsGravityField >(
                    4.0e5, 2.0e3, nominalCosine, nominalSine, variationSet, "BodyFixed" );

    const std::shared_ptr< Body > body = std::make_shared< Body >( );
    body->setBodyName( "TestBody" );
    body->setGravityFieldModel( gravityField );
    body->setMassProperties( std::make_shared< FromGravityFieldRigidBodyProperties >( gravityField, 0.4 ) );
    body->setGravityFieldVariationSet( variationSet );
    body->setIsBodyInPropagation( true );

    BOOST_CHECK( body->getMassProperties( )->isInertiaTensorDerivativeAvailable( ) );

    const Eigen::Vector3d angularVelocity = ( Eigen::Vector3d( ) << 1.0e-4, -2.0e-4, 3.0e-4 ).finished( );
    Eigen::Vector7d rotationalState = Eigen::Vector7d::Zero( );
    rotationalState( 0 ) = 1.0;
    rotationalState.tail( 3 ) = angularVelocity;
    body->setCurrentRotationalStateToLocalFrame( rotationalState );
    const auto inertialTorque = createInertialTorqueModel( body, "TestBody" );

    // Update gravity field, and rigid body properties derived from it, to given time
    std::function< void( const double ) > updateEnvironment = [ = ]( const double time ) {
        gravityField->update( time );
        body->getMassProperties( )->resetCurrentTime( );
        body->getMassProperties( )->update( time );
    };

    for( const double time : { 0.0, 10.0, 40.0 } )
    {
        const double step = 1.0;
        updateEnvironment( time - step );
        const Eigen::Matrix3d previousInertia = body->getBodyInertiaTensor( );
        updateEnvironment( time + step );
        const Eigen::Matrix3d nextInertia = body->getBodyInertiaTensor( );
        const Eigen::Matrix3d numericalDerivative = ( nextInertia - previousInertia ) / ( 2.0 * step );
        updateEnvironment( time );
        const Eigen::Matrix3d inertiaDerivative = body->getBodyInertiaTensorDerivative( );

        // Compare inertia tensor derivative with numerical derivative
        BOOST_REQUIRE_GT( inertiaDerivative.norm( ), 0.0 );
        BOOST_CHECK_SMALL( ( inertiaDerivative - numericalDerivative ).norm( ) / inertiaDerivative.norm( ), 1.0e-9 );

        inertialTorque->updateMembers( time );
        const Eigen::Vector3d derivativeTorque =
                inertialTorque->getTorque( ) + angularVelocity.cross( body->getBodyInertiaTensor( ) * angularVelocity );
        const Eigen::Vector3d expectedDerivativeTorque = -inertiaDerivative * angularVelocity;
        // Check that the inertial torque includes -dI/dt * omega
        BOOST_REQUIRE_GT( expectedDerivativeTorque.norm( ), 0.0 );
        BOOST_CHECK_SMALL( ( derivativeTorque - expectedDerivativeTorque ).norm( ) / expectedDerivativeTorque.norm( ), 1.0e-13 );
    }
}

BOOST_AUTO_TEST_CASE( testDirectRigidBodyProperties )
{
    using namespace tudat;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;
    using namespace orbital_element_conversions;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    for( unsigned int propagateMass = 0; propagateMass < 2; propagateMass++ )
    {
        unsigned int maximumNumberOfRuns = 2;
        if( propagateMass == 1 )
        {
            maximumNumberOfRuns = 3;
        }
        for( unsigned int i = 0; i < maximumNumberOfRuns; i++ )
        {
            double initialTime = 1.0E7;
            // Create vehicle objects.
            double vehicleMass = 5.0E3;
            double dryVehicleMass = 2.0E3;

            // Define simulation body settings.
            BodyListSettings bodySettings = getDefaultBodySettings( { "Earth", "Moon", "Sun" }, "Earth", "ECLIPJ2000" );

            // Create Earth object
            simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );

            bodies.createEmptyBody( "Vehicle" );

            std::shared_ptr< RigidBodyPropertiesSettings > rigidBodyProperties;
            if( i == 0 )
            {
                rigidBodyProperties = std::make_shared< ConstantRigidBodyPropertiesSettings >(
                        vehicleMass, ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( ), 3.0 * Eigen::Matrix3d::Identity( ) );
            }
            else if( i == 1 )
            {
                rigidBodyProperties = std::make_shared< FromFunctionRigidBodyPropertiesSettings >(
                        &dummyMassFunction, &dummyCenterOfMassFunction, &dummyInertiaTensorFunction );
            }
            else if( i == 2 )
            {
                rigidBodyProperties = std::make_shared< MassDependentMassDistributionSettings >(
                        vehicleMass, &dummyCenterOfMassFunction2, &dummyInertiaTensorFunction2 );
            }
            addRigidBodyProperties( bodies, "Vehicle", rigidBodyProperties );

            double thrustMagnitude1 = 1.0E3;
            double specificImpulse1 = 250.0;
            double massFlow1 = propulsion::computePropellantMassRateFromSpecificImpulse( thrustMagnitude1, specificImpulse1 );

            addEngineModel( "Vehicle",
                            "Engine",
                            std::make_shared< ConstantThrustMagnitudeSettings >( thrustMagnitude1, specificImpulse1 ),
                            bodies );

            bodies.at( "Vehicle" )
                    ->setRotationalEphemeris( createRotationModel(
                            orbitalStateBasedRotationSettings( "Earth", true, false, "ECLIPJ2000", "VehicleFixed" ), "Vehicle", bodies ) );

            // Define propagator settings variables.
            SelectedAccelerationMap accelerationMap;
            std::vector< std::string > bodiesToPropagate;
            std::vector< std::string > centralBodies;

            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
            accelerationsOfVehicle[ "Earth" ].push_back(
                    std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
            accelerationsOfVehicle[ "Moon" ].push_back(
                    std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
            accelerationsOfVehicle[ "Sun" ].push_back(
                    std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
            accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >( "Engine" ) );

            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

            bodiesToPropagate.push_back( "Vehicle" );
            centralBodies.push_back( "Earth" );

            // Create acceleration models and propagation settings.
            basic_astrodynamics::AccelerationMap accelerationModelMap =
                    createAccelerationModelsMap( bodies, accelerationMap, bodiesToPropagate, centralBodies );

            // Set Keplerian elements for Vehicle.
            Eigen::Vector6d vehicleInitialStateInKeplerianElements;
            vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 8000.0E3;
            vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
            vehicleInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
            vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
            vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
            vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

            double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
            const Eigen::Vector6d vehicleInitialState =
                    convertKeplerianToCartesianElements( vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave;
            dependentVariablesToSave.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >( current_body_mass_dependent_variable, "Vehicle" ) );
            dependentVariablesToSave.push_back( std::make_shared< SingleDependentVariableSaveSettings >( body_center_of_mass, "Vehicle" ) );
            dependentVariablesToSave.push_back( std::make_shared< SingleDependentVariableSaveSettings >( body_inertia_tensor, "Vehicle" ) );

            std::shared_ptr< IntegratorSettings<> > integratorSettings = std::make_shared< IntegratorSettings<> >( rungeKutta4, 1.0 );
            std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                    std::make_shared< propagators::PropagationTimeTerminationSettings >( initialTime + 1000.0 );
            std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >( centralBodies,
                                                                                        accelerationModelMap,
                                                                                        bodiesToPropagate,
                                                                                        vehicleInitialState,
                                                                                        initialTime,
                                                                                        integratorSettings,
                                                                                        terminationSettings,
                                                                                        cowell,
                                                                                        dependentVariablesToSave );

            std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > > > massRateModels;
            massRateModels[ "Vehicle" ].push_back(
                    createMassRateModel( "Vehicle", std::make_shared< FromThrustMassRateSettings >( 1 ), bodies, accelerationModelMap ) );

            std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings =
                    std::make_shared< MassPropagatorSettings< double > >( std::vector< std::string >{ "Vehicle" },
                                                                          massRateModels,
                                                                          ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ),
                                                                          initialTime,
                                                                          integratorSettings,
                                                                          terminationSettings,
                                                                          dependentVariablesToSave );

            std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
            propagatorSettingsVector.push_back( translationalPropagatorSettings );
            if( propagateMass == 1 )
            {
                propagatorSettingsVector.push_back( massPropagatorSettings );
            }

            std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
                    std::make_shared< MultiTypePropagatorSettings< double > >(
                            propagatorSettingsVector, integratorSettings, initialTime, terminationSettings, dependentVariablesToSave );

            // Create simulation object and propagate dynamics.
            auto dynamicsSimulator = std::dynamic_pointer_cast< SingleArcDynamicsSimulator<> >(
                    createDynamicsSimulator< double, double >( bodies, propagatorSettings ) );

            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableSolution =
                    dynamicsSimulator->getDependentVariableHistory( );
            for( auto it : dependentVariableSolution )
            {
                double currentMass = it.second( 0 );
                Eigen::Vector3d currentCenterOfMass = it.second.segment( 1, 3 );
                Eigen::Matrix3d currentInertiaTensor = getMatrixFromVectorRotationRepresentation( it.second.segment( 4, 9 ) );
                if( propagateMass == 0 )
                {
                    if( i == 0 )
                    {
                        BOOST_CHECK_CLOSE_FRACTION( currentMass, vehicleMass, std::numeric_limits< double >::epsilon( ) );
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentCenterOfMass,
                                                           ( ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( ) ),
                                                           std::numeric_limits< double >::epsilon( ) );
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                currentInertiaTensor, ( 3.0 * Eigen::Matrix3d::Identity( ) ), std::numeric_limits< double >::epsilon( ) );
                    }
                    else if( i == 1 )
                    {
                        BOOST_CHECK_CLOSE_FRACTION( currentMass, dummyMassFunction( it.first ), std::numeric_limits< double >::epsilon( ) );
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                currentCenterOfMass, ( dummyCenterOfMassFunction( it.first ) ), std::numeric_limits< double >::epsilon( ) );
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentInertiaTensor,
                                                           ( dummyInertiaTensorFunction( it.first ) ),
                                                           std::numeric_limits< double >::epsilon( ) );
                    }
                }
                else if( propagateMass == 1 )
                {
                    BOOST_CHECK_CLOSE_FRACTION( currentMass,
                                                vehicleMass - ( it.first - initialTime ) * massFlow1,
                                                1.0E3 * std::numeric_limits< double >::epsilon( ) );
                    if( i == 0 )
                    {
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentCenterOfMass,
                                                           ( ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( ) ),
                                                           std::numeric_limits< double >::epsilon( ) );
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                currentInertiaTensor, ( 3.0 * Eigen::Matrix3d::Identity( ) ), std::numeric_limits< double >::epsilon( ) );
                    }
                    else if( i == 1 )
                    {
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                currentCenterOfMass, ( dummyCenterOfMassFunction( it.first ) ), std::numeric_limits< double >::epsilon( ) );
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentInertiaTensor,
                                                           ( dummyInertiaTensorFunction( it.first ) ),
                                                           std::numeric_limits< double >::epsilon( ) );
                    }
                    else if( i == 2 )
                    {
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentCenterOfMass,
                                                           ( dummyCenterOfMassFunction2( currentMass ) ),
                                                           std::numeric_limits< double >::epsilon( ) );
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentInertiaTensor,
                                                           ( dummyInertiaTensorFunction2( currentMass ) ),
                                                           std::numeric_limits< double >::epsilon( ) );
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
