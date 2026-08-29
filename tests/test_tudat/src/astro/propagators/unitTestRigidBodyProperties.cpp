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

BOOST_AUTO_TEST_CASE( testGravityLinkedInertiaAvailabilityAndOwnership )
{
    using namespace gravitation;
    using namespace simulation_setup;

    const std::shared_ptr< Body > body = std::make_shared< Body >( );
    body->setBodyName( "TestBody" );

    // A gravity model alone is sufficient for mass, but not necessarily for inertia.
    body->setGravityFieldModel( std::make_shared< GravityFieldModel >( 4.0e5 ) );
    BOOST_REQUIRE( body->getMassProperties( ) != nullptr );
    BOOST_CHECK( !body->getMassProperties( )->isInertiaTensorAvailable( ) );
    BOOST_CHECK( !body->getMassProperties( )->isInertiaTensorDerivativeAvailable( ) );
    BOOST_CHECK_THROW( body->getBodyInertiaTensor( ), std::runtime_error );

    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients( 0, 0 ) = 1.0;
    cosineCoefficients( 2, 0 ) = -1.0e-3;
    cosineCoefficients( 2, 2 ) = 2.0e-4;
    const std::shared_ptr< SphericalHarmonicsGravityField > sphericalField =
            std::make_shared< SphericalHarmonicsGravityField >( 4.0e5, 2.0e3, cosineCoefficients, sineCoefficients, "BodyFixed" );
    body->setGravityFieldModel( sphericalField );
    body->setMassProperties( std::make_shared< FromGravityFieldRigidBodyProperties >( sphericalField, 0.4 ) );
    BOOST_CHECK( body->getMassProperties( )->isInertiaTensorAvailable( ) );
    BOOST_CHECK( body->getMassProperties( )->isInertiaTensorDerivativeAvailable( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( body->getBodyInertiaTensor( ), getInertiaTensorFromGravityField( sphericalField, 0.4 ), 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( body->getBodyInertiaTensorDerivative( ), Eigen::Matrix3d::Zero( ), 5.0e-15 );

    // An explicitly configured rigid-body object owns its tensor and is not replaced with the gravity field.
    const Eigen::Matrix3d explicitInertia = 7.0 * Eigen::Matrix3d::Identity( );
    const std::shared_ptr< TimeDependentRigidBodyProperties > explicitProperties =
            std::make_shared< TimeDependentRigidBodyProperties >( 12.0, Eigen::Vector3d::Zero( ), explicitInertia );
    body->setMassProperties( explicitProperties );
    body->setGravityFieldModel( std::make_shared< GravityFieldModel >( 9.0e5 ) );
    BOOST_CHECK( body->getMassProperties( ) == explicitProperties );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( body->getBodyInertiaTensor( ), explicitInertia, 5.0e-15 );
}

BOOST_AUTO_TEST_CASE( testGravityDerivedSettingsCompatibilityAndCallbacks )
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
    for( const std::string bodyName : { "Legacy", "Canonical", "NoInertia", "Explicit" } )
    {
        settings.addSettings( bodyName );
    }
    settings.at( "Legacy" )->gravityFieldSettings = makeGravitySettings( scaledMeanMomentOfInertia );
    settings.at( "Canonical" )->gravityFieldSettings = makeGravitySettings( );
    settings.at( "Canonical" )->rigidBodyPropertiesSettings = fromGravityFieldRigidBodyPropertiesSettings( scaledMeanMomentOfInertia );
    settings.at( "NoInertia" )->gravityFieldSettings = makeGravitySettings( );

    const double explicitMass = 123.0;
    const Eigen::Vector3d explicitCenterOfMass = ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( );
    const Eigen::Matrix3d explicitInertia = 7.0 * Eigen::Matrix3d::Identity( );
    settings.at( "Explicit" )->gravityFieldSettings = makeGravitySettings( scaledMeanMomentOfInertia );
    settings.at( "Explicit" )->rigidBodyPropertiesSettings =
            constantRigidBodyPropertiesSettings( explicitMass, explicitCenterOfMass, explicitInertia );

    SystemOfBodies bodies = createSystemOfBodies( settings );
    const std::shared_ptr< Body > legacyBody = bodies.at( "Legacy" );
    const std::shared_ptr< Body > canonicalBody = bodies.at( "Canonical" );
    const std::shared_ptr< Body > noInertiaBody = bodies.at( "NoInertia" );
    const std::shared_ptr< Body > explicitBody = bodies.at( "Explicit" );

    // The deprecated settings carrier and the canonical rigid-body settings create identical
    // runtime environments, with a single scaled-mean value owned by rigid-body properties.
    const std::shared_ptr< FromGravityFieldRigidBodyProperties > legacyProperties =
            std::dynamic_pointer_cast< FromGravityFieldRigidBodyProperties >( legacyBody->getMassProperties( ) );
    const std::shared_ptr< FromGravityFieldRigidBodyProperties > canonicalProperties =
            std::dynamic_pointer_cast< FromGravityFieldRigidBodyProperties >( canonicalBody->getMassProperties( ) );
    BOOST_REQUIRE( legacyProperties != nullptr );
    BOOST_REQUIRE( canonicalProperties != nullptr );
    BOOST_CHECK_EQUAL( legacyBody->getGravityFieldModel( )->getRigidBodyProperties( ), legacyBody->getMassProperties( ) );
    BOOST_CHECK_EQUAL( canonicalBody->getGravityFieldModel( )->getRigidBodyProperties( ), canonicalBody->getMassProperties( ) );
    BOOST_CHECK_CLOSE_FRACTION( legacyProperties->getScaledMeanMomentOfInertia( ), scaledMeanMomentOfInertia, 5.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( canonicalProperties->getScaledMeanMomentOfInertia( ), scaledMeanMomentOfInertia, 5.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( legacyBody->getBodyMass( ), canonicalBody->getBodyMass( ), 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( legacyBody->getBodyFixedCenterOfMass( ), canonicalBody->getBodyFixedCenterOfMass( ), 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( legacyBody->getBodyInertiaTensor( ), canonicalBody->getBodyInertiaTensor( ), 5.0e-15 );

    // Complete gravity coefficients do not imply an inertia tensor without a finite scaled mean.
    BOOST_REQUIRE( noInertiaBody->getMassProperties( ) != nullptr );
    BOOST_CHECK( !noInertiaBody->getMassProperties( )->isInertiaTensorAvailable( ) );
    BOOST_CHECK_THROW( noInertiaBody->getBodyInertiaTensor( ), std::runtime_error );

    // Explicit non-gravity rigid-body settings retain precedence over the compatibility carrier.
    BOOST_CHECK( std::dynamic_pointer_cast< FromGravityFieldRigidBodyProperties >( explicitBody->getMassProperties( ) ) == nullptr );
    BOOST_CHECK_EQUAL( explicitBody->getGravityFieldModel( )->getRigidBodyProperties( ), explicitBody->getMassProperties( ) );
    BOOST_CHECK_CLOSE_FRACTION( explicitBody->getBodyMass( ), explicitMass, 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( explicitBody->getBodyFixedCenterOfMass( ), explicitCenterOfMass, 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( explicitBody->getBodyInertiaTensor( ), explicitInertia, 5.0e-15 );

    // Runtime changes propagate from gravity data to gravity-derived properties. A mu change
    // updates both mass and the inertia normalization, while coefficient changes update COM and inertia.
    const double initialMass = canonicalBody->getBodyMass( );
    const Eigen::Matrix3d initialInertia = canonicalBody->getBodyInertiaTensor( );
    canonicalBody->getGravityFieldModel( )->resetGravitationalParameter( 2.0 * gravitationalParameter );
    BOOST_CHECK_CLOSE_FRACTION( canonicalBody->getBodyMass( ), 2.0 * initialMass, 5.0e-15 );
    const Eigen::Matrix3d expectedInertiaAfterMassUpdate = 2.0 * initialInertia;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( canonicalBody->getBodyInertiaTensor( ), expectedInertiaAfterMassUpdate, 5.0e-15 );

    const std::shared_ptr< SphericalHarmonicsGravityField > canonicalGravityField =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( canonicalBody->getGravityFieldModel( ) );
    BOOST_REQUIRE( canonicalGravityField != nullptr );
    const Eigen::Vector3d initialCenterOfMass = canonicalBody->getBodyFixedCenterOfMass( );
    Eigen::MatrixXd modifiedCosineCoefficients = cosineCoefficients;
    modifiedCosineCoefficients( 1, 1 ) *= 2.0;
    modifiedCosineCoefficients( 2, 0 ) *= 1.5;
    canonicalGravityField->setCosineCoefficients( modifiedCosineCoefficients );
    BOOST_CHECK_GT( ( canonicalBody->getBodyFixedCenterOfMass( ) - initialCenterOfMass ).norm( ), 0.0 );
    BOOST_CHECK_GT( ( canonicalBody->getBodyInertiaTensor( ) - expectedInertiaAfterMassUpdate ).norm( ), 0.0 );

    // The same gravity changes do not overwrite explicitly configured rigid-body properties.
    explicitBody->getGravityFieldModel( )->resetGravitationalParameter( 3.0 * gravitationalParameter );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( explicitBody->getBodyInertiaTensor( ), explicitInertia, 5.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( explicitBody->getBodyMass( ), explicitMass, 5.0e-15 );
}

BOOST_AUTO_TEST_CASE( testPrescribedGravityVariationUpdatesInertiaDuringPropagation )
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

    BOOST_CHECK( !body->getGravityFieldVariation( integrated_gravity_field_variation ).first );
    BOOST_CHECK( body->getMassProperties( )->isInertiaTensorAvailable( ) );
    BOOST_CHECK( !body->getMassProperties( )->isInertiaTensorDerivativeAvailable( ) );

    body->updateCurrentGravityField( 0.0 );
    const Eigen::Matrix3d initialInertia = body->getBodyInertiaTensor( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( initialInertia, getInertiaTensorFromGravityField( gravityField, 0.4 ), 5.0e-15 );

    body->updateCurrentGravityField( 10.0 );
    const Eigen::Matrix3d variedInertia = body->getBodyInertiaTensor( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( variedInertia, getInertiaTensorFromGravityField( gravityField, 0.4 ), 5.0e-15 );
    BOOST_CHECK_GT( ( variedInertia - initialInertia ).norm( ), 0.0 );
    BOOST_CHECK( !body->getMassProperties( )->isInertiaTensorDerivativeAvailable( ) );

    const Eigen::Vector3d angularVelocity = ( Eigen::Vector3d( ) << 1.0e-4, -2.0e-4, 3.0e-4 ).finished( );
    Eigen::Vector7d rotationalState = Eigen::Vector7d::Zero( );
    rotationalState( 0 ) = 1.0;
    rotationalState.tail( 3 ) = angularVelocity;
    body->setCurrentRotationalStateToLocalFrame( rotationalState );
    const std::shared_ptr< basic_astrodynamics::InertialTorqueModel > inertialTorque = createInertialTorqueModel( body, "TestBody" );
    inertialTorque->updateMembers( 10.0 );
    const Eigen::Vector3d expectedTorque = -angularVelocity.cross( variedInertia * angularVelocity );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( inertialTorque->getTorque( ), expectedTorque, 5.0e-15 );
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

            std::shared_ptr< IntegratorSettings<> > integratorSettings = std::make_shared< IntegratorSettings<> >( rungeKutta4, 0.0, 1.0 );
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
