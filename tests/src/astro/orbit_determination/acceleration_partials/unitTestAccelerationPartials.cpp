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

#include <limits>
#include <string>
#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include <boost/test/unit_test.hpp>

#include <boost/lambda/lambda.hpp>

#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/ephemerides/customRotationalEphemeris.h"
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantDragCoefficient.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/empiricalAccelerationCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/radiationPressureCoefficient.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/ppnParameters.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/directTidalTimeLag.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantThrust.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/relativity/relativisticAccelerationCorrection.h"
#include "tudat/simulation/estimation_setup/createAccelerationPartials.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/thrustSettings.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::relativity;
using namespace tudat::gravitation;
using namespace tudat::aerodynamics;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::unit_conversions;
using namespace tudat::orbit_determination;
using namespace tudat::acceleration_partials;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::electromagnetism;
using namespace tudat::basic_astrodynamics;
using namespace tudat::propulsion;

BOOST_AUTO_TEST_SUITE( test_acceleration_partials )

BOOST_AUTO_TEST_CASE( testPanelledRadiationPressureAccelerationPartials )
{
    tudat::spice_interface::loadStandardSpiceKernels( );

    for( int testIndex = 0; testIndex < 7; testIndex++ )
    {
        // Create empty bodies, earth and sun.
        std::shared_ptr< Body > vehicle = std::make_shared< Body >( );
        double vehicleMass = 400.0;
        vehicle->setConstantBodyMass( vehicleMass );
        std::shared_ptr< Body > sun = std::make_shared< Body >( );
        SystemOfBodies bodies;
        bodies.addBody( vehicle, "Vehicle" );
        bodies.addBody( sun, "Sun" );

        // Load spice kernels.

        // Set current state of sun and earth.
        sun->setState( Eigen::Vector6d::Zero( ) );//getBodyCartesianStateAtEpoch( "Sun", "SSB", "J2000", "NONE", 1.0E6 ) );
        sun->setRadiationSourceModel(
            createRadiationSourceModel(
                getDefaultRadiationSourceModelSettings( "Sun", TUDAT_NAN, TUDAT_NAN ),
                "Sun", bodies ) );
        vehicle->setState(
                    ( Eigen::Vector6d( ) << 1.4E11, 1.0E11, 1.1E11, 0.0, 0.0, 0.0 ).finished( ) );//getBodyCartesianStateAtEpoch(  "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );
        vehicle->setRotationalEphemeris(
    //        std::make_shared< ephemerides::CustomRotationalEphemeris >(
    //            [=](const double){return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); }, "ECLIPJ2000", "VehicleFixed" ) );
                    std::make_shared< tudat::ephemerides::SimpleRotationalEphemeris >( 0.2, 0.4, -0.2, 1.0E-5, 0.0, "ECLIPJ2000", "VehicleFixed" ) );
        vehicle->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );

        // Create links to set and get state functions of bodies.
        std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
                std::bind( &Body::setState, sun, std::placeholders::_1  );
        std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
                std::bind( &Body::setState, vehicle, std::placeholders::_1  );
        std::function< Eigen::Vector6d( ) > sunStateGetFunction =
                std::bind( &Body::getState, sun );
        std::function< Eigen::Vector6d( ) > vehicleStateGetFunction =
                std::bind( &Body::getState, vehicle );

        double specularReflectivity = 0.5;
        double diffuseReflectivity = 0.1;
        bool useInstantaneousReradiation = false;
        std::vector< std::shared_ptr< BodyPanelSettings > > panelSettingsList;

        if( testIndex < 6 )
        {
            if( testIndex % 2 == 1 )
            {
                useInstantaneousReradiation = true;
            }
            std::string panelName = "SolarPanel";
            double scalingValue = 1.0;
            if( testIndex > 3 )
            {
                scalingValue = 0.5;
                panelName = "PanelName";
            }
            panelSettingsList.push_back( std::make_shared< BodyPanelSettings >(
                std::make_shared< FrameFixedBodyPanelGeometrySettings >( -Eigen::Vector3d::UnitX( ), 1.0 ),
                std::make_shared< SpecularDiffuseBodyPanelReflectionLawSettings >( specularReflectivity, diffuseReflectivity, useInstantaneousReradiation ), "SolarPanel" ) );
            if( testIndex > 2 )
            {
                panelSettingsList.push_back( std::make_shared< BodyPanelSettings >(
                    std::make_shared< FrameFixedBodyPanelGeometrySettings >( -Eigen::Vector3d::UnitY( ), 3.254 ),
                    std::make_shared< SpecularDiffuseBodyPanelReflectionLawSettings >( specularReflectivity, diffuseReflectivity, useInstantaneousReradiation ), "SolarPanel" ) );
                panelSettingsList.push_back( std::make_shared< BodyPanelSettings >(
                    std::make_shared< FrameFixedBodyPanelGeometrySettings >( -Eigen::Vector3d::UnitZ( ), 8.654 ),
                    std::make_shared< SpecularDiffuseBodyPanelReflectionLawSettings >( specularReflectivity, diffuseReflectivity, useInstantaneousReradiation ), "SolarPanel" ) );
                panelSettingsList.push_back( std::make_shared< BodyPanelSettings >(
                    std::make_shared< FrameFixedBodyPanelGeometrySettings >( Eigen::Vector3d::UnitX( ), 1.346  ),
                    std::make_shared< SpecularDiffuseBodyPanelReflectionLawSettings >(
                        scalingValue * specularReflectivity, scalingValue * diffuseReflectivity, useInstantaneousReradiation ), panelName ) );
                panelSettingsList.push_back( std::make_shared< BodyPanelSettings >(
                    std::make_shared< FrameFixedBodyPanelGeometrySettings >( Eigen::Vector3d::UnitY( ), 10.4783 ),
                    std::make_shared< SpecularDiffuseBodyPanelReflectionLawSettings >(
                        scalingValue * specularReflectivity, scalingValue * diffuseReflectivity, useInstantaneousReradiation ), panelName ) );
                panelSettingsList.push_back( std::make_shared< BodyPanelSettings >(
                    std::make_shared< FrameFixedBodyPanelGeometrySettings >( Eigen::Vector3d::UnitZ( ), 6.4235 ),
                    std::make_shared< SpecularDiffuseBodyPanelReflectionLawSettings >(
                        scalingValue * specularReflectivity, scalingValue * diffuseReflectivity, useInstantaneousReradiation ), panelName ) );
            }
        }
        else
        {
            panelSettingsList.push_back( std::make_shared< BodyPanelSettings >(
                bodyTrackingPanelGeometry( "Sun", true , 1.0 ),
                std::make_shared< SpecularDiffuseBodyPanelReflectionLawSettings >( specularReflectivity, diffuseReflectivity, useInstantaneousReradiation ), "SolarPanel" ) );
            panelSettingsList.push_back( std::make_shared< BodyPanelSettings >(
                bodyTrackingPanelGeometry( "Sun", false , 1.346 ),
                std::make_shared< SpecularDiffuseBodyPanelReflectionLawSettings >( specularReflectivity, diffuseReflectivity, useInstantaneousReradiation ), "SolarPanel" ) );

        }


        addBodyExteriorPanelledShape(
            std::make_shared< FullPanelledBodySettings >( panelSettingsList ), "Vehicle", bodies );

        auto paneledRadiationPressureTargetSettings =
            std::make_shared< RadiationPressureTargetModelSettings >( paneled_target );

        std::shared_ptr<electromagnetism::PaneledRadiationPressureTargetModel> radiationPressureInterface =
                std::dynamic_pointer_cast<electromagnetism::PaneledRadiationPressureTargetModel>(
                        createRadiationPressureTargetModel(
                                paneledRadiationPressureTargetSettings, "Vehicle", bodies ).at( 0 ) );
        vehicle->addRadiationPressureTargetModel( radiationPressureInterface );

        // Create acceleration model.
        std::shared_ptr< RadiationPressureAcceleration > accelerationModel =
                std::dynamic_pointer_cast< RadiationPressureAcceleration >(
                    createRadiationPressureAccelerationModel(
                        vehicle, sun, "Vehicle", "Sun", bodies ) );
        accelerationModel->updateMembers( 0.0 );
        // Create partial-calculating object.
        std::shared_ptr< AccelerationPartial > accelerationPartial =
            createAnalyticalAccelerationPartial( accelerationModel, { "Vehicle", vehicle }, { "Sun", sun}, bodies );

        std::shared_ptr< EstimatableParameter< double > > parallelScalingFactor =
            std::make_shared< RadiationPressureScalingFactor >( accelerationModel, source_direction_radiation_pressure_scaling_factor, "Vehicle", "Sun" );
        std::shared_ptr< EstimatableParameter< double > > perpendicularScalingFactor =
            std::make_shared< RadiationPressureScalingFactor >( accelerationModel, source_perpendicular_direction_radiation_pressure_scaling_factor, "Vehicle", "Sun" );
        std::shared_ptr< EstimatableParameter< double > > diffuseReflectivityParameter =
            createDoubleParameterToEstimate< double, double >( std::make_shared< EstimatableParameterSettings >( "Vehicle", diffuse_reflectivity, "SolarPanel" ), bodies );
        std::shared_ptr< EstimatableParameter< double > > specularReflectivityParameter =
            createDoubleParameterToEstimate< double, double >( std::make_shared< EstimatableParameterSettings >( "Vehicle", specular_reflectivity, "SolarPanel" ), bodies );

        // Calculate analytical partials.
        accelerationPartial->update( 0.0 );
        Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );

        Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ) );

        Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );

        Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ) );

        //Eigen::MatrixXd partialWrtEmissivities = accelerationPartial->wrtParameter( panelEmissivitiesParameter );
        Eigen::MatrixXd partialWrtParallelScaling = accelerationPartial->wrtParameter( parallelScalingFactor );
        Eigen::MatrixXd partialWrtPerpendicularScaling = accelerationPartial->wrtParameter( perpendicularScalingFactor );
        Eigen::Vector3d partialWrtSpecularReflectivity = accelerationPartial->wrtParameter( specularReflectivityParameter );
        Eigen::Vector3d partialWrtDiffuseReflectivity = accelerationPartial->wrtParameter( diffuseReflectivityParameter );

        // Declare numerical partials.
        Eigen::Matrix3d testPartialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
        Eigen::Matrix3d testPartialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
        Eigen::Matrix3d testPartialWrtSunPosition = Eigen::Matrix3d::Zero( );
        Eigen::Matrix3d testPartialWrtSunVelocity = Eigen::Matrix3d::Zero( );
        Eigen::MatrixXd testPartialWrtEmissivities = Eigen::MatrixXd::Zero( 3, 2 );
        Eigen::MatrixXd testPartialWrtParallelScaling = Eigen::MatrixXd::Zero( 3, 1 );
        Eigen::MatrixXd testPartialWrtPerpendicularScaling = Eigen::MatrixXd::Zero( 3, 1 );
        Eigen::Vector3d testPartialWrtSpecularReflectivity = Eigen::Vector3d::Zero( );
        Eigen::Vector3d testPartialWrtDiffuseReflectivity = Eigen::Vector3d::Zero( );

        // Declare perturbations in position for numerical partial/
        Eigen::Vector3d positionPerturbation;
        positionPerturbation<< 10000.0, 10000.0, 10000.0;
        Eigen::Vector3d velocityPerturbation;
        velocityPerturbation<< 0.1, 0.1, 0.1;

        // Calculate numerical partials.
       std::function< void( const double ) > updateFunction1 = std::bind(
           &RadiationPressureTargetModel::updateMembers, vehicle->getRadiationPressureTargetModel( ), std::placeholders::_1 );
       std::function< void( const double ) > updateFunction2 = std::bind(
            &RadiationSourceModel::updateMembers, sun->getRadiationSourceModel( ), std::placeholders::_1  );
        std::function< void( ) > updateFunction = [=]( ){
            updateFunction2( 0.0 );
            updateFunction1( 0.0 );
        };

        testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
                    sunStateSetFunction, accelerationModel, sun->getState( ), positionPerturbation, 0, updateFunction );
        testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, accelerationModel, vehicle->getState( ), positionPerturbation, 0, updateFunction );
        testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
                    sunStateSetFunction, accelerationModel, sun->getState( ),velocityPerturbation, 3, updateFunction );
        testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, accelerationModel, vehicle->getState( ), velocityPerturbation, 3, updateFunction );
        //    testPartialWrtEmissivities = calculateAccelerationWrtParameterPartials(
        //                panelEmissivitiesParameter, accelerationModel, emissivityPerturbations );

        testPartialWrtParallelScaling = calculateAccelerationWrtParameterPartials(
            parallelScalingFactor, accelerationModel, 10.0, updateFunction );
        testPartialWrtPerpendicularScaling = calculateAccelerationWrtParameterPartials(
            perpendicularScalingFactor, accelerationModel, 10.0, updateFunction );

        testPartialWrtDiffuseReflectivity = calculateAccelerationWrtParameterPartials(
            diffuseReflectivityParameter, accelerationModel, 0.1, updateFunction );
        testPartialWrtSpecularReflectivity = calculateAccelerationWrtParameterPartials(
            specularReflectivityParameter, accelerationModel, 0.1, updateFunction );

        // Compare numerical and analytical results.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                           partialWrtSunPosition, 1.0e-6 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                           partialWrtSunVelocity, std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                           partialWrtVehiclePosition, 1.0e-6 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                           partialWrtVehicleVelocity, std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtParallelScaling,
                                           partialWrtParallelScaling, 1.0E-13 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPerpendicularScaling,
                                           partialWrtPerpendicularScaling, 1.0E-13 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtSpecularReflectivity,
                                           testPartialWrtSpecularReflectivity, 1.0E-8 );
        if( testIndex % 2 == 0 )
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtDiffuseReflectivity,
                                               testPartialWrtDiffuseReflectivity, 1.0E-13 );
        }
        else
        {
            BOOST_CHECK_EQUAL( testPartialWrtDiffuseReflectivity.norm( ), 0.0 );
            BOOST_CHECK_EQUAL( partialWrtDiffuseReflectivity.norm( ), 0.0 );
        }
    }

}


BOOST_AUTO_TEST_CASE( testCentralGravityPartials )
{
    // Create empty bodies, earth and sun.
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    std::shared_ptr< Body > sun = std::make_shared< Body >( );

    SystemOfBodies bodies;
    bodies.addBody( earth, "Earth" );
    bodies.addBody( sun, "Sun" );

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set current state of sun and earth.
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "Sun", "J2000", "NONE", 1.0E6 ) );
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "Sun", "J2000", "NONE", 1.0E6 ) );

    // Get sun gravitational parameter and set gravity field model.
    double sunsGravitationalParameter = getBodyGravitationalParameter( "Sun" );
    std::shared_ptr< GravityFieldModel > sunGravityFieldModel =
            std::make_shared< GravityFieldModel >( sunsGravitationalParameter );
    sun->setGravityFieldModel( sunGravityFieldModel );
    double earthGravitationalParameter = getBodyGravitationalParameter( "Earth" );
    std::shared_ptr< GravityFieldModel > earthGravityFieldModel =
            std::make_shared< GravityFieldModel >( earthGravitationalParameter );
    earth->setGravityFieldModel( earthGravityFieldModel );

    // Create acceleration due to sun on earth.
    std::shared_ptr< CentralGravitationalAccelerationModel3d > gravitationalAcceleration =
            createCentralGravityAcceleratioModel( earth, sun, "Earth", "Sun", 1 );

    // Create central gravity partial.
    std::shared_ptr< AccelerationPartial > centralGravitationPartial =
            createAnalyticalAccelerationPartial( gravitationalAcceleration, std::make_pair( "Earth", earth ),
                                                 std::make_pair( "Sun", sun ), bodies );

    // Create gravitational parameter object.
    std::shared_ptr< EstimatableParameter< double > > sunGravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( sunGravityFieldModel, "Sun" );
    std::shared_ptr< EstimatableParameter< double > > earthGravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( earthGravityFieldModel, "Earth" );

    // Calculate analytical partials.
    centralGravitationPartial->update( 0.0 );
    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    centralGravitationPartial->wrtPositionOfAcceleratedBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    centralGravitationPartial->wrtVelocityOfAcceleratedBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
    centralGravitationPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    centralGravitationPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::Vector3d partialWrtSunGravitationalParameter = centralGravitationPartial->wrtParameter(
                sunGravitationalParameterParameter );
    Eigen::Vector3d partialWrtEarthGravitationalParameter = centralGravitationPartial->wrtParameter(
                earthGravitationalParameterParameter );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10000.0, 10000.0, 10000.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &Body::setState, earth, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            std::bind( &Body::setState, sun, std::placeholders::_1 );
    std::function< Eigen::Vector6d ( ) > earthStateGetFunction =
            std::bind( &Body::getState, earth );
    std::function< Eigen::Vector6d ( ) > sunStateGetFunction =
            std::bind( &Body::getState, sun );

    // Calculate numerical partials.
    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), positionPerturbation, 0 );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), velocityPerturbation, 3 );
    testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, gravitationalAcceleration, sun->getState( ), positionPerturbation, 0 );
    testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, gravitationalAcceleration, sun->getState( ), velocityPerturbation, 3 );
    Eigen::Vector3d testPartialWrtSunGravitationalParameter = calculateAccelerationWrtParameterPartials(
                sunGravitationalParameterParameter, gravitationalAcceleration, 1.0E12 );
    Eigen::Vector3d testPartialWrtEarthGravitationalParameter = calculateAccelerationWrtParameterPartials(
                earthGravitationalParameterParameter, gravitationalAcceleration, 1.0E12 );

    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                       partialWrtEarthPosition, 1.0E-8 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                       partialWrtEarthVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                       partialWrtSunPosition, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                       partialWrtSunVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunGravitationalParameter,
                                       partialWrtSunGravitationalParameter, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtEarthGravitationalParameter,
                                       testPartialWrtEarthGravitationalParameter, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtEarthGravitationalParameter,
                                       partialWrtSunGravitationalParameter, std::numeric_limits< double >::epsilon(  ) );
}

BOOST_AUTO_TEST_CASE( testCannonballRadiationPressureAccelerationPartials )
{
    // Create empty bodies, earth and sun.
    std::shared_ptr< Body > vehicle = std::make_shared< Body >( );
    vehicle->setConstantBodyMass( 400.0 );
    std::shared_ptr< Body > sun = std::make_shared< Body >( );

    SystemOfBodies bodies;
    bodies.addBody( vehicle, "Vehicle" );
    bodies.addBody( sun, "Sun" );

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set current state of sun and earth.
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "SSB", "J2000", "NONE", 1.0E6 ) );
    vehicle->setState( getBodyCartesianStateAtEpoch(  "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );

    // Create links to set and get state functions of bodies.
    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            std::bind( &Body::setState, sun, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            std::bind( &Body::setState, vehicle, std::placeholders::_1 );
    std::function< Eigen::Vector6d( ) > sunStateGetFunction =
            std::bind( &Body::getState, sun );
    std::function< Eigen::Vector6d( ) > vehicleStateGetFunction =
            std::bind( &Body::getState, vehicle );
    bodies.at( "Sun" )->setRadiationSourceModel(
        createRadiationSourceModel(
            getDefaultRadiationSourceModelSettings( "Sun", TUDAT_NAN, TUDAT_NAN ), "Sun", bodies ) );

    // Create radiation pressure properties of vehicle
    std::shared_ptr<CannonballRadiationPressureTargetModel > radiationPressureInterface =
        std::make_shared< CannonballRadiationPressureTargetModel >( mathematical_constants::PI * 0.3 * 0.3, 1.2 );
    vehicle->addRadiationPressureTargetModel( radiationPressureInterface );

    // Create acceleration model.
    std::shared_ptr< RadiationPressureAcceleration > accelerationModel =
        createRadiationPressureAccelerationModel( bodies.at( "Vehicle" ), bodies.at( "Sun" ), "Vehicle", "Sun", bodies );

    // Create partial-calculating object.
    std::shared_ptr< AccelerationPartial > accelerationPartial =
            createAnalyticalAccelerationPartial( accelerationModel, std::make_pair( "Vehicle", vehicle ),
                                                 std::make_pair( "Sun", sun ), bodies );

    // Create parameter object
    std::string vehicleName = "Vehicle";
    std::shared_ptr< EstimatableParameter< double > > radiationPressureCoefficient =
            std::make_shared< RadiationPressureCoefficient >( radiationPressureInterface, vehicleName );

    std::shared_ptr< EstimatableParameter< double > > radiationPressureScalingSourceDirection =
        std::make_shared< RadiationPressureScalingFactor >( accelerationModel, source_direction_radiation_pressure_scaling_factor, vehicleName, "Sun" );

    std::vector< double > timeLimits;
    timeLimits.push_back( 0.0 );
    timeLimits.push_back( 3600.0 );
    timeLimits.push_back( 7200.0 );
    timeLimits.push_back( 10800.0 );

    std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > arcWiseRadiationPressureCoefficient =
            std::make_shared< ArcWiseRadiationPressureCoefficient >( radiationPressureInterface, timeLimits, vehicleName );



    // Calculate analytical partials.
    double currentTime = 0.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::Vector3d partialWrtRadiationPressureCoefficient = accelerationPartial->wrtParameter(
                radiationPressureCoefficient );
    Eigen::Vector3d partialWrtScalingFactor = accelerationPartial->wrtParameter(
        radiationPressureScalingSourceDirection );


    // Get arc-wise radiation pressure coefficient partials
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );
    currentTime = 1000.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise2 = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );
    currentTime = 4000.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise3 = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );
    currentTime = 7000.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise4 = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );
    currentTime = 10000.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise5 = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );
    currentTime = 12000.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise6 = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );

    // Check whether arc-wise radiation pressure partials are properly segmented
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 4; j++ )
        {
            if( j != 0 )
            {
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise2( i, j ), 0.0 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise2( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
            }

            if( j != 1 )
            {
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise3( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise4( i, j ), 0.0 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise3( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise4( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
            }

            if( j != 2 )
            {
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise5( i, j ), 0.0 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise5( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
            }

            if( j != 3 )
            {
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise6( i, j ), 0.0 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise6( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
            }
        }
    }

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Vector3d testPartialWrtRadiationPressureCoefficient = Eigen::Vector3d::Zero( );
    Eigen::Vector3d testPartialWrtScalingFactor = Eigen::Vector3d::Zero( );


    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10000.0, 10000.0, 10000.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Calculate numerical partials.
    std::function< void( ) > updateFunction =
            std::bind( &RadiationPressureTargetModel::updateMembers, radiationPressureInterface, 0.0 );
    testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, accelerationModel, sun->getState( ), positionPerturbation, 0, updateFunction );
    testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, vehicle->getState( ), positionPerturbation, 0, updateFunction );
    testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, accelerationModel, sun->getState( ),velocityPerturbation, 3, updateFunction );
    testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, vehicle->getState( ), velocityPerturbation, 3, updateFunction );
    testPartialWrtRadiationPressureCoefficient = calculateAccelerationWrtParameterPartials(
                radiationPressureCoefficient, accelerationModel, 1.0E-2, updateFunction );
    testPartialWrtScalingFactor = calculateAccelerationWrtParameterPartials(
        radiationPressureScalingSourceDirection, accelerationModel, 100.0, updateFunction );


    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                       partialWrtSunPosition, 1.0E-8 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                       partialWrtSunVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                       partialWrtVehiclePosition, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                       partialWrtVehicleVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtRadiationPressureCoefficient,
                                       partialWrtRadiationPressureCoefficient, 1.0E-12 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtScalingFactor,
                                       partialWrtScalingFactor, 1.0E-12 );
}

BOOST_AUTO_TEST_CASE( testThirdBodyGravityPartials )
{
    // Create empty bodies, earth and sun.
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    std::shared_ptr< Body > sun = std::make_shared< Body >( );
    std::shared_ptr< Body > moon = std::make_shared< Body >( );

    SystemOfBodies bodies;
    bodies.addBody( earth, "Earth" );
    bodies.addBody( sun, "Sun" );
    bodies.addBody( moon, "Moon" );

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set current state of sun and earth.
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "Sun", "J2000", "NONE", 1.0E6 ) );
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "Sun", "J2000", "NONE", 1.0E6 ) );
    moon->setState( getBodyCartesianStateAtEpoch(  "Moon", "Sun", "J2000", "NONE", 1.0E6 ) );

    // Get sun gravitational parameter and set gravity field model.
    double sunsGravitationalParameter = getBodyGravitationalParameter( "Sun" );
    std::shared_ptr< GravityFieldModel > sunGravityFieldModel =
            std::make_shared< GravityFieldModel >( sunsGravitationalParameter );
    sun->setGravityFieldModel( sunGravityFieldModel );

    double moonsGravitationalParameter = getBodyGravitationalParameter( "Moon" );
    std::shared_ptr< GravityFieldModel > moonGravityFieldModel =
            std::make_shared< GravityFieldModel >( moonsGravitationalParameter );
    moon->setGravityFieldModel( moonGravityFieldModel );

    double earthGravitationalParameter = getBodyGravitationalParameter( "Earth" );
    std::shared_ptr< GravityFieldModel > earthGravityFieldModel =
            std::make_shared< GravityFieldModel >( earthGravitationalParameter );
    earth->setGravityFieldModel( earthGravityFieldModel );

    // Create acceleration due to moon on earth.
    std::shared_ptr< ThirdBodyCentralGravityAcceleration > gravitationalAcceleration =
            createThirdBodyCentralGravityAccelerationModel(
                moon, sun, earth, "Moon", "Sun", "Earth" );

    // Create central gravity partial.
    std::shared_ptr< AccelerationPartial > thirdBodyGravitationPartial =
            createAnalyticalAccelerationPartial( gravitationalAcceleration, std::make_pair( "Moon", moon ),
                                                 std::make_pair( "Sun", sun ), bodies );

    // Create gravitational parameter object.
    std::shared_ptr< EstimatableParameter< double > > gravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( sunGravityFieldModel, "Sun" );
    std::shared_ptr< EstimatableParameter< double > > moonGravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( moonGravityFieldModel, "Moon" );
    std::shared_ptr< EstimatableParameter< double > > earthGravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( earthGravityFieldModel, "Earth" );

    // Calculate analytical partials.
    thirdBodyGravitationPartial->update( 1.0E6 );
    Eigen::MatrixXd partialWrtMoonPosition = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtPositionOfAcceleratedBody( partialWrtMoonPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtMoonVelocity = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtVelocityOfAcceleratedBody( partialWrtMoonVelocity.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtPositionOfAdditionalBody( "Earth", partialWrtEarthPosition.block( 0, 0, 3, 3 )  );
    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtVelocityOfAdditionalBody( "Earth", partialWrtEarthVelocity.block( 0, 0, 3, 3 )  );

    Eigen::MatrixXd partialWrtSunGravitationalParameter = thirdBodyGravitationPartial->wrtParameter(
                gravitationalParameterParameter );
    Eigen::MatrixXd partialWrtMoonGravitationalParameter = thirdBodyGravitationPartial->wrtParameter(
                moonGravitationalParameterParameter );
    Eigen::MatrixXd partialWrtEarthGravitationalParameter = thirdBodyGravitationPartial->wrtParameter(
                earthGravitationalParameterParameter );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtMoonPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMoonVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10000.0, 10000.0, 10000.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector6d ) > moonStateSetFunction =
            std::bind( &Body::setState, moon, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            std::bind( &Body::setState, sun, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &Body::setState, earth, std::placeholders::_1 );

    // Calculate numerical partials.
    testPartialWrtMoonPosition = calculateAccelerationWrtStatePartials(
                moonStateSetFunction, gravitationalAcceleration, moon->getState( ), positionPerturbation, 0 );
    testPartialWrtMoonVelocity = calculateAccelerationWrtStatePartials(
                moonStateSetFunction, gravitationalAcceleration, moon->getState( ), velocityPerturbation, 3 );
    testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, gravitationalAcceleration, sun->getState( ), positionPerturbation, 0 );
    testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, gravitationalAcceleration, sun->getState( ), velocityPerturbation, 3 );
    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), positionPerturbation, 0 );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), velocityPerturbation, 3 );
    Eigen::Vector3d testPartialWrtSunGravitationalParameter = calculateAccelerationWrtParameterPartials(
                gravitationalParameterParameter, gravitationalAcceleration, 1.0E16 );
    Eigen::Vector3d testPartialWrtEarthGravitationalParameter = calculateAccelerationWrtParameterPartials(
                earthGravitationalParameterParameter, gravitationalAcceleration, 1.0E16 );
    Eigen::Vector3d testPartialWrtMoonGravitationalParameter = calculateAccelerationWrtParameterPartials(
                moonGravitationalParameterParameter, gravitationalAcceleration, 1.0E16 );

    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMoonPosition,
                                       partialWrtMoonPosition, 1.0E-7 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMoonVelocity,
                                       partialWrtMoonVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                       partialWrtSunPosition, 1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                       partialWrtSunVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                       partialWrtEarthPosition, 1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                       partialWrtEarthVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunGravitationalParameter,
                                       partialWrtSunGravitationalParameter, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMoonGravitationalParameter,
                                       partialWrtMoonGravitationalParameter, std::numeric_limits< double >::epsilon(  ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthGravitationalParameter,
                                       partialWrtEarthGravitationalParameter, std::numeric_limits< double >::epsilon(  ) );
}


void updateFlightConditionsWithPerturbedState(
        const std::shared_ptr< aerodynamics::FlightConditions > flightConditions,
        const double timeToUpdate )
{
    flightConditions->resetCurrentTime( );
    flightConditions->updateConditions( timeToUpdate );
}

BOOST_AUTO_TEST_CASE( testAerodynamicAccelerationPartials )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    using namespace tudat;
    // Create Earth object
    BodyListSettings defaultBodySettings =
            getDefaultBodySettings( { "Earth" } );
    defaultBodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ) );
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Create vehicle objects.
    double vehicleMass = 5.0E3;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );


    bool areCoefficientsInAerodynamicFrame = 1;
    Eigen::Vector3d aerodynamicCoefficients = ( Eigen::Vector3d( ) << 2.5, -0.1, 0.5 ).finished( );

    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                2.0, 4.0, Eigen::Vector3d::Zero( ), aerodynamicCoefficients, Eigen::Vector3d::Zero( ),
                aerodynamics::getAerodynamicCoefficientFrame( areCoefficientsInAerodynamicFrame, 1 ),
                aerodynamics::getAerodynamicCoefficientFrame( areCoefficientsInAerodynamicFrame, 1 ) );
    bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ) );






    // Set spherical elements for vehicle.
    Eigen::Vector6d vehicleSphericalEntryState;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.0;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.7E3;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            -0.9 * mathematical_constants::PI / 180.0;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;

    // Convert vehicle state from spherical elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = tudat::ephemerides::transformStateToTargetFrame(
                convertSphericalOrbitalToCartesianState(
                    vehicleSphericalEntryState ), 0.0, bodies.at( "Earth" )->getRotationalEphemeris( ) );


    bodies.at( "Earth" )->setStateFromEphemeris( 0.0 );
    bodies.at( "Earth" )->setCurrentRotationToLocalFrameFromEphemeris( 0.0 );

    bodies.at( "Vehicle" )->setState( systemInitialState );


    std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationModel =
            simulation_setup::createAerodynamicAcceleratioModel(
                bodies.at( "Vehicle" ), bodies.at( "Earth" ), "Vehicle", "Earth" );
    bodies.at( "Vehicle" )->getFlightConditions( )->updateConditions( 0.0 );
    accelerationModel->updateMembers( 0.0 );

    std::shared_ptr< AccelerationPartial > aerodynamicAccelerationPartial =
            createAnalyticalAccelerationPartial(
                accelerationModel, std::make_pair( "Vehicle", bodies.at( "Vehicle" ) ),
                std::make_pair( "Earth", bodies.at( "Earth" ) ), bodies );

    // Create gravitational parameter object.
    std::shared_ptr< EstimatableParameter< double > > dragCoefficientParameter = std::make_shared<
            ConstantDragCoefficient >( std::dynamic_pointer_cast< aerodynamics::CustomAerodynamicCoefficientInterface >(
                                           bodies.at( "Vehicle" )->getAerodynamicCoefficientInterface( ) ), "Vehicle" );

    // Calculate analytical partials.
    aerodynamicAccelerationPartial->update( 0.0 );
    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    aerodynamicAccelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    aerodynamicAccelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    aerodynamicAccelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    aerodynamicAccelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::Vector3d partialWrtDragCoefficient = aerodynamicAccelerationPartial->wrtParameter(
                dragCoefficientParameter );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );

    std::function< void( ) > environmentUpdateFunction =
            std::bind( &updateFlightConditionsWithPerturbedState, bodies.at( "Vehicle" )->getFlightConditions( ), 0.0 );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1.0, 1.0, 1.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            std::bind( &Body::setState, bodies.at( "Vehicle" ), std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &Body::setState, bodies.at( "Earth" ), std::placeholders::_1 );
    std::function< Eigen::Vector6d ( ) > vehicleStateGetFunction =
            std::bind( &Body::getState, bodies.at( "Vehicle" ) );
    std::function< Eigen::Vector6d ( ) > earthStateGetFunction =
            std::bind( &Body::getState, bodies.at( "Earth" ) );

    // Calculate numerical partials.
    testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, bodies.at( "Vehicle" )->getState( ), positionPerturbation, 0,
                environmentUpdateFunction);
    testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, bodies.at( "Vehicle" )->getState( ), velocityPerturbation, 3,
                environmentUpdateFunction );
    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, accelerationModel, bodies.at( "Earth" )->getState( ), positionPerturbation, 0,
                environmentUpdateFunction );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, accelerationModel, bodies.at( "Earth" )->getState( ), velocityPerturbation, 3,
                environmentUpdateFunction );

    Eigen::Vector3d testPartialWrtDragCoefficient = calculateAccelerationWrtParameterPartials(
                dragCoefficientParameter, accelerationModel, 1.0E-4, environmentUpdateFunction );

    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                       partialWrtVehiclePosition, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                       partialWrtVehicleVelocity, 1.0E-6  );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                       partialWrtEarthPosition, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                       partialWrtEarthVelocity, 1.0E-6 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtDragCoefficient,
                                       partialWrtDragCoefficient, 1.0E-10 );
}


BOOST_AUTO_TEST_CASE( testRelativisticAccelerationPartial )
{
    // Create earth and vehicle bodies.
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    std::shared_ptr< Body > vehicle = std::make_shared< Body >( );

    // Create links to set and get state functions of bodies.
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &Body::setState, earth, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            std::bind( &Body::setState, vehicle, std::placeholders::_1  );
    std::function< Eigen::Vector6d( ) > earthStateGetFunction =
            std::bind( &Body::getState, earth );
    std::function< Eigen::Vector6d( ) > vehicleStateGetFunction =
            std::bind( &Body::getState, vehicle );

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    // Set vehicle and earth state.
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );
    Eigen::Vector6d vehicleKeplerElements;
    vehicleKeplerElements << 6378.0E3 + 249E3, 0.0004318, convertDegreesToRadians( 96.5975 ),
            convertDegreesToRadians( 217.6968 ), convertDegreesToRadians( 268.2663 ), convertDegreesToRadians( 142.3958 );
    vehicle->setState( earth->getState( ) + convertKeplerianToCartesianElements( vehicleKeplerElements,
                                                                                 getBodyGravitationalParameter( "Earth" ) ) );


    SystemOfBodies bodies;
    bodies.addBody( vehicle, "Vehicle" );;
    bodies.addBody( earth, "Earth" );

    // Create gravity field.
    std::shared_ptr< GravityFieldSettings > gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );
    std::shared_ptr< gravitation::GravityFieldModel > earthGravityField =
            createGravityFieldModel( gravityFieldSettings, "Earth", bodies );
    earth->setGravityFieldModel( earthGravityField );

    // Create acceleration model.
    std::function< double( ) > ppnParameterGammaFunction = std::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet );
    std::function< double( ) > ppnParameterBetaFunction = std::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet );
    std::shared_ptr< RelativisticAccelerationCorrection > accelerationModel =
            std::make_shared< RelativisticAccelerationCorrection >
            ( std::bind( &Body::getState, vehicle ),
              std::bind( &Body::getState, earth ),
              std::bind( &GravityFieldModel::getGravitationalParameter, earthGravityField ),
              ppnParameterGammaFunction, ppnParameterBetaFunction );

    // Create acceleration partial object.
    std::shared_ptr< RelativisticAccelerationPartial > accelerationPartial = std::make_shared< RelativisticAccelerationPartial >(
                accelerationModel, "Vehicle", "Earth" );

    // Create parameter objects.
    std::shared_ptr< EstimatableParameter< double > > gravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( earthGravityField, "Earth" );
    std::shared_ptr< EstimatableParameter< double > > ppnParameterGamma = std::make_shared<
            PPNParameterGamma >( ppnParameterSet );
    std::shared_ptr< EstimatableParameter< double > > ppnParameterBeta = std::make_shared<
            PPNParameterBeta >( ppnParameterSet );

    // Calculate analytical partials.
    accelerationModel->updateMembers( );
    accelerationPartial->update( );

    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ) );

    Eigen::Vector3d partialWrtEarthGravitationalParameter = accelerationPartial->wrtParameter(
                gravitationalParameterParameter );

    Eigen::Vector3d partialWrtGamma = accelerationPartial->wrtParameter( ppnParameterGamma );
    Eigen::Vector3d partialWrtBeta = accelerationPartial->wrtParameter( ppnParameterBeta );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10.0, 10.0, 10.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Calculate numerical partials.
    Eigen::Matrix3d testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, vehicle->getState( ), positionPerturbation, 0 );
    Eigen::Matrix3d testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, vehicle->getState( ), velocityPerturbation, 3 );
    Eigen::Matrix3d testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, accelerationModel, earth->getState( ), positionPerturbation, 0 );
    Eigen::Matrix3d testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, accelerationModel, earth->getState( ), velocityPerturbation, 3 );
    Eigen::Vector3d testPartialWrtEarthGravitationalParameter = calculateAccelerationWrtParameterPartials(
                gravitationalParameterParameter, accelerationModel, 1.0E10 );
    Eigen::Vector3d testPartialWrtPpnParameterGamma = calculateAccelerationWrtParameterPartials(
                ppnParameterGamma, accelerationModel, 100.0 );
    Eigen::Vector3d testPartialWrtPpnParameterBeta = calculateAccelerationWrtParameterPartials(
                ppnParameterBeta, accelerationModel, 100.0 );

    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                       partialWrtEarthPosition, 1.0e-7 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                       partialWrtEarthVelocity, 1.0e-7 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                       partialWrtVehiclePosition, 1.0e-7 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                       partialWrtVehicleVelocity, 1.0e-7 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthGravitationalParameter,
                                       partialWrtEarthGravitationalParameter, 1.0e-8 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterGamma, partialWrtGamma, 1.0e-8 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterBeta, partialWrtBeta, 1.0e-8 );
}


BOOST_AUTO_TEST_CASE( testEmpiricalAccelerationPartial )
{

    // Create earth and vehicle bodies.
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    std::shared_ptr< Body > vehicle = std::make_shared< Body >( );

    // Create links to set and get state functions of bodies.
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &Body::setState, earth, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            std::bind( &Body::setState, vehicle, std::placeholders::_1  );
    std::function< Eigen::Vector6d( ) > earthStateGetFunction =
            std::bind( &Body::getState, earth );
    std::function< Eigen::Vector6d( ) > vehicleStateGetFunction =
            std::bind( &Body::getState, vehicle );

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    // Set vehicle and earth state.
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );
    Eigen::Vector6d vehicleKeplerElements;
    vehicleKeplerElements << 6378.0E3 + 249E3, 0.0004318, convertDegreesToRadians( 96.5975 ),
            convertDegreesToRadians( 217.6968 ), convertDegreesToRadians( 268.2663 ), convertDegreesToRadians( 142.3958 );
    vehicle->setState( earth->getState( ) + convertKeplerianToCartesianElements( vehicleKeplerElements,
                                                                                 getBodyGravitationalParameter( "Earth" ) ) );


    SystemOfBodies bodies;
    bodies.addBody( vehicle, "Vehicle" );;
    bodies.addBody( earth, "Earth" );

    // Create gravity field.
    std::shared_ptr< GravityFieldSettings > gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );
    std::shared_ptr< gravitation::GravityFieldModel > earthGravityField =
            createGravityFieldModel( gravityFieldSettings, "Earth", bodies );
    earth->setGravityFieldModel( earthGravityField );

    // Create rotation model
    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > simpleRotationalEphemeris =
            std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000" , "IAU_Earth", 0.0 ),
                2.0 * mathematical_constants::PI / 86400.0,
                1.0E7,
                "ECLIPJ2000" , "IAU_Earth" );
    earth->setRotationalEphemeris( simpleRotationalEphemeris );
    earth->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );

    // Define empirical acceleration
    Eigen::Vector3d constantAcceleration = 0.0 * Eigen::Vector3d( 0.038, -0.7528, 0.00752 );
    Eigen::Vector3d sineAcceleration = Eigen::Vector3d( 0.984, 0.0427, -0.0764238 );
    Eigen::Vector3d cosineAcceleration = Eigen::Vector3d( -0.0024785, 1.839, -0.73288 );

    // Create acceleration model.
    std::shared_ptr< EmpiricalAcceleration > accelerationModel =
            std::make_shared< EmpiricalAcceleration >
            ( constantAcceleration, sineAcceleration, cosineAcceleration,
              std::bind( &Body::getState, vehicle ),
              std::bind( &GravityFieldModel::getGravitationalParameter, earthGravityField ),
              std::bind( &Body::getState, earth ) );

    // Create acceleration partial object.
    std::shared_ptr< EmpiricalAccelerationPartial > accelerationPartial = std::make_shared< EmpiricalAccelerationPartial >(
                accelerationModel, "Vehicle", "Earth" );

    // Define list of empirical accelerations w.r.t. which partials are to be computed
    std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > allEmpiricalShapesVector;
    allEmpiricalShapesVector.push_back( basic_astrodynamics::constant_empirical );
    allEmpiricalShapesVector.push_back( basic_astrodynamics::cosine_empirical );

    std::map< basic_astrodynamics::EmpiricalAccelerationComponents, std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > >
            empiricalComponentsToEstimate;
    empiricalComponentsToEstimate[ basic_astrodynamics::radial_empirical_acceleration_component ] = allEmpiricalShapesVector;
    empiricalComponentsToEstimate[ basic_astrodynamics::along_track_empirical_acceleration_component ] = allEmpiricalShapesVector;

    allEmpiricalShapesVector.push_back( basic_astrodynamics::sine_empirical );
    empiricalComponentsToEstimate[ basic_astrodynamics::across_track_empirical_acceleration_component ] = allEmpiricalShapesVector;

    // Create time-independent empirical acceleration object.
    std::shared_ptr< EmpiricalAccelerationCoefficientsParameter > empiricalAccelerationParameter = std::make_shared<
            EmpiricalAccelerationCoefficientsParameter >(
                std::vector< std::shared_ptr< EmpiricalAcceleration > >( { accelerationModel } ), "Vehicle", "Earth",
                empiricalComponentsToEstimate );

    {
        // Calculate analytical partials.
        accelerationModel->updateMembers( );
        accelerationPartial->update( );
        Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtEmpiricalCoefficients;
        accelerationPartial->wrtEmpiricalAccelerationCoefficient(
                    empiricalAccelerationParameter, partialWrtEmpiricalCoefficients );

        // Declare perturbations in position for numerical partial
        Eigen::Vector3d positionPerturbation;
        positionPerturbation << 1.0, 1.0, 1.0;
        Eigen::Vector3d velocityPerturbation;
        velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;
        int parameterSize = empiricalAccelerationParameter->getParameterSize( );
        Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Constant( parameterSize, 1.0E-5 );

        // Calculate numerical partials.
        Eigen::MatrixXd testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, accelerationModel, vehicle->getState( ), 10.0 * positionPerturbation, 0 );
        Eigen::MatrixXd testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, accelerationModel, vehicle->getState( ), velocityPerturbation, 3 );
        Eigen::MatrixXd testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                    earthStateSetFunction, accelerationModel, earth->getState( ), positionPerturbation, 0 );
        Eigen::MatrixXd testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                    earthStateSetFunction, accelerationModel, earth->getState( ), velocityPerturbation, 3 );
        Eigen::MatrixXd testPartialWrtEmpiricalCoefficients = calculateAccelerationWrtParameterPartials(
                    empiricalAccelerationParameter, accelerationModel, parameterPerturbation );

        // Compare numerical and analytical results.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                           partialWrtEarthPosition, 1.0e-3 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                           partialWrtEarthVelocity, 1.0e-3 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                           partialWrtVehiclePosition, 1.0e-3 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                           partialWrtVehicleVelocity, 1.0e-3 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEmpiricalCoefficients,
                                           partialWrtEmpiricalCoefficients, 1.0e-3 );
    }

    // Define arc split times of arc-wise empirical accelerations
    std::vector< double > arcStartTimes;
    arcStartTimes.push_back( 0.0 );
    arcStartTimes.push_back( 1.0E4 );
    arcStartTimes.push_back( 2.0E4 );
    arcStartTimes.push_back( 3.0E4 );
    arcStartTimes.push_back( 7.0E4 );
    std::shared_ptr< ArcWiseEmpiricalAccelerationCoefficientsParameter > arcWiseEmpiricalAccelerationParameter = std::make_shared<
            ArcWiseEmpiricalAccelerationCoefficientsParameter >(
                std::vector< std::shared_ptr< EmpiricalAcceleration > >( { accelerationModel } ), "Vehicle", "Earth",
                empiricalComponentsToEstimate, arcStartTimes );

    // Define list of times at which to test empirical acceleration
    std::vector< double > evaluationTimes;

    evaluationTimes.push_back( 1.0 );
    evaluationTimes.push_back( 0.5E4 );
    evaluationTimes.push_back( 1.0E4 - 1.0 );
    evaluationTimes.push_back( 1.0E4 + 1.0 );
    evaluationTimes.push_back( 1.2E4 );
    evaluationTimes.push_back( 3.5E4 );
    evaluationTimes.push_back( 1.0E5 );

    Eigen::VectorXd accelerationPerturbationVector = Eigen::VectorXd::Zero( arcWiseEmpiricalAccelerationParameter->getParameterSize( ) );
    for( int i = 0; i < arcWiseEmpiricalAccelerationParameter->getParameterSize( ); i++  )
    {
        accelerationPerturbationVector( i ) = 1.0E-6;
    }

    // Iterate over all test times.
    for( unsigned int i = 0; i < evaluationTimes.size( ); i++ )
    {
        // Update models to current time
        accelerationModel->resetCurrentTime( );
        accelerationModel->updateMembers( evaluationTimes.at( i ) );
        accelerationPartial->update( evaluationTimes.at( i ) );

        // Compute analytical partials
        Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtEmpiricalCoefficients = accelerationPartial->wrtParameter(
                    arcWiseEmpiricalAccelerationParameter );

        // Set numerical partial settings
        Eigen::Vector3d positionPerturbation;
        positionPerturbation << 1.0, 1.0, 1.0;
        Eigen::Vector3d velocityPerturbation;
        velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;
        int parameterSize = empiricalAccelerationParameter->getParameterSize( );
        Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Constant( parameterSize, 1.0E-5 );

        // Calculate numerical partials.
        Eigen::MatrixXd testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, accelerationModel, vehicle->getState( ), positionPerturbation, 0, emptyFunction,
                    evaluationTimes.at( i ) );
        Eigen::MatrixXd testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, accelerationModel, vehicle->getState( ), velocityPerturbation, 3, emptyFunction,
                    evaluationTimes.at( i ) );
        Eigen::MatrixXd testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                    earthStateSetFunction, accelerationModel, earth->getState( ), positionPerturbation, 0, emptyFunction,
                    evaluationTimes.at( i ) );
        Eigen::MatrixXd testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                    earthStateSetFunction, accelerationModel, earth->getState( ), velocityPerturbation, 3, emptyFunction,
                    evaluationTimes.at( i ) );
        Eigen::MatrixXd testPartialWrtEmpiricalCoefficients = calculateAccelerationWrtParameterPartials(
                    arcWiseEmpiricalAccelerationParameter, accelerationModel, accelerationPerturbationVector,
                    emptyFunction, evaluationTimes.at( i ) );


        //Compare numerical and analytical results.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                           partialWrtEarthPosition, 1.0e-5 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                           partialWrtEarthVelocity, 1.0e-6 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                           partialWrtVehiclePosition, 1.0e-5 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                           partialWrtVehicleVelocity, 1.0e-6 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEmpiricalCoefficients,
                                           partialWrtEmpiricalCoefficients, 1.0e-6 );
    }
}

BOOST_AUTO_TEST_CASE( testDirectDissipationAccelerationPartial )
{

    // Create bodies
    std::shared_ptr< Body > jupiter = std::make_shared< Body >( );
    std::shared_ptr< Body > io = std::make_shared< Body >( );

    // Create links to set and get state functions of bodies.
    std::function< void( Eigen::Vector6d ) > ioStateSetFunction =
            std::bind( &Body::setState, io, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > jupiterStateSetFunction =
            std::bind( &Body::setState, jupiter, std::placeholders::_1  );
    std::function< Eigen::Vector6d( ) > ioStateGetFunction =
            std::bind( &Body::getState, io );
    std::function< Eigen::Vector6d( ) > jupiterStateGetFunction =
            std::bind( &Body::getState, jupiter );

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    // Set state.
    jupiter->setState( Eigen::Vector6d::Zero( ) );
    Eigen::Vector6d ioKeplerElements =
            ( Eigen::Vector6d( ) << 1.0 * 421.8E6, 1.0 * 0.004, 0.001, 2.0, 3.0, 0.4 ).finished( );
    io->setState( convertKeplerianToCartesianElements(
                      ioKeplerElements,
                      getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ) ) );


    SystemOfBodies bodies;
    bodies.addBody( io, "Io" );
    bodies.addBody( jupiter, "Jupiter" );


    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients( 0, 0 ) = 1.0;
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );

    // Create jupiter gravity field.
    std::shared_ptr< GravityFieldSettings > jupiterGravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >
            ( getBodyGravitationalParameter( "Jupiter" ), getAverageRadius( "Jupiter" ),
              cosineCoefficients, sineCoefficients, "IAU_Jupiter" );
    std::shared_ptr< gravitation::GravityFieldModel > jupiterGravityField =
            createGravityFieldModel( jupiterGravityFieldSettings, "Jupiter", bodies );
    jupiter->setGravityFieldModel( jupiterGravityField );

    // Create io gravity field.
    std::shared_ptr< GravityFieldSettings > ioGravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >
            ( getBodyGravitationalParameter( "Io" ), getAverageRadius( "Io" ),
              cosineCoefficients, sineCoefficients, "IAU_Io" );
    std::shared_ptr< gravitation::GravityFieldModel > ioGravityField =
            createGravityFieldModel( ioGravityFieldSettings, "Io", bodies );
    io->setGravityFieldModel( ioGravityField );

    // Create rotation model
    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > simpleRotationalEphemeris =
            std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ),
                2.0 * mathematical_constants::PI / ( 9.925 * 3600.0 ), 0.0,
                "ECLIPJ2000" , "IAU_Jupiter" );
    jupiter->setRotationalEphemeris( simpleRotationalEphemeris );
    jupiter->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );

    double loveNumber = 0.1;
    double timeLag = 100.0;
    for( unsigned int useRadialTerm = 0; useRadialTerm < 2; useRadialTerm++ )
    {
        for( unsigned int usePlanetTide = 0; usePlanetTide < 2; usePlanetTide++ )
        {

            // Create acceleration model.
            std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > accelerationModel =
                    std::dynamic_pointer_cast<  gravitation::DirectTidalDissipationAcceleration >(
                        simulation_setup::createAccelerationModel(
                            io, jupiter, std::make_shared< simulation_setup::DirectTidalDissipationAccelerationSettings >(
                                loveNumber, timeLag, bool( useRadialTerm ), bool( usePlanetTide ) ) , "Io", "Jupiter" ) );

            // Create acceleration partial object.
            std::shared_ptr< acceleration_partials::DirectTidalDissipationAccelerationPartial > accelerationPartial =
                    std::make_shared< acceleration_partials::DirectTidalDissipationAccelerationPartial >(
                        accelerationModel, "Io", "Jupiter" );

            // Create gravitational parameter object.
            std::shared_ptr< EstimatableParameter< double > > jupiterGravitationalParameterParameter = std::make_shared<
                    GravitationalParameter >( jupiterGravityField, "Jupiter" );
            std::shared_ptr< EstimatableParameter< double > > ioGravitationalParameterParameter = std::make_shared<
                    GravitationalParameter >( ioGravityField, "Io" );
            std::shared_ptr< EstimatableParameter< double > > tidalTimeLagParameter = std::make_shared< DirectTidalTimeLag >(
                        std::vector< std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > >{ accelerationModel },
                        usePlanetTide ? "Jupiter" : "Io" );


            {
                // Calculate analytical partials.
                accelerationModel->updateMembers( );
                accelerationPartial->update( );
                Eigen::MatrixXd partialWrtJupiterPosition = Eigen::Matrix3d::Zero( );
                accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtJupiterPosition.block( 0, 0, 3, 3 ) );
                Eigen::MatrixXd partialWrtJupiterVelocity = Eigen::Matrix3d::Zero( );
                accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtJupiterVelocity.block( 0, 0, 3, 3 ) );
                Eigen::MatrixXd partialWrtIoPosition = Eigen::Matrix3d::Zero( );
                accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtIoPosition.block( 0, 0, 3, 3 ) );
                Eigen::MatrixXd partialWrtIoVelocity = Eigen::Matrix3d::Zero( );
                accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtIoVelocity.block( 0, 0, 3, 3 ) );
                Eigen::MatrixXd partialWrtJupiterGravitationalParameter = accelerationPartial->wrtParameter(
                            jupiterGravitationalParameterParameter );
                Eigen::MatrixXd partialWrtIoGravitationalParameter = accelerationPartial->wrtParameter(
                            ioGravitationalParameterParameter );
                Eigen::MatrixXd partialWrtTidalTimeLag = accelerationPartial->wrtParameter(
                            tidalTimeLagParameter );

                // Declare perturbations in position for numerical partial
                Eigen::Vector3d positionPerturbation;
                positionPerturbation << 10.0, 10.0, 10.0;
                Eigen::Vector3d velocityPerturbation;
                velocityPerturbation << 1.0E-1, 1.0E-1, 1.0E-1;
                double jupiterGravityFieldPerturbation = 1.0E8;
                double ioGravityFieldPerturbation = 1.0E8;
                double timelagPerturbation = 1.0;;

                // Calculate numerical partials.
                Eigen::MatrixXd testPartialWrtIoPosition = calculateAccelerationWrtStatePartials(
                            ioStateSetFunction, accelerationModel, io->getState( ), positionPerturbation, 0 );
                Eigen::MatrixXd testPartialWrtIoVelocity = calculateAccelerationWrtStatePartials(
                            ioStateSetFunction, accelerationModel, io->getState( ), velocityPerturbation, 3 );
                Eigen::MatrixXd testPartialWrtJupiterPosition = calculateAccelerationWrtStatePartials(
                            jupiterStateSetFunction, accelerationModel, jupiter->getState( ), positionPerturbation, 0 );
                Eigen::MatrixXd testPartialWrtJupiterVelocity = calculateAccelerationWrtStatePartials(
                            jupiterStateSetFunction, accelerationModel, jupiter->getState( ), velocityPerturbation, 3 );
                Eigen::MatrixXd testPartialWrtJupiterGravitationalParameter = calculateAccelerationWrtParameterPartials(
                            jupiterGravitationalParameterParameter, accelerationModel, jupiterGravityFieldPerturbation );
                Eigen::MatrixXd testPartialWrtIoGravitationalParameter = calculateAccelerationWrtParameterPartials(
                            ioGravitationalParameterParameter, accelerationModel, ioGravityFieldPerturbation );
                Eigen::MatrixXd testPartialWrtTidalTimeLag = calculateAccelerationWrtParameterPartials(
                            tidalTimeLagParameter, accelerationModel, timelagPerturbation );

                // Compare numerical and analytical results.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtJupiterPosition,
                                                   partialWrtJupiterPosition, 1.0e-5 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtJupiterVelocity,
                                                   partialWrtJupiterVelocity, 1.0e-5 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtIoPosition,
                                                   partialWrtIoPosition, 1.0e-5 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtIoVelocity,
                                                   partialWrtIoVelocity, 1.0e-5 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtJupiterGravitationalParameter,
                                                   partialWrtJupiterGravitationalParameter, 1.0e-6 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtIoGravitationalParameter,
                                                   partialWrtIoGravitationalParameter, 1.0e-6 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtTidalTimeLag,
                                                   partialWrtTidalTimeLag, 1.0e-6 );

            }
        }
    }
}




BOOST_AUTO_TEST_CASE( testPanelledSurfaceRadiationPressureAccelerationPartials )
{
    // Create empty bodies, earth and sun.
    std::shared_ptr< Body > vehicle = std::make_shared< Body >( );
    vehicle->setConstantBodyMass( 400.0 );
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    std::shared_ptr< Body > sun = std::make_shared< Body >( );

    SystemOfBodies bodies;
    bodies.addBody( vehicle, "Vehicle" );
    bodies.addBody( earth, "Earth" );
    bodies.addBody( sun, "Sun" );

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    earth->setShapeModel( createBodyShapeModel( getDefaultBodyShapeSettings( "Earth", TUDAT_NAN, TUDAT_NAN ), "Earth" ) );
    // Set current state of sun and earth.
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "SSB", "J2000", "NONE", 1.0E6 ) );
    earth->setState( getBodyCartesianStateAtEpoch( "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );

    const auto tle = std::make_shared<Tle>(
        "1 08820U 76039  A 77047.52561960  .00000002 +00000-0 +00000-0 0  9994\n"
        "2 08820 109.8332 127.3884 0044194 201.3006 158.6132 06.38663945018402");
    const auto tleEphemeris = std::make_shared<TleEphemeris>("Earth", "J2000", tle);
    vehicle->setState( tleEphemeris->getCartesianState(tle->getEpoch( ) ) +
                           getBodyCartesianStateAtEpoch( "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );

    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > simpleRotationalEphemeris =
        std::make_shared< ephemerides::SimpleRotationalEphemeris >(
            spice_interface::computeRotationQuaternionBetweenFrames( "J2000" , "IAU_Earth", 0.0 ),
            2.0 * mathematical_constants::PI / 86400.0, 0.0,
            "J2000" , "IAU_Earth" );
    earth->setRotationalEphemeris( simpleRotationalEphemeris );
    earth->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );

    // Create links to set and get state functions of bodies.
    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
        std::bind( &Body::setState, sun, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
        std::bind( &Body::setState, earth, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
        std::bind( &Body::setState, vehicle, std::placeholders::_1 );
    std::function< Eigen::Vector6d( ) > sunStateGetFunction =
        std::bind( &Body::getState, sun );
    std::function< Eigen::Vector6d( ) > earthStateGetFunction =
        std::bind( &Body::getState, earth );
    std::function< Eigen::Vector6d( ) > vehicleStateGetFunction =
        std::bind( &Body::getState, vehicle );
    bodies.at( "Sun" )->setRadiationSourceModel(
        createRadiationSourceModel(
            getDefaultRadiationSourceModelSettings( "Sun", TUDAT_NAN, TUDAT_NAN ), "Sun", bodies ) );
    bodies.at( "Sun" )->getRadiationSourceModel( )->updateMembers( 1.0E7 );

    bodies.at( "Earth" )->setRadiationSourceModel(
        createRadiationSourceModel(
            getKnockeEarthRadiationPressureSettings(  ), "Earth", bodies ) );
    bodies.at( "Earth" )->getRadiationSourceModel( )->updateMembers( 1.0E7 );

    // Create radiation pressure properties of vehicle
    std::shared_ptr<CannonballRadiationPressureTargetModel > radiationPressureInterface =
        std::make_shared< CannonballRadiationPressureTargetModel >( mathematical_constants::PI * 0.3 * 0.3, 1.2 );
    vehicle->addRadiationPressureTargetModel( radiationPressureInterface );
    bodies.at( "Vehicle" )->getRadiationPressureTargetModel( )->updateMembers( 1.0E7 );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parametersNames;
    parametersNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
        "Vehicle", getBodyCartesianStateAtEpoch( "Earth", "SSB", "J2000", "NONE", 1.0E6 ), "SSB", "J2000" ) );
    parametersNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
        "Earth", tleEphemeris->getCartesianState(tle->getEpoch( ) ), "Earth", "J2000" ) );
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
        createParametersToEstimate< double >( parametersNames, bodies, nullptr );

    // Create acceleration model.
    std::shared_ptr< RadiationPressureAcceleration > accelerationModel =
        createRadiationPressureAccelerationModel( bodies.at( "Vehicle" ), bodies.at( "Earth" ), "Vehicle", "Earth", bodies );
    accelerationModel->updateMembers( 1.0E7 );

    // Create partial-calculating object.
    std::shared_ptr< AccelerationPartial > accelerationPartial =
        createAnalyticalAccelerationPartial( accelerationModel, std::make_pair( "Vehicle", vehicle ),
                                             std::make_pair( "Earth", earth ), bodies, parametersToEstimate );

    std::string vehicleName = "Vehicle";
    std::shared_ptr< EstimatableParameter< double > > radiationPressureCoefficient =
        std::make_shared< RadiationPressureCoefficient >( radiationPressureInterface, vehicleName );
    std::shared_ptr< EstimatableParameter< double > > parallelScalingFactor =
        std::make_shared< RadiationPressureScalingFactor >( accelerationModel, source_direction_radiation_pressure_scaling_factor, "Vehicle", "Earth" );
    std::shared_ptr< EstimatableParameter< double > > perpendicularScalingFactor =
        std::make_shared< RadiationPressureScalingFactor >( accelerationModel, source_perpendicular_direction_radiation_pressure_scaling_factor, "Vehicle", "Earth" );

    // Calculate analytical partials.
    double currentTime = 0.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::Vector3d partialWrtRadiationPressureCoefficient = accelerationPartial->wrtParameter(
        radiationPressureCoefficient );
    Eigen::MatrixXd partialWrtParallelScaling = accelerationPartial->wrtParameter(
        parallelScalingFactor );
    Eigen::MatrixXd partialWrtPerpendicularScaling = accelerationPartial->wrtParameter(
        perpendicularScalingFactor );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Vector3d testPartialWrtRadiationPressureCoefficient = Eigen::Vector3d::Zero( );
    Eigen::MatrixXd testPartialWrtParallelScaling = Eigen::MatrixXd::Zero( 3, 1 );
    Eigen::MatrixXd testPartialWrtPerpendicularScaling = Eigen::MatrixXd::Zero( 3, 1 );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 100.0, 100.0, 100.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Calculate numerical partials.
    std::function< void( ) > updateFunction =
        std::bind( &RadiationPressureTargetModel::updateMembers, radiationPressureInterface, 0.0 );
    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
        earthStateSetFunction, accelerationModel, earth->getState( ), positionPerturbation, 0, updateFunction );
    testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
        vehicleStateSetFunction, accelerationModel, vehicle->getState( ), positionPerturbation, 0, updateFunction );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
        earthStateSetFunction, accelerationModel, earth->getState( ),velocityPerturbation, 3, updateFunction );
    testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
        vehicleStateSetFunction, accelerationModel, vehicle->getState( ), velocityPerturbation, 3, updateFunction );
    testPartialWrtRadiationPressureCoefficient = calculateAccelerationWrtParameterPartials(
        radiationPressureCoefficient, accelerationModel, 1.0E-2, updateFunction );
    testPartialWrtParallelScaling = calculateAccelerationWrtParameterPartials(
        parallelScalingFactor, accelerationModel, 10.0, updateFunction );
    testPartialWrtPerpendicularScaling = calculateAccelerationWrtParameterPartials(
        perpendicularScalingFactor, accelerationModel, 10.0, updateFunction );

    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                       partialWrtEarthPosition, 1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                       partialWrtEarthVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                       partialWrtVehiclePosition, 1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                       partialWrtVehicleVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtRadiationPressureCoefficient,
                                       partialWrtRadiationPressureCoefficient, 1.0E-12 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtParallelScaling,
                                       partialWrtParallelScaling, 1.0E-13 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPerpendicularScaling,
                                       partialWrtPerpendicularScaling, 1.0E-13 );
}


BOOST_AUTO_TEST_CASE( testThrustPartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    for( unsigned int   i = 0; i < 2; i++ )
    {
        // Create empty bodies, earth and sun.
        std::shared_ptr< Body > vehicle = std::make_shared< Body >( );

        SystemOfBodies bodies;
        bodies.addBody( vehicle, "Vehicle" );

        double vehicleMass = 5.0E3;
        bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );

        Eigen::Vector3d thrustDirection;
        thrustDirection << -1.4, 2.4, 5.6;

        std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
                [=](const double){ return thrustDirection.normalized( ); };
        bodies.at( "Vehicle" )->setRotationalEphemeris(
                    createRotationModel(
                        std::make_shared< BodyFixedDirectionBasedRotationSettings >(
                            thrustDirectionFunction, "ECLIPJ2000", "VehicleFixed" ),
                        "Vehicle", bodies ) );

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        double thrustMagnitude1 = 1.0E3;
        double specificImpulse1 = 250.0;
        addEngineModel( "Vehicle", "Engine1",
                        std::make_shared< ConstantThrustMagnitudeSettings >(
                            thrustMagnitude1, specificImpulse1 ), bodies );

        double thrustMagnitude2 = 20.0E3;
        double specificImpulse2 = 50.0;
        addEngineModel( "Vehicle", "Engine2",
                        std::make_shared< ConstantThrustMagnitudeSettings >(
                            thrustMagnitude2, specificImpulse2 ), bodies, Eigen::Vector3d::UnitY( ) );

        // Create acceleration due to sun on earth.
        std::shared_ptr< ThrustAcceleration > thrustAcceleration;
        if( i == 0 )
        {
            thrustAcceleration = std::dynamic_pointer_cast< ThrustAcceleration >(
                    tudat::simulation_setup::createThrustAcceleratioModel(
                        std::make_shared< ThrustAccelerationSettings >( "Engine1" ), bodies, "Vehicle" ) );
        }
        else if( i == 1 )
        {
            thrustAcceleration = std::dynamic_pointer_cast< ThrustAcceleration >(
                    tudat::simulation_setup::createThrustAcceleratioModel(
                        std::make_shared< ThrustAccelerationSettings >(
                                std::vector< std::string >( { "Engine1", "Engine2" } ) ), bodies, "Vehicle" ) );
        }

        std::shared_ptr< EstimatableParameter< double > > constantThrustParameter1 = std::make_shared<
                ConstantThrustMagnitudeParameter >(
                    std::dynamic_pointer_cast< propulsion::ConstantThrustMagnitudeWrapper >(
                        vehicle->getVehicleSystems( )->getEngineModels( ).at( "Engine1" )->getThrustMagnitudeWrapper( ) ),
                    "Vehicle", "Engine1" );
        std::shared_ptr< EstimatableParameter< double > > constantThrustParameter2 = std::make_shared<
                ConstantThrustMagnitudeParameter >(
                    std::dynamic_pointer_cast< propulsion::ConstantThrustMagnitudeWrapper >(
                        vehicle->getVehicleSystems( )->getEngineModels( ).at( "Engine2" )->getThrustMagnitudeWrapper( ) ),
                    "Vehicle", "Engine2" );

        // Create central gravity partial.
        std::shared_ptr< ThrustAccelerationPartial > thrustPartial =
                std::dynamic_pointer_cast< ThrustAccelerationPartial >(
                    createAnalyticalAccelerationPartial( thrustAcceleration, std::make_pair( "Vehicle", vehicle ),
                                                         std::make_pair( "Vehicle", vehicle ), bodies ) );

        // Calculate analytical partials.
        thrustPartial->update( 0.0 );

        Eigen::MatrixXd partialWrtMass = Eigen::Vector3d::Zero( );
        thrustPartial->wrtBodyMass( partialWrtMass.block( 0, 0, 3, 1 ) );
        Eigen::Vector3d partialWrtEngine1Thrust = thrustPartial->wrtParameter(
                    constantThrustParameter1 );
        Eigen::Vector3d partialWrtEngine2Thrust = thrustPartial->wrtParameter(
                    constantThrustParameter2 );

        // Declare numerical partials.
        Eigen::Vector3d testPartialWrtMass = Eigen::Vector3d::Zero( );
        Eigen::Vector3d testPartialWrtEngine1Thrust = Eigen::Vector3d::Zero( );
        Eigen::Vector3d testPartialWrtEngine2Thrust = Eigen::Vector3d::Zero( );

        // Declare perturbations in position for numerical partial/
        double massPerturbation = 0.01;

        // Create state access/modification functions for bodies.
        std::function< void( double ) > massSetFunction =
                std::bind( &Body::setConstantBodyMass, vehicle, std::placeholders::_1 );

        // Calculate numerical partials.
        testPartialWrtMass = calculateAccelerationWrtMassPartials(
                    massSetFunction, thrustAcceleration, vehicleMass, massPerturbation );
        testPartialWrtEngine1Thrust = calculateAccelerationWrtParameterPartials(
                    constantThrustParameter1, thrustAcceleration, 1.0 );
        testPartialWrtEngine2Thrust = calculateAccelerationWrtParameterPartials(
                    constantThrustParameter2, thrustAcceleration, 1.0 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtMass, testPartialWrtMass, 1.0E-9 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEngine1Thrust, partialWrtEngine1Thrust, 1.0E-9 );

        BOOST_CHECK_SMALL( std::fabs( partialWrtEngine2Thrust( 2 ) ), ( 1.0E-14 * partialWrtEngine2Thrust.norm( ) ) );
        testPartialWrtEngine2Thrust( 2 ) = 0.0;
        partialWrtEngine2Thrust( 2 ) = 0.0;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEngine2Thrust, partialWrtEngine2Thrust, 1.0E-9 );

        bool isPartialChecked = false;
        if( i == 0 )
        {
            for( int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_EQUAL( partialWrtEngine2Thrust( j ), 0.0 );
                isPartialChecked = true;
            }
        }
        if( !isPartialChecked )
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtMass, testPartialWrtMass, 1.0E-4 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtEngine1Thrust, testPartialWrtEngine1Thrust, 1.0E-4 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtEngine2Thrust, testPartialWrtEngine2Thrust, 1.0E-4 );

        }
    }
}


BOOST_AUTO_TEST_CASE( testYarkovskyPartials )
{
    // Create empty bodies, earth and sun.
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    std::shared_ptr< Body > sun = std::make_shared< Body >( );

    SystemOfBodies bodies;
    bodies.addBody( earth, "Earth" );
    bodies.addBody( sun, "Sun" );

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set current state of sun and earth.
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "Sun", "J2000", "NONE", 1.0E6 ) );
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "Sun", "J2000", "NONE", 1.0E6 ) );

    double yarkovskyParameter = 1.0E-6;
    // Create acceleration due to sun on earth.
    std::shared_ptr< YarkovskyAcceleration > yarkovskyAccelerationModel =
        std::dynamic_pointer_cast< YarkovskyAcceleration >(
            createAccelerationModel(
            earth, sun, yarkovskyAcceleration( yarkovskyParameter ), "Earth", "Sun" ) );

    // Create central gravity partial.
    std::shared_ptr< AccelerationPartial > yarkovskyPartial =
        createAnalyticalAccelerationPartial( yarkovskyAccelerationModel, std::make_pair( "Earth", earth ),
                                             std::make_pair( "Sun", sun ), bodies );

    // Create gravitational parameter object.
    std::shared_ptr< EstimatableParameter< double > > yarkovskyParameterParameter = std::make_shared<
        YarkovskyParameter >( yarkovskyAccelerationModel, "Earth", "Sun" );

    // Calculate analytical partials.
    yarkovskyPartial->update( 0.0 );
    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    yarkovskyPartial->wrtPositionOfAcceleratedBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    yarkovskyPartial->wrtVelocityOfAcceleratedBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
    yarkovskyPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    yarkovskyPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::Vector3d partialWrtSunYarkovskyParameter = yarkovskyPartial->wrtParameter(
        yarkovskyParameterParameter );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10000.0, 10000.0, 10000.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
        std::bind( &Body::setState, earth, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
        std::bind( &Body::setState, sun, std::placeholders::_1 );
    std::function< Eigen::Vector6d ( ) > earthStateGetFunction =
        std::bind( &Body::getState, earth );
    std::function< Eigen::Vector6d ( ) > sunStateGetFunction =
        std::bind( &Body::getState, sun );

    // Calculate numerical partials.
    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
        earthStateSetFunction, yarkovskyAccelerationModel, earth->getState( ), positionPerturbation, 0 );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
        earthStateSetFunction, yarkovskyAccelerationModel, earth->getState( ), velocityPerturbation, 3 );
    testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
        sunStateSetFunction, yarkovskyAccelerationModel, sun->getState( ), positionPerturbation, 0 );
    testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
        sunStateSetFunction, yarkovskyAccelerationModel, sun->getState( ), velocityPerturbation, 3 );
    Eigen::Vector3d testPartialWrtSunYarkovskyParameter = calculateAccelerationWrtParameterPartials(
        yarkovskyParameterParameter, yarkovskyAccelerationModel, 1.0E-6 );

    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                       partialWrtEarthPosition, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                       partialWrtEarthVelocity, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                       partialWrtSunPosition, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                       partialWrtSunVelocity, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunYarkovskyParameter,
                                       partialWrtSunYarkovskyParameter, 1.0E-8 );
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
