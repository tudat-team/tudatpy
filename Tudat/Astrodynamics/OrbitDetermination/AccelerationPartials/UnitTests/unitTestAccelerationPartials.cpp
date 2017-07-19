/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/exponentialAtmosphere.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantDragCoefficient.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/empiricalAccelerationCoefficients.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/gravitationalParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/radiationPressureCoefficient.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/ppnParameters.h"
#include "Tudat/Astrodynamics/Relativity/metric.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/numericalAccelerationPartial.h"
#include "Tudat/Astrodynamics/Relativity/relativisticAccelerationCorrection.h"
#include "Tudat/SimulationSetup/EstimationSetup/createAccelerationPartials.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"

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
using namespace tudat::electro_magnetism;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_acceleration_partials )

BOOST_AUTO_TEST_CASE( testCentralGravityPartials )
{
    // Create empty bodies, earth and sun.
    boost::shared_ptr< Body > earth = boost::make_shared< Body >( );
    boost::shared_ptr< Body > sun = boost::make_shared< Body >( );

    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = earth;
    bodyMap[ "Sun" ] = sun;

    // Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    // Set current state of sun and earth.
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "Sun", "J2000", "NONE", 1.0E6 ) );
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "Sun", "J2000", "NONE", 1.0E6 ) );

    // Get sun gravitational parameter and set gravity field model.
    double sunsGravitationalParameter = getBodyGravitationalParameter( "Sun" );
    boost::shared_ptr< GravityFieldModel > sunGravityFieldModel =
            boost::make_shared< GravityFieldModel >( sunsGravitationalParameter );
    sun->setGravityFieldModel( sunGravityFieldModel );
    double earthGravitationalParameter = getBodyGravitationalParameter( "Earth" );
    boost::shared_ptr< GravityFieldModel > earthGravityFieldModel =
            boost::make_shared< GravityFieldModel >( earthGravitationalParameter );
    earth->setGravityFieldModel( earthGravityFieldModel );

    // Create acceleration due to sun on earth.
    boost::shared_ptr< CentralGravitationalAccelerationModel3d > gravitationalAcceleration =\
            createCentralGravityAcceleratioModel( earth, sun, "Earth", "Sun", 1 );

    // Create central gravity partial.
    boost::shared_ptr< AccelerationPartial > centralGravitationPartial =
            createAnalyticalAccelerationPartial( gravitationalAcceleration, std::make_pair( "Earth", earth ),
                                                 std::make_pair( "Sun", sun ), bodyMap );

    // Create gravitational parameter object.
    boost::shared_ptr< EstimatableParameter< double > > sunGravitationalParameterParameter = boost::make_shared<
            GravitationalParameter >( sunGravityFieldModel, "Sun" );
    boost::shared_ptr< EstimatableParameter< double > > earthGravitationalParameterParameter = boost::make_shared<
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
    boost::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            boost::bind( &Body::setState, earth, _1 );
    boost::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            boost::bind( &Body::setState, sun, _1 );
    boost::function< Eigen::Vector6d ( ) > earthStateGetFunction =
            boost::bind( &Body::getState, earth );
    boost::function< Eigen::Vector6d ( ) > sunStateGetFunction =
            boost::bind( &Body::getState, sun );

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

BOOST_AUTO_TEST_CASE( testRadiationPressureAccelerationPartials )
{
    // Create empty bodies, earth and sun.
    boost::shared_ptr< Body > vehicle = boost::make_shared< Body >( );
    vehicle->setConstantBodyMass( 400.0 );
    boost::shared_ptr< Body > sun = boost::make_shared< Body >( );

    NamedBodyMap bodyMap;
    bodyMap[ "Vehicle" ] = vehicle;
    bodyMap[ "Sun" ] = sun;

    // Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    // Set current state of sun and earth.
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "SSB", "J2000", "NONE", 1.0E6 ) );
    vehicle->setState( getBodyCartesianStateAtEpoch(  "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );

    // Create links to set and get state functions of bodies.
    boost::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            boost::bind( &Body::setState, sun, _1 );
    boost::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            boost::bind( &Body::setState, vehicle, _1 );
    boost::function< Eigen::Vector6d( ) > sunStateGetFunction =
            boost::bind( &Body::getState, sun );
    boost::function< Eigen::Vector6d( ) > vehicleStateGetFunction =
            boost::bind( &Body::getState, vehicle );

    // Create radiation pressure properties of vehicle
    boost::shared_ptr< RadiationPressureInterface > radiationPressureInterface =
            createRadiationPressureInterface( boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                                                  "Sun", mathematical_constants::PI * 0.3 * 0.3, 1.2 ), "Vehicle", bodyMap );
    radiationPressureInterface->updateInterface( 0.0 );
    vehicle->setRadiationPressureInterface( "Sun", radiationPressureInterface );

    // Create acceleration model.
    boost::shared_ptr< CannonBallRadiationPressureAcceleration > accelerationModel =
            boost::make_shared< CannonBallRadiationPressureAcceleration >(
                boost::bind( &Body::getPosition, sun ),
                boost::bind( &Body::getPosition, vehicle ),
                boost::bind( &RadiationPressureInterface::getCurrentRadiationPressure,
                             radiationPressureInterface ),
                boost::bind( &RadiationPressureInterface::getRadiationPressureCoefficient, radiationPressureInterface ),
                boost::bind( &RadiationPressureInterface::getArea, radiationPressureInterface ),
                boost::bind( &Body::getBodyMass, vehicle ) );

    // Create partial-calculating object.
    boost::shared_ptr< AccelerationPartial > accelerationPartial =
            createAnalyticalAccelerationPartial( accelerationModel, std::make_pair( "Vehicle", vehicle ),
                                                 std::make_pair( "Sun", sun ), bodyMap );

    // Create parameter object
    std::string vehicleName = "Vehicle";
    boost::shared_ptr< EstimatableParameter< double > > radiationPressureCoefficient =
            boost::make_shared< RadiationPressureCoefficient >( radiationPressureInterface, vehicleName );

    // Calculate analytical partials.
    accelerationPartial->update( 0.0 );
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

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Vector3d testPartialWrtRadiationPressureCoefficient = Eigen::Vector3d::Zero( );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10000.0, 10000.0, 10000.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Calculate numerical partials.
    boost::function< void( ) > updateFunction =
            boost::bind( &RadiationPressureInterface::updateInterface, radiationPressureInterface, 0.0 );
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
}

BOOST_AUTO_TEST_CASE( testThirdBodyGravityPartials )
{
    // Create empty bodies, earth and sun.
    boost::shared_ptr< Body > earth = boost::make_shared< Body >( );
    boost::shared_ptr< Body > sun = boost::make_shared< Body >( );
    boost::shared_ptr< Body > moon = boost::make_shared< Body >( );

    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = earth;
    bodyMap[ "Sun" ] = sun;
    bodyMap[ "Moon" ] = moon;

    // Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    // Set current state of sun and earth.
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "Sun", "J2000", "NONE", 1.0E6 ) );
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "Sun", "J2000", "NONE", 1.0E6 ) );
    moon->setState( getBodyCartesianStateAtEpoch(  "Moon", "Sun", "J2000", "NONE", 1.0E6 ) );

    // Get sun gravitational parameter and set gravity field model.
    double sunsGravitationalParameter = getBodyGravitationalParameter( "Sun" );
    boost::shared_ptr< GravityFieldModel > sunGravityFieldModel =
            boost::make_shared< GravityFieldModel >( sunsGravitationalParameter );
    sun->setGravityFieldModel( sunGravityFieldModel );

    double moonsGravitationalParameter = getBodyGravitationalParameter( "Moon" );
    boost::shared_ptr< GravityFieldModel > moonGravityFieldModel =
            boost::make_shared< GravityFieldModel >( moonsGravitationalParameter );
    moon->setGravityFieldModel( moonGravityFieldModel );

    double earthGravitationalParameter = getBodyGravitationalParameter( "Earth" );
    boost::shared_ptr< GravityFieldModel > earthGravityFieldModel =
            boost::make_shared< GravityFieldModel >( earthGravitationalParameter );
    earth->setGravityFieldModel( earthGravityFieldModel );

    // Create acceleration due to moon on earth.
    boost::shared_ptr< ThirdBodyCentralGravityAcceleration > gravitationalAcceleration =
            createThirdBodyCentralGravityAccelerationModel(
                moon, sun, earth, "Moon", "Sun", "Earth" );

    // Create central gravity partial.
    boost::shared_ptr< AccelerationPartial > thirdBodyGravitationPartial =
            createAnalyticalAccelerationPartial( gravitationalAcceleration, std::make_pair( "Moon", moon ),
                                                 std::make_pair( "Sun", sun ), bodyMap );

    // Create gravitational parameter object.
    boost::shared_ptr< EstimatableParameter< double > > gravitationalParameterParameter = boost::make_shared<
            GravitationalParameter >( sunGravityFieldModel, "Sun" );
    boost::shared_ptr< EstimatableParameter< double > > moonGravitationalParameterParameter = boost::make_shared<
            GravitationalParameter >( moonGravityFieldModel, "Moon" );
    boost::shared_ptr< EstimatableParameter< double > > earthGravitationalParameterParameter = boost::make_shared<
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

    Eigen::Vector3d partialWrtSunGravitationalParameter = thirdBodyGravitationPartial->wrtParameter(
                gravitationalParameterParameter );
    Eigen::Vector3d partialWrtMoonGravitationalParameter = thirdBodyGravitationPartial->wrtParameter(
                moonGravitationalParameterParameter );
    Eigen::Vector3d partialWrtEarthGravitationalParameter = thirdBodyGravitationPartial->wrtParameter(
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
    boost::function< void( Eigen::Vector6d ) > moonStateSetFunction =
            boost::bind( &Body::setState, moon, _1 );
    boost::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            boost::bind( &Body::setState, sun, _1 );
    boost::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            boost::bind( &Body::setState, earth, _1 );

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
        const boost::shared_ptr< aerodynamics::FlightConditions > flightConditions,
        const double timeToUpdate )
{
    flightConditions->resetCurrentTime( TUDAT_NAN );
    flightConditions->updateConditions( timeToUpdate );
}

BOOST_AUTO_TEST_CASE( testAerodynamicAccelerationPartials )
{

    //Load spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    using namespace tudat;
    // Create Earth object
    std::map< std::string, boost::shared_ptr< BodySettings > > defaultBodySettings =
            getDefaultBodySettings( boost::assign::list_of( "Earth" ) );
    defaultBodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ) );
    NamedBodyMap bodyMap = createBodies( defaultBodySettings );

    // Create vehicle objects.
    double vehicleMass = 5.0E3;
    bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );


    bool areCoefficientsInAerodynamicFrame = 1;
    Eigen::Vector3d aerodynamicCoefficients = ( Eigen::Vector3d( )<< 2.5, -0.1, 0.5 ).finished( );

    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                2.0, 4.0, 1.5, Eigen::Vector3d::Zero( ), aerodynamicCoefficients, Eigen::Vector3d::Zero( ),
                areCoefficientsInAerodynamicFrame, 1 );
    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


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
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                vehicleSphericalEntryState );

    bodyMap.at( "Earth" )->setStateFromEphemeris( 0.0 );
    bodyMap.at( "Vehicle" )->setState( systemInitialState );


    boost::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationModel =
            simulation_setup::createAerodynamicAcceleratioModel(
                bodyMap[ "Vehicle" ], bodyMap[ "Earth" ], "Vehicle", "Earth" );
    bodyMap.at( "Vehicle" )->getFlightConditions( )->updateConditions( 0.0 );
    accelerationModel->updateMembers( 0.0 );

    boost::shared_ptr< AccelerationPartial > aerodynamicAccelerationPartial =
            createAnalyticalAccelerationPartial(
                accelerationModel, std::make_pair( "Vehicle", bodyMap[ "Vehicle" ] ),
            std::make_pair( "Earth", bodyMap[ "Earth" ] ), bodyMap );

    // Create gravitational parameter object.
    boost::shared_ptr< EstimatableParameter< double > > dragCoefficientParameter = boost::make_shared<
            ConstantDragCoefficient >( boost::dynamic_pointer_cast< aerodynamics::CustomAerodynamicCoefficientInterface >(
                                           bodyMap[ "Vehicle" ]->getAerodynamicCoefficientInterface( ) ), "Vehicle" );

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

    boost::function< void( ) > environmentUpdateFunction =
            boost::bind( &updateFlightConditionsWithPerturbedState, bodyMap.at( "Vehicle" )->getFlightConditions( ), 0.0 );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1.0, 1.0, 1.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    // Create state access/modification functions for bodies.
    boost::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            boost::bind( &Body::setState, bodyMap.at( "Vehicle" ), _1 );
    boost::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            boost::bind( &Body::setState, bodyMap.at( "Earth" ), _1 );
    boost::function< Eigen::Vector6d ( ) > vehicleStateGetFunction =
            boost::bind( &Body::getState, bodyMap.at( "Vehicle" ) );
    boost::function< Eigen::Vector6d ( ) > earthStateGetFunction =
            boost::bind( &Body::getState, bodyMap.at( "Earth" ) );

    // Calculate numerical partials.
    testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, bodyMap.at( "Vehicle" )->getState( ), positionPerturbation, 0,
                environmentUpdateFunction);
    testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, bodyMap.at( "Vehicle" )->getState( ), velocityPerturbation, 3,
                environmentUpdateFunction );
    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, accelerationModel, bodyMap.at( "Earth" )->getState( ), positionPerturbation, 0,
                environmentUpdateFunction );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, accelerationModel, bodyMap.at( "Earth" )->getState( ), velocityPerturbation, 3,
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
    boost::shared_ptr< Body > earth = boost::make_shared< Body >( );
    boost::shared_ptr< Body > vehicle = boost::make_shared< Body >( );

    // Create links to set and get state functions of bodies.
    boost::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            boost::bind( &Body::setState, earth, _1  );
    boost::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            boost::bind( &Body::setState, vehicle, _1  );
    boost::function< Eigen::Vector6d( ) > earthStateGetFunction =
            boost::bind( &Body::getState, earth );
    boost::function< Eigen::Vector6d( ) > vehicleStateGetFunction =
            boost::bind( &Body::getState, vehicle );

    // Load spice kernel.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    // Set vehicle and earth state.
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );
    Eigen::Vector6d vehicleKeplerElements;
    vehicleKeplerElements<<6378.0E3 + 249E3, 0.0004318, convertDegreesToRadians( 96.5975 ),
            convertDegreesToRadians( 217.6968 ), convertDegreesToRadians( 268.2663 ), convertDegreesToRadians( 142.3958 );
    vehicle->setState( earth->getState( ) + convertKeplerianToCartesianElements( vehicleKeplerElements,
                                                                                 getBodyGravitationalParameter( "Earth" ) ) );


    NamedBodyMap bodyMap;
    bodyMap[ "Vehicle" ] = vehicle;
    bodyMap[ "Earth" ] = earth;

    // Create gravity field.
    boost::shared_ptr< GravityFieldSettings > gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );
    boost::shared_ptr< gravitation::GravityFieldModel > earthGravityField =
            createGravityFieldModel( gravityFieldSettings, "Earth", bodyMap );
    earth->setGravityFieldModel( earthGravityField );

    // Create acceleration model.
    boost::function< double( ) > ppnParameterGammaFunction = boost::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet );
    boost::function< double( ) > ppnParameterBetaFunction = boost::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet );
    boost::shared_ptr< RelativisticAccelerationCorrection > accelerationModel =
            boost::make_shared< RelativisticAccelerationCorrection >
            ( boost::bind( &Body::getState, vehicle ),
              boost::bind( &Body::getState, earth ),
              boost::bind( &GravityFieldModel::getGravitationalParameter, earthGravityField ),
              ppnParameterGammaFunction, ppnParameterBetaFunction );

    // Create acceleration partial object.
    boost::shared_ptr< RelativisticAccelerationPartial > accelerationPartial = boost::make_shared< RelativisticAccelerationPartial >(
                accelerationModel, "Vehicle", "Earth" );

    // Create parameter objects.
    boost::shared_ptr< EstimatableParameter< double > > gravitationalParameterParameter = boost::make_shared<
            GravitationalParameter >( earthGravityField, "Earth" );
    boost::shared_ptr< EstimatableParameter< double > > ppnParameterGamma = boost::make_shared<
            PPNParameterGamma >( ppnParameterSet );
    boost::shared_ptr< EstimatableParameter< double > > ppnParameterBeta = boost::make_shared<
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
    positionPerturbation<< 10.0, 10.0, 10.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation<< 1.0, 1.0, 1.0;

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
    boost::shared_ptr< Body > earth = boost::make_shared< Body >( );
    boost::shared_ptr< Body > vehicle = boost::make_shared< Body >( );

    // Create links to set and get state functions of bodies.
    boost::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            boost::bind( &Body::setState, earth, _1  );
    boost::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            boost::bind( &Body::setState, vehicle, _1  );
    boost::function< Eigen::Vector6d( ) > earthStateGetFunction =
            boost::bind( &Body::getState, earth );
    boost::function< Eigen::Vector6d( ) > vehicleStateGetFunction =
            boost::bind( &Body::getState, vehicle );

    // Load spice kernel.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    // Set vehicle and earth state.
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );
    Eigen::Vector6d vehicleKeplerElements;
    vehicleKeplerElements<<6378.0E3 + 249E3, 0.0004318, convertDegreesToRadians( 96.5975 ),
            convertDegreesToRadians( 217.6968 ), convertDegreesToRadians( 268.2663 ), convertDegreesToRadians( 142.3958 );
    vehicle->setState( earth->getState( ) + convertKeplerianToCartesianElements( vehicleKeplerElements,
                                                                                 getBodyGravitationalParameter( "Earth" ) ) );


    NamedBodyMap bodyMap;
    bodyMap[ "Vehicle" ] = vehicle;
    bodyMap[ "Earth" ] = earth;

    // Create gravity field.
    boost::shared_ptr< GravityFieldSettings > gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );
    boost::shared_ptr< gravitation::GravityFieldModel > earthGravityField =
            createGravityFieldModel( gravityFieldSettings, "Earth", bodyMap );
    earth->setGravityFieldModel( earthGravityField );

    // Create rotation model
    boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > simpleRotationalEphemeris =
            boost::make_shared< ephemerides::SimpleRotationalEphemeris >(
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
    boost::shared_ptr< EmpiricalAcceleration > accelerationModel =
            boost::make_shared< EmpiricalAcceleration >
            ( constantAcceleration, sineAcceleration, cosineAcceleration,
              boost::bind( &Body::getState, vehicle ),
              boost::bind( &GravityFieldModel::getGravitationalParameter, earthGravityField ),
              boost::bind( &Body::getState, earth ) );

    // Create acceleration partial object.
    boost::shared_ptr< EmpiricalAccelerationPartial > accelerationPartial = boost::make_shared< EmpiricalAccelerationPartial >(
                accelerationModel, "Vehicle", "Earth" );

    // Define list of empirical accelerations w.r.t. which partials are to be computed
    std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > allEmpiricalShapesVector;
    allEmpiricalShapesVector.push_back( basic_astrodynamics::constant_empirical );
    allEmpiricalShapesVector.push_back( basic_astrodynamics::cosine_empirical );

    std::map< basic_astrodynamics::EmpiricalAccelerationComponents, std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > >
            empiricalComponentsToEstimate;
    empiricalComponentsToEstimate[ basic_astrodynamics::radial_empicial_acceleration_component ] = allEmpiricalShapesVector;
    empiricalComponentsToEstimate[ basic_astrodynamics::along_track_empicial_acceleration_component ] = allEmpiricalShapesVector;

    allEmpiricalShapesVector.push_back( basic_astrodynamics::sine_empirical );
    empiricalComponentsToEstimate[ basic_astrodynamics::across_track_empicial_acceleration_component ] = allEmpiricalShapesVector;

    // Create time-independent empirical acceleration object.
    boost::shared_ptr< EmpiricalAccelerationCoefficientsParameter > empiricalAccelerationParameter = boost::make_shared<
            EmpiricalAccelerationCoefficientsParameter >( accelerationModel, "Vehicle", empiricalComponentsToEstimate );

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
        positionPerturbation<< 1.0, 1.0, 1.0;
        Eigen::Vector3d velocityPerturbation;
        velocityPerturbation<< 1.0E-3, 1.0E-3, 1.0E-3;
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
    boost::shared_ptr< ArcWiseEmpiricalAccelerationCoefficientsParameter > arcWiseEmpiricalAccelerationParameter = boost::make_shared<
            ArcWiseEmpiricalAccelerationCoefficientsParameter >( accelerationModel, "Vehicle", empiricalComponentsToEstimate, arcStartTimes );

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
        accelerationModel->resetTime( TUDAT_NAN );
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
        positionPerturbation<< 1.0, 1.0, 1.0;
        Eigen::Vector3d velocityPerturbation;
        velocityPerturbation<< 1.0E-3, 1.0E-3, 1.0E-3;
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

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




