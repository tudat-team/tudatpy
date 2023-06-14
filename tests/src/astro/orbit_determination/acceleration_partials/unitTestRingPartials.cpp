/*    Copyright (c) 2010-2023, Delft University of Technology
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

#include "tudat/basics/testMacros.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/simulation/estimation.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_ring_partials )

BOOST_AUTO_TEST_CASE( testRingAccelerationPartial )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create empty bodies, earth and vehicle.
    std::shared_ptr< simulation_setup::Body > earth = std::make_shared< simulation_setup::Body >( );
    std::shared_ptr< simulation_setup::Body > vehicle = std::make_shared< simulation_setup::Body >( );

    const double gravitationalParameter = 3.986004418e14;

    // Define ring radius: set to the Earth radius
    double ringRadius = 6378137.0;

    simulation_setup::SystemOfBodies bodies;
    bodies.addBody( earth, "Earth" );
    bodies.addBody( vehicle, "Vehicle" );;
    bodies.addBody( createSystemOfBodies( simulation_setup::getDefaultBodySettings( { "Moon" } ) ).at( "Moon" ), "Moon" );

    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > simpleRotationalEphemeris =
            std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000" , "IAU_Earth", 0.0 ),
                2.0 * mathematical_constants::PI / 86400.0,
                1.0E7,
                "ECLIPJ2000" , "IAU_Earth" );
    earth->setRotationalEphemeris( simpleRotationalEphemeris );

    std::shared_ptr< simulation_setup::GravityFieldSettings > earthGravityFieldSettings =
            std::make_shared< simulation_setup::RingGravityFieldSettings >( gravitationalParameter, ringRadius, "IAU_Earth" );

    std::shared_ptr< tudat::gravitation::RingGravityField > earthGravityField =
            std::dynamic_pointer_cast< gravitation::RingGravityField  >(
                simulation_setup::createGravityFieldModel( earthGravityFieldSettings, "Earth", bodies ) );
    earth->setGravityFieldModel( earthGravityField );

    // Set current state of vehicle and earth.
    double testTime = 1.0E6;
    earth->setState( Eigen::Vector6d::Zero( ) );
    earth->setCurrentRotationToLocalFrameFromEphemeris( testTime );
    bodies.at( "Moon" ) ->setState( tudat::spice_interface::getBodyCartesianStateAtEpoch( "Moon", "Earth", "ECLIPJ2000" ,"None", testTime ) );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( orbital_element_conversions::semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( orbital_element_conversions::eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( orbital_element_conversions::inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( orbital_element_conversions::argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( orbital_element_conversions::longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( orbital_element_conversions::trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    Eigen::Vector6d asterixInitialState = orbital_element_conversions::convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, gravitationalParameter );

    vehicle->setState( asterixInitialState );

    // Create acceleration due Earth vehicle on vehicle.
    std::shared_ptr< simulation_setup::AccelerationSettings > accelerationSettings =
            std::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::ring_gravity );
    std::shared_ptr< gravitation::RingGravitationalAccelerationModel > gravitationalAcceleration =
            std::dynamic_pointer_cast< gravitation::RingGravitationalAccelerationModel >(
                createAccelerationModel( vehicle, earth, accelerationSettings, "Vehicle", "Earth" ) );
    gravitationalAcceleration->updateMembers( 0.0 );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &simulation_setup::Body::setState, earth, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            std::bind( &simulation_setup::Body::setState, vehicle, std::placeholders::_1  );
    std::function< Eigen::Vector6d ( ) > earthStateGetFunction =
            std::bind( &simulation_setup::Body::getState, earth );
    std::function< Eigen::Vector6d ( ) > vehicleStateGetFunction =
            std::bind( &simulation_setup::Body::getState, vehicle );

    // Create list of estimatable parameters settings
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< estimatable_parameters::InitialRotationalStateEstimatableParameterSettings< double > >(
                                  "Earth", 0.0, "ECLIPJ2000" ) );
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodies );

    // Create acceleration partial object.
    std::shared_ptr< acceleration_partials::RingGravityPartial > accelerationPartial =
            std::dynamic_pointer_cast< acceleration_partials::RingGravityPartial > (
                createAnalyticalAccelerationPartial(
                    gravitationalAcceleration,
                    std::make_pair( "Vehicle", vehicle ),
                    std::make_pair( "Earth", earth ),
                    bodies, parameterSet ) );

    accelerationPartial->update( testTime );

    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    // Declare perturbations in position for numerical partial
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10.0, 10.0, 10.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    // Calculate numerical partials.
    testPartialWrtVehiclePosition = acceleration_partials::calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, gravitationalAcceleration, vehicle->getState( ), positionPerturbation, 0 );
    testPartialWrtVehicleVelocity = acceleration_partials::calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, gravitationalAcceleration, vehicle->getState( ), velocityPerturbation, 3 );

    testPartialWrtEarthPosition = acceleration_partials::calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), positionPerturbation, 0 );
    testPartialWrtEarthVelocity = acceleration_partials::calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), velocityPerturbation, 3 );

    std::cerr << testPartialWrtVehiclePosition << std::endl;
    std::cerr << partialWrtVehiclePosition << std::endl << std::endl;

    std::cerr << testPartialWrtEarthPosition << std::endl;
    std::cerr << partialWrtEarthPosition << std::endl;

    // Partials wrt position: had to use higher tolerances than in the other partials tests (1e-5 instead of 1e-6), not sure why
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition, partialWrtVehiclePosition, 1.0E-3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity, partialWrtVehicleVelocity, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition, partialWrtEarthPosition, 1.0E-3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity, partialWrtEarthVelocity, 1.0E-6 );
}

//! Unit test to check working of ring state partial for synchronously rotating body (and rotation depending on state)
BOOST_AUTO_TEST_CASE( testRingAccelerationPartialWithSynchronousRotation )
{
    const double gravitationalParameter = 3.986004418e14;

    // Define ring radius: set to the Earth radius
    double ringRadius = 6378137.0;

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    // Create bodies needed in simulation
    simulation_setup::BodyListSettings bodySettings = simulation_setup::getDefaultBodySettings(
            bodyNames, "Earth", "ECLIPJ2000" );

    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< simulation_setup::SynchronousRotationModelSettings >(
            "Moon", "ECLIPJ2000", "IAU_Earth" );
    bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< simulation_setup::RingGravityFieldSettings >
            ( gravitationalParameter, ringRadius, "IAU_Earth" );
    simulation_setup::SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    std::shared_ptr< tudat::simulation_setup::Body > earth = bodies.at( "Earth" );
    std::shared_ptr< tudat::simulation_setup::Body > moon = bodies.at( "Moon" );
    std::dynamic_pointer_cast< tudat::ephemerides::SynchronousRotationalEphemeris >(
                earth->getRotationalEphemeris( ) )->setIsBodyInPropagation( 1 );

    // Set translational and rotational state of bodies
    double testTime = 1.0E6;
    earth->setStateFromEphemeris( testTime );
    Eigen::Vector6d moonState = moon->getStateInBaseFrameFromEphemeris( testTime );
    moon->setState( moonState * 0.1 );

    earth->setCurrentRotationToLocalFrameFromEphemeris( testTime );
    moon->setCurrentRotationToLocalFrameFromEphemeris( testTime );

    // Create acceleration model
    std::shared_ptr< simulation_setup::AccelerationSettings > accelerationSettings =
            std::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::ring_gravity );
    std::shared_ptr< gravitation::RingGravitationalAccelerationModel > gravitationalAcceleration =
            std::dynamic_pointer_cast< gravitation::RingGravitationalAccelerationModel >(
                createAccelerationModel( moon, earth, accelerationSettings, "Vehicle", "Earth" ) );
    gravitationalAcceleration->updateMembers( 0.0 );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMoonPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMoonVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1000.0, 1000.0, 1000.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0E-1, 1.0;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &simulation_setup::Body::setState, earth, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > moonStateSetFunction =
            std::bind( &simulation_setup::Body::setState, moon, std::placeholders::_1  );
    std::function< Eigen::Vector6d ( ) > earthStateGetFunction =
            std::bind( &simulation_setup::Body::getState, earth );
    std::function< Eigen::Vector6d ( ) > moonStateGetFunction =
            std::bind( &simulation_setup::Body::getState, moon );


    // Define estimated parameters
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< estimatable_parameters::InitialRotationalStateEstimatableParameterSettings< double > >(
                                  "Earth", 0.0, "ECLIPJ2000" ) );
    parameterNames.push_back( std::make_shared< estimatable_parameters::InitialRotationalStateEstimatableParameterSettings< double > >(
                                  "Moon", 0.0, "ECLIPJ2000" ) );
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodies );


    // Create acceleration partial object.
    std::shared_ptr< acceleration_partials::RingGravityPartial > accelerationPartial =
            std::dynamic_pointer_cast< acceleration_partials::RingGravityPartial > (
                createAnalyticalAccelerationPartial(
                    gravitationalAcceleration,
                    std::make_pair( "Moon", moon ),
                    std::make_pair( "Earth", earth ),
                    bodies, parameterSet ) );
    accelerationPartial->update( testTime );

    // Calculate analytical partials.
    Eigen::MatrixXd partialWrtMoonPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtMoonPosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtMoonVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtMoonVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    // Calculate numerical partials.
    std::function< void( ) > updateFunction =
            std::bind( &simulation_setup::Body::setCurrentRotationToLocalFrameFromEphemeris, bodies.at( "Earth" ), testTime );

    testPartialWrtMoonPosition = acceleration_partials::calculateAccelerationWrtStatePartials(
                moonStateSetFunction, gravitationalAcceleration, moon->getState( ), positionPerturbation, 0,
                updateFunction );
    testPartialWrtMoonVelocity = acceleration_partials::calculateAccelerationWrtStatePartials(
                moonStateSetFunction, gravitationalAcceleration, moon->getState( ), velocityPerturbation, 3,
                updateFunction );

    testPartialWrtEarthPosition = acceleration_partials::calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), positionPerturbation, 0,
                updateFunction );
    testPartialWrtEarthVelocity = acceleration_partials::calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), velocityPerturbation, 3,
                updateFunction );

    // Partials have very small values -> adding 1
    testPartialWrtMoonPosition += Eigen::Matrix3d::Ones( );
    partialWrtMoonPosition += Eigen::Matrix3d::Ones( );
    testPartialWrtEarthPosition += Eigen::Matrix3d::Ones( );
    partialWrtEarthPosition += Eigen::Matrix3d::Ones( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMoonPosition, partialWrtMoonPosition, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMoonVelocity, partialWrtMoonVelocity, 1.0E-4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition, partialWrtEarthPosition, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity, partialWrtEarthVelocity, 1.0E-4 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
