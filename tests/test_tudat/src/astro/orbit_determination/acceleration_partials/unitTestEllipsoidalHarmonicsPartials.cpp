/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          Garmier, R., Barriot, J.-P. (2001), Celestial Mechanics and Dynamical Astronomy, 79(4), 235-275.
 *          Hu, X., Jekeli, C. (2015), Journal of Geodesy, 89(2), 159-177.
 *          See tudat/math/basic/ellipsoidalHarmonics.h for the full citations.
 */

#define BOOST_TEST_MAIN

#include <string>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/estimation_setup/createAccelerationPartials.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_ellipsoidal_harmonics_partials )

std::string getDataDirectory( )
{
    return paths::getTudatTestDataPath( ) + "/ellipsoidal_harmonics/";
}

const double testTime = 1.0E3;

//! Comet with the degree 15 field of 67P and the given rotation, and a vehicle well outside the reference ellipsoid.
SystemOfBodies createTestBodies( const std::shared_ptr< RotationModelSettings > rotationModelSettings )
{
    std::shared_ptr< Body > comet = std::make_shared< Body >( );
    std::shared_ptr< Body > vehicle = std::make_shared< Body >( );
    vehicle->setConstantBodyMass( 400.0 );

    SystemOfBodies bodies;
    bodies.addBody( comet, "Comet" );
    bodies.addBody( vehicle, "Vehicle" );

    comet->setRotationalEphemeris( createRotationModel( rotationModelSettings, "Comet" ) );
    comet->setGravityFieldModel( createGravityFieldModel(
            ellipsoidalHarmonicsGravitySettingsFromFiles(
                    ( Eigen::Vector3d( ) << 2.802437330193671E+03, 2.090307360485069E+03, 1.867442428762624E+03 ).finished( ),
                    getDataDirectory( ) + "67p_ellipsoidal_harmonic_coefficients.txt",
                    getDataDirectory( ) + "67p_lame_coefficients.txt",
                    getDataDirectory( ) + "67p_transform_factors.txt",
                    15,
                    24,
                    3200,
                    0.0,
                    "Comet_fixed" ),
            "Comet",
            bodies ) );

    comet->setState( Eigen::Vector6d::Zero( ) );
    comet->setCurrentRotationToLocalFrameFromEphemeris( testTime );
    vehicle->setState( ( Eigen::Vector6d( ) << 5.0E3, 4.0E3, -3.0E3, 1.0E-2, -2.0E-3, 5.0E-3 ).finished( ) );

    return bodies;
}

//! Creates the acceleration of the comet on the vehicle and its analytical partial.
std::pair< std::shared_ptr< EllipsoidalHarmonicsGravitationalAccelerationModel >,
           std::shared_ptr< acceleration_partials::AccelerationPartial > >
createAccelerationAndPartial( const SystemOfBodies& bodies,
                              const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet )
{
    std::shared_ptr< EllipsoidalHarmonicsGravitationalAccelerationModel > acceleration =
            std::dynamic_pointer_cast< EllipsoidalHarmonicsGravitationalAccelerationModel >(
                    createAccelerationModel( bodies.at( "Vehicle" ),
                                             bodies.at( "Comet" ),
                                             std::make_shared< AccelerationSettings >( ellipsoidal_harmonic_gravity ),
                                             "Vehicle",
                                             "Comet" ) );
    acceleration->updateMembers( testTime );

    std::shared_ptr< acceleration_partials::AccelerationPartial > partial =
            createAnalyticalAccelerationPartial( acceleration,
                                                 std::make_pair( "Vehicle", bodies.at( "Vehicle" ) ),
                                                 std::make_pair( "Comet", bodies.at( "Comet" ) ),
                                                 bodies,
                                                 parameterSet );
    partial->update( testTime );

    return std::make_pair( acceleration, partial );
}

//! The partials of the acceleration w.r.t. the positions and velocities of both bodies must match numerical ones.
BOOST_AUTO_TEST_CASE( testStatePartials )
{
    // A fixed rotation that is not the identity, so the partial has to be rotated correctly.
    const Eigen::Quaterniond rotationToInertialFrame( Eigen::AngleAxisd( 0.6, Eigen::Vector3d( 0.2, 0.5, -0.8 ).normalized( ) ) );
    SystemOfBodies bodies = createTestBodies(
            constantRotationModelSettings( "ECLIPJ2000", "Comet_fixed", Eigen::Matrix3d( rotationToInertialFrame ).transpose( ) ) );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< estimatable_parameters::InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Vehicle", bodies.at( "Vehicle" )->getState( ), "Comet" ) );
    const auto accelerationAndPartial = createAccelerationAndPartial( bodies, createParametersToEstimate( parameterNames, bodies ) );
    const auto acceleration = accelerationAndPartial.first;
    const auto partial = accelerationAndPartial.second;

    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    partial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    partial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtCometPosition = Eigen::Matrix3d::Zero( );
    partial->wrtPositionOfAcceleratingBody( partialWrtCometPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtCometVelocity = Eigen::Matrix3d::Zero( );
    partial->wrtVelocityOfAcceleratingBody( partialWrtCometVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    // The body is only a few kilometres across, so the steps are much smaller than for a planet.
    const Eigen::Vector3d positionPerturbation = Eigen::Vector3d::Constant( 1.0 );
    const Eigen::Vector3d velocityPerturbation = Eigen::Vector3d::Constant( 1.0E-3 );
    std::function< void( Eigen::Vector6d ) > cometStateSetFunction =
            std::bind( &Body::setState, bodies.at( "Comet" ), std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            std::bind( &Body::setState, bodies.at( "Vehicle" ), std::placeholders::_1 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            acceleration_partials::calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, acceleration, bodies.at( "Vehicle" )->getState( ), positionPerturbation, 0 ),
            partialWrtVehiclePosition,
            1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            acceleration_partials::calculateAccelerationWrtStatePartials(
                    cometStateSetFunction, acceleration, bodies.at( "Comet" )->getState( ), positionPerturbation, 0 ),
            partialWrtCometPosition,
            1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            acceleration_partials::calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, acceleration, bodies.at( "Vehicle" )->getState( ), velocityPerturbation, 3 ),
            partialWrtVehicleVelocity,
            1.0E-12 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            acceleration_partials::calculateAccelerationWrtStatePartials(
                    cometStateSetFunction, acceleration, bodies.at( "Comet" )->getState( ), velocityPerturbation, 3 ),
            partialWrtCometVelocity,
            1.0E-12 );
}

//! The partials w.r.t. the comet's gravitational parameter, rotation rate and pole position must match numerical ones.
BOOST_AUTO_TEST_CASE( testParameterPartials )
{
    using namespace tudat::estimatable_parameters;

    // A rotation model with an estimatable pole and rate. The rate is that of 67P.
    const Eigen::Quaterniond initialOrientation( Eigen::AngleAxisd( 0.6, Eigen::Vector3d( 0.2, 0.5, -0.8 ).normalized( ) ) );
    SystemOfBodies bodies = createTestBodies( simpleRotationModelSettings(
            "ECLIPJ2000", "Comet_fixed", initialOrientation, 0.0, 2.0 * mathematical_constants::PI / 44676.0 ) );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( gravitationalParameter( "Comet" ) );
    parameterNames.push_back( constantRotationRate( "Comet" ) );
    parameterNames.push_back( rotationPolePosition( "Comet" ) );
    std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createParametersToEstimate( parameterNames, bodies );
    const auto accelerationAndPartial = createAccelerationAndPartial( bodies, parameterSet );
    const auto acceleration = accelerationAndPartial.first;
    const auto partial = accelerationAndPartial.second;

    const std::shared_ptr< EstimatableParameter< double > > gravitationalParameterObject = parameterSet->getDoubleParameters( ).at( 0 );
    const std::shared_ptr< EstimatableParameter< double > > rotationRateObject = parameterSet->getDoubleParameters( ).at( 1 );
    const std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > polePositionObject = parameterSet->getVectorParameters( ).at( 2 );
    BOOST_REQUIRE( gravitationalParameterObject->getParameterName( ).first == gravitational_parameter );
    BOOST_REQUIRE( rotationRateObject->getParameterName( ).first == constant_rotation_rate );
    BOOST_REQUIRE( polePositionObject->getParameterName( ).first == rotation_pole_position );

    const std::function< void( const double ) > rotationUpdate =
            std::bind( &Body::setCurrentRotationToLocalFrameFromEphemeris, bodies.at( "Comet" ), std::placeholders::_1 );

    // The acceleration is linear in the gravitational parameter, so a large step is exact.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( Eigen::Vector3d( partial->wrtParameter( gravitationalParameterObject ) ),
                                       Eigen::Vector3d( acceleration_partials::calculateAccelerationWrtParameterPartials(
                                               gravitationalParameterObject, acceleration, 1.0 ) ),
                                       1.0E-10 );

    // Orientation steps of about 1E-5 radians.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            Eigen::Vector3d( partial->wrtParameter( rotationRateObject ) ),
            Eigen::Vector3d( acceleration_partials::calculateAccelerationWrtParameterPartials(
                    rotationRateObject, acceleration, 1.0E-8, &acceleration_partials::emptyFunction, testTime, rotationUpdate ) ),
            1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            Eigen::MatrixXd( partial->wrtParameter( polePositionObject ) ),
            Eigen::MatrixXd( acceleration_partials::calculateAccelerationWrtParameterPartials( polePositionObject,
                                                                                               acceleration,
                                                                                               Eigen::Vector2d::Constant( 1.0E-5 ),
                                                                                               &acceleration_partials::emptyFunction,
                                                                                               testTime,
                                                                                               rotationUpdate ) ),
            1.0E-5 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
