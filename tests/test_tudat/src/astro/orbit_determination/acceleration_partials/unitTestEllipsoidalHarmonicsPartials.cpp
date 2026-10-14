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

#include "tudat/astro/orbit_determination/acceleration_partials/ellipsoidalHarmonicsAccelerationPartial.h"
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

/*! The partial of the acceleration w.r.t. position must match a numerical one.
 *
 * The analytic partial is the Hessian of the potential, rotated into the propagation frame. It
 * is checked here against the generic numerical partial Tudat uses for exactly this purpose.
 */
BOOST_AUTO_TEST_CASE( testEllipsoidalHarmonicsAccelerationPartial )
{
    const Eigen::Vector3d semiAxes =
            ( Eigen::Vector3d( ) << 2.802437330193671E+03, 2.090307360485069E+03, 1.867442428762624E+03 ).finished( );

    std::shared_ptr< Body > comet = std::make_shared< Body >( );
    std::shared_ptr< Body > vehicle = std::make_shared< Body >( );
    vehicle->setConstantBodyMass( 400.0 );

    SystemOfBodies bodies;
    bodies.addBody( comet, "Comet" );
    bodies.addBody( vehicle, "Vehicle" );

    // A fixed, non-identity rotation, so the partial has to be rotated correctly rather than
    // just passed through.
    const Eigen::Quaterniond rotationToInertialFrame( Eigen::AngleAxisd( 0.6, Eigen::Vector3d( 0.2, 0.5, -0.8 ).normalized( ) ) );
    comet->setRotationalEphemeris( createRotationModel(
            constantRotationModelSettings( "ECLIPJ2000", "Comet_fixed", Eigen::Matrix3d( rotationToInertialFrame ).transpose( ) ),
            "Comet" ) );

    std::shared_ptr< GravityFieldSettings > gravityFieldSettings =
            ellipsoidalHarmonicsGravitySettingsFromFiles( semiAxes,
                                                          getDataDirectory( ) + "67p_ellipsoidal_harmonic_coefficients.txt",
                                                          getDataDirectory( ) + "67p_lame_coefficients.txt",
                                                          getDataDirectory( ) + "67p_transform_factors.txt",
                                                          15,
                                                          24,
                                                          3200,
                                                          0.0,
                                                          "Comet_fixed" );
    comet->setGravityFieldModel( createGravityFieldModel( gravityFieldSettings, "Comet", bodies ) );

    const double testTime = 1.0E3;
    comet->setState( Eigen::Vector6d::Zero( ) );
    comet->setCurrentRotationToLocalFrameFromEphemeris( testTime );

    // Well outside the reference ellipsoid and away from the long axis.
    Eigen::Vector6d vehicleState = Eigen::Vector6d::Zero( );
    vehicleState << 5.0E3, 4.0E3, -3.0E3, 1.0E-2, -2.0E-3, 5.0E-3;
    vehicle->setState( vehicleState );

    std::shared_ptr< AccelerationSettings > accelerationSettings =
            std::make_shared< AccelerationSettings >( basic_astrodynamics::ellipsoidal_harmonic_gravity );
    std::shared_ptr< EllipsoidalHarmonicsGravitationalAccelerationModel > gravitationalAcceleration =
            std::dynamic_pointer_cast< EllipsoidalHarmonicsGravitationalAccelerationModel >(
                    createAccelerationModel( vehicle, comet, accelerationSettings, "Vehicle", "Comet" ) );
    BOOST_REQUIRE( gravitationalAcceleration != nullptr );
    gravitationalAcceleration->updateMembers( testTime );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< estimatable_parameters::InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Vehicle", vehicleState, "Comet" ) );
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodies );

    std::shared_ptr< acceleration_partials::EllipsoidalHarmonicsGravityPartial > accelerationPartial =
            std::dynamic_pointer_cast< acceleration_partials::EllipsoidalHarmonicsGravityPartial >(
                    createAnalyticalAccelerationPartial( gravitationalAcceleration,
                                                         std::make_pair( "Vehicle", vehicle ),
                                                         std::make_pair( "Comet", comet ),
                                                         bodies,
                                                         parameterSet ) );
    BOOST_REQUIRE( accelerationPartial != nullptr );

    accelerationPartial->update( testTime );

    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtCometPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtCometPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtCometVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtCometVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    std::function< void( Eigen::Vector6d ) > cometStateSetFunction = std::bind( &Body::setState, comet, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction = std::bind( &Body::setState, vehicle, std::placeholders::_1 );

    // The body is only a few kilometres across, so the steps have to be far smaller than the
    // ones a planet-sized test would use.
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1.0, 1.0, 1.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    const Eigen::Matrix3d testPartialWrtVehiclePosition = acceleration_partials::calculateAccelerationWrtStatePartials(
            vehicleStateSetFunction, gravitationalAcceleration, vehicle->getState( ), positionPerturbation, 0 );
    const Eigen::Matrix3d testPartialWrtVehicleVelocity = acceleration_partials::calculateAccelerationWrtStatePartials(
            vehicleStateSetFunction, gravitationalAcceleration, vehicle->getState( ), velocityPerturbation, 3 );
    const Eigen::Matrix3d testPartialWrtCometPosition = acceleration_partials::calculateAccelerationWrtStatePartials(
            cometStateSetFunction, gravitationalAcceleration, comet->getState( ), positionPerturbation, 0 );
    const Eigen::Matrix3d testPartialWrtCometVelocity = acceleration_partials::calculateAccelerationWrtStatePartials(
            cometStateSetFunction, gravitationalAcceleration, comet->getState( ), velocityPerturbation, 3 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition, partialWrtVehiclePosition, 1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtCometPosition, partialWrtCometPosition, 1.0E-5 );

    // A gravitational acceleration does not depend on either velocity.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity, partialWrtVehicleVelocity, 1.0E-12 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtCometVelocity, partialWrtCometVelocity, 1.0E-12 );

    // The partial w.r.t. the position of the accelerating body is minus the other one.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( Eigen::Matrix3d( partialWrtCometPosition ), Eigen::Matrix3d( -partialWrtVehiclePosition ), 1.0E-14 );

    // The partial is the Hessian of the potential, so in the body-fixed frame it is symmetric
    // and traceless outside the body.
    const Eigen::Matrix3d bodyFixedPartial = accelerationPartial->getCurrentBodyFixedPartialWrtPosition( );
    BOOST_CHECK_SMALL(
            ( bodyFixedPartial - bodyFixedPartial.transpose( ) ).cwiseAbs( ).maxCoeff( ) / bodyFixedPartial.cwiseAbs( ).maxCoeff( ),
            1.0E-9 );
    BOOST_CHECK_SMALL( std::fabs( bodyFixedPartial.trace( ) ) / bodyFixedPartial.cwiseAbs( ).maxCoeff( ), 1.0E-6 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
