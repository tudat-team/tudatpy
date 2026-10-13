/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <string>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/gravitation/ellipsoidalHarmonicsGravityModel.h"
#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/ellipsoidalHarmonicsFileReader.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::input_output;

BOOST_AUTO_TEST_SUITE( test_ellipsoidal_harmonics_gravity_model )

//! Reference ellipsoid of 67P/Churyumov-Gerasimenko, longest semi-axis first, in metres.
Eigen::Vector3d getReferenceSemiAxes( )
{
    return ( Eigen::Vector3d( ) << 2.802437330193671E+03, 2.090307360485069E+03, 1.867442428762624E+03 ).finished( );
}

std::string getDataDirectory( )
{
    return paths::getTudatTestDataPath( ) + "/ellipsoidal_harmonics/";
}

std::shared_ptr< GravityFieldSettings > getTestFieldSettings( const std::string& frameName )
{
    return ellipsoidalHarmonicsGravitySettingsFromFiles( getReferenceSemiAxes( ),
                                                         getDataDirectory( ) + "67p_ellipsoidal_harmonic_coefficients.txt",
                                                         getDataDirectory( ) + "67p_lame_coefficients.txt",
                                                         getDataDirectory( ) + "67p_transform_factors.txt",
                                                         15,
                                                         24,
                                                         3200,
                                                         0.0,
                                                         frameName );
}

/*! Builds a two-body system: a comet with the ellipsoidal harmonic field, and a spacecraft.
 *
 * The comet is given a fixed rotation so the body-fixed and propagation frames are related by a
 * known, non-identity rotation, which is what the acceleration model has to undo and redo.
 */
SystemOfBodies createTestBodies( const Eigen::Quaterniond& rotationToInertialFrame )
{
    BodyListSettings bodySettings( "SSB", "ECLIPJ2000" );
    bodySettings.addSettings( "67P" );
    bodySettings.at( "67P" )->gravityFieldSettings = getTestFieldSettings( "67P_fixed" );
    bodySettings.at( "67P" )->rotationModelSettings =
            constantRotationModelSettings( "ECLIPJ2000", "67P_fixed", Eigen::Matrix3d( rotationToInertialFrame ).transpose( ) );
    bodySettings.at( "67P" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "SSB", "ECLIPJ2000" );

    bodySettings.addSettings( "Spacecraft" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.at( "Spacecraft" )->setConstantBodyMass( 1.0E3 );

    return bodies;
}

//! The acceleration must be the gradient of the potential, rotated into the propagation frame.
BOOST_AUTO_TEST_CASE( test_acceleration_is_rotated_gradient_of_potential )
{
    // A deliberately non-trivial rotation between the body-fixed and propagation frames.
    const Eigen::Quaterniond rotationToInertialFrame( Eigen::AngleAxisd( 0.7, Eigen::Vector3d( 0.3, -0.5, 0.8 ).normalized( ) ) );

    SystemOfBodies bodies = createTestBodies( rotationToInertialFrame );

    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Spacecraft" ][ "67P" ].push_back( ellipsoidalHarmonicsAcceleration( ) );

    const AccelerationMap accelerationModels = createAccelerationModelsMap( bodies, accelerationSettings, { "Spacecraft" }, { "67P" } );

    const std::shared_ptr< AccelerationModel3d > accelerationModel = accelerationModels.at( "Spacecraft" ).at( "67P" ).at( 0 );

    // The settings carry no acceleration-specific options, so the factory must still have made
    // the ellipsoidal harmonics model rather than falling back on something else.
    BOOST_CHECK( std::dynamic_pointer_cast< EllipsoidalHarmonicsGravitationalAccelerationModel >( accelerationModel ) != nullptr );
    BOOST_CHECK_EQUAL( getAccelerationModelType( accelerationModel ), ellipsoidal_harmonic_gravity );

    const std::shared_ptr< EllipsoidalHarmonicsGravityField > gravityField =
            std::dynamic_pointer_cast< EllipsoidalHarmonicsGravityField >( bodies.at( "67P" )->getGravityFieldModel( ) );
    BOOST_REQUIRE( gravityField != nullptr );

    const Eigen::Vector3d inertialPosition = ( Eigen::Vector3d( ) << 4.0E3, -3.0E3, 2.5E3 ).finished( );
    Eigen::Vector6d state = Eigen::Vector6d::Zero( );
    state.segment( 0, 3 ) = inertialPosition;
    bodies.at( "Spacecraft" )->setState( state );
    bodies.at( "67P" )->setState( Eigen::Vector6d::Zero( ) );
    bodies.at( "67P" )->setCurrentRotationToLocalFrameFromEphemeris( 0.0 );

    accelerationModel->updateMembers( 0.0 );
    const Eigen::Vector3d acceleration = accelerationModel->getAcceleration( );

    // Same thing computed by hand: rotate into the body-fixed frame, take the gradient there,
    // and rotate the result back.
    const Eigen::Vector3d bodyFixedPosition = rotationToInertialFrame.inverse( ) * inertialPosition;
    const Eigen::Vector3d expected = rotationToInertialFrame * gravityField->getGradientOfPotential( bodyFixedPosition );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( acceleration, expected, 1.0E-13 );

    // The acceleration must point generally towards the body, and be of the size the body's
    // gravitational parameter implies at this distance.
    BOOST_CHECK( acceleration.dot( inertialPosition ) < 0.0 );
    const double pointMassMagnitude = gravityField->getGravitationalParameter( ) / inertialPosition.squaredNorm( );
    BOOST_CHECK_CLOSE_FRACTION( acceleration.norm( ), pointMassMagnitude, 0.3 );
}

//! Without a rotation model the acceleration cannot be formed, so it must be refused.
BOOST_AUTO_TEST_CASE( test_missing_rotation_model_is_rejected )
{
    BodyListSettings bodySettings( "SSB", "ECLIPJ2000" );
    bodySettings.addSettings( "67P" );
    bodySettings.at( "67P" )->gravityFieldSettings = getTestFieldSettings( "67P_fixed" );
    bodySettings.at( "67P" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "SSB", "ECLIPJ2000" );
    bodySettings.addSettings( "Spacecraft" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Spacecraft" ][ "67P" ].push_back( ellipsoidalHarmonicsAcceleration( ) );

    BOOST_CHECK_THROW( createAccelerationModelsMap( bodies, accelerationSettings, { "Spacecraft" }, { "67P" } ), std::runtime_error );
}

//! A rotation model whose frame does not match the field's must be refused.
BOOST_AUTO_TEST_CASE( test_mismatched_body_fixed_frame_is_rejected )
{
    BodyListSettings bodySettings( "SSB", "ECLIPJ2000" );
    bodySettings.addSettings( "67P" );
    // The field says one frame name, the rotation model another.
    bodySettings.at( "67P" )->gravityFieldSettings = getTestFieldSettings( "67P_fixed" );
    bodySettings.at( "67P" )->rotationModelSettings =
            constantRotationModelSettings( "ECLIPJ2000", "some_other_frame", Eigen::Matrix3d::Identity( ) );
    bodySettings.at( "67P" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "SSB", "ECLIPJ2000" );
    bodySettings.addSettings( "Spacecraft" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Spacecraft" ][ "67P" ].push_back( ellipsoidalHarmonicsAcceleration( ) );

    BOOST_CHECK_THROW( createAccelerationModelsMap( bodies, accelerationSettings, { "Spacecraft" }, { "67P" } ), std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
