/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/astro/gravitation/ellipsoidalHarmonicsGravityModel.h"
#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
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

//! The acceleration must be the gradient of the potential, rotated from the body-fixed frame into the inertial frame.
BOOST_AUTO_TEST_CASE( testAcceleration )
{
    // A non-trivial rotation between the body-fixed and inertial frames.
    const Eigen::Quaterniond rotationToInertialFrame( Eigen::AngleAxisd( 0.7, Eigen::Vector3d( 0.3, -0.5, 0.8 ).normalized( ) ) );

    BodyListSettings bodySettings( "SSB", "ECLIPJ2000" );
    bodySettings.addSettings( "67P" );
    bodySettings.at( "67P" )->gravityFieldSettings =
            ellipsoidalHarmonicsGravitySettingsFromFiles( getReferenceSemiAxes( ),
                                                          getDataDirectory( ) + "67p_ellipsoidal_harmonic_coefficients.txt",
                                                          getDataDirectory( ) + "67p_lame_coefficients.txt",
                                                          getDataDirectory( ) + "67p_transform_factors.txt",
                                                          15,
                                                          24,
                                                          3200,
                                                          0.0,
                                                          "67P_fixed" );
    bodySettings.at( "67P" )->rotationModelSettings =
            constantRotationModelSettings( "ECLIPJ2000", "67P_fixed", Eigen::Matrix3d( rotationToInertialFrame ).transpose( ) );
    bodySettings.at( "67P" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "SSB", "ECLIPJ2000" );
    bodySettings.addSettings( "Spacecraft" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Spacecraft" ][ "67P" ].push_back( ellipsoidalHarmonicsAcceleration( ) );
    const std::shared_ptr< AccelerationModel3d > accelerationModel =
            createAccelerationModelsMap( bodies, accelerationSettings, { "Spacecraft" }, { "67P" } ).at( "Spacecraft" ).at( "67P" ).at( 0 );
    BOOST_CHECK_EQUAL( getAccelerationModelType( accelerationModel ), ellipsoidal_harmonic_gravity );

    const Eigen::Vector3d inertialPosition = ( Eigen::Vector3d( ) << 4.0E3, -3.0E3, 2.5E3 ).finished( );
    Eigen::Vector6d state = Eigen::Vector6d::Zero( );
    state.segment( 0, 3 ) = inertialPosition;
    bodies.at( "Spacecraft" )->setState( state );
    bodies.at( "67P" )->setState( Eigen::Vector6d::Zero( ) );
    bodies.at( "67P" )->setCurrentRotationToLocalFrameFromEphemeris( 0.0 );
    accelerationModel->updateMembers( 0.0 );

    // Same result by hand: rotate the position into the body-fixed frame, take the gradient there and rotate it back.
    const Eigen::Vector3d expected = rotationToInertialFrame *
            bodies.at( "67P" )->getGravityFieldModel( )->getGradientOfPotential( rotationToInertialFrame.inverse( ) * inertialPosition );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( accelerationModel->getAcceleration( ), expected, 1.0E-13 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
