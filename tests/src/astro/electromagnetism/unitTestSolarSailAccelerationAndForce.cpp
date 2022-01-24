/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Ganeff, M.I. Solar radiation pressure benchmark data script,
 *          solarRadiationPressureBenchmarkData.m, available at http://tudat.tudelft.nl, 2012.
 *      Giancoli, D.C. Physics for Scientists and Engineers with Modern Physics 
 *          Fourth Edition. New Jersey: Prentice-Hall, Inc., 1985.
 *      Hirsh, S.M. Solar radiation pressure benchmark data script,
 *          solarRadiationUnitTests.cpp, available at 
 *          https://github.com/sethhirsh/solarRadiationUnitTests, 2013.
 *      Irizarry, V. The effect of solar radiation pressure on Ulysses orbit,
 *          http://ccar.colorado.edu/asen5050/projects/projects_2001/aponte/Ulysses.htm, 2001,
 *          last accessed: 7th October, 2013.
 *      Willmott, P. General astro Library, http://www.amsat-bda.org/GAL_Downloads.html,
 *          2011, last accessed: 7th October, 2013.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/electromagnetism/solarSailForce.h"
#include "tudat/astro/electromagnetism/solarSailAcceleration.h"
#include "tudat/astro/electromagnetism/cannonBallRadiationPressureForce.h"
#include "tudat/astro/electromagnetism/cannonBallRadiationPressureAcceleration.h"
#include "tudat/astro/basic_astro/unitConversions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_solar_sail_acceleration_and_force_models )

// Set radiation pressure at 1 AU [N/m^2].
const double radiationPressureAtOneAU = 4.56e-6;

// Set 1 AU in metres [m].
const double astronomicalUnitInMeters = 1.49598e11;

// Set gravitational parameter of the Earth [m^3 s^-2].
// Reference: Mass taken from Table 1.2, pg. 6, (de Pater, 2010), value
//            of gravitational constant taken from http://ssd.jpl.nasa.gov/?constants#ref.
const double earthGravitationalParameter = 3.9859383624e14;

//! Test solar sail force in the case where the non-ideal solar sail model is equivalent to the cannonball radiation pressure model.
BOOST_AUTO_TEST_CASE( testSolarSailModelVsCannonBall )
{
    // Data for the cannon ball model.

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.21;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 0.5;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource
            = Eigen::Vector3d( astronomicalUnitInMeters, astronomicalUnitInMeters, 0.0 );

    // Compute distance between spacecraft and source [m].
    const double positionVectorToSourceNorm = positionVectorToSource.norm();

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector from spacecraft to source [-].
    positionVectorToSource.normalize( );

    // Compute radiation pressure force from cannonball radiation pressure model [N].
    const Eigen::Vector3d computedRadiationPressureForce = electromagnetism::computeCannonBallRadiationPressureForce(
                radiationPressureAtTarget, positionVectorToSource, areaSubjectToRadiationPressure, radiationPressureCoefficient );

    // Define solar sail model settings.
    const double frontEmissivityCoefficient = 0.4;
    const double backEmissivityCoefficient = frontEmissivityCoefficient;
    const double frontLambertianCoefficient = 0.4;
    const double backLambertianCoefficient = 0.4;
    const double reflectivityCoefficient = 0.21;
    const double specularReflectionCoefficient = 1.0;
    const double coneAngle = 0.0;
    const double clockAngle = 0.0;

    // Compute velocity vector of the spacecraft w.r.t. central body [m/s].
    Eigen::Vector3d velocityVector = std::sqrt( earthGravitationalParameter / positionVectorToSourceNorm )
        * Eigen::Vector3d( std::sqrt( 2.0 ) / 2.0, std::sqrt( 2.0 ) / 2.0, 0.0 );

    // Compute force from solar sail model [N].
    const Eigen::Vector3d computedSolarSailForce = electromagnetism::computeSolarSailForce(
                frontEmissivityCoefficient, backEmissivityCoefficient, frontLambertianCoefficient, backLambertianCoefficient,
                reflectivityCoefficient, specularReflectionCoefficient, positionVectorToSource, velocityVector.normalized(),
                radiationPressureAtTarget, areaSubjectToRadiationPressure, coneAngle, clockAngle);

    // Compare computed and expected radiation pressure force vectors.
    BOOST_CHECK_EQUAL( computedRadiationPressureForce.z( ), computedSolarSailForce.z( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureForce.segment( 0, 2 ), computedSolarSailForce.segment( 0, 2 ),
                                       1.0e-15 );
}


//! Test solar sail acceleration model on spacecraft at approximately 1 AU away from the Sun, with a solar sail model
//! equivalent to the cannonball radiation pressure model.
BOOST_AUTO_TEST_CASE( testSolarSailAccelerationModelEarth )
{
    // Benchmark data is obtained using the General astro Library (Willmott, 2011).

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.3;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 2.0;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource = Eigen::Vector3d( astronomicalUnitInMeters, 0.0, 0.0 );

    // Compute distance between spacecraft and source [m].
    const double positionVectorToSourceNorm = positionVectorToSource.norm();

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units [-].
    positionVectorToSource.normalize( );

    // Set mass of accelerated body [kg].
    const double mass = 4.0;

    // Compute radiation pressure acceleration from cannonball radiation pressure model [m/s^2].
    const Eigen::Vector3d computedRadiationPressureAcceleration = electromagnetism::computeCannonBallRadiationPressureAcceleration(
                radiationPressureAtTarget, positionVectorToSource, areaSubjectToRadiationPressure, radiationPressureCoefficient, mass );

    // Define solar sail model settings.
    const double frontEmissivityCoefficient = 0.4;
    const double backEmissivityCoefficient = frontEmissivityCoefficient;
    const double frontLambertianCoefficient = 0.4;
    const double backLambertianCoefficient = 0.4;
    const double reflectivityCoefficient = 0.3;
    const double specularReflectionCoefficient = 1.0;
    const double coneAngle = 0.0;
    const double clockAngle = 0.0;

    // Compute velocity vector of the spacecraft w.r.t. central body [m/s].
    Eigen::Vector3d velocityVector = std::sqrt( earthGravitationalParameter / positionVectorToSourceNorm )
            * Eigen::Vector3d( std::sqrt( 2.0 ) / 2.0, std::sqrt( 2.0 ) / 2.0, 0.0 );

    // Compute solar sail acceleration [m/s^2].
    const Eigen::Vector3d computedSolarSailAcceleration = electromagnetism::computeSolarSailAcceleration(
                frontEmissivityCoefficient, backEmissivityCoefficient, frontLambertianCoefficient, backLambertianCoefficient,
                reflectivityCoefficient, specularReflectionCoefficient, positionVectorToSource, velocityVector.normalized(),
                radiationPressureAtTarget, areaSubjectToRadiationPressure, coneAngle, clockAngle, mass);

    // Compare computed and expected radiation pressure acceleration vectors
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureAcceleration, computedSolarSailAcceleration,
                                       1.0e-15 );
}


//! Test solar sail acceleration model on spacecraft at approximately the distance of Venus away from the Sun, with a solar sail model
//! equivalent to the cannonball radiation pressure model.
BOOST_AUTO_TEST_CASE( testSolarSailAccelerationModelVenus )
{
    // Benchmark data is obtained using the General astro Library (Willmott, 2011).

    // Set the distance from the Sun to Venus [m].
    const double distanceSunToVenus = 0.732 * astronomicalUnitInMeters;

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.5;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 0.005;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource
            = Eigen::Vector3d( distanceSunToVenus / std::sqrt( 2.0 ), distanceSunToVenus / std::sqrt( 2.0 ), 0.0 );

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units [-].
    positionVectorToSource.normalize( );

    // Set mass of accelerated body [kg].
    const double mass = 0.0022;

    // Compute radiation pressure acceleration from cannonball radiation pressure model [m/s^2].
    const Eigen::Vector3d computedRadiationPressureAcceleration = electromagnetism::computeCannonBallRadiationPressureAcceleration(
                radiationPressureAtTarget, positionVectorToSource, areaSubjectToRadiationPressure, radiationPressureCoefficient, mass );

    // Define solar sail model settings.
    const double frontEmissivityCoefficient = 0.4;
    const double backEmissivityCoefficient = frontEmissivityCoefficient;
    const double frontLambertianCoefficient = 0.4;
    const double backLambertianCoefficient = 0.4;
    const double reflectivityCoefficient = 0.5;
    const double specularReflectionCoefficient = 1.0;
    const double coneAngle = 0.0;
    const double clockAngle = 0.0;

    // Compute velocity vector of the spacecraft w.r.t. central body [m/s].
    Eigen::Vector3d velocityVector = std::sqrt( earthGravitationalParameter / distanceSunToVenus)
        * Eigen::Vector3d( std::sqrt( 2.0 ) / 2.0, std::sqrt( 2.0 ) / 2.0, 0.0 );

    // Compute solar sail acceleration [m/s^2].
    const Eigen::Vector3d computedSolarSailAcceleration
        = electromagnetism::computeSolarSailAcceleration(
            frontEmissivityCoefficient, backEmissivityCoefficient, frontLambertianCoefficient,
            backLambertianCoefficient, reflectivityCoefficient, specularReflectionCoefficient, positionVectorToSource,
            velocityVector.normalized(), radiationPressureAtTarget, areaSubjectToRadiationPressure, coneAngle, clockAngle, mass);


    // Check there is no acceleration generated along the z-axis.
    BOOST_CHECK_SMALL( computedRadiationPressureAcceleration.z( ), std::numeric_limits< double >::min( ) );

    // Compare computed and expected radiation pressure acceleration vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureAcceleration.segment( 0, 2 ),
                                       computedSolarSailAcceleration.segment( 0, 2 ),
                                       1.0e-14 );
}


//! Test solar sail force model on spacecraft at approximately the distance of Uranus away from the Sun, with a solar sail model
//! equivalent to the cannonball radiation pressure model.
BOOST_AUTO_TEST_CASE( testSolarSailForceModelUranus )
{
    // Benchmark data is obtained using the General astro Library (Willmott, 2011).

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.8;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 69939064094327.4;

    // Define distance between Uranus and the Sun, which is also the distance between the spacecraft and the Sun [m].
    const double distanceSunToUranus = 9.529 * astronomicalUnitInMeters;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource
            = Eigen::Vector3d( distanceSunToUranus / std::sqrt( 2.0 ), distanceSunToUranus / std::sqrt( 2.0 ), 0.0 );

    // Compute radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units [-].
    positionVectorToSource.normalize( );

    // Compute radiation pressure force from cannonball radiation pressure model [N].
    const Eigen::Vector3d computedRadiationPressureForce = electromagnetism::computeCannonBallRadiationPressureForce(
                radiationPressureAtTarget, positionVectorToSource, areaSubjectToRadiationPressure, radiationPressureCoefficient );

    // Define solar sail model settings.
    const double frontEmissivityCoefficient = 0.4;
    const double backEmissivityCoefficient = frontEmissivityCoefficient;
    const double frontLambertianCoefficient = 0.4;
    const double backLambertianCoefficient = 0.4;
    const double reflectivityCoefficient = 0.8;
    const double specularReflectionCoefficient = 1.0;
    const double coneAngle = 0.0;
    const double clockAngle = 0.0;

    // Compute velocity vector of the spacecraft w.r.t. the central body [m/s].
    Eigen::Vector3d velocityVector = std::sqrt( earthGravitationalParameter / distanceSunToUranus )
        * Eigen::Vector3d( std::sqrt( 2.0 ) / 2.0, std::sqrt( 2.0 ) / 2.0, 0.0 );

    // Compute force from the solar sail model [N].
    const Eigen::Vector3d computedSolarSailForce = electromagnetism::computeSolarSailForce(
                frontEmissivityCoefficient, backEmissivityCoefficient, frontLambertianCoefficient,
                backLambertianCoefficient, reflectivityCoefficient, specularReflectionCoefficient, positionVectorToSource,
                velocityVector.normalized(), radiationPressureAtTarget, areaSubjectToRadiationPressure, coneAngle, clockAngle);


    // Compare computed and expected radiation pressure force vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureForce,
                                       computedSolarSailForce,
                                       1.0e-14 );
}


//! Test solar sail force model on a spacecraft at a random distance from the Sun, with a solar sail model equivalent to the
//! cannonball radiation pressure model.
BOOST_AUTO_TEST_CASE( testSolarSailForceModelRandom )
{
    // Benchmark data is obtained using the General astro Library (Willmott, 2011).

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.4058;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 514701.9505;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource = Eigen::Vector3d( 94359740.25, 90831886.1, 14668782.92 );

    // Compute distance of the spacecraft w.r.t. the Sun [m].
    const double positionVectorToSourceNorm = positionVectorToSource.norm();

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units [-].
    positionVectorToSource.normalize( );

    // Compute radiation pressure force from cannonball radiation pressure model [N].
    const Eigen::Vector3d computedRadiationPressureForce = electromagnetism::computeCannonBallRadiationPressureForce(
                radiationPressureAtTarget, positionVectorToSource,
                areaSubjectToRadiationPressure, radiationPressureCoefficient );

    // Define solar sail model settings.
    const double frontEmissivityCoefficient = 0.4;
    const double backEmissivityCoefficient = frontEmissivityCoefficient;
    const double frontLambertianCoefficient = 0.4;
    const double backLambertianCoefficient = 0.4;
    const double reflectivityCoefficient = 0.4058;
    const double specularReflectionCoefficient = 1.0;
    const double coneAngle = 0.0;
    const double clockAngle = 0.0;

    // Compute velocity vector of the spacecraft w.r.t. the central body [m/s].
    Eigen::Vector3d velocityVector = std::sqrt( earthGravitationalParameter / positionVectorToSourceNorm )
        * Eigen::Vector3d( std::sqrt( 2.0 ) / 2.0, std::sqrt( 2.0 ) / 2.0, 0.0 );

    // Compute force from solar sail model [N].
    const Eigen::Vector3d computedSolarSailForce = electromagnetism::computeSolarSailForce(
            frontEmissivityCoefficient, backEmissivityCoefficient, frontLambertianCoefficient,
            backLambertianCoefficient, reflectivityCoefficient, specularReflectionCoefficient, positionVectorToSource,
            velocityVector.normalized(), radiationPressureAtTarget, areaSubjectToRadiationPressure, coneAngle, clockAngle);

    // Compare computed and expected radiation pressure force vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureForce, computedSolarSailForce, 1.0e-14 );
}

//! Test radiation force model on a hand at approximately 1 AU from the Sun, with a solar sail model equivalent to the
//! cannonball radiation pressure model.
BOOST_AUTO_TEST_CASE( testSolarSailForceModelGiancoliData )
{
    // Benchmark data is obtained from Physics for Scientists and Engineers
    // with Modern Physics Volume 2 (Ch. 31, Ex. 7) (Giancoli, 1985).

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.0;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 0.02;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource = Eigen::Vector3d( astronomicalUnitInMeters, 0.0, 0.0 );

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units [-].
    positionVectorToSource.normalize( );

    // Compute radiation pressure force from cannonball radiation pressure model [N].
    const Eigen::Vector3d computedRadiationPressureForce = electromagnetism::computeCannonBallRadiationPressureForce(
                radiationPressureAtTarget, positionVectorToSource, areaSubjectToRadiationPressure, radiationPressureCoefficient );

    // Define solar sail settings.
    const double frontEmissivityCoefficient = 0.4;
    const double backEmissivityCoefficient = frontEmissivityCoefficient;
    const double frontLambertianCoefficient = 0.4;
    const double backLambertianCoefficient = 0.4;
    const double reflectivityCoefficient = 0.0;
    const double specularReflectionCoefficient = 1.0;
    const double coneAngle = 0.0;
    const double clockAngle = 0.0;

    // Compute velocity vector of the spacecraft w.r.t. the central body [m/s].
    Eigen::Vector3d velocityVector = std::sqrt(earthGravitationalParameter/astronomicalUnitInMeters)
        * Eigen::Vector3d( std::sqrt( 2.0 ) / 2.0, std::sqrt( 2.0 ) / 2.0, 0.0 );

    // Compute force from solar sail model [N].
    const Eigen::Vector3d computedSolarSailForce = electromagnetism::computeSolarSailForce(
                frontEmissivityCoefficient, backEmissivityCoefficient, frontLambertianCoefficient,
                backLambertianCoefficient, reflectivityCoefficient, specularReflectionCoefficient, positionVectorToSource,
                velocityVector.normalized(), radiationPressureAtTarget, areaSubjectToRadiationPressure, coneAngle, clockAngle);


    // Check there is no force component along the x-axis.
    BOOST_CHECK_SMALL( std::fabs( computedRadiationPressureForce.x( ) - computedSolarSailForce.x( ) ), 1.0e-6 );

    // Compare computed and expected radiation pressure force vectors.
    BOOST_CHECK_SMALL( computedSolarSailForce.y( ), std::numeric_limits< double >::min( ) );
    BOOST_CHECK_SMALL( computedSolarSailForce.z( ), std::numeric_limits< double >::min( ) );
}


//! Test radiation acceleration model on Ulysses satellite at 1AU, with a solar sail model equivalent to the
//! cannonball radiation pressure model.
BOOST_AUTO_TEST_CASE( testSolarSailAccelerationModelUlysses )
{
    // Benchmark data obtained from (Irizarry, 2001).

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.327;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 10.59;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource = Eigen::Vector3d( astronomicalUnitInMeters, 0.0, 0.0 );

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units [-].
    positionVectorToSource.normalize( );

    // Set mass of accelerated body [kg].
    const double mass = 370.0;

    // Compute radiation pressure acceleration from cannonball radiation pressure model [m/s^2].
    const Eigen::Vector3d computedRadiationPressureAcceleration = electromagnetism::computeCannonBallRadiationPressureAcceleration(
                radiationPressureAtTarget, positionVectorToSource,
                areaSubjectToRadiationPressure, radiationPressureCoefficient, mass );

    // Define solar sail model settings.
    const double frontEmissivityCoefficient = 0.4;
    const double backEmissivityCoefficient = frontEmissivityCoefficient;
    const double frontLambertianCoefficient = 0.4;
    const double backLambertianCoefficient = 0.4;
    const double reflectivityCoefficient = 0.327;
    const double specularReflectionCoefficient = 1.0;
    const double coneAngle = 0.0;
    const double clockAngle = 0.0;

    // Compute velocity vector of the spacecraft w.r.t. the central body [m/s].
    Eigen::Vector3d velocityVector = std::sqrt( earthGravitationalParameter / astronomicalUnitInMeters)
        * Eigen::Vector3d( std::sqrt( 2.0 ) / 2.0, std::sqrt( 2.0 ) / 2.0, 0.0 );

    // Compute solar sail acceleration from solar sail model [m/s^2].
    const Eigen::Vector3d computedSolarSailAcceleration = electromagnetism::computeSolarSailAcceleration(
            frontEmissivityCoefficient, backEmissivityCoefficient, frontLambertianCoefficient,
            backLambertianCoefficient, reflectivityCoefficient, specularReflectionCoefficient, positionVectorToSource,
            velocityVector.normalized(), radiationPressureAtTarget, areaSubjectToRadiationPressure, coneAngle, clockAngle, mass);

    // Check that there is no acceleration component along the x-axis.
    BOOST_CHECK_SMALL( std::fabs( computedRadiationPressureAcceleration.x( )
                                  - computedSolarSailAcceleration.x( ) ), 1.0e-8 );

    // Compare computed and expected radiation pressure acceleration vectors.
    BOOST_CHECK_SMALL( computedSolarSailAcceleration.y( ), std::numeric_limits< double >::min( ) );
    BOOST_CHECK_SMALL( computedSolarSailAcceleration.z( ), std::numeric_limits< double >::min( ) );
}


//! Test class implementation of radiation pressure acceleration model.

// Set position of source of radiation pressure at origin [m].
static Eigen::Vector3d sourcePosition = Eigen::Vector3d::Zero( );

// Get position of source of radiation pressure [m].
Eigen::Vector3d getSourcePosition( ) { return sourcePosition; }

// Set position of accelerated body [m].
static Eigen::Vector3d acceleratedBodyPosition = Eigen::Vector3d( - astronomicalUnitInMeters, 0.0, 0.0 );

// Get position of accelerated body [m].
Eigen::Vector3d getAcceleratedBodyPosition( ) { return acceleratedBodyPosition; }

// Set velocity of accelerated body [m/s].
static Eigen::Vector3d acceleratedBodyVelocity
    = std::sqrt( earthGravitationalParameter / astronomicalUnitInMeters ) * Eigen::Vector3d( 1.0 / std::sqrt( 2.0 ), 1.0 / std::sqrt( 2.0 ), 0.0 );

// Get velocity of accelerated body [m/s].
Eigen::Vector3d getAcceleratedBodyVelocity( ) { return acceleratedBodyVelocity; }

// Set velocity of central body [m/s].
static Eigen::Vector3d centralBodyVelocity
    = Eigen::Vector3d( 0.0, 0.0, 0.0 );

// Get velocity of central body [m/s].
Eigen::Vector3d getCentralBodyVelocity( ) { return centralBodyVelocity; }

// Get vector from accelerated body to source [m].
Eigen::Vector3d getVectorToSource( )
{
    return getSourcePosition( ) - getAcceleratedBodyPosition( );
}

// Set radiation pressure at location of acceleration body [N/m^2].
static double radiationPressure = radiationPressureAtOneAU
    * astronomicalUnitInMeters * astronomicalUnitInMeters / getVectorToSource( ).squaredNorm( );

// Get radiation pressure at location of acceleration body [N/m^2].
double getRadiationPressure( ) { return radiationPressure; }

// Set cone angle of the solar sail [rad].
static double coneAngle = 0.0;

// Get cone angle of the solar sail [rad].
double getConeAngle( ) { return coneAngle; }

// Set clock angle of the solar sail [rad].
static double clockAngle = 0.0;

// Get clock angle of the solar sail [rad].
double getClockAngle( ) { return clockAngle; }

// Set front emissivity coefficient of the solar sail [-].
static double frontEmissivityCoefficient = 0.4;

// Get front emissivity coefficient of the solar sail [-].
double getFrontEmissivityCoefficient( ) { return frontEmissivityCoefficient; }

// Set back emissivity coefficient of the solar sail [-].
static double backEmissivityCoefficient = 0.4;

// Get back emissivity coefficient of the solar sail [-].
double getBackEmissivityCoefficient( ) { return backEmissivityCoefficient; }

// Set front Lambertian coefficient of the solar sail [-].
static double frontLambertianCoefficient = 0.4;

// Get front Lambertian coefficient of the solar sail [-].
double getFrontLambertianCoefficient( ) { return frontLambertianCoefficient; }

// Set back Lambertian coefficient of the solar sail [-].
static double backLambertianCoefficient = 0.4;

// Get back Lambertian coefficient of the solar sail [-].
double getBackLambertianCoefficient( ) { return backLambertianCoefficient; }

// Set reflectivity coefficient of the solar sail [-].
static double reflectivityCoefficient = 0.3;

// Get reflectivity coefficient of the solar sail [-].
double getReflectivityCoefficient( ) { return reflectivityCoefficient; }

// Set specular reflection coefficient of the solar sail [-].
static double specularReflectionCoefficient = 1.0;

// Get specular reflection coefficient of the solar sail [-].
double getSpecularReflectionCoefficient( ) { return specularReflectionCoefficient; }

// Set area subject to radiation pressure [m^2].
static double areaSubjectToRadiationPressure = 2.0;

// Get area subject to radiation pressure [m^2].
double getAreaSubjectToRadiationPressure( ) { return areaSubjectToRadiationPressure; }

// Set mass of accelerated body [kg].
static double massOfAcceleratedBody = 4.0;

// Get mass of acceleration body [kg].
double getMassOfAcceleratedBody( ) { return massOfAcceleratedBody; }

// Set radiation pressure coefficient [-].
static double radiationPressureCoefficient = 1.0 + 0.3;

// Get radiation pressure coefficient [-].
double getRadiationPressureCoefficient( ) { return radiationPressureCoefficient; }


//! Test solar sail acceleration model constructor, with a solar sail model equivalent to the
//! cannonball radiation pressure model.
BOOST_AUTO_TEST_CASE( testSolarSailAccelerationModelClassConstructor )
{

    // Declare and initialize cannon-ball radiation pressure acceleration model.
    electromagnetism::CannonBallRadiationPressurePointer radiationPressureModel
        = std::make_shared< electromagnetism::CannonBallRadiationPressureAcceleration >(
            &getSourcePosition, &getAcceleratedBodyPosition, &getRadiationPressure,
            &getRadiationPressureCoefficient, &getAreaSubjectToRadiationPressure,
            &getMassOfAcceleratedBody );
    radiationPressureModel->updateMembers( 0.0 );

    // Declare and initialize solar sail radiation acceleration model.
    electromagnetism::SolarSailAccelerationPointer solarSailModel
        = std::make_shared< electromagnetism::SolarSailAcceleration >(
            &getSourcePosition, &getAcceleratedBodyPosition, &getAcceleratedBodyVelocity,
            &getCentralBodyVelocity, &getRadiationPressure, &getConeAngle, &getClockAngle,
            &getFrontEmissivityCoefficient, &getBackEmissivityCoefficient,
            &getFrontLambertianCoefficient, &getBackLambertianCoefficient,
            &getReflectivityCoefficient, &getSpecularReflectionCoefficient,
            &getAreaSubjectToRadiationPressure, &getMassOfAcceleratedBody );
    solarSailModel->updateMembers( 0.0 );

    // Compute solar sail acceleration [m/s^2].
    const Eigen::Vector3d computedSolarSailAcceleration = solarSailModel->getAcceleration( );

    // Compute radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d computedRadiationPressureAcceleration = radiationPressureModel->getAcceleration( );

    // Compare computed and expected radiation pressure acceleration vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedSolarSailAcceleration,
                                       computedRadiationPressureAcceleration,
                                       1.0e-15 );
}

//! Test radiation pressure acceleration model update-members function, with a solar sail model equivalent to the
//! cannonball radiation pressure model.
BOOST_AUTO_TEST_CASE( testSolarSailAccelerationModelClassUpdateMembers )
{
    // Declare and initialize cannon-ball radiation pressure acceleration model.
    electromagnetism::CannonBallRadiationPressurePointer radiationPressureModel
        = std::make_shared< electromagnetism::CannonBallRadiationPressureAcceleration >(
            &getSourcePosition, &getAcceleratedBodyPosition, &getRadiationPressure,
            &getRadiationPressureCoefficient, &getAreaSubjectToRadiationPressure,
            &getMassOfAcceleratedBody );
    radiationPressureModel->updateMembers( 0.0 );

    // Declare and initialize solar sail radiation acceleration model.
    electromagnetism::SolarSailAccelerationPointer solarSailModel
        = std::make_shared< electromagnetism::SolarSailAcceleration >(
            &getSourcePosition, &getAcceleratedBodyPosition, &getAcceleratedBodyVelocity,
            &getCentralBodyVelocity, &getRadiationPressure, &getConeAngle, &getClockAngle,
            &getFrontEmissivityCoefficient, &getBackEmissivityCoefficient,
            &getFrontLambertianCoefficient, &getBackLambertianCoefficient,
            &getReflectivityCoefficient, &getSpecularReflectionCoefficient,
            &getAreaSubjectToRadiationPressure, &getMassOfAcceleratedBody );
    solarSailModel->updateMembers( 0.0 );

    // Set the distance from the Sun to Venus [m].
    const double distanceSunToVenus = 0.732 * astronomicalUnitInMeters;

    // Update position of accelerated body [m].
    acceleratedBodyPosition = Eigen::Vector3d( - distanceSunToVenus / std::sqrt( 2.0 ), - distanceSunToVenus / std::sqrt( 2.0 ), 0.0 );
    acceleratedBodyVelocity = std::sqrt( earthGravitationalParameter / distanceSunToVenus )
            * Eigen::Vector3d( 1.0 / std::sqrt( 2.0 ), 1.0 / std::sqrt( 2.0 ), 0.0 );

    // Update radiation pressure at location of accelerated body [N/m^2].
    radiationPressure = radiationPressureAtOneAU * astronomicalUnitInMeters * astronomicalUnitInMeters
        / getVectorToSource( ).squaredNorm( );

    // Update radiation pressure coefficient [-].
    radiationPressureCoefficient = 1.0 + 0.5;

    //Update reflectivity coefficient [-].
    reflectivityCoefficient = 0.5;

    // Update area subject to radiation pressure [m^2].
    areaSubjectToRadiationPressure = 0.005;

    // Update mass of accelerated body [kg].
    massOfAcceleratedBody = 0.0022;

    // Update class members.
    radiationPressureModel->updateMembers( 0.0 );
    solarSailModel->updateMembers( 0.0 );

    // Compute radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d computedRadiationPressureAcceleration = radiationPressureModel->getAcceleration( );
    const Eigen::Vector3d computedSolarSailAcceleration = solarSailModel->getAcceleration( );

    // Check that there is no acceleration component along z-axis.
    BOOST_CHECK_SMALL( computedSolarSailAcceleration.z( ), std::numeric_limits< double >::min( ) );

    // Compare computed and expected radiation pressure acceleration vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureAcceleration.segment( 0, 2 ),
                                       computedSolarSailAcceleration.segment( 0, 2 ),
                                       1.0e-14 );
}


//! Test implementation of solar sail force model for a set of different cone and clock angles.
BOOST_AUTO_TEST_CASE( testImplementedSolarSailModelVsTheoreticalFirst )
{
    // Theoretical non-ideal solar sail force vectors for various test cases [N].
    std::vector< Eigen::Vector3d > expectedSolarSailForceVector;
    for ( int testCase = 0 ; testCase < 8 ; testCase++ )
    {
        expectedSolarSailForceVector.push_back(
                    ( Eigen::Vector3d( ) << - 8.52932077149712e-24, - 8.52932077149712e-24, 0.0 ).finished() );
    }
    for ( int testCase = 8 ; testCase < 12 ; testCase++ )
    {
        expectedSolarSailForceVector.push_back(
                    ( Eigen::Vector3d( ) << - 1.46272332079017e-06, - 1.46272332079017e-06, 0.0 ).finished() );
    }
    for ( int testCase = 12 ; testCase < 16 ; testCase++ )
    {
        expectedSolarSailForceVector.push_back(
                    ( Eigen::Vector3d( ) << - 5.64907119362171e-07, - 5.64907119362171e-07, 0.0 ).finished() );
    }
    for ( int testCase = 16 ; testCase < 18 ; testCase++ )
    {
        expectedSolarSailForceVector.push_back(
                    ( Eigen::Vector3d( ) << - 3.0864061085007e-08, - 3.0864061085007e-08, 0.0 ).finished() );
    }
    expectedSolarSailForceVector.push_back(
                ( Eigen::Vector3d( ) << - 1.32692625650149e-06, - 1.32692625650149e-06, 0.0 ).finished() );
    expectedSolarSailForceVector.push_back(
                ( Eigen::Vector3d( ) << - 1.80571791929449e-07, - 1.80571791929449e-07, - 1.80571791929449e-07 ).finished() );
    expectedSolarSailForceVector.push_back(
                ( Eigen::Vector3d( ) << - 3.40314752162177e-06, - 4.53753002882903e-06, - 6.80629504324355e-06 ).finished() );
    expectedSolarSailForceVector.push_back(
                ( Eigen::Vector3d( ) << - 4.73337860380423e-06, - 6.29769698345645e-06, - 7.16893475205534e-06 ).finished() );


    // Define set of cone angles for various test cases [rad].
    std::vector< double > coneAngleVector;
    coneAngleVector.push_back( mathematical_constants::PI / 2.0 );
    coneAngleVector.push_back( mathematical_constants::PI / 2.0 );
    coneAngleVector.push_back( mathematical_constants::PI / 2.0 );
    coneAngleVector.push_back( mathematical_constants::PI / 2.0 );
    coneAngleVector.push_back( - mathematical_constants::PI / 2.0 );
    coneAngleVector.push_back( - mathematical_constants::PI / 2.0 );
    coneAngleVector.push_back( - mathematical_constants::PI / 2.0 );
    coneAngleVector.push_back( - mathematical_constants::PI / 2.0 );
    coneAngleVector.push_back( 0.0 );
    coneAngleVector.push_back( 0.0 );
    coneAngleVector.push_back( 0.0 );
    coneAngleVector.push_back( 0.0 );
    coneAngleVector.push_back( mathematical_constants::PI / 4.0 );
    coneAngleVector.push_back( - mathematical_constants::PI / 4.0 );
    coneAngleVector.push_back( mathematical_constants::PI / 4.0 );
    coneAngleVector.push_back( - mathematical_constants::PI / 4.0 );
    coneAngleVector.push_back( unit_conversions::convertDegreesToRadians( 80.0 ) );
    coneAngleVector.push_back( unit_conversions::convertDegreesToRadians( - 80.0 ) );
    coneAngleVector.push_back( unit_conversions::convertDegreesToRadians( 15.0 ) );
    coneAngleVector.push_back( unit_conversions::convertDegreesToRadians( 15.0 ) );
    coneAngleVector.push_back( unit_conversions::convertDegreesToRadians( 15.0 ) );
    coneAngleVector.push_back( unit_conversions::convertDegreesToRadians( 15.0 ) );


    // Define set of clock angles for various test cases [rad].
    std::vector< double  > clockAngleVector;
    clockAngleVector.push_back( 0.0 );
    clockAngleVector.push_back( mathematical_constants::PI / 2.0 );
    clockAngleVector.push_back( mathematical_constants::PI );
    clockAngleVector.push_back( 2.0 * mathematical_constants::PI );
    clockAngleVector.push_back( 0.0 );
    clockAngleVector.push_back( mathematical_constants::PI / 2.0 );
    clockAngleVector.push_back( mathematical_constants::PI );
    clockAngleVector.push_back( 2.0 * mathematical_constants::PI );
    clockAngleVector.push_back( 0.0 );
    clockAngleVector.push_back( mathematical_constants::PI / 2.0 );
    clockAngleVector.push_back( mathematical_constants::PI );
    clockAngleVector.push_back( 2.0 * mathematical_constants::PI );
    clockAngleVector.push_back( mathematical_constants::PI / 4.0 );
    clockAngleVector.push_back( - mathematical_constants::PI / 4.0 );
    clockAngleVector.push_back( - mathematical_constants::PI / 4.0 );
    clockAngleVector.push_back( mathematical_constants::PI / 4.0 );
    clockAngleVector.push_back( unit_conversions::convertDegreesToRadians( 80.0 ) );
    clockAngleVector.push_back( unit_conversions::convertDegreesToRadians( 80.0 ) );
    clockAngleVector.push_back( unit_conversions::convertDegreesToRadians( 25.0 ) );
    clockAngleVector.push_back( unit_conversions::convertDegreesToRadians( 25.0 )  );
    clockAngleVector.push_back( unit_conversions::convertDegreesToRadians( 25.0 ) );
    clockAngleVector.push_back( unit_conversions::convertDegreesToRadians( 25.0 ) );


    // Define position vector.
    Eigen::Vector3d positionVectorToSource;

    // Define velocity vector of the spacecraft w.r.t. central body.
    Eigen::Vector3d velocityVector;


    for ( int testCase = 0 ; testCase < 20 ; testCase++ ){

        // Define different position and velocity vectors for various test cases.
        if ( testCase < 19 )
        {
            // Set position vector [m].
            positionVectorToSource = Eigen::Vector3d( astronomicalUnitInMeters, astronomicalUnitInMeters, 0.0 );

            // Compute distance between spacecraft and radiation pressure source [m].
            double positionVectorToSourceNorm = positionVectorToSource.norm();

            // Set velocity vector of the spacecraft w.r.t. central body [m/s].
            velocityVector = std::sqrt( earthGravitationalParameter / positionVectorToSourceNorm )
                                       * ( Eigen::Vector3d() << std::sqrt( 2.0 ) / 2.0, std::sqrt( 2.0 ) / 2.0, 0.0 ).finished();
        }
        else if ( testCase == 19 )
        {
            // Set position vector [m].
            positionVectorToSource = 2.0 * Eigen::Vector3d( astronomicalUnitInMeters, astronomicalUnitInMeters, astronomicalUnitInMeters );

            // Compute distance between spacecraft and radiation pressure source [m].
            double positionVectorToSourceNorm = positionVectorToSource.norm();

            // Set velocity vector of the spacecraft w.r.t. central body [m/s].
            velocityVector = std::sqrt( earthGravitationalParameter / positionVectorToSourceNorm )
                                       * ( Eigen::Vector3d() << 1.0, 1.0, 1.0 ).finished();
        }
        else if ( testCase == 20 )
        {
            // Set position vector [m].
            positionVectorToSource = Eigen::Vector3d( astronomicalUnitInMeters / 4.0, astronomicalUnitInMeters / 3.0,
                                                      astronomicalUnitInMeters / 2.0 );

            // Compute distance between spacecraft and radiation pressure source [m].
            double positionVectorToSourceNorm = positionVectorToSource.norm();

            // Set velocity vector of the spacecraft w.r.t. central body [m/s].
            velocityVector = std::sqrt( earthGravitationalParameter / positionVectorToSourceNorm )
                                       * ( Eigen::Vector3d() << 1.0 / 4.0, 1.0/ 3.0, 1.0 / 2.0 ).finished();
        }
        else if ( testCase == 21 ){

            // Set position vector [m].
            positionVectorToSource = Eigen::Vector3d( astronomicalUnitInMeters / 3.5, astronomicalUnitInMeters / 3.0,
                                                      astronomicalUnitInMeters / 2.5 );

            // Compute distance between spacecraft and radiation pressure source [m].
            double positionVectorToSourceNorm = positionVectorToSource.norm();

            // Set velocity vector of the spacecraft w.r.t. central body [m/s].
            velocityVector = std::sqrt( earthGravitationalParameter / positionVectorToSourceNorm )
                                       * ( Eigen::Vector3d() << 1.0 / 2.5, 1.0 / 2.6, 1.0 / 2.9 ).finished();
        }


        // Set radiation pressure at target [N/m^2].
        const double radiationPressureAtTarget = radiationPressureAtOneAU * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

        // Normalize position vector to get vector pointing to source in non-dimensional units [-].
        positionVectorToSource.normalize( );

        // Set area on target that is subject to radiation pressure [m^2].
        const double areaSubjectToRadiationPressure = 0.5;

        // Define solar sail model settings.
        const double frontEmissivityCoefficient = 0.05;
        const double backEmissivityCoefficient = 0.64;
        const double frontLambertianCoefficient = 0.79;
        const double backLambertianCoefficient = 0.55;
        const double reflectivityCoefficient = 0.88;
        const double specularReflectionCoefficient = 0.94;

        // Compute force from solar sail model [N].
        Eigen::Vector3d computedSolarSailForce = electromagnetism::computeSolarSailForce(
                    frontEmissivityCoefficient, backEmissivityCoefficient, frontLambertianCoefficient,
                    backLambertianCoefficient, reflectivityCoefficient, specularReflectionCoefficient, positionVectorToSource,
                    velocityVector.normalized(), radiationPressureAtTarget, areaSubjectToRadiationPressure,
                    coneAngleVector[ testCase ], clockAngleVector[ testCase ] );

        // Compare computed and expected radiation pressure force vectors.
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( expectedSolarSailForceVector[ testCase ]( j ) - computedSolarSailForce( j ) ), 1.0E-20 );
        }

    }


}


//! Test radiation pressure acceleration model update-members function, for changes in the coefficients also.
BOOST_AUTO_TEST_CASE( testSolarSailAccelerationModelClassUpdateMembersIncludingCoefficients )
{

    // Declare and initialize solar sail radiation acceleration model.
    electromagnetism::SolarSailAccelerationPointer solarSailModel
        = std::make_shared< electromagnetism::SolarSailAcceleration >(
            &getSourcePosition, &getAcceleratedBodyPosition, &getAcceleratedBodyVelocity,
            &getCentralBodyVelocity, &getRadiationPressure, &getConeAngle, &getClockAngle,
            &getFrontEmissivityCoefficient, &getBackEmissivityCoefficient,
            &getFrontLambertianCoefficient, &getBackLambertianCoefficient,
            &getReflectivityCoefficient, &getSpecularReflectionCoefficient,
            &getAreaSubjectToRadiationPressure, &getMassOfAcceleratedBody );
    solarSailModel->updateMembers( 0.0 );

    // Set expected radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d expectedSolarSailAcceleration = Eigen::Vector3d( 8.59654650354177e-06, 7.44856903789465e-06, 6.83437564207585e-06 );

    // Set the distance from the Sun to Venus [m].
    const double distanceSunToVenus = 0.732 * astronomicalUnitInMeters;

    // Update position of accelerated body [m].
    acceleratedBodyPosition = Eigen::Vector3d( distanceSunToVenus * 0.9, distanceSunToVenus * 0.8, distanceSunToVenus * 0.7 );
    const double acceleratedBodyPositionNorm = acceleratedBodyPosition.norm();
    acceleratedBodyVelocity = std::sqrt( earthGravitationalParameter / acceleratedBodyPositionNorm )
            * Eigen::Vector3d( 1.0 / std::sqrt( 2.0 ), 1.0 / std::sqrt( 2.0 ), 1.0 / std::sqrt( 2.0 ) );

    // Update radiation pressure at location of accelerated body [N/m^2].
    radiationPressure = radiationPressureAtOneAU * astronomicalUnitInMeters * astronomicalUnitInMeters / getVectorToSource( ).squaredNorm( );

    //Update cone angle [rad].
    coneAngle = 15.0 * mathematical_constants::PI / 180.0;

    //Update clock angle [rad].
    clockAngle = 25.0 * mathematical_constants::PI / 180.0;

    //Update front emissivity coefficient [-].
    frontEmissivityCoefficient=0.03;

    //Update back emissivity coefficient [-].
    backEmissivityCoefficient = 0.5;

    //Update front Lambertian coefficient [-].
    frontLambertianCoefficient = 0.70;

    //Update front Lambertian coefficient [-].
    backLambertianCoefficient = 0.4;

    //Update specular reflection coefficient [-].
    specularReflectionCoefficient = 0.9;

    //Update reflectivity coefficient [-].
    reflectivityCoefficient = 0.6;

    // Update area subject to radiation pressure [m^2].
    areaSubjectToRadiationPressure = 0.005;

    // Update mass of accelerated body [kg].
    massOfAcceleratedBody = 0.0022;

    // Update class members.
    solarSailModel->updateMembers( 0.0 );

    // Compute radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d computedSolarSailAcceleration = solarSailModel->getAcceleration( );

    // Compare computed and expected radiation pressure acceleration vectors.
    BOOST_CHECK_SMALL( std::fabs( expectedSolarSailAcceleration.x() - computedSolarSailAcceleration.x()), 1.0e-19);
    BOOST_CHECK_SMALL( std::fabs( expectedSolarSailAcceleration.y() - computedSolarSailAcceleration.y()), 1.0e-19 );
    BOOST_CHECK_SMALL( std::fabs( expectedSolarSailAcceleration.z() - computedSolarSailAcceleration.z()), 1.0e-19 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
