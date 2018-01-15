/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Mooij, E. The Motion of a Vehicle in a Planetary Atmosphere, TU Delft, 1997.
 *
 *    Notes
 *      The reference frame definitions/abbreviations can be found in the file
 *      referenceFrameTransformations.h.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{
namespace unit_tests
{

using namespace unit_conversions;

BOOST_AUTO_TEST_SUITE( test_reference_frame_transformations )

// Summary of tests.
// Test 1: Test inertial to rotating planetocentric frame transformation.
// Test 2: Check whether the transformed matrix of Test 1 is also correct.
// Test 3: Same test as Test 1 for the transformation quaternion.
// Test 4: Same test as Test 2 for the transformation quaternion.
// Test 5: Test airspeed-based aerodynamic to body frame transformation.
// Test 6: Test planetocentric to local vertical frame transformation quaternion.
// Test 7: Check whether the transformed matrix of Test 6 is also correct.
// Test 8: Test transformations between trajectory and local vertical frame.
// Test 9: Test transformations between trajectory and aerodynamic frame.
// Test 10: Test transformations between body and aerodynamic reference frame.

// Test inertial to rotating planetocentric frame transformations.
BOOST_AUTO_TEST_CASE( testRotatingPlanetocentricFrameTransformations )
{
    // Using declarations.
    using std::atan2;
    using std::cos;
    using std::sin;
    using std::pow;
    using std::sqrt;


    // Test 1: Test Inertial to rotating planetocentric frame transformation.
    {
        // Initialize initial location vector in inertial frame.
        Eigen::Vector3d startLocation;
        startLocation( 0 ) = 4.0;
        startLocation( 1 ) = 1.0;
        startLocation( 2 ) = 5.0;

        double horizontalStartLocationSize = sqrt( pow( startLocation( 0 ), 2.0 )
                                                   + pow( startLocation( 1 ), 2.0 ) );

        // Declare and initialize angle between vector and XR-axis.
        double startAngle = atan2( startLocation( 1 ), startLocation( 0 ) );

        // Rotate by 10 degrees around the positive Z-axis
        double angleInTime = convertDegreesToRadians( 10.0 );

        // Declare and initialize the angle between the XI- and XR-axis.
        double endAngle = startAngle - angleInTime;

        // Declare the expected location of the point in the planetocentric reference frame.
        Eigen::Vector3d expectedLocation;
        expectedLocation( 0 ) = horizontalStartLocationSize * cos( endAngle );
        expectedLocation( 1 ) = horizontalStartLocationSize * sin( endAngle );
        expectedLocation( 2 ) = startLocation( 2 );

        // Compute location of the point in the rotating frame subject to the transformation
        // matrix.
        Eigen::Vector3d transformedLocation;
        transformedLocation = reference_frames::
                getInertialToPlanetocentricFrameTransformationMatrix( angleInTime )
                * startLocation;

        // Check whether both vectors are equal within tolerances.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transformedLocation, expectedLocation,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test 2: Check whether the transformed matrix of Test 1 is also correct.
    // Compute the error in the calculation.
    {
        // Initialize initial location vector in inertial frame.
        Eigen::Vector3d startLocation;
        startLocation( 0 ) = 4.0;
        startLocation( 1 ) = 1.0;
        startLocation( 2 ) = 5.0;

        // Rotate by 10 degrees around the positive Z-axis
        double angleInTime = convertDegreesToRadians( 10.0 );

        // Compute location of the point in the rotating frame subject to the transformation matrix.
        Eigen::Vector3d transformedLocation;
        transformedLocation = reference_frames::
                getRotatingPlanetocentricToInertialFrameTransformationMatrix( angleInTime ) *
                reference_frames::
                getInertialToPlanetocentricFrameTransformationMatrix( angleInTime )
                * startLocation;

        // Check whether both vectors are equal within tolerances.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transformedLocation, startLocation,
                                           std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_CASE( testRotatingPlanetocentricFrameTransformationQuaternion )
{
    // Using declarations.
    using std::atan2;
    using std::cos;
    using std::sin;
    using std::pow;
    using std::sqrt;


    // Test 3: Same test as Test 1 for the transformation quaternion.
    // Compute location of the point in the Rotating frame subject to the transformation matrix.
    {
        // Initialize initial location vector in inertial frame.
        Eigen::Vector3d startLocation;
        startLocation( 0 ) = 4.0;
        startLocation( 1 ) = 1.0;
        startLocation( 2 ) = 5.0;

        double horizontalStartLocationSize = sqrt( pow( startLocation( 0 ), 2.0 )
                                                   + pow( startLocation( 1 ), 2.0 ) );

        // Declare and initialize angle between vector and XR-axis.
        double startAngle = atan2( startLocation( 1 ), startLocation( 0 ) );

        // Rotate by 10 degrees around the positive Z-axis
        double angleInTime = convertDegreesToRadians( 10.0 );

        // Declare and initialize the angle between the XI- and XR-axis.
        double endAngle = startAngle - angleInTime;

        // Declare the expected location of the point in the planetocentric reference frame.
        Eigen::Vector3d expectedLocation;
        expectedLocation( 0 ) = horizontalStartLocationSize * cos( endAngle );
        expectedLocation( 1 ) = horizontalStartLocationSize * sin( endAngle );
        expectedLocation( 2 ) = startLocation( 2 );

        // Compute location of the point in the rotating frame subject to the transformation matrix.
        Eigen::Vector3d transformedLocation;
        transformedLocation = reference_frames::
                getInertialToPlanetocentricFrameTransformationQuaternion( angleInTime )
                * startLocation;

        // Check whether both vectors are equal within tolerances.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transformedLocation, expectedLocation,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test 4: Same test as Test 2 for the transformation quaternion.
    // Compute the error in the calculation.
    {
        // Initialize initial location vector in inertial frame.
        Eigen::Vector3d startLocation;
        startLocation( 0 ) = 4.0;
        startLocation( 1 ) = 1.0;
        startLocation( 2 ) = 5.0;

        // Rotate by 10 degrees around the positive Z-axis
        double angleInTime = convertDegreesToRadians( 10.0 );

        // Compute location of the point in the rotating frame subject to the transformation matrix.
        Eigen::Vector3d transformedLocation;
        transformedLocation = reference_frames::
                getRotatingPlanetocentricToInertialFrameTransformationQuaternion( angleInTime ) *
                reference_frames::
                getInertialToPlanetocentricFrameTransformationQuaternion( angleInTime )
                * startLocation;

        // Check whether both vectors are equal within tolerances.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transformedLocation, startLocation,
                                           std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_CASE( testAirspeedBasedAerodynamicToBodyFrameTransformation )
{
    // Using declarations.
    using std::atan2;
    using std::cos;


    // Test 5: Test airspeed-Based Aerodynamic to body frame transformation.
    // Declare and initialize the start location and angles.
    {
        // Initialize initial location vector in inertial frame.
        Eigen::Vector3d startLocation;
        startLocation( 0 ) = -10.0;
        startLocation( 1 ) = 0.0;
        startLocation( 2 ) = 0.0;

        double angleOfAttack = convertDegreesToRadians( 45.0 );
        double angleOfSideslip = convertDegreesToRadians( 60.0 );;

        // Compute expected location.
        // As there is only an angle of attack, the following simplified equations can be used.
        Eigen::Vector3d expectedLocation;
        expectedLocation( 0 ) = startLocation( 0 ) * cos( angleOfSideslip ) * cos( angleOfAttack );
        expectedLocation( 1 ) = startLocation( 0 ) * sin( angleOfSideslip );
        expectedLocation( 2 ) = startLocation( 0 ) * cos( angleOfSideslip ) * sin( angleOfAttack );

        // Compute location of the point in the rotating frame subject to the transformation matrix.
        Eigen::Vector3d transformedLocation;
        transformedLocation = reference_frames::
                getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
                    angleOfAttack, angleOfSideslip ) * startLocation;

        // Check whether both vectors are equal within tolerances.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transformedLocation, expectedLocation, 1.0e-14 );
    }
}

BOOST_AUTO_TEST_CASE( testRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion )
{
    // Using declarations.
    using std::atan2;
    using std::cos;
    using std::sin;


    // Test 6: Test Rotating planetocentric to local vertical frame transformation quaternion.
    {
        // Initialize initial location vector.
        Eigen::Vector3d startLocation;
        startLocation( 0 ) = 10.0;
        startLocation( 1 ) = 5.0;
        startLocation( 2 ) = 2.0;

        // Declare rotation angle for planet and set to 90 degrees.
        double longitude = convertDegreesToRadians( 60.0 );

        // Declare latitude and set to 45 degrees.
        double latitude = convertDegreesToRadians( 20.0 );

        // Declare the expected location of the point in the planet reference frame.
        Eigen::Vector3d expectedLocation;
        expectedLocation( 0 ) = -cos( longitude ) * sin( latitude ) * startLocation( 0 ) -
                sin( longitude ) * sin( latitude ) * startLocation( 1 )  +
                cos( latitude ) * startLocation( 2 );
        expectedLocation( 1 ) = -sin( longitude ) * startLocation( 0 ) +
                cos( longitude ) * startLocation( 1 ) +
                0.0;
        expectedLocation( 2 ) = -cos( longitude ) * cos( latitude ) * startLocation( 0 ) -
                sin( longitude ) * cos( latitude ) * startLocation( 1 )  -
                sin( latitude ) * startLocation( 2 );

        // Compute location of the point in the rotating frame subject to the transformation matrix.
        Eigen::Vector3d transformedLocation;
        transformedLocation = reference_frames::
                getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
                    longitude, latitude ) * startLocation;

        // Check whether both vectors are equal within tolerances.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transformedLocation, expectedLocation, 1.0e-14 );
    }

    // Test 7: Check whether the transformed matrix of Test 6 is also correct.
    // Compute the error in the calculation.
    {
        // Initialize initial location vector.
        Eigen::Vector3d startLocation;
        startLocation( 0 ) = 10.0;
        startLocation( 1 ) = 5.0;
        startLocation( 2 ) = 2.0;

        // Declare rotation angle for planet and set to 90 degrees.
        double longitude = convertDegreesToRadians( 60.0 );

        // Declare latitude and set to 45 degrees.
        double latitude = convertDegreesToRadians( 20.0 );

        // Compute location of the point in the rotating frame subject to the transformation matrix.
        Eigen::Vector3d transformedLocation;
        transformedLocation = reference_frames::
                getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion( longitude,
                                                                                       latitude )
                * reference_frames::
                getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion( longitude,
                                                                                       latitude )
                * startLocation;

        // Check whether both vectors are equal within tolerances.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transformedLocation, startLocation, 1.0e-14 );
    }
}

// Test inertial to rotating planetocentric frame transformations, without assumption of equal x-y
// planes.
BOOST_AUTO_TEST_CASE( testRotatingPlanetocentricWithEquatorChangeFrameTransformations )
{
    using namespace reference_frames;

    // Data from pck00010.tpc Spice kernel.
    const double venusPoleRightAscension = convertDegreesToRadians( 272.76 );
    const double venusPoleDeclination = convertDegreesToRadians( 67.16 );
    const double venusPrimeMeridianAtJ2000 = convertDegreesToRadians( 160.20 );

    // Set rotation, as calculated with Spice
    // (see Astrodynamics/Ephemerides/unitTestSimpleRotationalEphemeris.cpp).
    Eigen::Matrix3d spiceInitialRotationToTargetFrameMatrix;
    spiceInitialRotationToTargetFrameMatrix
            << -0.9548214974296336, 0.2665104385944917, 0.1314841974018291,
            -0.296591573568662, -0.882413772579987, -0.3652114078848295,
            0.01869081416890202, -0.3877088083617989, 0.9215923900425705;

    // Calculate given rotation to Venus-fixed frame.
    Eigen::Matrix3d calculatedTransformationMatrix = Eigen::Matrix3d(
                getInertialToPlanetocentricFrameTransformationQuaternion(
                    venusPoleDeclination, venusPoleRightAscension, venusPrimeMeridianAtJ2000 ) );

    // Compare calculated with Spice result.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                spiceInitialRotationToTargetFrameMatrix, calculatedTransformationMatrix, 1.0E-15 );

    // Calculate inverse rotation.
    calculatedTransformationMatrix = Eigen::Matrix3d(
                getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
                    venusPoleDeclination, venusPoleRightAscension,
                    venusPrimeMeridianAtJ2000 ) );

    // Multiply inverse rotation with Spice rotation to Venus-fixed frame. The two should be
    // orthonormal, so result should be the identity matrix.
    Eigen::Matrix3d matrixToCheck =
            spiceInitialRotationToTargetFrameMatrix * calculatedTransformationMatrix;

    // Check whether result is indeed (close to) identity matrix.
    for ( unsigned int i = 0; i < 3; i++ )
    {
        for ( unsigned int j = 0; j < 3; j++ )
        {
            if ( i == j )
            {
                BOOST_CHECK_CLOSE_FRACTION( matrixToCheck( i, j ), 1.0, 1.0E-15 );
            }
            else
            {
                BOOST_CHECK_SMALL( matrixToCheck( i, j ), 1.0E-15 );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testTrajectoryToLocalVerticalFrameTransformations )
{
    // Using declarations.
    using std::cos;
    using std::sin;


    // Test 8: Test trajectory to local vertical frame and inverse transformations.
    {
        // Initialize initial location vector.
        Eigen::Vector3d startLocation;
        startLocation( 0 ) = 10.0;
        startLocation( 1 ) = 5.0;
        startLocation( 2 ) = 2.0;

        // Declare rotation angles.
        for ( double headingAngle = convertDegreesToRadians( 0.0 );
              headingAngle < convertDegreesToRadians( 360.0 );
              headingAngle += convertDegreesToRadians( 30.0 ) )
        {
            for ( double flightPathAngle = convertDegreesToRadians( -90.0 );
                  flightPathAngle < convertDegreesToRadians( 90.0 );
                  flightPathAngle += convertDegreesToRadians( 30.0 ) )
            {
                // Declare the expected location of the point in the V-frame (per Mooij 1997).
                Eigen::Vector3d expectedLocation;
                Eigen::Matrix3d rotation;
                rotation << cos( headingAngle ) * cos( flightPathAngle ),
                        -sin( headingAngle ),
                        cos( headingAngle ) * sin( flightPathAngle ),
                        sin( headingAngle ) * cos( flightPathAngle ),
                        cos( headingAngle ),
                        sin( headingAngle ) * sin( flightPathAngle ),
                        -sin( flightPathAngle ),
                        0.0,
                        cos( flightPathAngle );
                expectedLocation = rotation * startLocation;

                // Compute location of the point in the V-frame using the tested function.
                Eigen::Vector3d transformedLocationQuat =
                        reference_frames::
                        getTrajectoryToLocalVerticalFrameTransformationQuaternion(
                            flightPathAngle, headingAngle ) * startLocation;
                Eigen::Vector3d transformedLocationMat =
                        reference_frames::getTrajectoryToLocalVerticalFrameTransformationMatrix(
                            flightPathAngle, headingAngle ) * startLocation;

                // Compute product of inverse rotations.
                Eigen::Vector3d inverseRotationProductMat =
                        reference_frames::getTrajectoryToLocalVerticalFrameTransformationMatrix(
                            flightPathAngle, headingAngle ) *
                        reference_frames::getLocalVerticalFrameToTrajectoryTransformationMatrix(
                            flightPathAngle, headingAngle ) * startLocation;
                Eigen::Vector3d inverseRotationProductQuat =
                        reference_frames::
                        getTrajectoryToLocalVerticalFrameTransformationQuaternion(
                            flightPathAngle, headingAngle ) *
                        reference_frames::
                        getLocalVerticalFrameToTrajectoryTransformationQuaternion(
                            flightPathAngle, headingAngle ) * startLocation;

                // Check whether transformed vectors match within tolerances.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            transformedLocationQuat, expectedLocation, 1.0e-14 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            transformedLocationMat, expectedLocation, 2.0e-14 );

                // Check whether product of inverse rotations equals unity.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            inverseRotationProductMat, startLocation, 1.0e-14 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            inverseRotationProductQuat, startLocation, 1.0e-14 );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testTrajectoryToAerodynamicFrameTransformations )
{
    // Using declarations.
    using std::cos;
    using std::sin;


    // Test 9: Test trajectory to aerodynamic and inverse transformations.
    {
        // Initialize initial location vector.
        Eigen::Vector3d startLocation;
        startLocation( 0 ) = 10.0;
        startLocation( 1 ) = 5.0;
        startLocation( 2 ) = 2.0;

        // Declare rotation angles.
        for ( double bankAngle = convertDegreesToRadians( 0.0 );
              bankAngle < convertDegreesToRadians( 360.0 );
              bankAngle += convertDegreesToRadians( 30.0 ) )
        {
            // Declare the expected location of the point in the aerodynamic reference frame
            // (per Mooij 1997).
            Eigen::Vector3d expectedLocation;
            Eigen::Matrix3d rotation;
            rotation << 1.0, 0.0, 0.0,
                    0.0, cos( bankAngle ), -sin( bankAngle ),
                    0.0, sin( bankAngle ), cos( bankAngle );
            expectedLocation = rotation * startLocation;

            // Compute location of the point in the rotating frame subject to the
            // transformation matrix.
            Eigen::Vector3d transformedLocationQuat =
                    reference_frames::getTrajectoryToAerodynamicFrameTransformationQuaternion(
                        bankAngle ) * startLocation;
            Eigen::Vector3d transformedLocationMat =
                    reference_frames::getTrajectoryToAerodynamicFrameTransformationMatrix(
                        bankAngle ) * startLocation;

            // Compute product of inverse rotations.
            Eigen::Vector3d inverseRotationProductMat =
                    reference_frames::getAerodynamicToTrajectoryFrameTransformationMatrix(
                        bankAngle ) *
                    reference_frames::getTrajectoryToAerodynamicFrameTransformationMatrix(
                        bankAngle ) * startLocation;
            Eigen::Vector3d inverseRotationProductQuat =
                    reference_frames::getAerodynamicToTrajectoryFrameTransformationQuaternion(
                        bankAngle ) *
                    reference_frames::getTrajectoryToAerodynamicFrameTransformationQuaternion(
                        bankAngle ) * startLocation;

            // Check whether transformed vectors match within tolerances.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        transformedLocationQuat, expectedLocation, 1.0e-14 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transformedLocationMat, expectedLocation, 1.0e-14 );

            // Check whether product of inverse rotations equals unity.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( inverseRotationProductMat, startLocation, 1.0e-14 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        inverseRotationProductQuat, startLocation, 1.0e-14 );
        }
    }
}

BOOST_AUTO_TEST_CASE( testAerodynamicToBodyFrameTransformations )
{
    // Using declarations.
    using std::cos;
    using std::sin;


    // Test 10: Test body to airspeed-based aerodynamic frame and inverse transformations.
    {
        // Initialize initial location vector.
        Eigen::Vector3d startLocation;
        startLocation( 0 ) = 10.0;
        startLocation( 1 ) = 5.0;
        startLocation( 2 ) = 2.0;

        // Declare rotation angles.
        for ( double angleOfAttack = convertDegreesToRadians( -180.0 );
              angleOfAttack < convertDegreesToRadians( 180.0 );
              angleOfAttack += convertDegreesToRadians( 30.0 ) )
        {
            for ( double angleOfSideslip = convertDegreesToRadians( -90.0 );
                  angleOfSideslip < convertDegreesToRadians( 90.0 );
                  angleOfSideslip += convertDegreesToRadians( 30.0 ) )
            {
                // Declare the expected location of the point in the AA-frame.
                // (per Mooij 1997).
                Eigen::Vector3d expectedLocation;
                Eigen::Matrix3d rotation;
                rotation << cos( angleOfAttack ) * cos( angleOfSideslip ),
                        sin( angleOfSideslip ),
                        sin( angleOfAttack ) * cos( angleOfSideslip ),
                        -cos( angleOfAttack ) * sin( angleOfSideslip ),
                        cos( angleOfSideslip ),
                        -sin( angleOfAttack ) * sin( angleOfSideslip ),
                        -sin( angleOfAttack ),
                        0.0,
                        cos( angleOfAttack );
                expectedLocation = rotation * startLocation;

                // Compute location of the point in the AA-frame using the tested function.
                Eigen::Vector3d transformedLocationQuat =
                        reference_frames::
                        getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
                            angleOfAttack, angleOfSideslip ) * startLocation;
                Eigen::Vector3d transformedLocationMat =
                        reference_frames::
                        getBodyToAirspeedBasedAerodynamicFrameTransformationMatrix(
                            angleOfAttack, angleOfSideslip ) * startLocation;

                // Compute product of inverse rotations
                Eigen::Vector3d inverseRotationProductMat =
                        reference_frames::
                        getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
                            angleOfAttack, angleOfSideslip ) *
                        reference_frames::
                        getBodyToAirspeedBasedAerodynamicFrameTransformationMatrix(
                            angleOfAttack, angleOfSideslip ) * startLocation;
                Eigen::Vector3d inverseRotationProductQuat =
                        reference_frames::
                        getAirspeedBasedAerodynamicToBodyFrameTransformationQuaternion(
                            angleOfAttack, angleOfSideslip ) *
                        reference_frames::
                        getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
                            angleOfAttack, angleOfSideslip ) * startLocation;

                // Check whether transformed vectors match within tolerances.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            transformedLocationQuat, expectedLocation, 1.0e-14 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            transformedLocationMat, expectedLocation, 1.0e-14 );

                // Check whether product of inverse rotations equals unity.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            inverseRotationProductMat, startLocation, 1.0e-14 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            inverseRotationProductQuat, startLocation, 1.0e-14 );
            }
        }
    }
}

//! Function to test if a matrix is an identity matrix.
void testIsIdentityMatrix( const Eigen::Matrix3d matrixToTest )
{
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            if( i == j )
            {
                BOOST_CHECK_CLOSE_FRACTION( matrixToTest( i, j ), 1.0, std::numeric_limits< double >::epsilon( ) );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( matrixToTest( i, j ) ), std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}

//! Function to test if a rotation matrix defines a right-handed frame.
void testIsRotationMatrixRightHanded( const Eigen::Matrix3d matrixToTest )
{
    Eigen::Vector3d expectedUnitVectorZ =
           ( Eigen::Vector3d( matrixToTest.block( 0, 0, 3, 1 ) ) ).cross(
                Eigen::Vector3d( matrixToTest.block( 0, 1, 3, 1 ) ) );
    Eigen::Vector3d expectedUnitVectorX =
           ( Eigen::Vector3d( matrixToTest.block( 0, 1, 3, 1 ) ) ).cross(
                Eigen::Vector3d( matrixToTest.block( 0, 2, 3, 1 ) ) );
    Eigen::Vector3d expectedUnitVectorY  =
           ( Eigen::Vector3d( matrixToTest.block( 0, 2, 3, 1 ) ) ).cross(
                Eigen::Vector3d( matrixToTest.block( 0, 0, 3, 1 ) ) );

    for( unsigned int i = 0; i < 3; i++ )\
    {
        BOOST_CHECK_CLOSE_FRACTION(
                    expectedUnitVectorX( i ), matrixToTest( i, 0 ), 2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION(
                    expectedUnitVectorY( i ), matrixToTest( i, 1 ), 2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION(
                    expectedUnitVectorZ( i ), matrixToTest( i, 2 ), 2.0 * std::numeric_limits< double >::epsilon( ) );
    }
}

// Test velocity based LVLH frame transformation.
BOOST_AUTO_TEST_CASE( testVelocityBasedLvlhFrameTransformations )
{
    // Using declarations.
    using std::atan2;
    using std::cos;
    using std::sin;
    using std::pow;
    using std::sqrt;

    // Test 11: Test velocity=based LVLH frame to inertial (I) frame transformation: for arbitrary state vectors,
    // test all expected properties of rotation matrices
    {
        // Position vectors to semi-random values
        Eigen::Vector6d vehicleStateCartesian, centralBodyStateCartesian, relativeState;
        vehicleStateCartesian << 3.2, -1.7, 8.2, -1.4E-3, -5.6E-4, 9.2E-4;
        centralBodyStateCartesian << 7.4, 6.3, -3.6, 5.3E-4, 7.64E3, -4.4E-4;

        Eigen::Matrix3d nAxisAwayFromBodyMatrix = reference_frames::getVelocityBasedLvlhToInertialRotation(
                    vehicleStateCartesian, centralBodyStateCartesian, true );
        Eigen::Matrix3d nAxisAwayTowardsBodyMatrix = reference_frames::getVelocityBasedLvlhToInertialRotation(
                    vehicleStateCartesian, centralBodyStateCartesian, false );

        // Test if central body is properly processed.
        {
            Eigen::Matrix3d nAxisAwayFromBodyMatrixCentral = reference_frames::getVelocityBasedLvlhToInertialRotation(
                        vehicleStateCartesian - centralBodyStateCartesian, Eigen::Vector6d::Zero( ), true );
            Eigen::Matrix3d nAxisAwayTowardsBodyMatrixCentral = reference_frames::getVelocityBasedLvlhToInertialRotation(
                        vehicleStateCartesian - centralBodyStateCartesian, Eigen::Vector6d::Zero( ), false );

            for( unsigned int i = 0; i < 3; i++ )
            {
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_CLOSE_FRACTION(
                                nAxisAwayFromBodyMatrixCentral( i, j ), nAxisAwayFromBodyMatrix( i, j ),
                                2.0 * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_CLOSE_FRACTION(
                                nAxisAwayTowardsBodyMatrix( i, j ), nAxisAwayTowardsBodyMatrixCentral( i, j ),
                                2.0 * std::numeric_limits< double >::epsilon( ) );
                }

            }
        }

        // Compute Relative states
        relativeState = vehicleStateCartesian - centralBodyStateCartesian;

        // Test if n axis indeed points away from body
        double positionDotProductWithNAwayFromBody =
                Eigen::Vector3d( nAxisAwayFromBodyMatrix.block( 0, 1, 3, 1 ) ).dot(
                    relativeState.segment( 0, 3 ).normalized( ) );
        BOOST_CHECK_EQUAL( positionDotProductWithNAwayFromBody > 0, 1 );

        // Test if n axis indeed points towards body
        double positionDotProuctWithNTowardsBody =
                Eigen::Vector3d( nAxisAwayTowardsBodyMatrix.block( 0, 1, 3, 1 ) ).dot(
                    relativeState.segment( 0, 3 ).normalized( ) );
        BOOST_CHECK_EQUAL( positionDotProuctWithNTowardsBody < 0, 1 );

        // Check if axes point exactly in opposite directions
        BOOST_CHECK_CLOSE_FRACTION( positionDotProductWithNAwayFromBody, -positionDotProuctWithNTowardsBody,
                                    std::numeric_limits< double >::epsilon( ) );

        // Test if two matrices compare as they should
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_CLOSE_FRACTION(
                            nAxisAwayFromBodyMatrix( i, j ), ( j == 0 ? 1.0 : -1.0 ) * nAxisAwayTowardsBodyMatrix( i, j ),
                            2.0 * std::numeric_limits< double >::epsilon( ) );
            }
        }

        // Test orthonormality of matrices
        Eigen::Matrix3d matrixProduct;
        matrixProduct = nAxisAwayFromBodyMatrix * nAxisAwayFromBodyMatrix.transpose( );
        testIsIdentityMatrix( matrixProduct );

        matrixProduct = nAxisAwayTowardsBodyMatrix * nAxisAwayTowardsBodyMatrix.transpose( );
        testIsIdentityMatrix( matrixProduct );

        matrixProduct = nAxisAwayFromBodyMatrix * nAxisAwayFromBodyMatrix.inverse( );
        testIsIdentityMatrix( matrixProduct );

        matrixProduct = nAxisAwayTowardsBodyMatrix * nAxisAwayTowardsBodyMatrix.inverse( );
        testIsIdentityMatrix( matrixProduct );

        // Test if matrices define right-handed frames
        testIsRotationMatrixRightHanded( nAxisAwayFromBodyMatrix );
        testIsRotationMatrixRightHanded( nAxisAwayTowardsBodyMatrix );

        Eigen::Vector3d lvlhFrameThrust1, lvlhFrameThrust2;

        // Define semi-random thrust vector along velocity vector
        double thrustScale = 25.2;
        Eigen::Vector3d inertialThrustVectorAlongVelocity = thrustScale * relativeState.segment( 3, 3 );
        {
            // Compute local thrust vectors
            lvlhFrameThrust1 =
                    nAxisAwayFromBodyMatrix.transpose( ) * inertialThrustVectorAlongVelocity;
            lvlhFrameThrust2 =
                    nAxisAwayTowardsBodyMatrix.transpose( ) * inertialThrustVectorAlongVelocity;

            // Compute thrust vector magnitudes
            double thrustMagnitudeAlongVelocityVector = inertialThrustVectorAlongVelocity.norm( );

            // Test if entry 0 of two vectors are the same
            BOOST_CHECK_CLOSE_FRACTION( lvlhFrameThrust1( 0 ), lvlhFrameThrust2( 0 ),
                                        std::numeric_limits< double >::epsilon( ) );

            // Test thrust vector norms.
            BOOST_CHECK_CLOSE_FRACTION( thrustMagnitudeAlongVelocityVector, lvlhFrameThrust1.norm( ),
                                        5.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( thrustMagnitudeAlongVelocityVector, lvlhFrameThrust2.norm( ),
                                        5.0 * std::numeric_limits< double >::epsilon( ) );

            // Test if entry 0 magnitude is equal to norm
            BOOST_CHECK_CLOSE_FRACTION( std::fabs( lvlhFrameThrust1( 0 ) ), thrustMagnitudeAlongVelocityVector,
                                        std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( std::fabs( lvlhFrameThrust2( 0 ) ), thrustMagnitudeAlongVelocityVector,
                                        std::numeric_limits< double >::epsilon( ) );

            // Test if other entries are 0 (to numerical precision).
            BOOST_CHECK_SMALL( lvlhFrameThrust1( 1 ),
                               thrustMagnitudeAlongVelocityVector * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( lvlhFrameThrust2( 1 ),
                               thrustMagnitudeAlongVelocityVector * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( lvlhFrameThrust1( 2 ),
                               thrustMagnitudeAlongVelocityVector * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( lvlhFrameThrust2( 2 ),
                               thrustMagnitudeAlongVelocityVector * std::numeric_limits< double >::epsilon( ) );
        }

        // Define semi-random thrust vector perpendicular to plane
        Eigen::Vector3d inertialThrustPerpendicularToPlane = thrustScale *  (
                    Eigen::Vector3d( relativeState.segment( 0, 3 ) ).cross(
                        Eigen::Vector3d( relativeState.segment( 3, 3 ) ) ) );
        {
            // Compute local thrust vectors
            lvlhFrameThrust1 =
                    nAxisAwayFromBodyMatrix.transpose( ) * inertialThrustPerpendicularToPlane;
            lvlhFrameThrust2 =
                    nAxisAwayTowardsBodyMatrix.transpose( ) * inertialThrustPerpendicularToPlane;

            // Compute thrust vector magnitudes
            double thrustMagnitudePerpendicularToPlane = inertialThrustPerpendicularToPlane.norm( );

            // Test if entry 2 of two vectors are the same but opposite
            BOOST_CHECK_CLOSE_FRACTION( lvlhFrameThrust1( 2 ), -lvlhFrameThrust2( 2 ),
                                        std::numeric_limits< double >::epsilon( ) );

            // Test thrust vector norms.
            BOOST_CHECK_CLOSE_FRACTION( thrustMagnitudePerpendicularToPlane, lvlhFrameThrust1.norm( ),
                                        5.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( thrustMagnitudePerpendicularToPlane, lvlhFrameThrust2.norm( ),
                                        5.0 * std::numeric_limits< double >::epsilon( ) );


            // Test if entry 2 magnitude is equal to norm
            BOOST_CHECK_CLOSE_FRACTION( std::fabs( lvlhFrameThrust1( 2 ) ), thrustMagnitudePerpendicularToPlane,
                                        std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( std::fabs( lvlhFrameThrust2( 2 ) ), thrustMagnitudePerpendicularToPlane,
                                        std::numeric_limits< double >::epsilon( ) );

            // Test if other entries are 0 (to numerical precision).
            BOOST_CHECK_SMALL( lvlhFrameThrust1( 0 ),
                               thrustMagnitudePerpendicularToPlane * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( lvlhFrameThrust2( 0 ),
                               thrustMagnitudePerpendicularToPlane * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( lvlhFrameThrust1( 1 ),
                               thrustMagnitudePerpendicularToPlane * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( lvlhFrameThrust2( 1 ),
                               thrustMagnitudePerpendicularToPlane * std::numeric_limits< double >::epsilon( ) );
        }

        // Define semi-random thrust vector along local y-axis
        Eigen::Vector3d inertialThrustAlongLocalYAxis = inertialThrustPerpendicularToPlane.cross(
                    inertialThrustVectorAlongVelocity );
        {
            // Compute local thrust vectors
            lvlhFrameThrust1 =
                    nAxisAwayFromBodyMatrix.transpose( ) * inertialThrustAlongLocalYAxis;
            lvlhFrameThrust2 =
                    nAxisAwayTowardsBodyMatrix.transpose( ) * inertialThrustAlongLocalYAxis;

            // Compute thrust vector magnitudes
            double thrustMagnitudeAlongYAxis = inertialThrustAlongLocalYAxis.norm( );

            // Test if entry 1 of two vectors are the same
            BOOST_CHECK_CLOSE_FRACTION( lvlhFrameThrust1( 1 ), -lvlhFrameThrust2( 1 ),
                                        std::numeric_limits< double >::epsilon( ) );

            // Test thrust vector norms.
            BOOST_CHECK_CLOSE_FRACTION( thrustMagnitudeAlongYAxis, lvlhFrameThrust1.norm( ),
                                        5.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( thrustMagnitudeAlongYAxis, lvlhFrameThrust2.norm( ),
                                        5.0 * std::numeric_limits< double >::epsilon( ) );

            // Test if entry 1 magnitude is equal to norm
            BOOST_CHECK_CLOSE_FRACTION( std::fabs( lvlhFrameThrust1( 1 ) ), thrustMagnitudeAlongYAxis,
                                        std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( std::fabs( lvlhFrameThrust2( 1 ) ), thrustMagnitudeAlongYAxis,
                                        std::numeric_limits< double >::epsilon( ) );

            // Test if other entries are 0 (to numerical precision).
            BOOST_CHECK_SMALL( lvlhFrameThrust1( 0 ), thrustMagnitudeAlongYAxis * std::numeric_limits< double >::epsilon( ));
            BOOST_CHECK_SMALL( lvlhFrameThrust2( 0 ), thrustMagnitudeAlongYAxis * std::numeric_limits< double >::epsilon( ));
            BOOST_CHECK_SMALL( lvlhFrameThrust1( 2 ), thrustMagnitudeAlongYAxis * std::numeric_limits< double >::epsilon( ));
            BOOST_CHECK_SMALL( lvlhFrameThrust2( 2 ), thrustMagnitudeAlongYAxis * std::numeric_limits< double >::epsilon( ));
        }

    }

    // Test 11: Test velocity based LVLH frame to inertial (I) frame transformation.
    {
        // Initialize initial thrust vector in velocity based LVLH frame.
        Eigen::Vector3d startThrustVector;
        startThrustVector( 0 ) = 2.0; // T direction
        startThrustVector( 1 ) = 4.0; // N direction
        startThrustVector( 2 ) = 1.0; // W direction
        bool doesNaxisPointAwayFromCentralBody = false;

        Eigen::Vector6d vehicleStateCartesian, centralBodyStateCartesian;
        vehicleStateCartesian << 0.0, -1.0, 0.0, -1.0, -0.5, 0.0;
        centralBodyStateCartesian << 0.0, 2.0, 0.0, -0.5, 0.0, 0.0;

        // Compute angle between positive T axis and negative y axis
        double angle = atan2( ( vehicleStateCartesian( 3 ) - centralBodyStateCartesian( 3 ) ),
                              ( vehicleStateCartesian( 4 ) - centralBodyStateCartesian( 4 ) ) );

        // Declare the expected thrust vector in the inertial (I) reference frame.
        Eigen::Vector3d expectedThrustVector;
        expectedThrustVector( 0 ) = startThrustVector( 0 ) * sin( angle ) + startThrustVector( 1 ) * cos( angle );
        expectedThrustVector( 1 ) = startThrustVector( 0 ) * cos( angle ) - startThrustVector( 1 ) * sin( angle );
        expectedThrustVector( 2 ) = -startThrustVector( 2 ); // z direction

        // Compute location of the point in the rotating frame subject to the transformation
        // matrix.
        Eigen::Vector3d transformedThrustVector;
        transformedThrustVector = reference_frames::getVelocityBasedLvlhToInertialRotation(
                    vehicleStateCartesian, centralBodyStateCartesian, doesNaxisPointAwayFromCentralBody ) * startThrustVector;

        // Check whether both vectors are equal within tolerances.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedThrustVector, transformedThrustVector, 1.0e-15 );
    }

    // Test 12: Test velocity based LVLH frame to planetocentric frame transformation with Keplerian elements.
    {
        // Initialize initial thrust vector in velocity based LVLH frame.
        Eigen::Vector3d startThrustVector;
        startThrustVector( 0 ) = 15.0; // T direction
        startThrustVector( 1 ) = -1.0; // N direction
        startThrustVector( 2 ) = 2.0; // W direction

        Eigen::Vector6d vehicleStateKeplerian;
        vehicleStateKeplerian << 1.0, 0.5, convertDegreesToRadians( 60.0 ),
                convertDegreesToRadians( 180.0 ), convertDegreesToRadians( 15.0 ), convertDegreesToRadians( 90.0 );

        // Declare the expected thrust vector in the inertial (I) reference frame. Reference data comes from Matlab file
        Eigen::Vector3d expectedThrustVector;
        expectedThrustVector( 0 ) = 13.959420290978478079; // x direction
        expectedThrustVector( 1 ) = -1.988147005863377672; // y direction
        expectedThrustVector( 2 ) = -5.5840716885526107944; // z direction

        // Compute location of the point in the rotating frame subject to the transformation
        // matrix.
        Eigen::Matrix3d rotationMatrixFromKeplerElements =
                reference_frames::getVelocityBasedLvlhToPlanetocentricRotationKeplerian(
                    vehicleStateKeplerian ).toRotationMatrix( );
        Eigen::Vector3d transformedThrustVector;
        transformedThrustVector = rotationMatrixFromKeplerElements * startThrustVector;

        // Check whether both vectors are equal within tolerances.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedThrustVector, transformedThrustVector,
                                           ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );

        Eigen::Vector6d vehicleStateCartesian =
                orbital_element_conversions::convertKeplerianToCartesianElements(
                    vehicleStateKeplerian, 3.986E14 );

        Eigen::Matrix3d rotationMatrixFromCartesianElements =
                reference_frames::getVelocityBasedLvlhToInertialRotation(
                                    vehicleStateCartesian, Eigen::Vector6d::Zero( ), false );

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL(
                            std::fabs( rotationMatrixFromCartesianElements( i, j ) -
                                       rotationMatrixFromKeplerElements( i, j ) ),
                            5.0 * std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}


BOOST_AUTO_TEST_CASE( testEulerAngleRetrieval )
{
    const double angleX = 2.1;
    const double angleZ = -0.3;
    const double angleY = 1.45;

    Eigen::Matrix3d rotationMatrix =
            ( Eigen::AngleAxisd( -angleX, Eigen::Vector3d::UnitX( ) ) *
              Eigen::AngleAxisd( -angleZ, Eigen::Vector3d::UnitZ( ) ) *
              Eigen::AngleAxisd( -angleY, Eigen::Vector3d::UnitY( ) ) ).toRotationMatrix( );

    Eigen::Vector3d eulerAngles = reference_frames::get132EulerAnglesFromRotationMatrix(
                rotationMatrix );

    BOOST_CHECK_CLOSE_FRACTION( angleX, eulerAngles( 0 ), 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( angleZ, eulerAngles( 1 ), 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( angleY, eulerAngles( 2 ), 1.0E-15 );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
