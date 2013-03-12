/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110519    F.M. Engelen      File created.
 *      110628    K. Kumar          Minor comment and layout modifications.
 *      110701    K. Kumar          Updated unit tests to check relative error;
 *                                  Updated file path.
 *      110718    F.M. Engelen      Took out falacy in test (only checking the norm) and
 *                                  added ItoE and EtoI transformation.
 *      110726    K. Kumar          Minor modifications; updated relative error
 *                                  wrt to norm.
 *      110808    F.M. Engelen      Updated with better tests, changed test for vertical frame.
 *      120529    E.A.G. Heeren     Boostified unit tests.
 *      120614    P. Musegaas       Removed unneccessary using statements and normalizations.
 *      130312    D. Dirkx          Added unit test for planet-fixed <-> inertial without equal
 *                                  equatorial frame.
 *
 *    References
 *
 *    Notes
 *      The reference frame definitions/abbreviations can be found in the file
 *      referenceFrameTransformations.h.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_reference_frame_transformations )

// Summary of tests.
// Test 1: Test inertial to rotating planetocentric frame transformation.
// Test 2: Check whether the transformed matrix of Test 1 is also correct.
// Test 3: Same test as Test 1 for the transformation quaternion.
// Test 4: Same test as Test 2 for the transformation quaternion.
// Test 5: Test airspeed-based aerodynamic to body frame transformation.
// Test 6: Test planetocentric to local vertical frame transformation quaternion.
// Test 7: Check whether the transformed matrix of Test 6 is also correct.

// Test inertial to rotating planetocentric frame transformations.
BOOST_AUTO_TEST_CASE( testRotatingPlanetocentricFrameTransformations )
{
    // Using declarations.
    using std::atan2;
    using std::cos;
    using std::sin;
    using std::pow;
    using std::sqrt;
    using tudat::unit_conversions::convertDegreesToRadians;

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

        // Compute location of the point in the rotating frame subject to the transformation matrix.
        Eigen::Vector3d transformedLocation;
        transformedLocation = reference_frames::
                getInertialToPlanetocentricFrameTransformationMatrix( angleInTime ) * startLocation;

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
        transformedLocation = tudat::reference_frames::
                getRotatingPlanetocentricToInertialFrameTransformationMatrix( angleInTime ) *
                tudat::reference_frames::
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
    using tudat::unit_conversions::convertDegreesToRadians;

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
    using tudat::unit_conversions::convertDegreesToRadians;

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
                getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix( angleOfAttack,
                                                                            angleOfSideslip )
                * startLocation;

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
    using tudat::unit_conversions::convertDegreesToRadians;

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
                getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion( longitude,
                                                                                       latitude )
                * startLocation;

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
    using tudat::basic_astrodynamics::unit_conversions::convertDegreesToRadians;
    using namespace tudat::reference_frames;

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

BOOST_AUTO_TEST_SUITE_END( )

} //namespace unit_tests
} //namespace tudat

