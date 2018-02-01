/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedRotationalEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;

BOOST_AUTO_TEST_SUITE( test_rotational_ephemeris )

// Test functions to calculate rotation matrix derivative from angular velocity vector and vice
// versa.
BOOST_AUTO_TEST_CASE( testRotationalEphemeris )
{
    // Define names of frames.
    const std::string baseFrame = "J2000";
    const std::string targetFrame = "IAU_VENUS";

    // Set time at which rotational ephemeris it to be called for subsequent tests.

    // Test rotation to target frame at specified time.
    {

        // The following code block can be used to retrieve the benchmark data from Spice.
        //        spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) +
        // "pck00010.tpc" );
        //        const double secondsSinceJ2000 = 1.0E6;
        //        Eigen::Quaterniond spiceRotationMatrixToFrame =
        //                spice_interface::computeRotationQuaternionBetweenFrames(
        // baseFrame, targetFrame, secondsSinceJ2000 );
        //        Eigen::Matrix3d spiceRotationMatrixDerivativeToFrame =
        //                spice_interface::computeRotationMatrixDerivativeBetweenFrames(
        // baseFrame, targetFrame, secondsSinceJ2000 );
        //        Eigen::Vector3d spiceRotationalVelocityVectorOfTargetFrameExpressedInBaseFrame =
        //                spice_interface::getAngularVelocityVectorOfFrameInOriginalFrame(
        // baseFrame, targetFrame, secondsSinceJ2000 );

        // Set rotational characteristics at given time, as calculated with Spice
        // (see above commented lines).
        Eigen::Matrix3d spiceRotationMatrixDerivativeToFrame;
        spiceRotationMatrixDerivativeToFrame <<
                                                1.690407961416589e-07, 2.288121921543265e-07, 9.283170431475241e-08,
                -2.468632444964533e-07, 1.540516111965609e-07, 6.981529179974795e-08,
                0.0,           0.0,          0.0;

        Eigen::Matrix3d spiceRotationMatrixToFrame;
        spiceRotationMatrixToFrame << -0.8249537745726603, 0.5148010526833556, 0.2333048348715243,
                -0.5648910720519699, -0.7646317780963481, -0.3102197940834743,
                0.01869081416890206, -0.3877088083617987, 0.9215923900425707;

        Eigen::Vector3d spiceRotationalVelocityVectorOfTargetFrameExpressedInBaseFrame;
        spiceRotationalVelocityVectorOfTargetFrameExpressedInBaseFrame << -5.593131603532092e-09,
                1.160198999048488e-07,
                -2.75781861386115e-07;

        // Calculate rotational velocity from SPICE rotation matrix derivative.
        Eigen::Vector3d manualRotationalVelocityVectorOfTargetFrameExpressedInBaseFrame =
                getRotationalVelocityVectorInBaseFrameFromMatrices(
                    spiceRotationMatrixToFrame, spiceRotationMatrixDerivativeToFrame.transpose( ) );

        // Calculate rotation matrix derivative from SPICE rotational velocity vector.
        Eigen::Matrix3d rotationMatrixDerivative = getDerivativeOfRotationMatrixToFrame(
                    spiceRotationMatrixToFrame, spiceRotationalVelocityVectorOfTargetFrameExpressedInBaseFrame );

        // Calculate rotational velocity from previously calculated rotation matrix derivative.
        Eigen::Matrix3d backCalculatedRotationMatrixDerivative =
                getDerivativeOfRotationMatrixToFrame(
                    spiceRotationMatrixToFrame,
                    manualRotationalVelocityVectorOfTargetFrameExpressedInBaseFrame );

        // Calculate rotation matrix derivative from previously calculated rotational velocity
        // vector.
        Eigen::Vector3d backCalculatedRotationalVelocityVector =
                getRotationalVelocityVectorInBaseFrameFromMatrices(
                    spiceRotationMatrixToFrame, rotationMatrixDerivative.transpose( ) );

        // Check equivalence of results.
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( manualRotationalVelocityVectorOfTargetFrameExpressedInBaseFrame( i ) -
                               backCalculatedRotationalVelocityVector( i ), 2.0E-22 );
            BOOST_CHECK_SMALL( manualRotationalVelocityVectorOfTargetFrameExpressedInBaseFrame( i ) -
                               spiceRotationalVelocityVectorOfTargetFrameExpressedInBaseFrame( i ), 2.0E-22 );

            for( int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( rotationMatrixDerivative( i, j ) -
                                   backCalculatedRotationMatrixDerivative( i, j ), 2.0E-22 );
                BOOST_CHECK_SMALL( rotationMatrixDerivative( i, j ) -
                                   spiceRotationMatrixDerivativeToFrame( i, j ), 2.0E-22 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
