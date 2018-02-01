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

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
//#include "Tudat/External/SpiceInterface/spiceInterface.h"
//#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{
namespace unit_tests
{

using unit_conversions::convertDegreesToRadians;
using ephemerides::SimpleRotationalEphemeris;
using basic_astrodynamics::JULIAN_DAY_ON_J2000;
using basic_astrodynamics::JULIAN_DAY_AT_0_MJD;
using reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion;

BOOST_AUTO_TEST_SUITE( test_simple_rotational_ephemeris )

// Test simple rotational ephemeris class by compariing results to Venus-fixed frame with Spice.
BOOST_AUTO_TEST_CASE( testSimpleRotationalEphemeris )
{    
    // Data from pck00010.tpc spice kernel.
    const double venusPoleRightAscension = convertDegreesToRadians( 272.76 );
    const double venusPoleDeclination = convertDegreesToRadians( 67.16 );
    const double venusPrimeMeridianAtJ2000 = convertDegreesToRadians( 160.20 );
    const double venusRotationRate = convertDegreesToRadians( -1.4813688 ) /
            physical_constants::JULIAN_DAY;
    
    // Define names of frames.
    const std::string baseFrame = "J2000";
    const std::string targetFrame = "IAU_VENUS";

    // Calculate initial rotation quaternion to frame.
    const Eigen::Quaterniond initialRotationToTargetFrame =
            getInertialToPlanetocentricFrameTransformationQuaternion(
                venusPoleDeclination, venusPoleRightAscension, venusPrimeMeridianAtJ2000 );

    // Create rotational ephemeris objects from (one for each type of constructor)
    SimpleRotationalEphemeris venusRotationalEphemerisFromAngles(
                venusPoleRightAscension, venusPoleDeclination, venusPrimeMeridianAtJ2000,
                venusRotationRate, 0.0, baseFrame, targetFrame );
    tudat::ephemerides::SimpleRotationalEphemeris venusRotationalEphemerisFromInitialState(
                initialRotationToTargetFrame,
                venusRotationRate, 0.0, baseFrame, targetFrame );


    // Test initial rotation matrix, both from the direct constructed object, and the one from
    // constituent angles.
    {
        // The following code block can be used to retrieve the benchmark data from Spice.
        /*
        loadSpiceKernelInTudat( getTudatRootPath( ) +
              "Astrodynamics/Ephemerides/UnitTests/pck00010.tpc" );
        Eigen::Matrix3d spiceInitialRotationToTargetFrameMatrix;
        Eigen::Quaterniond spiceInitialRotationToTargetFrame =
                computeRotationQuaternionBetweenFrames( baseFrame, targetFrame, 0.0 );
       */

        // Set initial rotation, as calculated with Spice (see above commented lines)
        Eigen::Matrix3d spiceInitialRotationToTargetFrameMatrix;
        spiceInitialRotationToTargetFrameMatrix
                << -0.9548214974296336, 0.2665104385944917, 0.1314841974018291,
                -0.296591573568662, -0.882413772579987, -0.3652114078848295,
                0.01869081416890202, -0.3877088083617989, 0.9215923900425705;
        
        // Check calculated with spice initial state.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( spiceInitialRotationToTargetFrameMatrix,
                                           Eigen::Matrix3d( initialRotationToTargetFrame ),
                                           1.0E-15 );
        
        // Check ephemeris calculations with Spice initial state.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    Eigen::Matrix3d(
                        venusRotationalEphemerisFromAngles.getRotationToTargetFrame( 0.0 ) ),
                    Eigen::Matrix3d( initialRotationToTargetFrame ),
                    1.0E-15 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    Eigen::Matrix3d(
                        venusRotationalEphemerisFromInitialState.getRotationToTargetFrame( 0.0 ) ),
                    Eigen::Matrix3d( initialRotationToTargetFrame ),
                    1.0E-15 );
    }

    // Set time at which rotational ephemeris it to be called for subsequent tests.
    const double secondsSinceJ2000 = 1.0E6;

    // Test rotation to target frame at specified time.
    {
        /*
        spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00010.tpc" );
        // The following code block can be used to retrieve the benchmark data from Spice.

        Eigen::Quaterniond spiceRotationMatrix =
               computeRotationQuaternionBetweenFrames( baseFrame, targetFrame, secondsSinceJ2000 );\
        spice_interface::computeRotationMatrixDerivativeBetweenFrames( baseFrame, targetFrame, secondsSinceJ2000 );

        */
        Eigen::Matrix3d spiceRotationMatrixDerivative;
        spiceRotationMatrixDerivative << 1.690407961416589e-07, 2.288121921543265e-07, 9.283170431475241e-08,
                -2.468632444964533e-07, 1.540516111965609e-07, 6.981529179974795e-08,
                0.0,           0.0,          0.0;


        // Set rotation at given time, as calculated with Spice (see above commented lines)
        Eigen::Matrix3d spiceRotationMatrix;
        spiceRotationMatrix << -0.8249537745726603, 0.5148010526833556, 0.2333048348715243,
                -0.5648910720519699, -0.7646317780963481, -0.3102197940834743,
                0.01869081416890206, -0.3877088083617987, 0.9215923900425707;

        // Calculate rotation to frame, and its time derivative, at certain time;
        Eigen::Quaterniond ephemerisRotation =
                venusRotationalEphemerisFromAngles.getRotationToTargetFrame(
                    secondsSinceJ2000 );
        Eigen::Matrix3d ephemerisRotationDerivative =
                venusRotationalEphemerisFromAngles.getDerivativeOfRotationToTargetFrame(
                    secondsSinceJ2000 );;

        // Check Spice result with ephemerides results.
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL(
                            Eigen::Matrix3d( ephemerisRotation )( i, j ) -
                            spiceRotationMatrix( i, j ), 2.0E-15 );
                BOOST_CHECK_SMALL(
                            spiceRotationMatrixDerivative( i, j ) -
                            ephemerisRotationDerivative( i, j ), 2.0E-22 );
            }

        }

        // Compare inverse rotation and its time derivative.
        Eigen::Quaterniond inverseEphemerisRotation =
                venusRotationalEphemerisFromInitialState.getRotationToBaseFrame(
                    secondsSinceJ2000 );
        Eigen::Matrix3d inverseEphemerisRotationDerivative =
                venusRotationalEphemerisFromAngles.getDerivativeOfRotationToBaseFrame(
                    secondsSinceJ2000 );

        Eigen::Matrix3d inverseSpiceRotationMatrix = spiceRotationMatrix.inverse( );
        Eigen::Matrix3d inverseSpiceRotationMatrixDerivative = spiceRotationMatrixDerivative.transpose( );

        // Check Spice result with ephemerides results.
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL(
                            Eigen::Matrix3d( inverseEphemerisRotation )( i, j ) -
                            inverseSpiceRotationMatrix( i, j ), 2.0E-15 );
                BOOST_CHECK_SMALL(
                            inverseEphemerisRotationDerivative( i, j ) -
                            inverseSpiceRotationMatrixDerivative( i, j ), 2.0E-22 );
            }

        }
    }

    // Test rotation matrix derivative by means of finite differences
    {
        double timeStep = 10.0;

        // Calculate rotations to frame at two times.
        Eigen::Quaterniond ephemerisRotation1 =
                venusRotationalEphemerisFromAngles.getRotationToTargetFrame(
                    secondsSinceJ2000 );
        Eigen::Quaterniond ephemerisRotation2 =
                venusRotationalEphemerisFromAngles.getRotationToTargetFrame(
                    secondsSinceJ2000 + timeStep );

        // Numerically calculate matrix derivative.
        Eigen::Matrix3d numericalEphemerisRotationDerivative =
                ( Eigen::Matrix3d( ephemerisRotation2 ) - Eigen::Matrix3d( ephemerisRotation1 ) ) /
                timeStep;

        // Calculate matrix derivative directly and compare.
        Eigen::Matrix3d ephemerisRotationDerivative =
                venusRotationalEphemerisFromAngles.getDerivativeOfRotationToTargetFrame(
                    secondsSinceJ2000 + timeStep / 2.0 );
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL(
                            numericalEphemerisRotationDerivative( i, j ) -
                            ephemerisRotationDerivative( i, j ), 2.0E-16 );
            }

        }

        // Calculate rotations from frame at two times;
        ephemerisRotation1 =
                venusRotationalEphemerisFromAngles.getRotationToBaseFrame(
                    secondsSinceJ2000 );
        ephemerisRotation2 =
                venusRotationalEphemerisFromAngles.getRotationToBaseFrame(
                    secondsSinceJ2000 + timeStep );

        // Numerically calculate matrix derivative.
        ephemerisRotationDerivative =
                venusRotationalEphemerisFromAngles.getDerivativeOfRotationToBaseFrame(
                    secondsSinceJ2000 + timeStep / 2.0 );

        // Calculate matrix derivative directly and compare.
        numericalEphemerisRotationDerivative =
                ( Eigen::Matrix3d( ephemerisRotation2 ) - Eigen::Matrix3d( ephemerisRotation1 ) ) /
                timeStep;
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL(
                            numericalEphemerisRotationDerivative( i, j ) -
                            ephemerisRotationDerivative( i, j ), 2.0E-16 );
            }

        }
    }

    // Test computation of rotational velocity vector
    {
        Eigen::Vector3d spiceRotationalVelocityVectorInBaseFrame;
        spiceRotationalVelocityVectorInBaseFrame << -5.593131603532092e-09, 1.160198999048488e-07, -2.75781861386115e-07;

        // Set rotation at given time, as calculated with Spice (see above commented lines)
        Eigen::Matrix3d spiceRotationMatrixToTargetFrame;
        spiceRotationMatrixToTargetFrame << -0.8249537745726603, 0.5148010526833556, 0.2333048348715243,
                -0.5648910720519699, -0.7646317780963481, -0.3102197940834743,
                0.01869081416890206, -0.3877088083617987, 0.9215923900425707;

        Eigen::Vector3d spiceRotationalVelocityVectorInTargetFrame =
                spiceRotationMatrixToTargetFrame * spiceRotationalVelocityVectorInBaseFrame;

        Eigen::Vector3d rotationalVelocityVectorInBaseFrame =
                venusRotationalEphemerisFromAngles.getRotationalVelocityVectorInBaseFrame( secondsSinceJ2000 );
        Eigen::Vector3d rotationalVelocityVectorInBaseFrame2 =
                venusRotationalEphemerisFromInitialState.getRotationalVelocityVectorInBaseFrame( secondsSinceJ2000 );

        Eigen::Vector3d rotationalVelocityVectorInTargetFrame =
                venusRotationalEphemerisFromAngles.getRotationalVelocityVectorInTargetFrame( secondsSinceJ2000 );
        Eigen::Vector3d rotationalVelocityVectorInTargetFrame2 =
                venusRotationalEphemerisFromInitialState.getRotationalVelocityVectorInTargetFrame( secondsSinceJ2000 );

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( rotationalVelocityVectorInBaseFrame( i ) -
                                   spiceRotationalVelocityVectorInBaseFrame( i ) ), 1.0E-21 );
            BOOST_CHECK_SMALL(
                        std::fabs( rotationalVelocityVectorInBaseFrame2( i ) -
                                   spiceRotationalVelocityVectorInBaseFrame( i ) ), 1.0E-21 );

            BOOST_CHECK_SMALL(
                        std::fabs( rotationalVelocityVectorInTargetFrame( i ) -
                                   spiceRotationalVelocityVectorInTargetFrame( i ) ), 1.0E-21 );
            BOOST_CHECK_SMALL(
                        std::fabs( rotationalVelocityVectorInTargetFrame2( i ) -
                                   spiceRotationalVelocityVectorInTargetFrame( i ) ), 1.0E-21 );
        }
    }

    // Test rotation from target frame at specified time (is checked by checking if it is inverse
    // of opposite rotation).
    {
        // Test orthonormality of matrices from object created from initial angles.
        Eigen::Matrix3d productOfOppositeRotations =
                Eigen::Matrix3d( venusRotationalEphemerisFromAngles.getRotationToTargetFrame(
                                     secondsSinceJ2000 ) ) *
                Eigen::Matrix3d( venusRotationalEphemerisFromAngles.getRotationToBaseFrame(
                                     secondsSinceJ2000 ) );

        for ( int i = 0; i < 3; i++ )
        {
            for ( int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( productOfOppositeRotations( i, j ) -
                                   Eigen::Matrix3d::Identity( )( i, j ), 1.0E-15 );
            }
        }

        // Test orthonormality of matrices from object created from initial state.
        productOfOppositeRotations =
                Eigen::Matrix3d( venusRotationalEphemerisFromInitialState.getRotationToTargetFrame(
                                     secondsSinceJ2000 ) ) *
                Eigen::Matrix3d( venusRotationalEphemerisFromInitialState.getRotationToBaseFrame(
                                     secondsSinceJ2000 ) );

        for( int i = 0; i < 3; i++ )
        {
            for( int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( productOfOppositeRotations( i, j ) -
                                   Eigen::Matrix3d::Identity( )( i, j ), 1.0E-15 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
