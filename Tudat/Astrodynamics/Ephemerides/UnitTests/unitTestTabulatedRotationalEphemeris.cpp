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
#include "Tudat/Astrodynamics/Ephemerides/tabulatedRotationalEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;

BOOST_AUTO_TEST_SUITE( test_tabulated_rotational_ephemeris )

//! Test to see if rotational ephemeris provides correct output, by comparing output from class functions to direct Spice calls.
BOOST_AUTO_TEST_CASE( testTabulatedRotationalEphemeris )
{
    // Define frames between which rotation will be computed
    const std::string baseFrame = "J2000";
    const std::string targetFrame = "IAU_Mars";

    // Load Spice kernel
    spice_interface::loadStandardSpiceKernels( );

    // Define time interval over which rotations are to be generated
    double startTime = 0.0;
    double finalTime = 2.0 * physical_constants::JULIAN_DAY;
    double timeStep = 10.0;

    // Create a list of concatenated quaternions and angular velocity vectors for rotation.
    std::map< double, Eigen::Matrix< double, 7, 1 > > rotationMap;
    double currentTime = startTime;
    while( currentTime < finalTime )
    {
        Eigen::Matrix< double, 7, 1 > currentState;
        Eigen::Quaterniond currentRotationToBaseFrame = spice_interface::computeRotationQuaternionBetweenFrames(
                    targetFrame, baseFrame, currentTime );
        currentState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat(
                    currentRotationToBaseFrame );
        currentState.segment( 4, 3 ) = currentRotationToBaseFrame.inverse( ) *
                spice_interface::getAngularVelocityVectorOfFrameInOriginalFrame(
                    baseFrame, targetFrame, currentTime );
        rotationMap[ currentTime ] = currentState;
        currentTime += timeStep;
    }

    // Create interpolator for rotation
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > rotationInterpolator =
            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( rotationMap, 6 );

    // Create tabulated rotational model
    std::shared_ptr< TabulatedRotationalEphemeris< double, double > > tabulatedEphemeris =
            std::make_shared< TabulatedRotationalEphemeris<  double, double > >( rotationInterpolator );

    //  Declare variables for rotation properties retrieved from class/spice
    Eigen::Matrix3d currentRotationMatrixToTargetFrame, currentRotationMatrixToBaseFrame;
    Eigen::Matrix3d expectedRotationMatrixToTargetFrame, expectedRotationMatrixToBaseFrame;

    Eigen::Matrix3d currentRotationMatrixDerivativeToTargetFrame, currentRotationMatrixDerivativeToBaseFrame;
    Eigen::Matrix3d expectedRotationMatrixDerivativeToTargetFrame, expectedRotationMatrixDerivativeToBaseFrame;

    Eigen::Vector3d bodyAngularVelocityVectorInBaseFrame, bodyAngularVelocityVectorInTargetFrame;
    Eigen::Vector3d expectedBodyAngularVelocityVectorInBaseFrame, expectedBodyAngularVelocityVectorInTargetFrame;

    // Set initial time slightly away from start of interpolator
    currentTime = startTime + 3600.0;

    // Set interpolation time step such that interpolation is not on the nodes
    double interpolationTimeStep = 100.0 * mathematical_constants::PI;

    // Compare rotational properties from interpolator/Spice
    while( currentTime < finalTime - 3600.0 )
    {

        // Compare rotation matrices
        currentRotationMatrixToTargetFrame = tabulatedEphemeris->getRotationToTargetFrame( currentTime );
        expectedRotationMatrixToTargetFrame = spice_interface::computeRotationQuaternionBetweenFrames(
                    baseFrame, targetFrame, currentTime );

        currentRotationMatrixToBaseFrame = tabulatedEphemeris->getRotationToBaseFrame( currentTime );
        expectedRotationMatrixToBaseFrame = spice_interface::computeRotationQuaternionBetweenFrames(
                    targetFrame, baseFrame, currentTime );

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL(
                            std::fabs( currentRotationMatrixToTargetFrame( i, j ) -
                                       expectedRotationMatrixToTargetFrame( i, j ) ), 1.0E-10 );
                BOOST_CHECK_SMALL(
                            std::fabs( currentRotationMatrixToBaseFrame( i, j ) -
                                       expectedRotationMatrixToBaseFrame( i, j ) ), 1.0E-10 );
            }
        }

        // Compare derivatives of rotation matrices
        currentRotationMatrixDerivativeToTargetFrame = tabulatedEphemeris->getDerivativeOfRotationToTargetFrame( currentTime );
        expectedRotationMatrixDerivativeToTargetFrame = spice_interface::computeRotationMatrixDerivativeBetweenFrames(
                    baseFrame, targetFrame, currentTime );

        currentRotationMatrixDerivativeToBaseFrame = tabulatedEphemeris->getDerivativeOfRotationToBaseFrame( currentTime );
        expectedRotationMatrixDerivativeToBaseFrame = spice_interface::computeRotationMatrixDerivativeBetweenFrames(
                    targetFrame, baseFrame, currentTime );

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL(
                            std::fabs( currentRotationMatrixDerivativeToTargetFrame( i, j ) -
                                       expectedRotationMatrixDerivativeToTargetFrame( i, j ) ), 1.0E-15 );
                BOOST_CHECK_SMALL(
                            std::fabs( currentRotationMatrixDerivativeToBaseFrame( i, j ) -
                                       expectedRotationMatrixDerivativeToBaseFrame( i, j ) ), 1.0E-15 );
            }
        }

        // Compare angular velocity vectors
        bodyAngularVelocityVectorInBaseFrame = tabulatedEphemeris->getRotationalVelocityVectorInBaseFrame( currentTime );
        bodyAngularVelocityVectorInTargetFrame = tabulatedEphemeris->getRotationalVelocityVectorInTargetFrame( currentTime );

        expectedBodyAngularVelocityVectorInBaseFrame = spice_interface::getAngularVelocityVectorOfFrameInOriginalFrame(
                    baseFrame, targetFrame, currentTime );
        expectedBodyAngularVelocityVectorInTargetFrame =
                expectedRotationMatrixToTargetFrame * expectedBodyAngularVelocityVectorInBaseFrame;

        for( unsigned int i = 0; i < 3; i++ )
        {

            BOOST_CHECK_SMALL(
                        std::fabs( bodyAngularVelocityVectorInBaseFrame( i ) -
                                   expectedBodyAngularVelocityVectorInBaseFrame( i ) ), 1.0E-10 );
            BOOST_CHECK_SMALL(
                        std::fabs( bodyAngularVelocityVectorInTargetFrame( i ) -
                                   expectedBodyAngularVelocityVectorInTargetFrame( i ) ), 1.0E-10 );
        }

        currentTime += interpolationTimeStep;
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
