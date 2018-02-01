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
 *      [1] J.R. Clynch, Introduction to Earth Gravity, March 2006,
 *          http://www.gmat.unsw.edu.au/snap/gps/clynch_pdfs/intgrav.pdf
 *      [2] WolframAlpha: computational knowledge engine, coriolis calculator
 *          http://www.wolframalpha.com/input/?i=coriolis+effect&a=*C.coriolis+effect-_*Formula-
 *
 */

#define BOOST_TEST_MAIN

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/ReferenceFrames/apparentAccelerationModel.h"

namespace tudat
{
namespace unit_tests
{

using reference_frames::computeApparentAcceleration;
using reference_frames::ApparentAccelerationModel;

BOOST_AUTO_TEST_SUITE( test_reference_frame_transformations )

// Summary of tests.
// Test 1: Check apparent acceleration of static object at Earth equator (=only centrifugal).
// Test 2: Check apparent acceleration of moving object at Earth pole (=only Coriolis).
// Test 3: Check apparent acceleration of moving object in an accelerating, rotational reference
//         frame.
// Test 4: Check apparent acceleration of moving object in inertial reference frame (=0).
// Test 5: Test class implementation.

//! Test apparent acceleration of static object at Earth equator (=only centrifugal).
// Reference data from [1] (page 5), recomputed using Excel for higher precision.
BOOST_AUTO_TEST_CASE( testCentrifugalAccelerationEarth )
{
    // Specify reference frame dynamics.
    const Eigen::Vector3d referenceFrameAcceleration = Eigen::Vector3d::Zero( );
    const Eigen::Vector3d rotationRate( 0.0, 0.0, 7.29212351699e-5 );
    const Eigen::Vector3d rotationAcceleration = Eigen::Vector3d::Zero( );

    // Specify state.
    const Eigen::Vector3d position( 6378000.0, 0.0, 0.0 );
    const Eigen::Vector3d velocity = Eigen::Vector3d::Zero( );

    // Compute apparent acceleration.
    const Eigen::Vector3d apparentAcceleration = computeApparentAcceleration(
            referenceFrameAcceleration, rotationRate, rotationAcceleration, position, velocity );
    const Eigen::Vector3d expectedAcceleration( 0.0339150567038567, 0.0, 0.0 );

    // Check that resulting acceleration matches expected one.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( apparentAcceleration, expectedAcceleration, 1.0e-12 );
}

//! Test apparent acceleration of moving object at Earth pole (=only Coriolis).
// Reference data from pen/paper analysis, checked using [2].
BOOST_AUTO_TEST_CASE( testCoriolisAccelerationAtEarthPoles )
{
    // Specify reference frame dynamics.
    const Eigen::Vector3d referenceFrameAcceleration = Eigen::Vector3d::Zero( );
    const Eigen::Vector3d rotationRate( 0.0, 0.0, 7.29212351699e-5 );
    const Eigen::Vector3d rotationAcceleration = Eigen::Vector3d::Zero( );

    // Specify state.
    const Eigen::Vector3d position( 0.0, 0.0, 6378000.0 );
    const Eigen::Vector3d velocity( 10.0, 0.0, 0.0 );

    // Compute apparent acceleration.
    const Eigen::Vector3d apparentAcceleration = computeApparentAcceleration(
            referenceFrameAcceleration, rotationRate, rotationAcceleration, position, velocity );
    const Eigen::Vector3d expectedAcceleration( 0.0, -0.0014584247033981, 0.0 );

    // Check that resulting acceleration matches expected one.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( apparentAcceleration, expectedAcceleration, 1.0e-13 );
}

//! Test apparent acceleration of moving object in a general non-inertial reference frame.
// Checked with computations performed in MATLAB.
BOOST_AUTO_TEST_CASE( testApparentAccelerationGeneral )
{
    // Specify reference frame dynamics.
    const Eigen::Vector3d referenceFrameAcceleration( 1.0, 2.0, 3.0 );
    const Eigen::Vector3d rotationRate( 0.1, 0.2, 0.3 );
    const Eigen::Vector3d rotationAcceleration( 0.03, 0.02, 0.01 );

    // Specify state.
    const Eigen::Vector3d position( 100.0, 300.0, 500.0 );
    const Eigen::Vector3d velocity( 10.0, 40.0, 25.0 );

    // Compute apparent acceleration.
    const Eigen::Vector3d apparentAcceleration = computeApparentAcceleration(
            referenceFrameAcceleration, rotationRate, rotationAcceleration, position, velocity );
    const Eigen::Vector3d expectedAcceleration( -2.0, 9.0, -10.0 );

    // Check that resulting acceleration matches expected one.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( apparentAcceleration, expectedAcceleration, 1.0e-14 );
}

//! Test apparent acceleration of moving object in an inertial reference frame (=0).
BOOST_AUTO_TEST_CASE( testApparentAccelerationStatic )
{
    // Specify reference frame dynamics.
    const Eigen::Vector3d referenceFrameAcceleration = Eigen::Vector3d::Zero( );
    const Eigen::Vector3d rotationRate = Eigen::Vector3d::Zero( );
    const Eigen::Vector3d rotationAcceleration = Eigen::Vector3d::Zero( );

    // Specify state.
    const Eigen::Vector3d position( 100.0, 300.0, 500.0 );
    const Eigen::Vector3d velocity( 10.0, 40.0, 25.0 );

    // Compute apparent acceleration.
    const Eigen::Vector3d apparentAcceleration = computeApparentAcceleration(
            referenceFrameAcceleration, rotationRate, rotationAcceleration, position, velocity );
    const Eigen::Vector3d expectedAcceleration( 0.0, 0.0, 0.0 );

    // Check that resulting acceleration matches expected one.
    TUDAT_CHECK_MATRIX_BASE( apparentAcceleration, expectedAcceleration );
}

//! Test class implementation of apparent acceleration model.
BOOST_AUTO_TEST_CASE( testApparentAccelerationClass )
{
    // Specify reference frame dynamics.
    const Eigen::Vector3d referenceFrameAcceleration( 1.0, 2.0, 3.0 );
    const Eigen::Vector3d rotationRate( 0.1, 0.2, 0.3 );
    const Eigen::Vector3d rotationAcceleration( 0.03, 0.02, 0.01 );
    const Eigen::Vector3d position( 100.0, 300.0, 500.0 );
    const Eigen::Vector3d velocity( 10.0, 40.0, 25.0 );

    // Create functions returning the parameters defined above.
    boost::function< Eigen::Vector3d( ) > accelerationFunction =
            boost::lambda::constant( referenceFrameAcceleration );
    boost::function< Eigen::Vector3d( ) > rotationRateFunction =
            boost::lambda::constant( rotationRate );
    boost::function< Eigen::Vector3d( ) > rotationAccelerationFunction =
            boost::lambda::constant( rotationAcceleration );
    boost::function< Eigen::Vector3d( ) > positionFunction =
            boost::lambda::constant( position );
    boost::function< Eigen::Vector3d( ) > velocityFunction =
            boost::lambda::constant( velocity );

    // Create object of apparentAcceleraitonModel class; pass functions in constructor
    ApparentAccelerationModel apparentAccelerationModel(
           accelerationFunction, rotationRateFunction, rotationAccelerationFunction,
           positionFunction, velocityFunction );
    const Eigen::Vector3d expectedAcceleration( -2.0, 9.0, -10.0 );

    // Compute apparent acceleration.
    const Eigen::Vector3d computedAcceleration = apparentAccelerationModel.getAcceleration( );

    // Check that resulting acceleration matches expected one.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedAcceleration,
                                       expectedAcceleration, 1.0e-14 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
