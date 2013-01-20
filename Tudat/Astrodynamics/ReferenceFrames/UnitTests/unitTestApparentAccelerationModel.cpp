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
 *      120716    A. Ronse          Creation of code.
 *
 *    References
 *      [1] J.R. Clynch, Introduction to Earth Gravity, March 2006,
 *          http://www.gmat.unsw.edu.au/snap/gps/clynch_pdfs/intgrav.pdf
 *      [2] WolframAlpha: computational knowledge engine, coriolis calculator
 *          http://www.wolframalpha.com/input/?i=coriolis+effect&a=*C.coriolis+effect-_*Formula-
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/ReferenceFrames/apparentAccelerationModel.h"

namespace tudat
{
namespace unit_tests
{

using tudat::reference_frames::computeApparentAcceleration;
using tudat::reference_frames::ApparentAccelerationModel;

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
