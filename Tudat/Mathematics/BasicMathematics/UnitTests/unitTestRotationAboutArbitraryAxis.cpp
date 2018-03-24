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
 *    Hameduddin, I. Rotate vector(s) about axis, rodrigues_rot.m, available at
 *        http://www.mathworks.com/matlabcentral/fileexchange/34426-rotate-vectors-about-axis,
 *        2012, last accessed: 11th January, 2014.
 *    Murray, G. Rotation Matrices and Formulas java script, RotationMatrix.java available at
 *        https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas, 2011.
 *        last accessed: 20th January, 2014.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Mathematics/BasicMathematics/rotationAboutArbitraryAxis.h"

namespace tudat
{
namespace unit_tests
{

//! Test suite for rotations about about arbitrary axes.
BOOST_AUTO_TEST_SUITE( test_RotationAboutArbitraryAxis )

BOOST_AUTO_TEST_CASE( test_RotationAboutArbitraryAxis_PointRotationWithCommonOrigin )
{

    //Benchmark data is obtained using Matlab Script (Hameduddin, 2012).

    //Set expected rotated vector.
    const Eigen::Vector3d expectedRotatedPosition = Eigen::Vector3d( 0.0, 1.414213562373095, 1.0 );

    //Set origin of rotation.
    const Eigen::Vector3d originOfRotation = Eigen::Vector3d( 0.0, 0.0, 0.0 );

    //Set angle of rotation [rad].
    const double angleOfRotation = mathematical_constants::PI / 4.0;

    //Set axis of rotation.
    const Eigen::Vector3d axisOfRotation = Eigen::Vector3d( 0.0, 0.0, 1.0 );

    //Set initial position of point.
    const Eigen::Vector3d initialPositionOfPoint = Eigen::Vector3d( 1.0, 1.0, 1.0 );

    //Compute rotated position.
    const Eigen::Vector3d computedRotatedPosition = basic_mathematics::
        computeRotationOfPointAboutArbitraryAxis( originOfRotation, angleOfRotation,
                                                  axisOfRotation, initialPositionOfPoint );

    // Compare computed and expected vectors.
    BOOST_CHECK_SMALL( computedRotatedPosition.x( ), std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRotatedPosition.segment( 1, 2 ),
                                       expectedRotatedPosition.segment( 1, 2),
                                       std::numeric_limits< double >::epsilon( ) );

}

BOOST_AUTO_TEST_CASE( test_RotationAboutArbitraryAxis_PointRotationWithDifferentOrigins )
{

    //Benchmark data is obtained using java script (Murray, 2013).

    //Set expected rotated vector.
    const Eigen::Vector3d expectedRotatedPosition = Eigen::Vector3d( 3.156561876696307,
                                                                     -5.97145870839039,
                                                                     -4.418680390057799 );
    //Set origin of rotation.
    const Eigen::Vector3d originOfRotation = Eigen::Vector3d( 4.0, 1.0, -1.0 );

    //Set angle of rotation [rad]
    const double angleOfRotation = 7.0;

    //Set axis of rotation.
    const Eigen::Vector3d axisOfRotation = Eigen::Vector3d( 2.0, -2.0, 3.0 );

    //Set initial position of point.
    const Eigen::Vector3d initialPositionOfPoint = Eigen::Vector3d( -1.0, -5.0, -1.0 );

    //Compute rotated position.
    const Eigen::Vector3d computedRotatedPosition = basic_mathematics::
        computeRotationOfPointAboutArbitraryAxis( originOfRotation,  angleOfRotation,
                                                  axisOfRotation, initialPositionOfPoint );

    // Compare computed and expected radiation pressure acceleration vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRotatedPosition,
                                       expectedRotatedPosition,
                                       std::numeric_limits<double>::epsilon( ) );

}

BOOST_AUTO_TEST_CASE( test_RotationAboutArbitraryAxis_PointRotationWithDifferentOrigins2 )
{

    //Benchmark data is obtained using java script (Murray, 2013).

    //Set expected rotated vector.
    const Eigen::Vector3d expectedRotatedPosition = Eigen::Vector3d( 15.71267522236938,
                                                                     -0.9723168350942417,
                                                                      9.449538267504582 );
    //Set origin of rotation.
    const Eigen::Vector3d originOfRotation = Eigen::Vector3d( 1.5, 3.8, 12.0 );

    //Set angle of rotation [rad]
    const double angleOfRotation = 3.0;

    //Set axis of rotation.
    const Eigen::Vector3d axisOfRotation = Eigen::Vector3d( -4.6, 6.75, 7.7 );

    //Set initial position of point.
    const Eigen::Vector3d initialPositionOfPoint = Eigen::Vector3d( -4.3, -5.2, 1.2 );

    //Compute rotated position.
    const Eigen::Vector3d computedRotatedPosition = basic_mathematics::
        computeRotationOfPointAboutArbitraryAxis( originOfRotation, angleOfRotation,
                                                  axisOfRotation, initialPositionOfPoint );

    // Compare computed and expected radiation pressure acceleration vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRotatedPosition, expectedRotatedPosition, 1.0e-14 );

}

BOOST_AUTO_TEST_CASE( test_RotationAboutArbitraryAxis_VectorRotationWithCommonOrigin )
{

    //Benchmark data is obtained using Matlab Script (Hameduddin, 2012).

    //Set expected rotated vector.
    const Eigen::Vector3d expectedRotatedVector = Eigen::Vector3d( -5.480598288892924,
                                                                    0.532794739754852,
                                                                    6.138336269794405 );
    //Set origin of rotation.
    const Eigen::Vector3d originOfRotation = Eigen::Vector3d( 0.0, 0.0, 0.0 );

    //Set angle of rotation [rad].
    const double angleOfRotation = 12.0;

    //Set axis of rotation.
    const Eigen::Vector3d axisOfRotation = Eigen::Vector3d( -1.0, -2.0, -3.0 );

    //Set initial position of vector tail.
    const Eigen::Vector3d initialPositionOfVectorTail = Eigen::Vector3d( 1.0, 3.0, 5.0 );

    //Set initial vector.
    const Eigen::Vector3d initialVector = Eigen::Vector3d( -6.0, 4.0, 4.0 );

    //Compute rotated position.
    const Eigen::Vector3d computedRotatedVector = basic_mathematics::
        computeRotationOfVectorAboutArbitraryAxis( originOfRotation, angleOfRotation,
                                                   axisOfRotation, initialPositionOfVectorTail,
                                                   initialVector );

    // Compare computed and expected radiation pressure acceleration vectors
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRotatedVector, expectedRotatedVector, 1.0e-14 );

}

BOOST_AUTO_TEST_CASE( test_RotationAboutArbitraryAxis_VectorRotationWithDifferentOrigins )
{

    //Benchmark data is obtained using java script (Murray, 2013).

    //Set expected rotated vector.
    const Eigen::Vector3d expectedRotatedVector = Eigen::Vector3d( -7.15233485964669,
                                                                   -2.4444137866784175,
                                                                    8.308967883857736 );
    //Set origin of rotation.
    const Eigen::Vector3d originOfRotation = Eigen::Vector3d( 1.5, 3.8, 12.0 );

    //Set angle of rotation [rad].
    const double angleOfRotation = 3.0;

    //Set axis of rotation.
    const Eigen::Vector3d axisOfRotation = Eigen::Vector3d( -4.6, 6.75, 7.7 );

    //Set initial position of vector tail.
    const Eigen::Vector3d initialPositionOfVectorTail = Eigen::Vector3d( -4.3, -5.2, 1.2 );

    //Set initial vector.
    const Eigen::Vector3d initialVector = Eigen::Vector3d( 0.3, 11.2, 0.8 );

    //Compute rotated position.
    const Eigen::Vector3d computedRotatedVector = basic_mathematics::
        computeRotationOfVectorAboutArbitraryAxis( originOfRotation, angleOfRotation,
                                                   axisOfRotation, initialPositionOfVectorTail,
                                                   initialVector );

  // Compare computed and expected rotated vectors.
  TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRotatedVector, expectedRotatedVector, 1.0e-14 );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
