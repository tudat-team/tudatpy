/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/gravitation/ringGravityField.h"
#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{
namespace unit_tests
{

//! Test the functionality of the polyhedron gravity field class.
BOOST_AUTO_TEST_SUITE( test_ring_gravity_field )

//! Test getters
BOOST_AUTO_TEST_CASE( testSettingAndGettingParameters )
{
    const double tolerance = 1.0E-14;

    // Define ring gravity parameters
    const double gravitationalParameter = 2.39e21 * physical_constants::GRAVITATIONAL_CONSTANT;
    const double ringRadius = 2.7 * physical_constants::ASTRONOMICAL_UNIT;
    const bool ellipticIntegralSFromDAndB = false;

    gravitation::RingGravityField gravityField = gravitation::RingGravityField(
            gravitationalParameter, ringRadius, ellipticIntegralSFromDAndB);

    // Gravitational parameter
    double retrievedGravitationalParameter = gravityField.getGravitationalParameter();
    BOOST_CHECK_CLOSE_FRACTION( gravitationalParameter, retrievedGravitationalParameter, tolerance );

    // Ring radius
    double retrievedRadius = gravityField.getRingRadius();
    BOOST_CHECK_CLOSE_FRACTION( ringRadius, retrievedRadius, tolerance );

    // Elliptic integrals flag
    bool retrievedEllipticIntegralSFromDAndB = gravityField.getEllipticIntegralSFromDAndB();
    BOOST_CHECK_EQUAL( ellipticIntegralSFromDAndB, retrievedEllipticIntegralSFromDAndB );

}

//! Test computation of potential, gradient of potential, hessian of potential
BOOST_AUTO_TEST_CASE( testGravityComputation )
{
    const double tolerance = 1.0E-14;

    // Define ring gravity parameters
    const double gravitationalParameter = 2.39e21 * physical_constants::GRAVITATIONAL_CONSTANT;
    const double ringRadius = 2.7 * physical_constants::ASTRONOMICAL_UNIT;

    for ( bool ellipticIntegralSFromDAndB : { true } )
    {

        gravitation::RingGravityField gravityField = gravitation::RingGravityField(
            gravitationalParameter, ringRadius, ellipticIntegralSFromDAndB);

        Eigen::Vector3d bodyFixedPosition;
        double expectedPotential;
        Eigen::Vector3d expectedGradient;

        (bodyFixedPosition << 0.0, 0.0, 0.0).finished();
        expectedPotential = gravitationalParameter / ringRadius;
        (expectedGradient << 0.0, 0.0, 0.0).finished();

        BOOST_CHECK_CLOSE_FRACTION(
                expectedPotential, gravityField.getGravitationalPotential( bodyFixedPosition ), tolerance );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                expectedGradient, gravityField.getGradientOfPotential( bodyFixedPosition ), tolerance );

//        for ( unsigned int positionId : { 0 } )
//        {
//            // Select expected values
//            if ( positionId == 0 )
//            {
//                (bodyFixedPosition << 0.0, 0.0, 0.0).finished();
//                expectedPotential = 3.19403761604211e-5;
//                (expectedGradient << 2.31329148957265e-6, 1.91973919943187e-6, 1.91973919943187e-6).finished();
//            }
//        }
    }

    // Compute potential, gradient, laplacian and hessian, and compare with results from D'Urso (2014)

//    Eigen::Vector3d bodyFixedPosition;
//    double expectedPotential, computedPotential;
//    Eigen::Vector3d expectedGradient, computedGradient;
//    double expectedLaplacian, computedLaplacian;
//    Eigen::Matrix3d expectedHessian, computedHessian;
//
//    for (unsigned int positionId : {0,1,2,3,4,5})
//    {
//        bool testPotential = true;
//        bool testGradient = true;
//        bool testLaplacian = true;
//        bool testHessian = true;
//
//        // Select expected values
//        if ( positionId == 0 )
//        {
//            (bodyFixedPosition << 0.0, 0.0, 0.0).finished();
//            expectedPotential = 3.19403761604211e-5;
//            (expectedGradient << 2.31329148957265e-6, 1.91973919943187e-6, 1.91973919943187e-6).finished();
//            expectedLaplacian = - 0.5 * mathematical_constants::PI * gravitationalConstant * density;
//            testHessian = false;
//        }
//        else if ( positionId == 1 )
//        {
//            (bodyFixedPosition << 5.0, 0.0, 0.0).finished();
//            expectedPotential = 3.99993558939122e-5;
//            (expectedGradient << 9.90115534890074e-7, 3.24128042248715e-6, 3.24128042248715e-6).finished();
//            expectedLaplacian = - 1.0 * mathematical_constants::PI * gravitationalConstant * density;
//            testHessian = false;
//        }
//        else if ( positionId == 2 )
//        {
//            (bodyFixedPosition << 0.0, 3.0, 2.0).finished();
//            expectedPotential = 4.03528375471853e-5;
//            (expectedGradient << 4.73368592565013e-6, 9.68164362892554e-7, 1.59674500375495e-6).finished();
//            expectedLaplacian = - 2.0 * mathematical_constants::PI * gravitationalConstant * density;
//            (expectedHessian <<
//                -4.02204713784183E-08, 1.87140408935899E-07, 3.51261972418670E-07,
//                1.87140408935899E-07, -5.01781942367494E-07, 8.58712984897779E-08,
//                3.51261972418670E-07, 8.58712984897779E-08, -5.77398275537941E-07).finished();
//        }
//        else if ( positionId == 3 )
//        {
//            (bodyFixedPosition << -5.0, 5.0, 5.0).finished();
//            expectedLaplacian = 0.0;
//            testPotential = false;
//            testGradient = false;
//            testHessian = false;
//        }
//        else if ( positionId == 4 )
//        {
//            (bodyFixedPosition << 10.0, 5.0, 5.0).finished();
//            expectedLaplacian = - 4.0 * mathematical_constants::PI * gravitationalConstant * density;
//            testPotential = false;
//            testGradient = false;
//            testHessian = false;
//        }
//        else
//        {
//            (bodyFixedPosition << 10.0, 5.0, 10.0 + 1e-10).finished();
//            expectedLaplacian = 0.0;
//            testPotential = false;
//            testGradient = false;
//            testHessian = false;
//        }
//
//        // Compute values and compare with expected ones
//        if ( testPotential )
//        {
//            computedPotential = gravityField.getGravitationalPotential( bodyFixedPosition );
//            BOOST_CHECK_CLOSE_FRACTION( expectedPotential, computedPotential, tolerance );
//        }
//        if ( testGradient )
//        {
//            computedGradient = gravityField.getGradientOfPotential( bodyFixedPosition );
//            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedGradient, computedGradient, tolerance );
//        }
//        if ( testLaplacian )
//        {
//            computedLaplacian = gravityField.getLaplacianOfPotential( bodyFixedPosition );
//            // Expected value is 0 for point 3, so add a constant to the laplacian
//            if ( positionId == 3 || positionId == 5 )
//            {
//                computedLaplacian += 0.1;
//                expectedLaplacian += 0.1;
//            }
//            BOOST_CHECK_CLOSE_FRACTION( expectedLaplacian, computedLaplacian, tolerance );
//
//        }
//        if ( testHessian )
//        {
//            computedHessian = gravityField.getHessianOfPotential( bodyFixedPosition );
//            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedHessian, computedHessian, tolerance );
//        }
//
//    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace tudat
} // namespace unit_tests
