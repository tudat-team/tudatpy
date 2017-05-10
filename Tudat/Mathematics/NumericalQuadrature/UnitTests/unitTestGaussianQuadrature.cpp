/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/NumericalQuadrature/gaussianQuadrature.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{


// Functions to be tested

double sinFunction( const double x ) {
    return std::sin( x );
}

double expFunction( const double x ) {
    return std::exp( x );
}

double polyFunction( const double x ) {
    return std::pow( x, 9 ) + 2 * std::pow( x, 7 ) - std::pow( x, 4 ) + 8 * std::pow( x, 2 ) - 11;
}


BOOST_AUTO_TEST_SUITE( test_gaussian_quadrature )

//! Test if quadrature is computed correctly (sine function with 1E4 data points).
BOOST_AUTO_TEST_CASE( testIntegralSineFunction )
{
    using namespace mathematical_constants;
    using namespace numerical_quadrature;

    const unsigned int steps = 8;
    GaussianQuadrature< double, double > integrator( sinFunction, 0, PI, steps );

    const double expectedIntegral = 2.0;
    const double computedIntegral = integrator.getQuadrature();

    // Check if computed sample mean matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedIntegral, expectedIntegral, 1E-10 );
}


//! Test if quadrature is computed correctly (exponential function).
BOOST_AUTO_TEST_CASE( testIntegralExpFunction )
{
    using namespace mathematical_constants;
    using namespace numerical_quadrature;

    const unsigned int steps = 7;
    GaussianQuadrature< double, double > integrator( expFunction, -2, 2, steps );

    const double expectedIntegral = 7.25372081569404;
    const double computedIntegral = integrator.getQuadrature();

    // Check if computed sample mean matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedIntegral, expectedIntegral, 1E-10 );
}


//! Test if quadrature is computed correctly (polynomial function).
BOOST_AUTO_TEST_CASE( testIntegralPolyFunction )
{
    using namespace mathematical_constants;
    using namespace numerical_quadrature;

    // A polynomial of order 2*steps - 1 is computed exactly. Order of tested polynomial is 9, thus steps = 5.
    const unsigned int steps = 5;
    GaussianQuadrature< double, double > integrator( polyFunction, -2, 4, steps );

    const double expectedIntegral = 120990;
    const double computedIntegral = integrator.getQuadrature();

    // Check if computed sample mean matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedIntegral, expectedIntegral, 1E-12 );
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
