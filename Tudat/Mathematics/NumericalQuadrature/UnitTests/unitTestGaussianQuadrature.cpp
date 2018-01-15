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

#include "Tudat/Mathematics/NumericalQuadrature/gaussianQuadrature.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{


// Functions to be tested

double sinFunction( const double x )
{
    return std::sin( x );
}

double expFunction( const double x )
{
    return std::exp( x );
}

double polyFunction( const double x )
{
    return std::pow( x, 9 ) + 2 * std::pow( x, 7 ) - std::pow( x, 4 ) + 8 * std::pow( x, 2 ) - 11;
}


//! Even-order derivatives for error assessment
double minSinFunction( const double x )
{
    return -std::sin( x );
}

boost::function< double( const double ) > dSinFunction( unsigned int n )
{
    if ( n % 2 == 0 )
    {
        if ( n % 4 == 0 )
        {
            return sinFunction;
        }
        else
        {
            return minSinFunction;
        }
    }
    else
    {
        throw std::runtime_error( "Unknown derivative" );
    }
}

boost::function< double( const double ) > dExpFunction( unsigned int n )
{
    return expFunction;
}

double d4PolyFunction( const double x )
{
    return 24 * ( 126 * std::pow( x, 5 ) + 70 * std::pow( x, 3 ) - 1 );
}

double d6PolyFunction( const double x )
{
    return 10080 * ( 6 * std::pow( x, 3 ) + x );
}

double d8PolyFunction( const double x )
{
    return 362800 * x;
}

boost::function< double( const double ) > dPolyFunction( unsigned int n )
{
    switch( n )
    {
    case 4: return d4PolyFunction;
    case 6: return d6PolyFunction;
    case 8: return d8PolyFunction;
    default: throw std::runtime_error( "Unknown derivative" );
    }
}


// Factorial
unsigned int factorial( const unsigned int x )
{
    if ( x > 0 )
    {
        return x * factorial( x - 1 );
    }
    else
    {
        return 1;
    }
}

// Error function for Gaussian quadrature [ from https://en.wikipedia.org/wiki/Gaussian_quadrature#Error_estimates ]
double gaussianQuadratureError( const unsigned int n, boost::function< double( const double ) > derivative2nth,
                                const double abscissa, const double lowerLimit, const double upperLimit )
{
    return std::pow( upperLimit - lowerLimit, 2 * n + 1 ) * std::pow( factorial( n ), 4 )
            / double( ( 2 * n + 1 ) * std::pow( factorial( 2 * n ), 3 ) ) * derivative2nth( abscissa );
}

// Check error is within bounds and decreasing for increasing number of nodes.
// Number of nodes ranges from minOrder to maxOrder (both included).
// It must hold that lowerLimit ≤ abscissa ≤ upperLimit.
// The derivative is a function taking as input an `int n` (the nth derivative) that returns a function taking as input
// a `double x` (the absicssa at which it will be evaluated).
void checkErrorWithinBounds( const unsigned int minOrder, const unsigned int maxOrder,
                             boost::function< double( const double ) > function,
                             const double lowerLimit, const double upperLimit,
                             boost::function< boost::function< double( const double ) >( const unsigned int ) > derivative,
                             const double abscissa, const double expectedSolution )
{
    double obtainedError = TUDAT_NAN;
    double errorBound = TUDAT_NAN;
    for ( unsigned int n = minOrder; n <= maxOrder; n++ )
    {
        const double previousObtainedError = obtainedError;
        numerical_quadrature::GaussianQuadrature< double, double > integrator( function, lowerLimit, upperLimit, n );
        obtainedError = std::fabs( expectedSolution - integrator.getQuadrature() );
        errorBound = std::fabs( gaussianQuadratureError( n, derivative( 2 * n ), abscissa, lowerLimit, upperLimit ) );
        BOOST_CHECK( obtainedError < errorBound );
        if ( n > minOrder )
        {
            BOOST_CHECK( obtainedError < previousObtainedError );
        }
    }
}


BOOST_AUTO_TEST_SUITE( test_gaussian_quadrature )

//! Test if quadrature is computed correctly (sine function with 1E4 data points).
BOOST_AUTO_TEST_CASE( testIntegralSineFunction )
{
    using namespace mathematical_constants;
    using namespace numerical_quadrature;

    const unsigned int order = 8;
    const double lowerLimit = 0.0;
    const double upperLimit = PI;
    GaussianQuadrature< double, double > integrator( sinFunction, lowerLimit, upperLimit, order );
    double computedSolution = integrator.getQuadrature();

    // Expected solution from Wolfram Alpha
    double expectedSolution = 2.0;

    // Check if computed solution matches expected value for a high order (8).
    BOOST_CHECK_CLOSE_FRACTION( computedSolution, expectedSolution, 1E-10 );


    /// Check that it is a Gaussian quadrature and not just any quadrature

    // Check if error is within bounds for order 2...4
    checkErrorWithinBounds( 2, 4, sinFunction, lowerLimit, upperLimit, dSinFunction, PI / 2, expectedSolution );

    // Set to order 2
    integrator.reset( sinFunction, lowerLimit, upperLimit, 2 );
    computedSolution = integrator.getQuadrature();

    // Expected solution from http://keisan.casio.com/exec/system/1330940731
    expectedSolution = 1.9358195746511370184019497173109914411780661179299;

    // Check if computed solution matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSolution, expectedSolution, 1E-12 );
}


//! Test if quadrature is computed correctly (exponential function).
BOOST_AUTO_TEST_CASE( testIntegralExpFunction )
{
    using namespace mathematical_constants;
    using namespace numerical_quadrature;

    const unsigned int order = 7;
    const double lowerLimit = -2.0;
    const double upperLimit = 2.0;
    GaussianQuadrature< double, double > integrator( expFunction, lowerLimit, upperLimit, order );
    double computedSolution = integrator.getQuadrature();

    // Expected solution from Wolfram Alpha
    double expectedSolution = 7.25372081569404;

    // Check if computed solution matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSolution, expectedSolution, 1E-10 );


    /// Check that it is a Gaussian quadrature and not just any quadrature

    // Check if error is within bounds for order 2...6
    checkErrorWithinBounds( 2, 6, expFunction, lowerLimit, upperLimit, dExpFunction, 1.0, expectedSolution );

    // Set to order 2
    integrator.reset( expFunction, lowerLimit, upperLimit, 2 );
    computedSolution = integrator.getQuadrature();

    // Expected solution from http://keisan.casio.com/exec/system/1330940731
    expectedSolution = 6.9764499206151121988504219845072175547665355367817;

    // Check if computed solution matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSolution, expectedSolution, 1E-12 );
}


//! Test if quadrature is computed correctly (polynomial function).
BOOST_AUTO_TEST_CASE( testIntegralPolyFunction )
{
    using namespace mathematical_constants;
    using namespace numerical_quadrature;

    // A polynomial of order 2*steps - 1 is computed exactly. Order of tested polynomial is 9, thus steps = 5.
    const unsigned int order = 5;
    const double lowerLimit = -2.0;
    const double upperLimit = 4.0;
    GaussianQuadrature< double, double > integrator( polyFunction, lowerLimit, upperLimit, order );
    double computedSolution = integrator.getQuadrature();

    // Expected solution from Wolfram Alpha
    double expectedSolution = 120990;

    // Check if computed solution matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSolution, expectedSolution, 1E-12 );

    // Check that this isn't the case for 4 nodes
    integrator.reset( polyFunction, lowerLimit, upperLimit, 4 );
    computedSolution = integrator.getQuadrature();
    BOOST_CHECK( std::fabs( computedSolution - expectedSolution) > 1 );


    /// Check that it is a Gaussian quadrature and not just any quadrature

    // Check if error is within bounds for order 2...4
    checkErrorWithinBounds( 2, 4, polyFunction, lowerLimit, upperLimit, dPolyFunction, 2.0, expectedSolution );

    // Set to order 2
    integrator.reset( polyFunction, lowerLimit, upperLimit, 2 );
    computedSolution = integrator.getQuadrature();

    // Expected solution from http://keisan.casio.com/exec/system/1330940731
    expectedSolution = 32214;

    // Check if computed solution matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSolution, expectedSolution, 1E-12 );
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
