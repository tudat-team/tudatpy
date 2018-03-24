/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/numericalDerivative.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_basic_functions )

typedef Eigen::MatrixXd( *DerivativeCallback )( const Eigen::VectorXd& );
typedef Eigen::VectorXd( *FunctionEvaluationCallback )( const Eigen::VectorXd& );

//! Run a test case with the exponential density.
void compareAnalyticalAndNumericalDerivative( const Eigen::VectorXd& input,
                                              DerivativeCallback analyticalCallback,
                                              FunctionEvaluationCallback numericalCallback )
{
    std::cout << std::numeric_limits< double >::epsilon( ) << std::endl;
    using numerical_derivatives::computeCentralDifference;

    // Compute the expected partial.
    Eigen::MatrixXd analyticalDerivative = analyticalCallback( input );

    // Test 2nd-order numerical derivative.
    {
        Eigen::MatrixXd numericalDerivative = computeCentralDifference( input, numericalCallback );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalDerivative, numericalDerivative,
                                           1.0e-7 );
    }

    // Test 4th-order numerical derivative.
    {
        Eigen::MatrixXd numericalDerivative = computeCentralDifference(
                    input, numericalCallback, 0.0, 0.0,
                    numerical_derivatives::order4 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalDerivative, numericalDerivative,
                                           1.0e-9 );

    }

    // Test 8th-order numerical derivative.
    {
        Eigen::MatrixXd numericalDerivative = computeCentralDifference(
                    input, numericalCallback, 0.0, 0.0,
                    numerical_derivatives::order8 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalDerivative, numericalDerivative,
                                           1.0e-9 );
    }
}

//! Callback that returns the analytical derivative of exp( ||r|| ).
Eigen::MatrixXd exponentialDensityAnalyticalDerivative( const Eigen::VectorXd& position )
{
    return position.transpose( ) / position.norm( ) * std::exp( position.norm( ) );
}

//! Callback to test numerical derivative calculation, returns exp( ||r|| ).
Eigen::VectorXd exponentialDensity( const Eigen::VectorXd& position )
{
    return Eigen::VectorXd::Constant( 1, std::exp( position.norm( ) ) );
}

//! Callback that returns the analytical derivative of -r/||r||.
Eigen::MatrixXd constantGravityAnalyticalDerivative( const Eigen::VectorXd& position )
{
    Eigen::MatrixXd partial = position * position.transpose( ) / std::pow( position.norm( ), 3.0 );
    partial.diagonal( ) = position.array( ).pow( 2.0 ) / std::pow( position.norm( ), 3.0 )
            - 1.0 / position.norm( );
    return partial;
}

//! Callback to test numerical derivative calculation, returns -r/||r||.
Eigen::VectorXd constantGravity( const Eigen::VectorXd& position )
{
    return -position / position.norm( );
}

BOOST_AUTO_TEST_CASE( testNumericalDerivatives )
{
    // Test numerical derivatives.
    Eigen::MatrixXd positions( 10, 3 );
    positions << 0.0416284088706, 0.365492068944, 0.805197604602,
                 1.63074391170,   8.04179355586,  6.74984731916,
                 55.4620045731,   86.8094364606,  95.4087064974,
                 41.3971344853,   80.6456253401,  359.560049206,
                 6389.36995846,   1891.72249537,  3768.41346114,
                 18357.5991764,   5355.13286809,  24582.5658116,
                 19887.9880951,   390769.463405,  949457.32454,
                 3634565.52581,   5564841.99331,  3208769.36002,
                 45156443.1799,   3463879.96686,  97241290.6455,
                 478348640.774,   6705325.08872,  335953979.068;

    // Loop through the input data and pass position data to check function.
    for ( int i = 0; i < positions.cols( ); i++ )
    {
        Eigen::VectorXd position = positions.block< 1, 3 >( i, 0 ).transpose( );

        compareAnalyticalAndNumericalDerivative( position,
                                                 &exponentialDensityAnalyticalDerivative,
                                                 &exponentialDensity );

        compareAnalyticalAndNumericalDerivative( position,
                                                 &constantGravityAnalyticalDerivative,
                                                 &constantGravity );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

