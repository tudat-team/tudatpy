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

#include <boost/test/unit_test.hpp>

#include <vector>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"

#include "Tudat/Mathematics/Interpolators/hermiteCubicSplineInterpolator.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{


BOOST_AUTO_TEST_SUITE( test_hermite_cubic_spline_interpolator )

// Test implementation of cubic spline class.
BOOST_AUTO_TEST_CASE( testHermiteCubicSplineInterpolator )
{
    // Test 1: Compare with analytical function 2 + 3x + 5x^2.
    // dy/dx = 10x + 3
    {
        // Declare and initialize independent variable values.
        std::vector< double > independentVariables;
        independentVariables.resize( 6 );
        independentVariables[ 0 ] = 1.0;
        independentVariables[ 1 ] = 3.0;
        independentVariables[ 2 ] = 5.0;
        independentVariables[ 3 ] = 7.0;
        independentVariables[ 4 ] = 9.0;
        independentVariables[ 5 ] = 11.0;

        // Declare and initialize dependent variable values.
        std::vector< double > dependentVariables;
        dependentVariables.resize( 6 );
        dependentVariables[ 0 ] = 10.0;
        dependentVariables[ 1 ] = 56.0;
        dependentVariables[ 2 ] = 142.0;
        dependentVariables[ 3 ] = 268.0;
        dependentVariables[ 4 ] = 434.0;
        dependentVariables[ 5 ] = 640.0;

        // Declare and initialize dependent variable values.
        std::vector< double > derivatives;
        derivatives.resize( 6 );
        derivatives[ 0 ] = 13.0;
        derivatives[ 1 ] = 33.0;
        derivatives[ 2 ] = 53.0;
        derivatives[ 3 ] = 73.0;
        derivatives[ 4 ] = 93.0;
        derivatives[ 5 ] = 113.0;

        // Declare and initialize target independent variable value.
        double targetIndependentVariableValue = 6.0;

        // Declare and initialize expected result of interpolation from analytical equation.
        double analyticalValue = 200.0;

        // Declare cubic spline object and initialize with input data.
        interpolators::HermiteCubicSplineInterpolatorDouble interpolator(
                    independentVariables, dependentVariables , derivatives );

        // Declare interpolated dependent variable value and execute interpolation.
        double interpolatedDependentVariableValue = interpolator.interpolate(
                    targetIndependentVariableValue );

        // Check if test result match analytical result.
        BOOST_CHECK_SMALL(  std::fabs( analyticalValue - interpolatedDependentVariableValue )
                            / analyticalValue, 1.0e-15 );

        // Second test
        targetIndependentVariableValue = 2.0;
        analyticalValue = 28.0;

        interpolatedDependentVariableValue = interpolator.interpolate(
                    targetIndependentVariableValue );

        // Check if test result match analytical result.
        BOOST_CHECK_SMALL(  std::fabs( analyticalValue - interpolatedDependentVariableValue )
                            / analyticalValue, 1.0e-15 );
    }
}

// Test implementation of cubic Hermite spline class against Matlab (pchipd from Mathworks central).
BOOST_AUTO_TEST_CASE( testHermiteCubicSplineInterpolatorWithMatlab )
{
    // Test 1: Compare with equispace solution for sin(2x)
    {
        // Declare and initialize independent variable values.
        std::vector< double > independentVariables;
        std::vector< double > dependentVariables;
        std::vector< double > derivatives;

        independentVariables.resize( 10 );
        dependentVariables.resize( 10 );
        derivatives.resize( 10 );

        for( unsigned int i = 0; i < 10; i++ )
        {
            independentVariables[ i ] = static_cast< double >( i ) * 2.0 / 3.0;
            dependentVariables[ i ] = std::sin( 2.0 * independentVariables[ i ] );
            derivatives[ i ] = 2.0 * std::cos( 2.0 * independentVariables[ i ] );
        }

        // Declare cubic spline object and initialize with input data.
        interpolators::HermiteCubicSplineInterpolatorDouble interpolator(
                    independentVariables, dependentVariables , derivatives );

        std::vector< double > matlabResults = { 0.0,        -0.500127989105607,        0.860566028073899, -0.995666072885886 };

        for( unsigned int i = 0; i < 4; i++ )
        {
            double targetIndependentVariableValue =  static_cast< double >( i ) * 11.0 / 6.0;
            // Declare interpolated dependent variable value and execute interpolation.
            double interpolatedDependentVariableValue = interpolator.interpolate(
                        targetIndependentVariableValue );

            // Check if test result match analytical result.
            BOOST_CHECK_SMALL( std::fabs( matlabResults.at( i ) - interpolatedDependentVariableValue ), 1.0e-15 );

        }
    }

    // Test 2: Compare with non-equispaced solution for sin(2x)
    {
        // Declare and initialize independent variable values.
        std::vector< double > independentVariables =
        { 0.100,        0.250,        0.600,        0.620,        1.0,        1.40,        1.51,        2.0,        3.0};
        std::vector< double > dependentVariables;
        std::vector< double > derivatives;

        for( unsigned int i = 0; i < independentVariables.size( ); i++ )
        {
            dependentVariables.push_back( std::sin( 2.0 * independentVariables[ i ] ) );
            derivatives.push_back( 2.0 * std::cos( 2.0 * independentVariables[ i ] ) );
        }

        // Declare cubic spline object and initialize with input data.
        interpolators::HermiteCubicSplineInterpolatorDouble interpolator(
                    independentVariables, dependentVariables , derivatives );

        std::vector< double > matlabResults =
        { -0.000114788941438321,        0.956262264186293,        -0.557211110925702,        -0.616700531016027 };

        for( unsigned int i = 0; i < 4; i++ )
        {
            double targetIndependentVariableValue =  static_cast< double >( i ) * 2.8 / 3.0;
            // Declare interpolated dependent variable value and execute interpolation.
            double interpolatedDependentVariableValue = interpolator.interpolate(
                        targetIndependentVariableValue );

            // Check if test result match analytical result.
            BOOST_CHECK_SMALL( std::fabs( matlabResults.at( i ) - interpolatedDependentVariableValue ), 2.0e-15 );

        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
