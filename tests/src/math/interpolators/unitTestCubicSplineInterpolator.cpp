/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <boost/test/unit_test.hpp>

#include <vector>

#include "tudat/basics/testMacros.h"
#include "tudat/io/matrixTextFileReader.h"

#include "tudat/math/interpolators/cubicSplineInterpolator.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

// NOTE: No benchmark data from established software package has been used. Matlab implementation
// is slightly different than that implemented here, yielding error of at worst 10e-7 compared
// when using error function, which is used here as test function.
// However, a Matlab implementation of the same algorithm was found on te Matlab forum,
// which is used here as benchmark.
// In addition, the code has been used rather extensively and has been checked for basic
// characteristics such as continuity.
BOOST_AUTO_TEST_SUITE( test_cubic_spline_interpolator )

// Test implementation of cubic spline class.
BOOST_AUTO_TEST_CASE( testCubicSplineInterpolator )
{
    // Test 1: Compare with analytical function 2 + 3x + 5x^2.
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

        // Declare and initialize target independent variable value.
        const double targetIndependentVariableValue = 6.0;

        // Declare and initialize expected result of interpolation from analytical equation.
        const double analyticalValue = 200.0;

        // Declare cubic spline object and initialize with input data.
        interpolators::CubicSplineInterpolatorDouble cubicSplineInterpolation( independentVariables, dependentVariables );

        // Declare interpolated dependent variable value and execute interpolation.
        const double interpolatedDependentVariableValue = cubicSplineInterpolation.interpolate( targetIndependentVariableValue );

        // Check if test result match analytical result.
        BOOST_CHECK_SMALL( std::fabs( analyticalValue - interpolatedDependentVariableValue ) / analyticalValue, 5.0e-3 );
    }
}

// Test exception handling implementation of cubic spline class.
BOOST_AUTO_TEST_CASE( testCubicSplineInterpolation_exception_empty_vectors )
{
    // Test 2: Interpolate with empty vectors.
    // Declare independent and dependent variable vectors.
    std::vector< double > independentVariables, dependentVariables;

    // Declare and initialize flag.
    bool areDependentAndIndependentVariablesInitialized = true;

    // Try to initialize with empty vectors.
    try
    {
        // Declare cubic spline object and initialize with input data.
        interpolators::CubicSplineInterpolatorDouble cubicSplineInterpolation( independentVariables, dependentVariables );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch( std::runtime_error const& )

    {
        areDependentAndIndependentVariablesInitialized = false;
    }

    // Check value of flag.
    BOOST_CHECK( !areDependentAndIndependentVariablesInitialized );
}

// Test cubic spline interpolator by comparing to Matlab code posted at
// http://www.mathworks.com/matlabcentral/newsreader/view_thread/173708.
BOOST_AUTO_TEST_CASE( test_cubicSplineInterpolator_matlab_forum_compare )
{
    using namespace interpolators;

    // Load input data used for generating matlab interpolation.
    Eigen::MatrixXd inputData =
            input_output::readMatrixFromFile( paths::getTudatTestDataPath( ) + "/interpolator_test_input_data.dat", "," );

    // Put data in STL vectors.
    std::vector< double > independentVariableValues;
    std::vector< double > dependentVariableValues;
    for( int i = 0; i < inputData.rows( ); i++ )
    {
        independentVariableValues.push_back( inputData( i, 0 ) );
        dependentVariableValues.push_back( inputData( i, 1 ) );
    }

    // Create cubic spline interpolator using hunting algorithm.
    CubicSplineInterpolatorDouble cubicSplineInterpolator( independentVariableValues, dependentVariableValues, huntingAlgorithm );

    // Load points at which interpolator is to be evaluated and data generated by Matlab.
    Eigen::MatrixXd benchmarkData =
            input_output::readMatrixFromFile( paths::getTudatTestDataPath( ) + "/cubic_spline_interpolator_test_output_data.dat", "," );

    // Perform interpolation for required data points.
    Eigen::VectorXd outputData = Eigen::VectorXd( benchmarkData.rows( ) );
    for( int i = 0; i < outputData.rows( ); i++ )
    {
        outputData[ i ] = cubicSplineInterpolator.interpolate( benchmarkData( i, 0 ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( benchmarkData.block( 0, 1, benchmarkData.rows( ), 1 ), outputData, 1.0e-13 );

    // Create cubic spline interpolator, now using binary search algorithm.
    cubicSplineInterpolator = CubicSplineInterpolatorDouble( independentVariableValues, dependentVariableValues, binarySearch );

    // Perform interpolation for required data points.
    outputData = Eigen::VectorXd( benchmarkData.rows( ) );
    for( int i = 0; i < outputData.rows( ); i++ )
    {
        outputData[ i ] = cubicSplineInterpolator.interpolate( benchmarkData( i, 0 ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( benchmarkData.block( 0, 1, benchmarkData.rows( ), 1 ), outputData, 1.0e-13 );
}

// Test cubic spline interpolator by comparing to Matlab implementation. Note that the two
// implementations are not identical, since the Matlab implementation imposes zero first
// derivatives whereas the present implementation imposes zero second derivatives at endpoints.
BOOST_AUTO_TEST_CASE( test_cubicSplineInterpolator_matlab_compare )
{
    using namespace interpolators;

    // Load input data used for generating matlab interpolation.
    Eigen::MatrixXd inputData =
            input_output::readMatrixFromFile( paths::getTudatTestDataPath( ) + "/interpolator_test_input_data.dat", "," );

    // Put data in STL vectors.
    std::vector< double > independentVariableValues;
    std::vector< double > dependentVariableValues;
    for( int i = 0; i < inputData.rows( ); i++ )
    {
        independentVariableValues.push_back( inputData( i, 0 ) );
        dependentVariableValues.push_back( inputData( i, 1 ) );
    }

    // Create cubic spline interpolator.
    CubicSplineInterpolatorDouble linearInterpolator( independentVariableValues, dependentVariableValues );

    // Load points at which interpolator is to be evaluated and data generated by Matlab.
    Eigen::MatrixXd benchmarkData = input_output::readMatrixFromFile(
            paths::getTudatTestDataPath( ) + "/cubic_spline_interpolator_approximate_test_output_data.dat", "," );

    // Perform interpolation for required data points.
    Eigen::VectorXd outputData = Eigen::VectorXd( benchmarkData.rows( ) );
    for( int i = 0; i < outputData.rows( ); i++ )
    {
        outputData[ i ] = linearInterpolator.interpolate( benchmarkData( i, 0 ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( benchmarkData.block( 0, 1, benchmarkData.rows( ), 1 ), outputData, 1.0e-5 );
}

// Test linear interpolation outside of independent variable range
BOOST_AUTO_TEST_CASE( test_cubicSplineInterpolator_boundary_case )
{
    using namespace interpolators;

    // Load input data used for generating matlab interpolation.
    Eigen::MatrixXd inputData =
            input_output::readMatrixFromFile( paths::getTudatTestDataPath( ) + "/interpolator_test_input_data.dat", "," );

    // Put data in STL vectors.
    std::vector< double > independentVariableValues;
    std::vector< double > dependentVariableValues;

    for( int i = 0; i < inputData.rows( ); i++ )
    {
        independentVariableValues.push_back( inputData( i, 0 ) );
        dependentVariableValues.push_back( inputData( i, 1 ) );
    }

    // Create linear interpolator using hunting algorithm.
    double valueOffset = 2.0;
    double valueBelowMinimumValue = independentVariableValues[ 0 ] - valueOffset;
    double valueAboveMaximumValue = independentVariableValues[ inputData.rows( ) - 1 ] + valueOffset;
    double interpolatedValue = TUDAT_NAN, expectedValue = TUDAT_NAN;
    bool exceptionIsCaught = false;

    for( unsigned int i = 0; i < 5; i++ )
    {
        CubicSplineInterpolatorDouble linearInterpolator(
                independentVariableValues, dependentVariableValues, huntingAlgorithm, static_cast< BoundaryInterpolationType >( i ) );

        if( static_cast< BoundaryInterpolationType >( i ) == throw_exception_at_boundary )
        {
            try
            {
                linearInterpolator.interpolate( valueBelowMinimumValue );
            }
            catch( tudat::exceptions::InterpolationOutOfBoundsError< double > const& )
            {
                exceptionIsCaught = true;
            }
            BOOST_CHECK_EQUAL( exceptionIsCaught, true );

            exceptionIsCaught = false;
            try
            {
                linearInterpolator.interpolate( valueAboveMaximumValue );
            }
            catch( tudat::exceptions::InterpolationOutOfBoundsError< double > const& )
            {
                exceptionIsCaught = true;
            }
            BOOST_CHECK_EQUAL( exceptionIsCaught, true );
        }
        else if( ( static_cast< BoundaryInterpolationType >( i ) == use_boundary_value ) ||
                 ( static_cast< BoundaryInterpolationType >( i ) == use_boundary_value_with_warning ) )
        {
            interpolatedValue = linearInterpolator.interpolate( valueBelowMinimumValue );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, dependentVariableValues.at( 0 ), 1.0E-15 );

            interpolatedValue = linearInterpolator.interpolate( valueAboveMaximumValue );
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, dependentVariableValues.at( inputData.rows( ) - 1 ), 1.0E-15 );
        }
        else if( ( static_cast< BoundaryInterpolationType >( i ) == extrapolate_at_boundary ) ||
                 ( static_cast< BoundaryInterpolationType >( i ) == extrapolate_at_boundary_with_warning ) )
        {
            interpolatedValue = linearInterpolator.interpolate( valueBelowMinimumValue );
            expectedValue = -1.0171383008483266;  // computed with MATLAB
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, expectedValue, 1.0E-15 );

            interpolatedValue = linearInterpolator.interpolate( valueAboveMaximumValue );
            expectedValue = 1.0171383008483266;  // computed with MATLAB
            BOOST_CHECK_CLOSE_FRACTION( interpolatedValue, expectedValue, 1.0E-15 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
