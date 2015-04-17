/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      110207    B. Romgens        File created.
 *      110215    K. Kumar          Minor modifications to layout, comments
 *                                  and variable-naming.
 *      110411    K. Kumar          Added unit test for
 *                                  convertCartesianToSpherical( ) function.
 *      110701    K. Kumar          Updated failing tests with relative errors.
 *      110708    K. Kumar          Added unit tests for computeSampleMean( )
 *                                  and computeSampleVariance( ) functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      111111    K. Kumar          Strange error with convertCylindricalToCartesian function;
 *                                  achieved precision of results is less than machine precision,
 *                                  fixed by using slightly larger precision tolerance.
 *      120202    K. Kumar          Separated from unitTestBasicMathematics.cpp into new
 *                                  Interpolators sub-directory.
 *      120529    E.A.G. Heeren     Boostified unit test.
 *      120615    T. Secretin       Minor layout changes.
 *      120716    D. Dirkx          Updated with interpolator architecture.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"

#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_linear_interpolation )

// Test implementation of linear interpolation function for data vectors.
BOOST_AUTO_TEST_CASE( test_linearInterpolation_vector )
{
    // Set vectors of data.
    const Eigen::Vector3d sortedIndependentVariables( 0.0, 1.0, 3.0 );
    const Eigen::Vector3d associatedDependentVariables( -20.0, 20.0, 21.0 );

    // Test 1: Test vector data with expected result of 0.0.
    {
        // Set target independent value in vector data.
        const double targetIndependentVariableValue = 0.5;

        // Compute interpolation.
        const double interpolatedValue = interpolators::computeLinearInterpolation(
                    sortedIndependentVariables, associatedDependentVariables,
                    targetIndependentVariableValue );

        // Verify that interpolated value corresponds to expected value (0.0).
        BOOST_CHECK_SMALL( std::fabs( interpolatedValue - 0.0 ),
                           std::numeric_limits< double >::min( ) );
    }

    // Test 2: Test vector data with expected result of 20.5.

    // Set target independent value in vector data.
    const double targetIndependentVariableValue = 2.0;

    // Compute interpolation.
    const double interpolatedValue = interpolators::computeLinearInterpolation(
                sortedIndependentVariables, associatedDependentVariables,
                targetIndependentVariableValue );

    // Verify that interpolated value corresponds to expected value (20.5).
    BOOST_CHECK_SMALL( std::fabs( interpolatedValue - 20.5 ),
                       std::numeric_limits< double >::min( ) );
}

// Test linear interpolation with map of vectors with keys as independent variable.
BOOST_AUTO_TEST_CASE( test_linearInterpolation_map )
{
    // Declare map of data and vectors for map value.
    std::map < double, Eigen::VectorXd > sortedIndepedentAndDependentVariables;
    const Eigen::Vector3d vectorOne( 10.0, -10.0, 70.0 );
    const Eigen::Vector3d vectorTwo( 20.0, -5.0, 80.0 );
    const Eigen::Vector3d vectorThree( 30.0, 60.0, 90.0 );

    // Set map values in map using vector data.
    sortedIndepedentAndDependentVariables[ 0.0 ] = vectorOne;
    sortedIndepedentAndDependentVariables[ 1.0 ] = vectorTwo;
    sortedIndepedentAndDependentVariables[ 2.0 ] = vectorThree;

    // Set target independent variable value for interpolation.
    const double targetIndependentVariableValue = 1.5;

    // Compute interpolation.
    Eigen::Vector3d interpolatedVector = interpolators::computeLinearInterpolation(
                sortedIndepedentAndDependentVariables,
                targetIndependentVariableValue );

    // Check that interpolated values correspond to expected elements [25, 27.5, 85].
    BOOST_CHECK_SMALL( std::fabs( interpolatedVector( 0 ) - 25.0 ),
                       std::numeric_limits< double >::min( ) );

    BOOST_CHECK_SMALL( std::fabs( interpolatedVector( 1 ) - 27.5 ),
                       std::numeric_limits< double >::min( ) );

    BOOST_CHECK_SMALL( std::fabs( interpolatedVector( 2 ) - 85.0 ),
                       std::numeric_limits< double >::min( ) );
}

// Test linear interpolation from benchmark values generetd by Matlab, interpolating the
// error function.
BOOST_AUTO_TEST_CASE( test_linearInterpolation_matlab_compare )
{
    using namespace interpolators;

    // Load input data used for generating matlab interpolation.
    Eigen::MatrixXd inputData = input_output::readMatrixFromFile(
                input_output::getTudatRootPath( ) +
                "Mathematics/Interpolators/UnitTests/interpolator_test_input_data.dat","," );

    // Put data in STL vectors.
    std::vector< double > independentVariableValues;
    std::vector< double > dependentVariableValues;

    for ( int i = 0; i < inputData.rows( ); i++ )
    {
        independentVariableValues.push_back( inputData( i, 0 ) );
        dependentVariableValues.push_back( inputData( i, 1 ) );
    }

    // Create linear interpolator using hunting algorithm.
    LinearInterpolatorDouble linearInterpolator(
                independentVariableValues, dependentVariableValues, huntingAlgorithm );

    // Load points at which interpolator is to be evaluated and data generated by Matlab.
    Eigen::MatrixXd benchmarkData = input_output::readMatrixFromFile(
                input_output::getTudatRootPath( ) +
                "Mathematics/Interpolators/UnitTests/linear_interpolator_test_output_data.dat",
                "," );

    // Perform interpolation for required data points.
    Eigen::VectorXd outputData = Eigen::VectorXd( benchmarkData.rows( ) );
    for ( int i = 0; i < outputData.rows( ); i++ )
    {
        outputData[ i ] = linearInterpolator.interpolate( benchmarkData( i, 0 ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( benchmarkData.block( 0, 1, benchmarkData.rows( ), 1 ),
                                       outputData, 1.0E-13 );

    // Create linear interpolator, now with nearest neighbvour search.
    linearInterpolator = LinearInterpolatorDouble(
                independentVariableValues, dependentVariableValues, binarySearch );

    // Perform interpolation for required data points.
    for ( int i = 0; i < outputData.rows( ); i++ )
    {
        outputData[ i ] = linearInterpolator.interpolate( benchmarkData( i, 0 ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( benchmarkData.block( 0, 1, benchmarkData.rows( ), 1 ),
                                       outputData, 1.0E-13 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
