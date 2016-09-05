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
 *      110621    F.M. Engelen      File created.
 *      110707    E.A.G. Heeren     Minor spelling/lay-out corrections.
 *      110714    E.A.G. Heeren     Minor spelling/lay-out corrections.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120529    E.A.G. Heeren     Boostified unit test.
 *      120615    T. Secretin       Added check for exception handling.
 *      120716    D. Dirkx          Updated with interpolator architecture.
 *
 *    References
 *
 *    Notes
 *
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


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
