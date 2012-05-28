/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *
 *    References
 *
 */

#include <cmath>
#include <iostream>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolation.h"

//! Test implementation of cubic spline class.
int main( )
{
    // Using declarations.
    using namespace tudat;

    // Summary of tests.
    // Test 1: Compare with analytical function 2 + 3x + 5x^2.

    // Declare and initilize test result to false.
    bool isCubicSplineInterpolationBad = false;

    // Declare and initialize independent variable values.
    Eigen::VectorXd independentVariables = Eigen::VectorXd( 6 );
    independentVariables( 0 ) = 1.0;
    independentVariables( 1 ) = 3.0;
    independentVariables( 2 ) = 5.0;
    independentVariables( 3 ) = 7.0;
    independentVariables( 4 ) = 9.0;
    independentVariables( 5 ) = 11.0;

    // Declare and initialize dependent variable values.
    Eigen::VectorXd dependentVariables = Eigen::VectorXd( 6 );
    dependentVariables( 0 ) = 10.0;
    dependentVariables( 1 ) = 56.0;
    dependentVariables( 2 ) = 142.0;
    dependentVariables( 3 ) = 268.0;
    dependentVariables( 4 ) = 434.0;
    dependentVariables( 5 ) = 640.0;

    // Declare and initialize target independent variable value.
    double targetIndependentVariableValue = 6.0;

    // Declare and initialize expected result of interpolation from analytical equation.
    double analyticalValue = 200.0;

    // Declare interpolated dependent variable value.
    double interpolatedDependentVariableValue;

    // Declare cubic spline object.
    tudat::mathematics::interpolators::CubicSplineInterpolation cubicSplineInterpolation;

    // Initialize cubic spline interpolation with input data.
    cubicSplineInterpolation.initializeCubicSplineInterpolation(
            independentVariables,dependentVariables );

    // Execute interpolation.
    interpolatedDependentVariableValue = cubicSplineInterpolation.interpolate(
                targetIndependentVariableValue );

    // Check if test result match analytical result or output cerr statements.
    if ( std::fabs( analyticalValue - interpolatedDependentVariableValue )
         / analyticalValue > 1.0e-2 )
    {
        isCubicSplineInterpolationBad = true;
        std::cerr << "Cubic spline interpolation is malfunctioning." << std::endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isCubicSplineInterpolationBad )
    {
        std::cerr << "cubicSplineInterpolation failed! " << std::endl;
    }

    return isCubicSplineInterpolationBad;
}
