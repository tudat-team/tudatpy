/*! \file unitTestCubicSplineInterpolation.cpp
 *    Source file that defines the unit test for the cubic spline interplation
 *    included in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : e.a.g.heeren@student.tudelft.nl
 *
 *    Date created      : 21 June, 2011
 *    Last modified     : 5 September, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110621    F.M. Engelen      File created.
 *      110707    E.A.G. Heeren     Minor spelling/lay-out corrections.
 *      110714    E.A.G. Heeren     Minor spelling/lay-out corrections.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Mathematics/cubicSplineInterpolation.h"
#include "Mathematics/LinearAlgebra/linearAlgebra.h"

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
    VectorXd independentVariables = VectorXd( 6 );
    independentVariables( 0 ) = 1.0;
    independentVariables( 1 ) = 3.0;
    independentVariables( 2 ) = 5.0;
    independentVariables( 3 ) = 7.0;
    independentVariables( 4 ) = 9.0;
    independentVariables( 5 ) = 11.0;

    // Declare and initialize dependent variable values.
    VectorXd dependentVariables = VectorXd( 6 );
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
    CubicSplineInterpolation cubicSplineInterpolation;

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

// End of file.
