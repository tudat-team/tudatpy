/*! \file unitTestSimpleLinearRegression.cpp
 *    Source file that defines the unitTestSimpleLinearRegression unit test,
 *    which tests the simple linear regression method implemented in Tudat.
 *
 *    Path              : /Mathematics/Statistics/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 2 July, 2011
 *    Last modified     : 2 August, 2011
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 *    Notes
 *      The chi-squared value and standard deviations of the coefficients
 *      of the linear fit are not tested against benchmark data at present;
 *      they are tested against the output of the SimpleLinearRegression class.
 *      Once benchmark data that includes these values is available, this test
 *      must, be modified.
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
 *      110702    K. Kumar          First creation of code.
 *      110726    K. Kumar          Changed filename and class name.
 *      110802    K. Kumar          Added standard deviation and chi-squared
 *                                  test; added note; renamed filename.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include "Mathematics/Statistics/unitTestSimpleLinearRegression.h"

// Using declarations.
using std::cerr;
using std::endl;
using mathematics::MACHINE_PRECISION_DOUBLES;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test implementation of simple linear regression method.
bool testSimpleLinearRegression( )
{
    // Test implementation of simple linear regression method against benchmark
    // data from pg. 487, example 1 of (Burden and Faires, 2001).

    // Test result initialised to false.
    bool isSimpleLinearRegressionErroneous = false;

    // Benchmark data.
    std::map< double, double > benchmarkInputData;
    benchmarkInputData[ 1.0 ] = 1.3;
    benchmarkInputData[ 2.0 ] = 3.5;
    benchmarkInputData[ 3.0 ] = 4.2;
    benchmarkInputData[ 4.0 ] = 5.0;
    benchmarkInputData[ 5.0 ] = 7.0;
    benchmarkInputData[ 6.0 ] = 8.8;
    benchmarkInputData[ 7.0 ] = 10.1;
    benchmarkInputData[ 8.0 ] = 12.5;
    benchmarkInputData[ 9.0 ] = 13.0;
    benchmarkInputData[ 10.0 ] = 15.6;

    // Expected coefficients of linear fit.
    double expectedCoefficientOfConstantTerm = -0.360;
    double expectedCoefficientOfLinearTerm = 1.538;

    // Expected standard deviations of fit coefficients.
    double expectedStandardDeviationOfCoefficientOfConstantTerm
            = 1.224508762402819;
    double expectedStandardDeviationOfCoefficientOfLinearTerm
            = 0.6514750878514816;

    // Expected chi-squared value.
    double expectedChiSquared = 2.344727272727272;

    // Declare simple linear regression object.
    SimpleLinearRegression simpleLinearRegression;

    // Set input data.
    simpleLinearRegression.setInputData( benchmarkInputData );

    // Compute linear fit.
    simpleLinearRegression.computeFit( );

    // Check if coefficients of fit are incorrect.
    if ( abs( simpleLinearRegression.getCoefficientOfConstantTerm( )
              - expectedCoefficientOfConstantTerm )
         / expectedCoefficientOfConstantTerm > MACHINE_PRECISION_DOUBLES
         || abs( simpleLinearRegression.getCoefficientOfLinearTerm( )
                 - expectedCoefficientOfLinearTerm )
            / expectedCoefficientOfLinearTerm > MACHINE_PRECISION_DOUBLES )
    {
        // Set test result to true.
        isSimpleLinearRegressionErroneous = true;

        // Cerr statements.
        cerr << "Computed linear fit coefficient values ( "
             << simpleLinearRegression.getCoefficientOfConstantTerm( ) << ", "
             << simpleLinearRegression.getCoefficientOfLinearTerm( ) << " ) "
             << " are not equal to the expected values ( "
             << expectedCoefficientOfConstantTerm << ", "
             << expectedCoefficientOfLinearTerm << " )." << endl;
        cerr << "Simple linear regression method failed." << endl;
    }

    // Compute linear fit errors.
    simpleLinearRegression.computeFitErrors( );

    // Check if standard deviations of fit coefficients are incorrect.
    if ( abs( simpleLinearRegression
              .getStandardDeviationOfCoefficientOfConstantTerm( )
              - expectedStandardDeviationOfCoefficientOfConstantTerm )
         / expectedStandardDeviationOfCoefficientOfConstantTerm
         > MACHINE_PRECISION_DOUBLES
         || abs( simpleLinearRegression
                 .getStandardDeviationOfCoefficientOfLinearTerm( )
                 - expectedStandardDeviationOfCoefficientOfLinearTerm )
            / expectedStandardDeviationOfCoefficientOfLinearTerm
         > MACHINE_PRECISION_DOUBLES )
    {
        // Set test result to true.
        isSimpleLinearRegressionErroneous = true;

        // Cerr statements.
        cerr << "Computed standard deviations of linear fit coefficient "
             << "values ( " << simpleLinearRegression
                .getStandardDeviationOfCoefficientOfConstantTerm( )
             << ", " << simpleLinearRegression
                .getStandardDeviationOfCoefficientOfLinearTerm( ) << " ) "
             << " are not equal to the expected values ( "
             << expectedStandardDeviationOfCoefficientOfConstantTerm
             << ", " << expectedStandardDeviationOfCoefficientOfLinearTerm
             << " )." << endl;
        cerr << "Simple linear regression method failed." << endl;
    }

    // Check if the chi-squared value is incorrect.
    if ( abs( simpleLinearRegression.getChiSquared( ) - expectedChiSquared )
         / expectedChiSquared > MACHINE_PRECISION_DOUBLES )
    {
        // Set test result to true.
        isSimpleLinearRegressionErroneous = true;

        // Cerr statements.
        cerr << "Computed chi-squared of linear fit ( "
             << simpleLinearRegression.getChiSquared( ) << " ) "
             << " is not equal to the expected value ( "
             << expectedChiSquared << " )." << endl;
        cerr << "Simple linear regression method failed." << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isSimpleLinearRegressionErroneous;
}

}

// End of file.
