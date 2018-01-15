/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>
#include <map>

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/Statistics/simpleLinearRegression.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_simple_linear_regression )

//! Test if simple linear regression method computes fit correctly.
BOOST_AUTO_TEST_CASE( testSimpleLinearRegressionBF )
{
    // Test 1: Test implementation of simple linear regression method against benchmark data from
    //         pg. 487, example 1 of (Burden and Faires, 2001). Standard deviations were benchmarked
    //         using MATLAB's lscov() function.

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
    const double expectedCoefficientOfConstantTerm = -0.359999999999999999;
    const double expectedCoefficientOfLinearTerm = 1.5381818181818181818;

    // Expected standard deviations of fit coefficients.
    const double expectedStandardDeviationOfCoefficientOfConstantTerm = 0.369832066721825;
    const double expectedStandardDeviationOfCoefficientOfLinearTerm = 0.059603834439483;

    // Expected chi-squared value.
    const double expectedChiSquared = 2.344727272727272;

    // Declare simple linear regression object and set input data.
    statistics::SimpleLinearRegression simpleLinearRegression( benchmarkInputData );

    // Compute linear fit.
    simpleLinearRegression.computeFit( );

    // Check that computed coefficient of constant term matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedCoefficientOfConstantTerm,
                                simpleLinearRegression.getCoefficientOfConstantTerm( ),
                                1.0e-14 );

    // Check that computed coefficient of linear term matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedCoefficientOfLinearTerm,
                                simpleLinearRegression.getCoefficientOfLinearTerm( ),
                                1.0e-15 );

    // Compute linear fit errors.
    simpleLinearRegression.computeFitErrors( );

    // Check that computed standard deviation of coefficient of constant term matches expected
    // value.
    BOOST_CHECK_CLOSE_FRACTION( expectedStandardDeviationOfCoefficientOfConstantTerm,
                                simpleLinearRegression
                                .getStandardDeviationOfCoefficientOfConstantTerm( ),
                                1.0e-13 );

    // Check that computed standard deviation of coefficient of linear term matches expected
    // value.
    BOOST_CHECK_CLOSE_FRACTION( expectedStandardDeviationOfCoefficientOfLinearTerm,
                                simpleLinearRegression
                                .getStandardDeviationOfCoefficientOfLinearTerm( ),
                                1.0e-13 );

    // Check that computed chi-squared fit matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedChiSquared,
                                simpleLinearRegression.getChiSquared( ),
                                1.0e-15 );
}

BOOST_AUTO_TEST_CASE( testSimpleLinearRegressionHorizontal )
{
    // Test 2: Test implementation of simple linear regression method in case of sample points
    // coinciding with the x-axis.

    std::map< double, double > benchmarkInputData;
    benchmarkInputData[ 1.0 ] = 0.0;
    benchmarkInputData[ 2.0 ] = 0.0;
    benchmarkInputData[ 3.0 ] = 0.0;
    benchmarkInputData[ 4.0 ] = 0.0;
    benchmarkInputData[ 5.0 ] = 0.0;
    benchmarkInputData[ 6.0 ] = 0.0;
    benchmarkInputData[ 7.0 ] = 0.0;
    benchmarkInputData[ 8.0 ] = 0.0;
    benchmarkInputData[ 9.0 ] = 0.0;
    benchmarkInputData[ 10.0 ] = 0.0;

    // Declare simple linear regression object and set input data.
    statistics::SimpleLinearRegression simpleLinearRegression( benchmarkInputData );

    // Compute linear fit.
    simpleLinearRegression.computeFit( );

    // Check that computed coefficient of constant term is zero.
    BOOST_CHECK_SMALL( simpleLinearRegression.getCoefficientOfConstantTerm( ),
                       std::numeric_limits< double >::min( ) );

    // Check that computed coefficient of linear term matches is zero.
    BOOST_CHECK_SMALL( simpleLinearRegression.getCoefficientOfLinearTerm( ),
                       std::numeric_limits< double >::min( ) );

    // Compute linear fit errors.
    simpleLinearRegression.computeFitErrors( );

    // Check that computed standard deviation of coefficient of constant term is zero.
    BOOST_CHECK_SMALL( simpleLinearRegression.getStandardDeviationOfCoefficientOfConstantTerm( ),
                       std::numeric_limits< double >::min( ) );


    // Check that computed standard deviation of coefficient of linear term is zero.
    BOOST_CHECK_SMALL( simpleLinearRegression.getStandardDeviationOfCoefficientOfLinearTerm( ),
                       std::numeric_limits< double >::min( ) );

    // Check that computed chi-squared fit is zero.
    BOOST_CHECK_SMALL( simpleLinearRegression.getChiSquared( ),
                       std::numeric_limits< double >::min( ) );
}

BOOST_AUTO_TEST_CASE( testSimpleLinearRegressionVertical )
{
    // Test 3: Test implementation of simple linear regression method in case of sample points
    // coinciding with the y-axis.

    std::map< double, double > benchmarkInputData;
    benchmarkInputData[ 0.0 ] = 1.3;
    benchmarkInputData[ 0.0 ] = 3.5;
    benchmarkInputData[ 0.0 ] = 4.2;
    benchmarkInputData[ 0.0 ] = 5.0;
    benchmarkInputData[ 0.0 ] = 7.0;
    benchmarkInputData[ 0.0 ] = 8.8;
    benchmarkInputData[ 0.0 ] = 10.1;
    benchmarkInputData[ 0.0 ] = 12.5;
    benchmarkInputData[ 0.0 ] = 13.0;
    benchmarkInputData[ 0.0 ] = 15.6;

    // Declare simple linear regression object and set input data.
    statistics::SimpleLinearRegression simpleLinearRegression( benchmarkInputData );

    // Compute linear fit.
    simpleLinearRegression.computeFit( );

    // Check that computed coefficient of constant term matches expected value.
    BOOST_CHECK( boost::math::isnan( simpleLinearRegression.getCoefficientOfConstantTerm( ) ) );

    // Check that computed coefficient of linear term matches expected value.
    BOOST_CHECK( boost::math::isnan( simpleLinearRegression.getCoefficientOfLinearTerm( ) ) );

    // Compute linear fit errors.
    simpleLinearRegression.computeFitErrors( );

    // Check that computed standard deviation of coefficient of constant term matches expected
    // value.
    BOOST_CHECK( boost::math::isnan( simpleLinearRegression
                             .getStandardDeviationOfCoefficientOfConstantTerm( ) ) );

    // Check that computed standard deviation of coefficient of linear term matches expected
    // value.
    BOOST_CHECK( boost::math::isnan(simpleLinearRegression
                            .getStandardDeviationOfCoefficientOfLinearTerm( ) ) );

    // Check that computed chi-squared fit matches expected value.
    BOOST_CHECK( boost::math::isnan(simpleLinearRegression.getChiSquared( ) ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
