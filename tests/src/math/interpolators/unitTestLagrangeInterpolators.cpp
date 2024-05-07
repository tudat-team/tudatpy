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

#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/math/interpolators/lagrangeInterpolator.h"

namespace tudat
{
namespace unit_tests
{

//! Function to evaluate polynomial
/*!
 *  Function to evaluate polynomial with coefficients and independent variable as input.
 *  \param coefficients Polynomial coefficients with the coefficient as map value and order as key.
 *  \param evaluationPoint Independent variable at which polynomial is to be evaluated.
 *  \return Polynomial value.
 */
double evaluatePolynomial( const std::map< int, double >& coefficients,
                           const double evaluationPoint )
{
    double polynomialValue = 0.0;
    for( std::map< int, double >::const_iterator it = coefficients.begin( );
         it != coefficients.end( ) ; it++ )
    {
        polynomialValue += it->second * std::pow( evaluationPoint, it->first );
    }
    return polynomialValue;
}

//! Function to retrieve polynomial coefficients
/*!
 *  Function to retrieve quasi-random polynomial coefficients, up to a given maximum order.
 *  \param polynomialOrder Order of polynomial.
 *  \return Polynomial coefficients with the coefficient as map value and order as key.
 */
std::map< int, double > getPolynomialCoefficients( const int polynomialOrder)
{
    std::map< int, double > allCoefficients;
    allCoefficients[ 0 ] = 8.05425;
    allCoefficients[ 1 ] = 2.540;
    allCoefficients[ 2 ] = -0.454;
    allCoefficients[ 3 ] = 1.1224;
    allCoefficients[ 4 ] = 0.03545;
    allCoefficients[ 5 ] = -0.004;
    allCoefficients[ 6 ] = 0.0784;
    allCoefficients[ 7 ] = -0.000334;
    allCoefficients[ 8 ] = -0.00004743;
    allCoefficients[ 9 ] = 0.000007284;
    allCoefficients[ 10 ] = 0.00000134;
    allCoefficients[ 11 ] = -0.000000324;

    // Copy subset of allCoefficients map into currentCoefficients.
    std::map< int, double > currentCoefficients( allCoefficients.begin(),
        std::next( allCoefficients.begin(), polynomialOrder + 1 ) );
    return currentCoefficients;

}


//! Create quasi-random vector of non-uniform independent variables
/*!
 *  Create quasi-random vector of non-uniform independent variables
 *  \return Non-uniform, but continuously increasing, set of independent variables.
 */
std::vector< double > getIndependentVariableVector( )
{
    std::vector< double > independentVariableVector;
    independentVariableVector.push_back( 0.0 );
    independentVariableVector.push_back( 0.1 );
    independentVariableVector.push_back( 0.2 );
    independentVariableVector.push_back( 0.3 );
    independentVariableVector.push_back( 0.45 );
    independentVariableVector.push_back( 0.7 );
    independentVariableVector.push_back( 1.0 );
    independentVariableVector.push_back( 1.4 );
    independentVariableVector.push_back( 2.0 );
    independentVariableVector.push_back( 2.1 );
    independentVariableVector.push_back( 2.5 );
    independentVariableVector.push_back( 4.1 );
    independentVariableVector.push_back( 5.7 );
    independentVariableVector.push_back( 6.3 );
    independentVariableVector.push_back( 8.9 );
    independentVariableVector.push_back( 10.2 );
    independentVariableVector.push_back( 11.8 );
    independentVariableVector.push_back( 12.4 );
    independentVariableVector.push_back( 15.5 );
    independentVariableVector.push_back( 16.4 );
    independentVariableVector.push_back( 22.0 );
    independentVariableVector.push_back( 25.0 );
    independentVariableVector.push_back( 30.89 );
    independentVariableVector.push_back( 35.21 );
    independentVariableVector.push_back( 40.38 );
    independentVariableVector.push_back( 43.23 );
    independentVariableVector.push_back( 52.3 );
    independentVariableVector.push_back( 72.0 );
    independentVariableVector.push_back( 89.0 );
    independentVariableVector.push_back( 104.0 );
    return independentVariableVector;
}


BOOST_AUTO_TEST_SUITE( test_lagrange_interpolation )

// Test whetehr Lagrange interpolator can properly reproduce polynomial interpolation
// Since Lagrange interpolation uses a unique (n-1)th order polynomial to fit n data points,
// using an (n-1)th order polynomial as depedent variables should yield an exact reporduction
// of the original polynomial (barring numerical losses).
BOOST_AUTO_TEST_CASE( test_lagrange_interpolation_polynomials )
{
    std::map< double, double > dataMap;
    std::map< int, double > coefficients;
    std::vector< double > independentVariableVector = getIndependentVariableVector( );

    // Test interpolator for 4;6;8;10 data points per interpolant
    // (i.e. 3rd, 5th, 7th and 9th order polynomial)
    for( unsigned int stages = 4; stages < 11; stages += 2 )
    {
        dataMap.clear( );

        // Get polynomial coefficients for current number of points
        coefficients = getPolynomialCoefficients( stages - 1 );

        // Generate dependent variables
        for( unsigned int i = 0; i < independentVariableVector.size( ); i++ )
        {
            dataMap[ independentVariableVector.at( i ) ] =
                    evaluatePolynomial( coefficients, independentVariableVector.at( i ) );
        }

        // Create interpolator
        interpolators::LagrangeInterpolator< double, double > interpolator =
                interpolators::LagrangeInterpolator< double, double >(
                    dataMap, stages, interpolators::huntingAlgorithm,
                    interpolators::lagrange_no_boundary_interpolation );

        // Iterate over all data points inside allowed (i.e. non-boundary) range
        int offsetEntries = stages / 2 - 1;
        for( unsigned int i = offsetEntries;
             i < independentVariableVector.size( ) - ( offsetEntries + 2 ); i++ )
        {
            // Test current interval at 10 equispaced points
            double currentStepSize =  ( independentVariableVector.at( i + 1 ) -
                                        independentVariableVector.at( i ) ) / 10.0;
            for( unsigned j = 0; j < 10; j ++ )
            {
                double currentDataPoint = independentVariableVector.at( i ) +
                        static_cast< double >( j ) * currentStepSize;

                // Check interpolated value against theoretical polynomial
                if( j < 9 )
                {
                    BOOST_CHECK_CLOSE_FRACTION( interpolator.interpolate( currentDataPoint ),
                                                evaluatePolynomial( coefficients, currentDataPoint ),
                                                5.0E-15 );
                }
                else
                {
                    BOOST_CHECK_CLOSE_FRACTION( interpolator.interpolate( currentDataPoint ),
                                                evaluatePolynomial( coefficients, currentDataPoint ),
                                                2.0E-14 );
                }

            }
        }
    }
}

// Test to check whether the various boundary handling methopds are properly implemented
BOOST_AUTO_TEST_CASE( test_lagrange_interpolation_boundary )
{
    std::vector< double > dataVector;
    std::map< int, double > coefficients;
    std::vector< double > independentVariableVector = getIndependentVariableVector( );

    unsigned int independentVariableVectorSize = independentVariableVector.size( );

    for( int testCase = 0; testCase < 4; testCase++ )
    {
        // Test interpolator for 4;6;8;10 data points per interpolant
        // (i.e. 3rd, 5th, 7th and 9th order polynomial)
        for( unsigned int stages = 4; stages < 11; stages += 2 )
        {
            dataVector.clear( );

            // Get polynomial coefficients for current number of points
            coefficients = getPolynomialCoefficients( stages - 1 );

            // Generate dependent variables
            for( unsigned int i = 0; i < independentVariableVectorSize; i++ )
            {
                dataVector.push_back( evaluatePolynomial(
                                          coefficients, independentVariableVector.at( i ) ) );
            }

            // Create interpolator with cubic spline interpolation at boundaries
            interpolators::LagrangeInterpolator< double, double > lagrangeInterpolator =
                    interpolators::LagrangeInterpolator< double, double >(
                        independentVariableVector, dataVector, stages,
                        interpolators::huntingAlgorithm, static_cast< interpolators::LagrangeInterpolatorBoundaryHandling >( testCase ) );

            // Create spline interpolator from edge points at lower bound.
            int dataPointsForSpline = ( stages / 2 > 4 ) ? ( stages / 2 ) : 4;
            std::map< double, double > boundaryMap;
            for( int i = 0; i < dataPointsForSpline; i++ )
            {
                boundaryMap[ independentVariableVector.at( i ) ] = dataVector.at( i );
            }
            interpolators::CubicSplineInterpolator< double, double > lowerBoundInterpolator =
                    interpolators::CubicSplineInterpolator< double, double >( boundaryMap );


            // Test whether the Lagrange interpolators correctly evaluate the cubic spline
            // polynomials at the lower edges.
            for( unsigned int i = 0; i < stages / 2 - 1; i++ )
            {
                double currentStepSize = ( independentVariableVector.at( i + 1 ) -
                                           independentVariableVector.at( i ) ) / 10.0;
                for( unsigned int j = 0; j < 10; j ++ )
                {
                    double currentTestIndependentVariable = independentVariableVector.at( i ) +
                            static_cast< double >( i ) * currentStepSize;
                    double interpolatedValue = lagrangeInterpolator.interpolate( currentTestIndependentVariable );
                    if( testCase < 2 )
                    {
                        BOOST_CHECK_EQUAL(
                            interpolatedValue,
                            lowerBoundInterpolator.interpolate(
                                currentTestIndependentVariable ));
                    }
                    else
                    {
                        BOOST_CHECK_EQUAL( !( interpolatedValue == interpolatedValue ), true );
                    }
                }
            }

            // Create spline interpolator from edge points at upper bound.
            boundaryMap.clear( );
            for( unsigned int i = independentVariableVectorSize - dataPointsForSpline;
                 i < independentVariableVectorSize; i++ )
            {
                boundaryMap[ independentVariableVector.at( i ) ] = dataVector.at( i );
            }
            interpolators::CubicSplineInterpolator< double, double > upperBoundInterpolator =
                    interpolators::CubicSplineInterpolator< double, double >( boundaryMap );


            // Test whether the Lagrange interpolators correctly evaluate the cubic spline
            // polynomials at the uper edges.
            for( unsigned int i = 0; i < stages / 2 - 1; i++ )
            {
                double currentStepSize =
                        ( independentVariableVector.at( independentVariableVectorSize - i - 1 ) -
                          independentVariableVector.at( independentVariableVectorSize - i - 2 ) ) /
                        10.0;
                for( unsigned int j = 0; j < 10; j ++ )
                {
                    double currentTestIndependentVariable = independentVariableVector.at(
                                independentVariableVectorSize - i - 2 ) +
                            static_cast< double >( i ) * currentStepSize;
                    double interpolatedValue = lagrangeInterpolator.interpolate( currentTestIndependentVariable );

                    if( testCase < 2 )
                    {
                        BOOST_CHECK_EQUAL( interpolatedValue,
                                           upperBoundInterpolator.interpolate(
                                               currentTestIndependentVariable ));
                    }
                    else
                    {
                        BOOST_CHECK_EQUAL( !( interpolatedValue == interpolatedValue ), true );
                    }
                }
            }
        }
    }

    // Test whether an error is thrown if lagrange_no_boundary_interpolation is
    // selected an interpolation at the boundaries is requested.
    bool runtimeErrorOccurred;
    {
        // Test interpolators with various number of stages
        for( unsigned int stages = 4; stages < 11; stages += 2 )
        {
            dataVector.clear( );

            // Get polynomial coefficients for current number of points
            coefficients = getPolynomialCoefficients( stages - 1 );

            // Generate dependent variables
            for( unsigned int i = 0; i < independentVariableVectorSize; i++ )
            {
                dataVector.push_back( evaluatePolynomial(
                                          coefficients, independentVariableVector.at( i ) ) );
            }

            // Create interpolator lagrange_no_boundary_interpolation
            interpolators::LagrangeInterpolator< double, double > lagrangeInterpolator =
                    interpolators::LagrangeInterpolator< double, double >(
                        independentVariableVector, dataVector, stages,
                        interpolators::huntingAlgorithm,
                        interpolators::lagrange_no_boundary_interpolation );

            // Test for each whether the interpolator correctly throws an exception if
            // interpolation at the lower boundary is requested.
            for( unsigned int i = 0; i < stages / 2 - 1; i++ )
            {
                double currentStepSize =
                        ( independentVariableVector.at( i + 1 ) -
                          independentVariableVector.at( i ) ) / 3.0;

                for( unsigned int j = 0; j < 3; j ++ )
                {
                    try
                    {
                        double currentTestIndependentVariable =
                                independentVariableVector.at( i ) +
                                static_cast< double >( j ) * currentStepSize;
                        lagrangeInterpolator.interpolate( currentTestIndependentVariable );

                    }
                    catch( std::runtime_error const& )
                    {
                        runtimeErrorOccurred = 1;
                    }
                    BOOST_CHECK_EQUAL( runtimeErrorOccurred, 1 );
                    runtimeErrorOccurred = 0;
                }
            }

            // Test for each whether the interpolator correctly throws an exception if
            // interpolation at the upper boundary is requested.
            for( unsigned int i = 0; i < stages / 2 - 1; i++ )
            {
                double currentStepSize = ( independentVariableVector.at(
                                               independentVariableVectorSize - i - 1 ) -
                                           independentVariableVector.at(
                                               independentVariableVectorSize - i - 2 ) ) / 3.0;
                for( unsigned int j = 0; j < 3; j ++ )
                {
                    try
                    {
                        double currentTestIndependentVariable =
                                independentVariableVector.at(
                                    independentVariableVectorSize - i - 2 ) +
                                static_cast< double >( j ) * currentStepSize;
                        lagrangeInterpolator.interpolate( currentTestIndependentVariable );

                    }
                    catch( std::runtime_error const& )
                    {
                        runtimeErrorOccurred = true;
                    }
                    BOOST_CHECK_EQUAL( runtimeErrorOccurred, 1 );
                    runtimeErrorOccurred = false;
                }
            }
        }
    }
}


// Test to check whether the various error handling methods are correctly implemented
BOOST_AUTO_TEST_CASE( test_lagrange_error_checks )
{
    std::map< double, double > dataMap;
    std::vector< double > dataVector;
    std::map< int, double > coefficients;
    std::vector< double > independentVariableVector = getIndependentVariableVector( );

    bool runtimeErrorOccurred;

    {
        // Create interpolator with empty data map
        runtimeErrorOccurred = false;
        try
        {
            interpolators::LagrangeInterpolator< double, double > interpolator =
                    interpolators::LagrangeInterpolator< double, double >(
                        dataMap, 8, interpolators::huntingAlgorithm,
                        interpolators::lagrange_no_boundary_interpolation );
        }
        catch( std::runtime_error const& )
        {
            runtimeErrorOccurred = true;
        }
        BOOST_CHECK_EQUAL( runtimeErrorOccurred, 1 );
        runtimeErrorOccurred = false;

        // Create interpolator with empty independent variable vector
        runtimeErrorOccurred = false;
        try
        {
            dataVector.push_back( 1.0 );
            interpolators::LagrangeInterpolator< double, double > interpolator =
                    interpolators::LagrangeInterpolator< double, double >(
                        independentVariableVector, dataVector, 8,
                        interpolators::huntingAlgorithm,
                        interpolators::lagrange_no_boundary_interpolation );
        }
        catch( std::runtime_error const& )
        {
            runtimeErrorOccurred = true;
        }
        BOOST_CHECK_EQUAL( runtimeErrorOccurred, 1 );
        runtimeErrorOccurred = false;
        dataVector.clear( );

        // Create interpolator with empty dependent variable vector
        bool runtimeErrorOccurred = false;
        try
        {
            independentVariableVector.push_back( 1.0 );
            interpolators::LagrangeInterpolator< double, double > interpolator =
                    interpolators::LagrangeInterpolator< double, double >(
                        independentVariableVector, dataVector, 8,
                        interpolators::huntingAlgorithm,
                        interpolators::lagrange_no_boundary_interpolation );
        }
        catch( std::runtime_error const& )
        {
            runtimeErrorOccurred = true;
        }
        BOOST_CHECK_EQUAL( runtimeErrorOccurred, 1 );
        runtimeErrorOccurred = false;
        independentVariableVector = getIndependentVariableVector( );
    }

    // Create interpolator with NaN first entry (cannot make zero value)
    {
        // Create interpolator with NaN first entry from map constructor
        runtimeErrorOccurred = false;
        try
        {
            coefficients = getPolynomialCoefficients( 7 );
            dataMap[ independentVariableVector.at( 0 ) ] = TUDAT_NAN;
            for( unsigned int i = 1; i < independentVariableVector.size( ); i++ )
            {
                dataMap[ independentVariableVector.at( i ) ] = evaluatePolynomial(
                            coefficients, independentVariableVector.at( i ) );
            }
            interpolators::LagrangeInterpolator< double, double > interpolator =
                    interpolators::LagrangeInterpolator< double, double >(
                        dataMap, 8, interpolators::huntingAlgorithm,
                        interpolators::lagrange_no_boundary_interpolation );
        }
        catch( std::runtime_error const& )
        {
            runtimeErrorOccurred = true;
        }
        BOOST_CHECK_EQUAL( runtimeErrorOccurred, 1 );
        runtimeErrorOccurred = false;
        coefficients.clear( );
        dataMap.clear( );

        // Create interpolator with NaN first entry from vectors constructor
        runtimeErrorOccurred = false;
        try
        {
            coefficients = getPolynomialCoefficients( 7 );
            dataMap[ independentVariableVector.at( 0 ) ] = TUDAT_NAN;
            for( unsigned int i = 1; i < independentVariableVector.size( ); i++ )
            {
                dataVector.push_back( evaluatePolynomial(
                                          coefficients, independentVariableVector.at( i ) ) );
            }
            interpolators::LagrangeInterpolator< double, double > interpolator =
                    interpolators::LagrangeInterpolator< double, double >(
                        independentVariableVector, dataVector, 8,
                        interpolators::huntingAlgorithm,
                        interpolators::lagrange_no_boundary_interpolation );
        }
        catch( std::runtime_error const& )
        {
            runtimeErrorOccurred = true;
        }
        BOOST_CHECK_EQUAL( runtimeErrorOccurred, 1 );
        runtimeErrorOccurred = false;
        coefficients.clear( );
        dataVector.clear( );
    }

    // Test error throwing when making Lagrange Interpolator with odd number of stages
    {
        // Test for a range of odd number of stages
        for( unsigned int numberOfStages = 1; numberOfStages < 12; numberOfStages+= 2 )
        {
            // Create interpolator with odd number of stages for map constructor
            runtimeErrorOccurred = false;
            try
            {
                coefficients = getPolynomialCoefficients( numberOfStages - 1 );
                for( unsigned int i = 0; i < independentVariableVector.size( ); i++ )
                {
                    dataMap[ independentVariableVector.at( i ) ] = evaluatePolynomial(
                                coefficients, independentVariableVector.at( i ) );
                }
                interpolators::LagrangeInterpolator< double, double > interpolator =
                        interpolators::LagrangeInterpolator< double, double >(
                            dataMap, numberOfStages, interpolators::huntingAlgorithm,
                            interpolators::lagrange_no_boundary_interpolation );
            }
            catch( std::runtime_error const& )
            {
                runtimeErrorOccurred = true;
            }
            BOOST_CHECK_EQUAL( runtimeErrorOccurred, 1 );
            runtimeErrorOccurred = false;
            coefficients.clear( );
            dataMap.clear( );

            // Create interpolator with odd number of stages for vectors constructor
            runtimeErrorOccurred = false;
            try
            {
                coefficients = getPolynomialCoefficients( numberOfStages - 1 );
                for( unsigned int i = 0; i < independentVariableVector.size( ); i++ )
                {
                    dataVector.push_back( evaluatePolynomial(
                                              coefficients, independentVariableVector.at( i ) ) );
                }
                interpolators::LagrangeInterpolator< double, double > interpolator =
                        interpolators::LagrangeInterpolator< double, double >(
                            independentVariableVector, dataVector, numberOfStages,
                            interpolators::huntingAlgorithm,
                            interpolators::lagrange_no_boundary_interpolation );
            }
            catch( std::runtime_error const& )
            {
                runtimeErrorOccurred = true;
            }
            BOOST_CHECK_EQUAL( runtimeErrorOccurred, 1 );
            runtimeErrorOccurred = false;
            coefficients.clear( );
            dataVector.clear( );
        }
    }

    // Test error throwing when making Lagrange Interpolator from vectors constructor with
    // differently sized (in)dependent variable vectors
    {
        int numberOfStages = 8;
        runtimeErrorOccurred = false;
        try
        {
            coefficients = getPolynomialCoefficients( numberOfStages - 1 );
            for( unsigned int i = 0; i < independentVariableVector.size( ) - 1; i++ )
            {
                dataVector.push_back( evaluatePolynomial(
                                          coefficients, independentVariableVector.at( i ) ) );
            }
            interpolators::LagrangeInterpolator< double, double > interpolator =
                    interpolators::LagrangeInterpolator< double, double >(
                        independentVariableVector, dataVector, numberOfStages,
                        interpolators::huntingAlgorithm,
                        interpolators::lagrange_no_boundary_interpolation );
        }
        catch( std::runtime_error const& )
        {
            runtimeErrorOccurred = true;
        }
        BOOST_CHECK_EQUAL( runtimeErrorOccurred, 1 );
        runtimeErrorOccurred = false;
        coefficients.clear( );
        dataVector.clear( );
    }
}



BOOST_AUTO_TEST_SUITE_END( )

}

}
