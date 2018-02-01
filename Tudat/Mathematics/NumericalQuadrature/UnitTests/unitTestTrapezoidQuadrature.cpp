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

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{


// Function to create linearly spaced data.
std::vector< double > linspace( double start, double end, int numberOfSamples )
{
    double spacing = ( end - start ) / ( static_cast< double >( numberOfSamples - 1 ) );

    std::vector< double > vector(0);
    for( int i = 0 ; i < numberOfSamples ; i++ )
    {
        vector.push_back( start + static_cast< double >( i ) * spacing );
    }

    return vector;
}

BOOST_AUTO_TEST_SUITE( test_trapezoid_integrator )

//! Test if qudrature is computed correctly (sine function with 1E4 data points).
BOOST_AUTO_TEST_CASE( testIntegralSineFunction )
{
    std::vector< int > numberOfSamplesList;
    numberOfSamplesList.push_back( 1E2 );
    numberOfSamplesList.push_back( 1E4 );
    numberOfSamplesList.push_back( 1E6 );

    std::vector< double > tolerances = { 2E-4, 2E-8, 2E-12 };

    double previousError = TUDAT_NAN;
    double currentError = TUDAT_NAN;

    for( unsigned int test = 0; test < 3; test++ )
    {
        // Generate independent variables
        int numberOfSamples = numberOfSamplesList.at( test );

        std::vector< double > bounds( 2 );
        bounds[ 0 ] = 0.0 ;
        bounds[ 1 ] = mathematical_constants::PI;
        std::vector< double > independentVariables = linspace( bounds[ 0 ], bounds[ 1 ], numberOfSamples );

        std::vector< double > dependentVariables;
        for( unsigned int i = 0 ; i < independentVariables.size( ) ; i++ )
        {
            dependentVariables.push_back( std::sin( independentVariables[ i ] ) );
        }

        tudat::numerical_quadrature::TrapezoidNumericalQuadrature< double, double > integrator(
                    independentVariables, dependentVariables );

        double expectedIntegral = 2.0;
        double computedIntegral = integrator.getQuadrature();

        // Check if computed sample mean matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( computedIntegral, expectedIntegral, tolerances.at( test ) );

        currentError = std::fabs( computedIntegral - expectedIntegral );

        // Test order of quadrature
        if( test > 0 )
        {
            BOOST_CHECK_CLOSE_FRACTION( previousError / currentError, 1.0E4, 0.1 );
        }

        previousError = currentError;

        // Reset data and recompute quadrature
        bounds[ 1 ] = 2.0 * mathematical_constants::PI;
        independentVariables = linspace( bounds[ 0 ], bounds[ 1 ], numberOfSamples );

        dependentVariables.clear( );
        for( unsigned int i = 0 ; i < independentVariables.size( ) ; i++ )
        {
            dependentVariables.push_back( std::sin( independentVariables[ i ] ) );
        }

        // Error should be close to zero, as integration (and its errors) are symmetrical
        integrator.resetData( independentVariables, dependentVariables );
        expectedIntegral = 0.0;
        computedIntegral = integrator.getQuadrature( );
        BOOST_CHECK_SMALL( computedIntegral - expectedIntegral, std::max( 5.0E-6, numberOfSamples * 5.0E-20 ) );
    }
}


//! Test if qudrature is computed correctly (exponential function).
BOOST_AUTO_TEST_CASE( testIntegralExpFunction )
{
    // Generate independent variables
    int numberOfSamples = 1E4;

    std::vector< double > bounds( 2 );
    bounds[ 0 ] = 0.0;
    bounds[ 1 ] = 2.0;
    std::vector< double > independentVariables = linspace( bounds[ 0 ], bounds[ 1 ], numberOfSamples );

    std::vector< double > dependentVariables(0);
    for( unsigned int i = 0 ; i < independentVariables.size() ; i++ )
    {
        dependentVariables.push_back( std::exp( independentVariables[ i ] ) );
    }

    // Create integrator
    tudat::numerical_quadrature::TrapezoidNumericalQuadrature< double, double > integrator(
                independentVariables, dependentVariables );

    double expectedIntegral = std::exp( 2.0 ) - std::exp( 0.0 );
    double computedIntegral = integrator.getQuadrature( );

    // Check if computed sample mean matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedIntegral, expectedIntegral, 1E-8 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
