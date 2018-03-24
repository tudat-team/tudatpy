/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/multi_array.hpp>

#include <limits>
#include <vector>
#include <cmath>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"

#include "Tudat/Mathematics/Interpolators/multiLinearInterpolator.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_multi_linear_interpolation )

// Test 1: Comparison to MATLAB solution of the example provided in matlab's interp2 function
// documentation.
BOOST_AUTO_TEST_CASE( test2Dimensions )
{
    // Create independent variable vector.
    std::vector< std::vector< double > > independentValues;

    // Resize independent variable vector to two dimensions and assign values.
    independentValues.resize( 2 );
    for ( int i = 0; i < 5; i++ )
    {
        independentValues[ 0 ].push_back( 1950.0 + static_cast< double >( i ) * 10.0 );
    }

    for ( int i = 0; i < 3; i++ )
    {
        independentValues[ 1 ].push_back( 10.0 + static_cast< double >( i ) * 10.0 );
    }

    // Create two-dimensional array for dependent values.
    boost::multi_array< double, 2 > dependentValues;
    dependentValues.resize( boost::extents[ 5 ][ 3 ] );

    dependentValues[ 0 ][ 0 ] = 150.697;
    dependentValues[ 0 ][ 1 ] = 199.592;
    dependentValues[ 0 ][ 2 ] = 187.625;
    dependentValues[ 1 ][ 0 ] = 179.323;
    dependentValues[ 1 ][ 1 ] = 195.072;
    dependentValues[ 1 ][ 2 ] = 250.287;
    dependentValues[ 2 ][ 0 ] = 203.212;
    dependentValues[ 2 ][ 1 ] = 179.092;
    dependentValues[ 2 ][ 2 ] = 322.767;
    dependentValues[ 3 ][ 0 ] = 226.505;
    dependentValues[ 3 ][ 1 ] = 153.706;
    dependentValues[ 3 ][ 2 ] = 426.730;
    dependentValues[ 4 ][ 0 ] = 249.633;
    dependentValues[ 4 ][ 1 ] = 120.281;
    dependentValues[ 4 ][ 2 ] = 598.243;

    // Initialize interpolator.
    interpolators::MultiLinearInterpolator< double, double, 2 > twoDimensionalInterpolator(
            independentValues, dependentValues );

    // Set interpolation target's independent values.
    std::vector< double > targetValue( 2 );
    targetValue[ 0 ] = 1975.0;
    targetValue[ 1 ] = 15.0;

    // Declare interpolated dependent variable value and execute interpolation.
    const double interpolationResult = twoDimensionalInterpolator.interpolate( targetValue );

    // Check if equal to MATLAB result up to machine precision.
    BOOST_CHECK_CLOSE_FRACTION( 190.62875, interpolationResult,
                                std::numeric_limits< double >::epsilon( ) );
}

// Test 2: 4-dimensional test. Comparison to matlab interpolated solution of analytical function.
BOOST_AUTO_TEST_CASE( test4Dimensions )
{
    // Create independent variable vector.
    std::vector< std::vector< double > > independentValues;

    // Resize independent variable vector to two dimensions and assign values.
    independentValues.resize( 4 );
    for ( int i = 0; i < 11; i++ )
    {
        independentValues[ 0 ].push_back( -1.0 + static_cast< double >( i ) * 0.2 );
        independentValues[ 1 ].push_back( -1.0 + static_cast< double >( i ) * 0.2 );
        independentValues[ 2 ].push_back( -1.0 + static_cast< double >( i ) * 0.2 );
    }

    for ( int i = 0; i < 6; i++ )
    {
        independentValues[ 3 ].push_back( static_cast< double >( i ) * 2.0 );
    }

    // Create four-dimensional array for dependent values based on analytical function
    // f = t*e^(-x^2 - y^2 - z^2 ).
    boost::multi_array< double, 4 > dependentValues;
    dependentValues.resize( boost::extents[ 11 ][ 11 ][ 11 ][ 6 ] );

    for ( int i = 0; i < 11; i++ )
    {
        for ( int j = 0; j < 11; j++ )
        {
            for ( int k = 0; k < 11; k++ )
            {
                for ( int l = 0; l < 6; l++ )
                {
                    dependentValues[ i ][ j ][ k ][ l ] =
                            independentValues[ 3 ][ l ] *
                            std::exp( -std::pow( independentValues[ 0 ][ i ], 2 ) -
                                      std::pow( independentValues[ 1 ][ j ], 2 ) -
                                      std::pow( independentValues[ 2 ][ k ], 2 ) );
                }
            }
        }
    }

    // Initialize interpolator.
    interpolators::MultiLinearInterpolator< double, double, 4 > twoDimensionalInterpolator(
            independentValues, dependentValues );

    // Set interpolation target's independent values.
    std::vector< double > targetValue( 4 );
    targetValue[ 0 ] = -1.0;
    targetValue[ 1 ] = 0.1;
    targetValue[ 2 ] = 0.5;
    targetValue[ 3 ] = 7.0;

    // Declare interpolated dependent variable value and execute interpolation.
    const double interpolationResult = twoDimensionalInterpolator.interpolate( targetValue );

    // Check if equal to MATLAB result up to machine precision
    BOOST_CHECK_CLOSE_FRACTION( 1.956391733957447, interpolationResult,
                                std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
