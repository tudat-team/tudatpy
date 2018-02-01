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

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/streamFilters.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/aerodynamicCoefficientReader.h"

namespace tudat
{
namespace unit_tests
{

void reduceTestCoefficients( Eigen::Vector3d& testCoefficients, const int testCase )
{
    switch( testCase )
    {
    case 0:
        break;
    case 1:
        testCoefficients( 0 ) = 0.0;
        break;
    case 2:
        testCoefficients( 1 ) = 0.0;
        break;
    case 3:
        testCoefficients( 0 ) = 0.0;
        testCoefficients( 2 ) = 0.0;
        break;
    default:
        break;
    }
}

BOOST_AUTO_TEST_SUITE( test_multi_array_reader )

// Test if multi-array file reader is working correctly
BOOST_AUTO_TEST_CASE( testMultiArrayReader )
{
    // Test functionality of 3-dimensional multi-array reader
    for( unsigned testCase = 0; testCase < 4; testCase++ )
    {
        std::string fileName = tudat::input_output::getTudatRootPath( )
                + "/Astrodynamics/Aerodynamics/UnitTests/dCDwTest.txt";

        std::map< int, std::string > files;
        if( ( testCase == 0 ) || ( testCase == 2 ) )
        {
            files[ 0 ] = tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Aerodynamics/UnitTests/aurora_CD.txt";
        }
        if( testCase != 2 )
        {
            files[ 1 ] = tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Aerodynamics/UnitTests/aurora_Cm.txt";
        }
        if( testCase != 3 )
        {
            files[ 2 ] = tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Aerodynamics/UnitTests/aurora_CL.txt";
        }

        std::pair< boost::multi_array< Eigen::Vector3d, 2 >, std::vector< std::vector< double > > > outputArray =
                input_output::readAerodynamicCoefficients< 2 >( files );

        // Check block size.
        BOOST_CHECK_EQUAL( outputArray.first.shape( )[ 0 ], 26 );
        BOOST_CHECK_EQUAL( outputArray.second.at( 0 ).size( ), 26 );

        BOOST_CHECK_EQUAL( outputArray.first.shape( )[ 1 ], 21 );
        BOOST_CHECK_EQUAL( outputArray.second.at( 1 ).size( ), 21 );

        // Check coefficients against selected manually read coefficients.
        Eigen::Vector3d testCoefficients;
        testCoefficients << 0.0545832, -0.0318930, -0.3213905;
        reduceTestCoefficients( testCoefficients, testCase );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testCoefficients, outputArray.first[ 0 ][ 0 ],
                std::numeric_limits< double >::epsilon( ) );

        testCoefficients << 0.0778386, -0.0348528,  -0.3352446;
        reduceTestCoefficients( testCoefficients, testCase );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testCoefficients, outputArray.first[ 3 ][ 0 ],
                std::numeric_limits< double >::epsilon( ) );

        testCoefficients << 0.0263660, -0.0093861,  -0.0625756;
        reduceTestCoefficients( testCoefficients, testCase );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testCoefficients, outputArray.first[ 3 ][ 4 ],
                std::numeric_limits< double >::epsilon( ) );

        // Check values of first independent variable
        for( unsigned int i = 0; i < 26; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( outputArray.second.at( 0 ).at( i ), static_cast< double >( i ),
                                        std::numeric_limits< double >::epsilon( ) );
        }

        // Check selected of second independent variable
        BOOST_CHECK_CLOSE_FRACTION( outputArray.second.at( 1 ).at( 0 ), -0.174532925199433,
                                    std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( outputArray.second.at( 1 ).at( 5 ), 0.0,
                                    std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( outputArray.second.at( 1 ).at( 8 ), 0.104719755119660,
                                    std::numeric_limits< double >::epsilon( ) );

    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat


