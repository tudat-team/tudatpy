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

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/streamFilters.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/multiDimensionalArrayReader.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_multi_array_reader )

// Test if multi-array file reader is working correctly
BOOST_AUTO_TEST_CASE( testMultiArrayReader )
{
    // Test functionality of 3-dimensional multi-array reader
    {
        std::string fileName = tudat::input_output::getTudatRootPath( )
                + "/Astrodynamics/Aerodynamics/UnitTests/dCDwTest.txt";

        for( unsigned int i = 0; i < 2; i++)
        {
            boost::multi_array< double, 3 > multiArrayFromFile =
                    tudat::input_output::MultiArrayFileReader< 3 >::readMultiArray( fileName );

            // Read only multi-array from file
            if( i == 0 )
            {
                multiArrayFromFile = tudat::input_output::MultiArrayFileReader< 3 >::readMultiArray( fileName );
            }
            // Read multi-array and independent variable values
            else
            {
                std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > > fileContents =
                        tudat::input_output::MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables( fileName );
                multiArrayFromFile = fileContents.first;
                std::vector< std::vector< double > > independentVariables = fileContents.second;

                // Test independent variable sizes
                BOOST_CHECK_EQUAL( independentVariables.size( ), 3 );
                BOOST_CHECK_EQUAL( independentVariables.at( 0 ).size( ), 11 );
                BOOST_CHECK_EQUAL( independentVariables.at( 1 ).size( ), 9 );
                BOOST_CHECK_EQUAL( independentVariables.at( 2 ).size( ), 5 );


            }

            // Test multi-array size
            BOOST_CHECK_EQUAL( multiArrayFromFile.shape( )[ 0 ], 11 );
            BOOST_CHECK_EQUAL( multiArrayFromFile.shape( )[ 1 ], 9 );
            BOOST_CHECK_EQUAL( multiArrayFromFile.shape( )[ 2 ], 5 );

            // Test selected multi-array values
            BOOST_CHECK_SMALL( multiArrayFromFile[ 1 ][ 3 ][ 1 ], std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( std::fabs( multiArrayFromFile[ 2 ][ 4 ][ 1 ] + 0.002 ),
                    std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( std::fabs( multiArrayFromFile[ 2 ][ 6 ][ 1 ] + 0.008 ),
                    std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( std::fabs( multiArrayFromFile[ 1 ][ 3 ][ 3 ] - 0.0028 ),
                    std::numeric_limits< double >::epsilon( ) );
        }
    }

    // Test functionality of 2-dimensional multi-array reader
    {
        std::string fileName = tudat::input_output::getTudatRootPath( )
                + "Astrodynamics/Propulsion/UnitTests/Isp_test.txt";

        for( unsigned int i = 0; i < 2; i++)
        {
            boost::multi_array< double, 2 > multiArrayFromFile =
                    tudat::input_output::MultiArrayFileReader< 2 >::readMultiArray( fileName );

            // Read only multi-array from file
            if( i == 0 )
            {
                multiArrayFromFile = tudat::input_output::MultiArrayFileReader< 2 >::readMultiArray( fileName );
            }
            // Read multi-array and independent variable values
            else
            {
                std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > fileContents =
                        tudat::input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables( fileName );
                multiArrayFromFile = fileContents.first;
                std::vector< std::vector< double > > independentVariables = fileContents.second;

                // Test independent variable sizes
                BOOST_CHECK_EQUAL( independentVariables.size( ), 2 );
                BOOST_CHECK_EQUAL( independentVariables.at( 0 ).size( ), 17 );
                BOOST_CHECK_EQUAL( independentVariables.at( 1 ).size( ), 5 );


            }

            // Test multi-array size
            BOOST_CHECK_EQUAL( multiArrayFromFile.shape( )[ 0 ], 17 );
            BOOST_CHECK_EQUAL( multiArrayFromFile.shape( )[ 1 ], 5 );

            // Test selected multi-array values
            BOOST_CHECK_SMALL( std::fabs( multiArrayFromFile[ 0 ][ 0 ] - 1900.0 ),
                    std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( std::fabs( multiArrayFromFile[ 2 ][ 4 ] - 2750.0 ),
                    std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( std::fabs( multiArrayFromFile[ 5 ][ 2 ] - 3100.0 ),
                    std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( std::fabs( multiArrayFromFile[ 16 ][ 0 ] - 590.0 ),
                    std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( std::fabs( multiArrayFromFile[ 16 ][ 4 ] - 888.0 ),
                    std::numeric_limits< double >::epsilon( ) );

        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

