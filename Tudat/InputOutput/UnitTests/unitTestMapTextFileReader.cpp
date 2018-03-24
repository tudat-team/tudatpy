/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/mapTextFileReader.h"

namespace tudat
{

namespace unit_tests
{

// Test if matrix text file reader is working correctly.
BOOST_AUTO_TEST_CASE( testMatrixTextFileReader )
{

    // Test 1
    {
        // Set expected matrix.
        std::map< int, std::vector< int > > expectedMap = {
            { 1, { 1, 2 } },
            { 2, { 1, 2, 4 } },
            { 3, { 1, 2, 4, 8 } },
            { 4, { 1, 2, 4, 8, 16 } }
        };

        // Read input file and store data in matrix.
        auto inputFileMap = input_output::readStlVectorMapFromFile< int, int >(
                    input_output::getTudatRootPath( ) + "/InputOutput/UnitTests/testMap1.txt" );

        bool allEqual = true;
        // Check if data input file matrix matches expected matrix.
        try
        {
            for ( auto ent: inputFileMap )
            {
                auto key = ent.first;
                auto values = ent.second;
                for ( unsigned int i = 0; i < values.size(); i++ )
                {
                    if ( inputFileMap[ key ][ i ] != expectedMap[ key ][ i ] )
                    {
                        std::cout << "Element at " << i << " for key " << key << " ("
                                  << inputFileMap[ key ][ i ] << ") does not match expected value of "
                                  << "(" << expectedMap[ key ][ i ] << ")." << std::endl;
                        allEqual = false;
                        break;
                    }
                }
                if ( ! allEqual )
                {
                    break;
                }
            }
        }
        catch( std::runtime_error &inconsistentSizes )
        {
            std::cout << "Maps have different keys or inconsistent sizes." << std::endl;
            allEqual = false;
        }

        BOOST_CHECK( allEqual );

    }


    // Test 2
    {
        // Set expected matrix.
        std::map< double, std::vector< float > > expectedMap = {
            { 0.0, {  0.5, 0.8 } },
            { 0.1, { -0.5, 1.0 } },
            { 0.2, {  0.5 } },
            { 0.3, {  1e6, 1.0, 1.0 } }
        };

        // Read input file and store data in matrix.
        auto inputFileMap = input_output::readStlVectorMapFromFile< double, float >(
                    input_output::getTudatRootPath( ) + "/InputOutput/UnitTests/testMap2.txt" );

        bool allEqual = true;
        // Check if data input file matrix matches expected matrix.
        try
        {
            for ( auto ent: inputFileMap )
            {
                auto key = ent.first;
                auto values = ent.second;
                for ( unsigned int i = 0; i < values.size(); i++ )
                {
                    if ( inputFileMap[ key ][ i ] != expectedMap[ key ][ i ] )
                    {
                        std::cout << "Element at " << i << " for key " << key << " ("
                                  << inputFileMap[ key ][ i ] << ") does not match expected value of "
                                  << "(" << expectedMap[ key ][ i ] << ")." << std::endl;
                        allEqual = false;
                        break;
                    }
                }
                if ( ! allEqual )
                {
                    break;
                }
            }
        }
        catch( std::runtime_error &inconsistentSizes )
        {
            std::cout << "Maps have different keys or inconsistent sizes." << std::endl;
            allEqual = false;
        }

        BOOST_CHECK( allEqual );

    }


    // Test 1
    {
        // Set expected matrix.
        std::map< std::string, std::vector< std::string > > expectedMap = {
            { "x", { "one", "two" } },
            { "y", { "three", "four", "five" } },
            { "z", { "six" } }
        };

        // Read input file and store data in matrix.
        auto inputFileMap = input_output::readStlVectorMapFromFile< std::string, std::string >(
                    input_output::getTudatRootPath( ) + "/InputOutput/UnitTests/testMap3.txt" );

        bool allEqual = true;
        // Check if data input file matrix matches expected matrix.
        try
        {
            for ( auto ent: inputFileMap )
            {
                auto key = ent.first;
                auto values = ent.second;
                for ( unsigned int i = 0; i < values.size(); i++ )
                {
                    if ( inputFileMap[ key ][ i ] != expectedMap[ key ][ i ] )
                    {
                        std::cout << "Element at " << i << " for key " << key << " ("
                                  << inputFileMap[ key ][ i ] << ") does not match expected value of "
                                  << "(" << expectedMap[ key ][ i ] << ")." << std::endl;
                        allEqual = false;
                        break;
                    }
                }
                if ( ! allEqual )
                {
                    break;
                }
            }
        }
        catch( std::runtime_error &inconsistentSizes )
        {
            std::cout << "Maps have different keys or inconsistent sizes." << std::endl;
            allEqual = false;
        }

        BOOST_CHECK( allEqual );

    }



}

} // namespace unit_tests

} // namespace tudat
