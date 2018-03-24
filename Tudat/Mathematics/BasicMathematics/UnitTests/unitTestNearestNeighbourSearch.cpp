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
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *
 */

#define BOOST_TEST_MAIN

#include <map>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"

namespace tudat
{
namespace unit_tests
{

//! Define Boost test suite.
BOOST_AUTO_TEST_SUITE( test_basic_functions )

//! Test if search for nearest left neighbor using binary search works correctly.
BOOST_AUTO_TEST_CASE( testNearestLeftNeighborUsingBinarySearch )
{
    using namespace basic_mathematics;

    // Case 1: test Eigen-interface.
    {
        // Populate vector of 10 sorted elements.
        Eigen::VectorXd vectorOfSortedData( 10 );
        vectorOfSortedData << 1.0, 4.5, 10.6, 14.98, 54.65, 88.9, 101.31, 144.63, 180.01, 201.94;

        // Declare vector of target values.
        Eigen::VectorXd vectorOfTargetValues( 5 );
        vectorOfTargetValues << 1.1, 4.6, 10.5, 54.55, 181.63;

        // Declare vector of expected indices.
        Eigen::VectorXi vectorOfExpectedIndices( 5 );
        vectorOfExpectedIndices << 0, 1, 1, 3, 8;

        // Compute nearest left neighbors and check if they match expectations.
        for ( int i = 0; i < vectorOfTargetValues.rows( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        vectorOfExpectedIndices[ i ],
                        computeNearestLeftNeighborUsingBinarySearch(
                            vectorOfSortedData, vectorOfTargetValues[ i ] ) );
        }
    }

    // Case 2: test map-interface with VectorXd.
    {
        // Populate map of 10 sorted elements.
        std::map< double, Eigen::VectorXd > mapOfSortedData;

        Eigen::VectorXd vectorOfData( 1 );
        vectorOfData << 1.0;

        mapOfSortedData[ 0.3 ] = vectorOfData;
        mapOfSortedData[ 3.65 ] = vectorOfData;
        mapOfSortedData[ 43.12 ] = vectorOfData;
        mapOfSortedData[ 2.23 ] = vectorOfData;
        mapOfSortedData[ 1.233 ] = vectorOfData;
        mapOfSortedData[ 6.78 ] = vectorOfData;
        mapOfSortedData[ 0.21 ] = vectorOfData;
        mapOfSortedData[ -1.23 ] = vectorOfData;
        mapOfSortedData[ -931.12 ] = vectorOfData;
        mapOfSortedData[ 124.52 ] = vectorOfData;

        // Declare vector of target values.
        Eigen::VectorXd vectorOfTargetValues( 5 );
        vectorOfTargetValues << -1.22, 3.66, -931.11, 43.12, 0.4;

        // Declare vector of expected indices.
        Eigen::VectorXi vectorOfExpectedIndices( 5 );
        vectorOfExpectedIndices << 1, 6, 0, 8, 3;

        // Compute nearest left neighbors and check if they match expectations.
        for ( int i = 0; i < vectorOfTargetValues.rows( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        vectorOfExpectedIndices[ i ],
                        computeNearestLeftNeighborUsingBinarySearch(
                            mapOfSortedData, vectorOfTargetValues[ i ] ) );
        }
    }

    // Case 3: test templated STL vector-interface.
    {
        // Populate vector of 10 sorted elements.
        std::vector< double > vectorOfSortedData( 10 );
        vectorOfSortedData[ 0 ] = 1.0;
        vectorOfSortedData[ 1 ] = 4.5;
        vectorOfSortedData[ 2 ] = 10.6;
        vectorOfSortedData[ 3 ] = 14.98;
        vectorOfSortedData[ 4 ] = 54.65;
        vectorOfSortedData[ 5 ] = 88.9;
        vectorOfSortedData[ 6 ] = 101.31;
        vectorOfSortedData[ 7 ] = 144.63;
        vectorOfSortedData[ 8 ] = 180.01;
        vectorOfSortedData[ 9 ] = 201.94;

        // Declare vector of target values.
        std::vector< double > vectorOfTargetValues( 5 );
        vectorOfTargetValues[ 0 ] = 1.1;
        vectorOfTargetValues[ 1 ] = 4.6;
        vectorOfTargetValues[ 2 ] = 10.5;
        vectorOfTargetValues[ 3 ] = 54.55;
        vectorOfTargetValues[ 4 ] = 181.63;

        // Declare vector of expected indices.
        std::vector< double > vectorOfExpectedIndices( 5 );
        vectorOfExpectedIndices[ 0 ] = 0;
        vectorOfExpectedIndices[ 1 ] = 1;
        vectorOfExpectedIndices[ 2 ] = 1;
        vectorOfExpectedIndices[ 3 ] = 3;
        vectorOfExpectedIndices[ 4 ] = 8;

        // Compute nearest left neighbors and check if they match expectations.
        for ( int i = 0; i < 5; i++ )
        {
            BOOST_CHECK_EQUAL( vectorOfExpectedIndices[ i ],
             computeNearestLeftNeighborUsingBinarySearch< double >(
                            vectorOfSortedData, vectorOfTargetValues[ i ] ) );
        }
    }

    // Case 4: test templated hunting algorithm ( with STL vector-interface ).
    {
        // Populate vector of 10 sorted elements.
        std::vector< double > vectorOfSortedData( 10 );
        vectorOfSortedData[ 0 ] = 1.0;
        vectorOfSortedData[ 1 ] = 4.5;
        vectorOfSortedData[ 2 ] = 10.6;
        vectorOfSortedData[ 3 ] = 14.98;
        vectorOfSortedData[ 4 ] = 54.65;
        vectorOfSortedData[ 5 ] = 88.9;
        vectorOfSortedData[ 6 ] = 101.31;
        vectorOfSortedData[ 7 ] = 144.63;
        vectorOfSortedData[ 8 ] = 180.01;
        vectorOfSortedData[ 9 ] = 201.94;

        // Declare vector of target values.
        std::vector< double > vectorOfTargetValues( 5 );
        vectorOfTargetValues[ 0 ] = 1.1;
        vectorOfTargetValues[ 1 ] = 4.6;
        vectorOfTargetValues[ 2 ] = 10.5;
        vectorOfTargetValues[ 3 ] = 54.55;
        vectorOfTargetValues[ 4 ] = 181.63;

        // Declare vector of expected indices.
        std::vector< double > vectorOfExpectedIndices( 5 );
        vectorOfExpectedIndices[ 0 ] = 0;
        vectorOfExpectedIndices[ 1 ] = 1;
        vectorOfExpectedIndices[ 2 ] = 1;
        vectorOfExpectedIndices[ 3 ] = 3;
        vectorOfExpectedIndices[ 4 ] = 8;

        // Compute nearest left neighbors and check if they match expectations.
        for ( int i = 0; i < 5; i++ )
        {
            // Check whether each initial guess yields correct result.
            for( int j = 0; j < 9; j++ )
            {
                BOOST_CHECK_EQUAL( vectorOfExpectedIndices[ i ],
                    findNearestLeftNeighbourUsingHuntingAlgorithm< double >(
                    vectorOfTargetValues[ i ], j, vectorOfSortedData ) );
            }
        }
    }

    // Case 5: test Eigen-interface of NearestNeighbourSearch.
    {
        // Populate vector of 10 sorted elements.
        Eigen::VectorXd vectorOfSortedData( 10 );
        vectorOfSortedData << 1.0, 4.5, 10.6, 14.98, 54.65, 88.9, 101.31, 144.63, 180.01, 201.94;

        // Declare vector of target values.
        Eigen::VectorXd vectorOfTargetValues( 10 );
        vectorOfTargetValues << 1.1, 2.74, 2.75, 2.76, 4.6, 10.5, 54.55, 181.63, 200.0, 205.0;

        // Declare vector of expected indices.
        Eigen::VectorXi vectorOfExpectedIndices( 10 );
        vectorOfExpectedIndices << 0, 0, 0, 1, 1, 2, 4, 8, 9, 9;

        // Compute nearest left neighbors and check if they match expectations.
        for ( int i = 0; i < vectorOfTargetValues.rows( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        vectorOfExpectedIndices[ i ],
                        computeNearestNeighborUsingBinarySearch(
                            vectorOfSortedData, vectorOfTargetValues[ i ] ) );
        }
    }

}

//! Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
