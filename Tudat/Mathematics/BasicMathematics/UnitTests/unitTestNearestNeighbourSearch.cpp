/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120207    K. Kumar          File created.
 *      160930    M. Van den Broeck Added unit test for int computeNearestNeighborUsingBinarySearch
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <map>

#include <boost/filesystem.hpp>
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
        Eigen::VectorXd vectorOfTargetValues( 8 );
        vectorOfTargetValues << 1.1, 2.74, 2.75, 2.76, 4.6, 10.5, 54.55, 181.63;

        // Declare vector of expected indices.
        Eigen::VectorXi vectorOfExpectedIndices( 8 );
        vectorOfExpectedIndices << 0, 0, 0, 1, 1, 2, 4, 8;

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
