/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      120207    K. Kumar          File created.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
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
                        tudat::mathematics::computeNearestLeftNeighborUsingBinarySearch(
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
                        tudat::mathematics::computeNearestLeftNeighborUsingBinarySearch(
                            mapOfSortedData, vectorOfTargetValues[ i ] ) );
        }
    }
}

//! Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
