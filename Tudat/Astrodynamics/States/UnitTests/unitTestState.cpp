/*    Copyright (c) 2010 - 2012 Delft University of Technology.
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
 *      111123    B. Tong Minh      First creation of the code.
 *      111201    K. Kumar          Minor corrections; fixed error in initialize length unit test.
 *      120508    K. Kumar          Boostified unit test.
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <Eigen/Core>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/States/state.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_state )

//! Test if state constructor is working correctly.
BOOST_AUTO_TEST_CASE( testStateConstructor )
{
    // Test default constructor.
    {
        // Declare uninitialized state.
        astrodynamics::states::State uninitializedState;

        // Check that size of uninitiatialized state is zero.
        BOOST_CHECK_EQUAL( uninitializedState.state.rows( ), 0 );
    }

    // Test constructor that initializes state, given an Eigen::VectorXd vector.
    {
        // Set size of state.
        const double stateSize = 10;

        // Create vector.
        const Eigen::VectorXd stateVector = ( Eigen::VectorXd( stateSize )
                                              << 1.0, 3.0, 4.6, 9.0, -0.3,
                                              10.6, 9.4, 0.56, 8.4, 10.3 ).finished( );

        // Create state from vector.
        astrodynamics::states::State initializedState( stateVector );

        // Check that state is as expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateVector, initializedState.state,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test constructor that initializes size of state (elements set to zero).
    {
        // Set size of state.
        const double stateSize = 10;

        // Create state from vector.
        astrodynamics::states::State sizeInitializedState( stateSize );

        // Check that state is as expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( Eigen::VectorXd::Zero( stateSize ),
                                           sizeInitializedState.state,
                                           std::numeric_limits< double >::epsilon( ) );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
