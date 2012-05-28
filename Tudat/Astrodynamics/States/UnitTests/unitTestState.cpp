/*    Copyright (c) 2010-2012, Delft University of Technology
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
