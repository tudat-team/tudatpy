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
 *      110110    K. Kumar          File created.
 *      110121    K. Kumar          Updated to comply with new protocol.
 *      120511    K. Kumar          Boostified unit test.
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/States/cartesianElements.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_cartesian_state )

//! Test if individual set- and get-functions are working correctly.
BOOST_AUTO_TEST_CASE( testCartesianStateIndividualSetFunctions )
{
    using astrodynamics::states::CartesianElements;
    using astrodynamics::states::xPositionIndex;
    using astrodynamics::states::yPositionIndex;
    using astrodynamics::states::zPositionIndex;
    using astrodynamics::states::xVelocityIndex;
    using astrodynamics::states::yVelocityIndex;
    using astrodynamics::states::zVelocityIndex;

    // Create Cartesian elements state.
    CartesianElements cartesianElementsState;

    // Create vector of Cartesian elements: x, y, z, xdot, ydot, zdot.
    Eigen::VectorXd cartesianElements( 6 );
    cartesianElements( xPositionIndex ) = 2.5e6;
    cartesianElements( yPositionIndex ) = 0.5e6;
    cartesianElements( zPositionIndex ) = 0.1e6;
    cartesianElements( xVelocityIndex ) = 125.0;
    cartesianElements( yVelocityIndex ) = 2000.0;
    cartesianElements( zVelocityIndex ) = 50.0;

    // Test 1: Set Cartesian elements individually in state object.
    cartesianElementsState.setCartesianElementX( cartesianElements( xPositionIndex ) );
    cartesianElementsState.setCartesianElementY( cartesianElements( yPositionIndex ) );
    cartesianElementsState.setCartesianElementZ( cartesianElements( zPositionIndex ) );
    cartesianElementsState.setCartesianElementXDot( cartesianElements( xVelocityIndex ) );
    cartesianElementsState.setCartesianElementYDot( cartesianElements( yVelocityIndex ) );
    cartesianElementsState.setCartesianElementZDot( cartesianElements( zVelocityIndex ) );

    // Check that result from using individual get-functions gives expected results.
    BOOST_CHECK_EQUAL( cartesianElements( xPositionIndex ),
                       cartesianElementsState.getCartesianElementX( ) );

    BOOST_CHECK_EQUAL( cartesianElements( yPositionIndex ),
                       cartesianElementsState.getCartesianElementY( ) );

    BOOST_CHECK_EQUAL( cartesianElements( zPositionIndex ),
                       cartesianElementsState.getCartesianElementZ( ) );

    BOOST_CHECK_EQUAL( cartesianElements( xVelocityIndex ),
                       cartesianElementsState.getCartesianElementXDot( ) );

    BOOST_CHECK_EQUAL( cartesianElements( yVelocityIndex ),
                       cartesianElementsState.getCartesianElementYDot( ) );

    BOOST_CHECK_EQUAL( cartesianElements( zVelocityIndex ),
                       cartesianElementsState.getCartesianElementZDot( ) );
}

//! Test if setting state vector is working correctly.
BOOST_AUTO_TEST_CASE( testCartesianStateSetFunction )
{
    using astrodynamics::states::CartesianElements;
    using astrodynamics::states::xPositionIndex;
    using astrodynamics::states::yPositionIndex;
    using astrodynamics::states::zPositionIndex;
    using astrodynamics::states::xVelocityIndex;
    using astrodynamics::states::yVelocityIndex;
    using astrodynamics::states::zVelocityIndex;

    // Create Cartesian elements state.
    CartesianElements cartesianElementsState;

    // Create vector of Cartesian elements: x, y, z, xdot, ydot, zdot.
    Eigen::VectorXd cartesianElements( 6 );
    cartesianElements( xPositionIndex ) = 2.5e6;
    cartesianElements( yPositionIndex ) = 0.5e6;
    cartesianElements( zPositionIndex ) = 0.1e6;
    cartesianElements( xVelocityIndex ) = 125.0;
    cartesianElements( yVelocityIndex ) = 2000.0;
    cartesianElements( zVelocityIndex ) = 50.0;

    // Test 2: Set Cartesian elements as vector in state object.
    cartesianElementsState.state = cartesianElements;

    // Check that result from using individual get-functions gives expected results.
    TUDAT_CHECK_MATRIX_BASE( cartesianElements, cartesianElementsState.state )
            BOOST_CHECK_EQUAL( cartesianElements, cartesianElementsState.state );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
