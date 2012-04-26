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
 *      110110    K. Kumar          File created.
 *      110121    K. Kumar          Updated to comply with new protocol.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// Test runs code and verifies result against expected value.
// If the tested code is erroneous, the test function returns a boolean
// true; if the code is correct, the function returns a boolean false.
// 

#include <iostream>
#include <limits>
#include "Tudat/Astrodynamics/States/cartesianElements.h"

//! Test Cartesian elements state class.
int main( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;
    using namespace tudat;

    // Two tests.
    // Test 1: Set Cartesian elements using individual set functions; and get
    //         Cartesian elements using individual get functions.
    // Test 2: Set Cartesian elements using setState( ) function; and get
    //         Cartesian elements using getState( ) function.

    // Initialize unit test result to false.
    bool isCartesianElementsErroneous = false;

    // Create Cartesian elements state objects.
    using tudat::astrodynamics::states::CartesianElements;
    CartesianElements cartesianElementsStateTest1;
    CartesianElements cartesianElementsStateTest2;

    // Expected test result.
    // Create vector of Cartesian elements: x, y, z, xdot, ydot, zdot.
    Eigen::VectorXd vectorOfCartesianElements( 6 );
    vectorOfCartesianElements( 0 ) = 2.5e6;
    vectorOfCartesianElements( 1 ) = 0.5e6;
    vectorOfCartesianElements( 2 ) = 0.1e6;
    vectorOfCartesianElements( 3 ) = 125.0;
    vectorOfCartesianElements( 4 ) = 2000.0;
    vectorOfCartesianElements( 5 ) = 50.0;

    // Test 1: Set Cartesian elements in state object.
    cartesianElementsStateTest1.setCartesianElementX( vectorOfCartesianElements( 0 ) );
    cartesianElementsStateTest1.setCartesianElementY( vectorOfCartesianElements( 1 ) );
    cartesianElementsStateTest1.setCartesianElementZ( vectorOfCartesianElements( 2 ) );
    cartesianElementsStateTest1.setCartesianElementXDot( vectorOfCartesianElements( 3 ) );
    cartesianElementsStateTest1.setCartesianElementYDot( vectorOfCartesianElements( 4 ) );
    cartesianElementsStateTest1.setCartesianElementZDot( vectorOfCartesianElements( 5 ) );

    // Test 1: Get Cartesian elements and store in a state vector.
    Eigen::VectorXd cartesianElementsStateVectorTest1( 6 );
    cartesianElementsStateVectorTest1( 0 ) = cartesianElementsStateTest1.getCartesianElementX( );
    cartesianElementsStateVectorTest1( 1 ) = cartesianElementsStateTest1.getCartesianElementY( );
    cartesianElementsStateVectorTest1( 2 ) = cartesianElementsStateTest1.getCartesianElementZ( );
    cartesianElementsStateVectorTest1( 3 )
            = cartesianElementsStateTest1.getCartesianElementXDot( );
    cartesianElementsStateVectorTest1( 4 )
            = cartesianElementsStateTest1.getCartesianElementYDot( );
    cartesianElementsStateVectorTest1( 5 )
            = cartesianElementsStateTest1.getCartesianElementZDot( );

    // Test 1: Difference between setting each Cartesian element and the
    // expected values.
    Eigen::VectorXd differenceBetweenResultsTest1( 6 );
    differenceBetweenResultsTest1 = cartesianElementsStateVectorTest1 - vectorOfCartesianElements;

    // Test 2: Set Cartesian elements using setState( ) function.
    cartesianElementsStateTest2.state = vectorOfCartesianElements;

    // Test 2: Difference between setting the Cartesian state as a whole and
    // the expected values.
    Eigen::VectorXd differenceBetweenResultsTest2;
    differenceBetweenResultsTest2 = cartesianElementsStateTest2.state - vectorOfCartesianElements;

    // Set test result to true if the test does not match the expected result.

    if ( differenceBetweenResultsTest1.norm( ) > std::numeric_limits< double >::epsilon( ) )
    {
        isCartesianElementsErroneous = true;

        // Generate error statements.
        cerr << "The computed value " << endl;
        cerr << "( " << cartesianElementsStateVectorTest1 << " )" << endl;
        cerr << "using individual set/get functions for the "
             << "CartesianElements class does not match the expected solution " << endl;
        cerr << "( " << vectorOfCartesianElements << " )." << endl;
        cerr << "The difference is: "
             << vectorOfCartesianElements - cartesianElementsStateVectorTest1 << endl;
    }

    if ( differenceBetweenResultsTest2.norm( ) > std::numeric_limits< double >::epsilon( ) )
    {
        isCartesianElementsErroneous = true;

        // Generate error statements.
        cerr << "The computed value " << endl;
        cerr << "( " << cartesianElementsStateTest2.state << " )" << endl;
        cerr << "using set/get functions for the whole state for the "
             << "CartesianElements class does not match the expected solution " << endl;
        cerr << "( " << vectorOfCartesianElements << " )." << endl;
        cerr << "The difference is: "
             << vectorOfCartesianElements - cartesianElementsStateTest2.state << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isCartesianElementsErroneous )
    {
        cerr << "testCartesianElements failed!" << endl;
    }

    return isCartesianElementsErroneous;
}
