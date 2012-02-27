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
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// Test runs code and verifies result against expected value.
// If the tested code is erroneous, the test function returns a boolean
// true; if the code is correct, the function returns a boolean false.
// 

#include <Eigen/Core>
#include "Tudat/Astrodynamics/States/state.h"

//! Test orbital element conversion code.
int main( )
{
    // Using declarations.
    using tudat::State;

    // Initialize unit test result to false.
    bool isStateTestErroneous = false;

    // Test the constructors.
    // Test default constructor.
    State uninitializedState;
    if ( uninitializedState.state.rows( ) != 0 )
    {
        std::cerr << "Uninitialized state vector length was " << uninitializedState.state.rows( )
                  << "; expected 0" << std::endl;
        isStateTestErroneous = true;
    }

    // Test constructor to initialize with a given Eigen VectorXd state.
    Eigen::VectorXd stateVector = Eigen::VectorXd::Random( 10 );
    State randomState( stateVector );
    if ( randomState.state.rows( ) != 10 )
    {
        std::cerr << "Random state vector length was " << randomState.state.rows( )
                  << "; expected 10" << std::endl;
        isStateTestErroneous = true;
    }

    if ( randomState.state != stateVector )
    {
        std::cerr << "Random state vector was not equal to the one passed to the constructor!"
                  << std::endl;
        isStateTestErroneous = true;
    }

    // Test constructor to initialize length of state.
    State lengthInitializedState( 13 );
    if ( lengthInitializedState.state.rows( ) != 13 )
    {
        std::cerr << "Length initialized state vector length was "
                  << lengthInitializedState.state.rows( ) << "; expected 13." << std::endl;
        isStateTestErroneous = true;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isStateTestErroneous )
    {
        std::cerr << "testState failed!" << std::endl;
    }

    return isStateTestErroneous;
}
