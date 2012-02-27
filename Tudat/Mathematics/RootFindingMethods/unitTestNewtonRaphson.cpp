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
 *      110111    E. Iorfida        First creation of the code.
 *                                  The code is tested with the function: f(x)= x^2 - 3.
 *      110111    K. Kumar          Updated to use address of global  functions instead of
 *                                  pointers for set functions; aligned code as required for
 *                                  namespaces; minor comment changes.
 *      110119    K. Kumar          Updated code to work with adaptor and abstract base
 *                                  implementation so that pointer-to-member functions are not
 *                                  required; filename changed; added cerr statements.
 *      110120    E. Iorfida        Added necessary class that contains functions, and related
 *                                  code, to allow a directly test with adaptor.
 *      110120    K. Kumar          Added global functions test; updated comments; modified layout.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// Test runs code and verifies result against expected value.
// If the tested code is erroneous, the test function returns a boolean
// true; if the code is correct, the function returns a boolean false.
// 
// Need to update test to also test class containing the mathematical
// functions.
// 

#include <cmath>
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

//! Struct for NewtonRaphson unit test code.
/*!
 * This struct contains functions, necessary to test NewtonRaphson method.
 */
struct NewtonRaphsonTest
{
public:

    //! Mathematical test function.
    /*!
     * Mathematical test function used by the Newton-Raphson algorithm.
     * \param inputValue Input value.
     */
    double computeTestFunction( double& inputValue ) { return std::pow( inputValue, 2.0 ) - 3.0; }

    //! First-derivative of mathematical test function.
    /*!
     * First-derivative of mathematical test function used by the
     * Newton-Raphson algorithm.
     * \param inputValue Input value.
     */
    double computeFirstDerivativeTestFunction( double& inputValue ) { return 2.0 * inputValue; }

protected:

private:
};

//! Global mathematical test function.
double computeGlobalTestFunction( double& inputValue )
{ return std::pow( inputValue, 2.0 ) - 3.0; }

//! Global first-derivative mathematical test function.
double computeGlobalFirstDerivativeTestFunction( double& inputValue ) { return 2.0 * inputValue; }

//! Test of Newton-Raphson method code.
int main( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;
    using namespace tudat;

    // Two tests.
    // Test 1: Test of Newton-Raphson code using global functions.
    // Test 2: Test of Newton-Raphson code using member-functions.

    // Test result initialised to false.
    bool isNewtonRaphsonMethodErroneous = false;

    // Expected test result.
    double expectedResult = std::sqrt( 3.0 );

    // Create pointers to new NewtonRaphson objects.
    NewtonRaphson* pointerToMyNewtonRaphsonTest1 = new NewtonRaphson;
    NewtonRaphson* pointerToMyNewtonRaphsonTest2 = new NewtonRaphson;

    // Test 1: Test of Newton-Raphson code using global functions.
    // Set values for the implementation of the code.
    pointerToMyNewtonRaphsonTest1->setTolerance( 1.0e-15 );
    pointerToMyNewtonRaphsonTest1->setInitialGuessOfRoot( 5.0 );

    // Set mathematical functions.
    pointerToMyNewtonRaphsonTest1->setMathematicalFunction( &computeGlobalTestFunction );
    pointerToMyNewtonRaphsonTest1->setFirstDerivativeMathematicalFunction(
                &computeGlobalFirstDerivativeTestFunction );

    // Compute Newton-Raphson method.
    pointerToMyNewtonRaphsonTest1->execute( );

    // Test 2: Test of Newton-Raphson code using member-functions.
    // Set values for the implementation of the code.
    pointerToMyNewtonRaphsonTest2->setTolerance( 1.0e-15 );
    pointerToMyNewtonRaphsonTest2->setInitialGuessOfRoot( 5.0 );

    // Create pointer to adaptor object of NewtonRaphsonAdaptor class.
    NewtonRaphsonAdaptor< NewtonRaphsonTest > newtonRaphsonAdaptorForNewtonRaphsonTest_;

    // Run code using NewtonRapshonTest class and adaptor class object.
    pointerToMyNewtonRaphsonTest2->setNewtonRaphsonAdaptor(
                &newtonRaphsonAdaptorForNewtonRaphsonTest_ );
    newtonRaphsonAdaptorForNewtonRaphsonTest_.setPointerToFunction(
                &NewtonRaphsonTest::computeTestFunction );
    newtonRaphsonAdaptorForNewtonRaphsonTest_.setPointerToFirstDerivativeFunction(
                &NewtonRaphsonTest::computeFirstDerivativeTestFunction );

    // Compute Newton-Raphson method.
    pointerToMyNewtonRaphsonTest2->execute( );

    // Set test result to true if either test does not match the expected
    // result.
    if ( std::fabs( pointerToMyNewtonRaphsonTest1->getComputedRootOfFunction( )
                    - expectedResult ) >= pointerToMyNewtonRaphsonTest1->getTolerance( )
         || std::fabs( pointerToMyNewtonRaphsonTest2->getComputedRootOfFunction( )
                       - expectedResult ) >= pointerToMyNewtonRaphsonTest2->getTolerance( ) )
    {
        // Set error flag to true.
        isNewtonRaphsonMethodErroneous = true;

        if ( std::fabs( pointerToMyNewtonRaphsonTest1->getComputedRootOfFunction( )
                        - expectedResult ) >= pointerToMyNewtonRaphsonTest2->getTolerance( ) )
        {
            // Generate error statements.
            cerr << "The computed value ( "
                 << pointerToMyNewtonRaphsonTest1->getComputedRootOfFunction( )
                 << " ) using the Newton-Raphson method and global functions "
                 << "does not match the expected solution (" << expectedResult << " )." << endl;
            cerr << "The difference is: "
                 << std::fabs( expectedResult
                               - pointerToMyNewtonRaphsonTest1->getComputedRootOfFunction( ) )
                 << endl;
        }

        if ( std::fabs( pointerToMyNewtonRaphsonTest2->getComputedRootOfFunction( )
                        - expectedResult ) )
        {
            isNewtonRaphsonMethodErroneous = true;

            // Generate error statements.
            cerr << "The computed value ( "
                 <<  pointerToMyNewtonRaphsonTest2->getComputedRootOfFunction( )
                 << " ) using the Newton-Raphson method and member-functions "
                 << "does not match the expected solution (" << expectedResult << " )." << endl;
            cerr << "The difference is: "
                 << std::fabs( expectedResult -
                               pointerToMyNewtonRaphsonTest2->getComputedRootOfFunction( ) )
                 << endl;
        }
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isNewtonRaphsonMethodErroneous )
    {
        cerr << "testNewtonRaphsonMethod failed!" << std::endl;
    }

    return isNewtonRaphsonMethodErroneous;
}
