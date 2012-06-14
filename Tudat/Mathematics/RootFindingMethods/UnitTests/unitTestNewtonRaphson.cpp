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

#define BOOST_TEST_MAIN

#include <cmath>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

namespace tudat
{
namespace unit_tests
{

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
{
    return std::pow( inputValue, 2.0 ) - 3.0;
}

//! Global first-derivative mathematical test function.
double computeGlobalFirstDerivativeTestFunction( double& inputValue ) { return 2.0 * inputValue; }

BOOST_AUTO_TEST_SUITE( test_newton_raphson )

//! Test if Newton-Raphson root-finder works correctly using global functions.
BOOST_AUTO_TEST_CASE( testNewtonRaphsonWithGlobalFunctions )
{
    // Set expected root.
    const double expectedRoot = std::sqrt( 3.0 );

    // Declare new Newton-Raphson object.
    boost::shared_ptr< tudat::NewtonRaphson > newtonRaphson
            = boost::make_shared< tudat::NewtonRaphson >( );

    // Set values for the implementation of the code.
    newtonRaphson->setTolerance( 1.0e-15 );
    newtonRaphson->setInitialGuessOfRoot( 5.0 );

    // Set mathematical functions.
    newtonRaphson->setMathematicalFunction( &computeGlobalTestFunction );
    newtonRaphson->setFirstDerivativeMathematicalFunction(
                &computeGlobalFirstDerivativeTestFunction );

    // Compute root.
    newtonRaphson->execute( );

    // Check if computed root matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedRoot, newtonRaphson->getComputedRootOfFunction( ),
                                newtonRaphson->getTolerance( ) );

}

//! Test if Newton-Raphson root-finder works correctly using member functions.
BOOST_AUTO_TEST_CASE( testNewtonRaphsonWithMemberFunctions )
{
    // Set expected root.
    const double expectedRoot = std::sqrt( 3.0 );

    // Declare new Newton-Raphson object.
    boost::shared_ptr< tudat::NewtonRaphson > newtonRaphson
            = boost::make_shared< tudat::NewtonRaphson >( );

    // Set values for the implementation of the code.
    newtonRaphson->setTolerance( 1.0e-15 );
    newtonRaphson->setInitialGuessOfRoot( 5.0 );

    // Declare NewtonRaphsonAdaptor object.
    tudat::NewtonRaphsonAdaptor< NewtonRaphsonTest > newtonRaphsonAdaptor_;

    // Set adaptor class object and member functions.
    newtonRaphson->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptor_ );
    newtonRaphsonAdaptor_.setPointerToFunction(
                &NewtonRaphsonTest::computeTestFunction );
    newtonRaphsonAdaptor_.setPointerToFirstDerivativeFunction(
                &NewtonRaphsonTest::computeFirstDerivativeTestFunction );

    // Compute root.
    newtonRaphson->execute( );

    // Check if computed root matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedRoot, newtonRaphson->getComputedRootOfFunction( ),
                                newtonRaphson->getTolerance( ) );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
