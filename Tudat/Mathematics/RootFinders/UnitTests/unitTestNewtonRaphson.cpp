/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      110111    E. Iorfida        File created.
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
 *      120712    P. Musegaas       Changed absolute tolerance into a safe variant of relative
 *                                  tolerance. Added new unit tests for this.
 *      120318    S. Billemont      Move to new root_finders codebase.
 *      120402    T. Secretin       Code-check.
 *      120810    P. Musegaas       Code-check. Merged two branches. Various edits.
 *
 *    References
 *
 *    Notes
 *      The test cases in this unit test should be more extensive.
 *
 */

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunction1.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunction2.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunction3.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunctionWithLargeRootDifference.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunctionWithZeroRoot.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( testsuite_rootfinders )

using namespace tudat;
using namespace root_finders;
using namespace root_finders::termination_conditions;

//! Check if Newton-Raphson converges on test function #1 (TestFunction1).
BOOST_AUTO_TEST_CASE( test_newtonRaphson_testFunction1 )
{
    // Create object containing the test functions.
    boost::shared_ptr< TestFunction1 > testFunction = boost::make_shared< TestFunction1 >( 1 );

    // The termination condition.
    NewtonRaphson::TerminationFunction terminationConditionFunction =
            boost::bind( &RootAbsoluteToleranceTerminationCondition::checkTerminationCondition,
                         boost::make_shared< RootAbsoluteToleranceTerminationCondition >(
                             testFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );

    // Test Newton-Raphson object.
    NewtonRaphson newtonRaphson( terminationConditionFunction );

    // Let Newton-Raphson search for the root.
    const double root = newtonRaphson.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Newton-Raphson converges on test function #2 (TestFunction2).
BOOST_AUTO_TEST_CASE( test_newtonRaphson_testFunction2 )
{
    // Create object containing the test functions.
    boost::shared_ptr< TestFunction2 > testFunction = boost::make_shared< TestFunction2 >( 1 );

    // The termination condition.
    NewtonRaphson::TerminationFunction terminationConditionFunction =
            boost::bind( &RootAbsoluteToleranceTerminationCondition::checkTerminationCondition,
                         boost::make_shared< RootAbsoluteToleranceTerminationCondition >(
                             testFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );
    
    // Test Newton-Raphson object.
    NewtonRaphson newtonRaphson( terminationConditionFunction );

    // Let Newton-Raphson search for the root.
    const double root = newtonRaphson.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Newton-Raphson converges on test function #3 (TestFunction3).
BOOST_AUTO_TEST_CASE( test_newtonRaphson_testFunction3 )
{
    // Create object containing the test functions.
    boost::shared_ptr< TestFunction3 > testFunction = boost::make_shared< TestFunction3 >( 1 );

    // The termination condition.
    NewtonRaphson::TerminationFunction terminationConditionFunction =
            boost::bind( &RootAbsoluteToleranceTerminationCondition::checkTerminationCondition,
                         boost::make_shared< RootAbsoluteToleranceTerminationCondition >(
                             testFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );

    // Test Newton-Raphson object.
    NewtonRaphson newtonRaphson( terminationConditionFunction );

    // Let Newton-Raphson search for the root.
    const double root = newtonRaphson.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Newton-Raphson converges on function with large root difference
//! (testFunctionWithLargeRootDifference).
// Not the best test case. Inheritance from old code.
BOOST_AUTO_TEST_CASE( test_newtonRaphson_testFunctionWithLargeRootDifference )
{
    // Declare tolerance.
    const double tolerance = 1.0e-10;

    // Declare expected roots.
    const double expectedRootLowCase = 1.00000000793634;
    const double expectedRootHighCase = 7937.3386333591;

    // Create objects containing the test functions. Values were obtained during a limit case
    // gravity assist calculation (while evaluating Cassini-1 trajectory).
    boost::shared_ptr< TestFunctionWithLargeRootDifference > testFunctionLowCase =
            boost::make_shared< TestFunctionWithLargeRootDifference >
            ( 1, -3.24859999867635e18, -3248600.0, 1.5707963267949 );
    boost::shared_ptr< TestFunctionWithLargeRootDifference > testFunctionHighCase =
            boost::make_shared< TestFunctionWithLargeRootDifference >
            ( 1, -3248600.0, -3.24859999867635e18, 1.5707963267949 );

    // The termination condition.
    NewtonRaphson::TerminationFunction terminationConditionFunction
            = boost::bind( &RootRelativeToleranceTerminationCondition::checkTerminationCondition,
                           boost::make_shared< RootRelativeToleranceTerminationCondition >(
                               1.0e-10 ), _1, _2, _3, _4, _5 );

    // Make the Newton-Raphson object.
    NewtonRaphson newtonRaphson( terminationConditionFunction );

    // Let Newton-Raphson search for the root for both cases.
    double rootLowCase = newtonRaphson.execute( testFunctionLowCase, 1.0 + 1.0e-10 );
    double rootHighCase = newtonRaphson.execute( testFunctionHighCase, 1.0 + 1.0e-2 );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( rootLowCase, expectedRootLowCase, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( rootHighCase, expectedRootHighCase, tolerance );
}

//! Check if Newton-Raphson converges on function with zero root (testFunctionWithZeroRoot).
// Not the best test case. Inheritance from old code. Not really relevant anymore. The basic
// idea is that Newton-Raphson should work for both a function that becomes zero, as well as for a
// function that does not become zero. A better case should be written.
BOOST_AUTO_TEST_CASE( test_newtonRaphson_testFunctionWithZeroRoot )
{
    // Create object containing the test functions.
    boost::shared_ptr< TestFunctionWithZeroRoot > testFunction =
            boost::make_shared< TestFunctionWithZeroRoot >( 1 );

    // The termination condition.
    NewtonRaphson::TerminationFunction terminationConditionFunction
            = boost::bind( &RootAbsoluteOrRelativeToleranceTerminationCondition::
                           checkTerminationCondition,
                           boost::make_shared< RootAbsoluteOrRelativeToleranceTerminationCondition >(
                               1.0e-308, 1.0e-15 ), _1, _2, _3, _4, _5 );

    // Test Newton-Raphson object.
    NewtonRaphson newtonRaphson( terminationConditionFunction );

    // Let Newton-Raphson search for the root.
    const double root = newtonRaphson.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_SMALL( root, 1.0e-150 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

BOOST_AUTO_TEST_SUITE_END( ) // testsuite_rootfinders

} // namespace unit_tests
} // namespace tudat
