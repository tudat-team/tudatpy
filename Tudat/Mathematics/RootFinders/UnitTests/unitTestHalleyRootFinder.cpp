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
 *      140219    E. Brandon        File copied from Newton-Raphson unit test.
 *
 *    References
 *
 *    Notes
 *      This unit test does not include testing with the function
 *      testFunctionWithLargeRootDifferences, because this test function does not provide a second
 *      derivative.
 *
 */

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/RootFinders/halleyRootFinder.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunction1.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunction2.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunction3.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunctionWithZeroRoot.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( testsuite_rootfinders )

using namespace tudat;
using namespace root_finders;
using namespace root_finders::termination_conditions;

//! Check if Halley method converges on test function #1 (TestFunction1).
BOOST_AUTO_TEST_CASE( test_halleyRootFinder_testFunction1 )
{
    // Create object containing the test functions.
    boost::shared_ptr< TestFunction1 > testFunction = boost::make_shared< TestFunction1 >( 2 );

    // The termination condition.
    HalleyRootFinder::TerminationFunction terminationConditionFunction =
            boost::bind( &RootAbsoluteToleranceTerminationCondition::checkTerminationCondition,
                         boost::make_shared< RootAbsoluteToleranceTerminationCondition >(
                             testFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );

    // Test Halley object.
    HalleyRootFinder halleyRootFinder( terminationConditionFunction );

    // Let the Halley method search for the root.
    const double root = halleyRootFinder.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Halley method converges on test function #2 (TestFunction2).
BOOST_AUTO_TEST_CASE( test_halleyRootFinder_testFunction2 )
{
    // Create object containing the test functions.
    boost::shared_ptr< TestFunction2 > testFunction = boost::make_shared< TestFunction2 >( 2 );

    // The termination condition.
    HalleyRootFinder::TerminationFunction terminationConditionFunction =
            boost::bind( &RootAbsoluteToleranceTerminationCondition::checkTerminationCondition,
                         boost::make_shared< RootAbsoluteToleranceTerminationCondition >(
                             testFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );

    // Test Halley object.
    HalleyRootFinder halleyRootFinder( terminationConditionFunction );

    // Let the Halley method search for the root.
    const double root = halleyRootFinder.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Halley method converges on test function #3 (TestFunction3).
BOOST_AUTO_TEST_CASE( test_halleyRootFinder_testFunction3 )
{
    // Create object containing the test functions.
    boost::shared_ptr< TestFunction3 > testFunction = boost::make_shared< TestFunction3 >( 2 );

    // The termination condition.
    HalleyRootFinder::TerminationFunction terminationConditionFunction =
            boost::bind( &RootAbsoluteToleranceTerminationCondition::checkTerminationCondition,
                         boost::make_shared< RootAbsoluteToleranceTerminationCondition >(
                             testFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );

    // Test Halley object.
    HalleyRootFinder halleyRootFinder( terminationConditionFunction );

    // Let the Halley method search for the root.
    const double root = halleyRootFinder.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Halley method converges on function with zero root (testFunctionWithZeroRoot).
// Not the best test case. Inheritance from old code. Not really relevant anymore. The basic
// idea is that Halley's method should work for both a function that becomes zero, as well as for a
// function that does not become zero. A better case should be written.
BOOST_AUTO_TEST_CASE( test_halleyRootFinder_testFunctionWithZeroRoot )
{
    // Create object containing the test functions.
    boost::shared_ptr< TestFunctionWithZeroRoot > testFunction =
            boost::make_shared< TestFunctionWithZeroRoot >( 2 );

    // The termination condition.
    HalleyRootFinder::TerminationFunction terminationConditionFunction
            = boost::bind( &RootAbsoluteToleranceTerminationCondition::
                           checkTerminationCondition,
                           boost::make_shared< RootAbsoluteToleranceTerminationCondition >(
                               1.0e-150 ), _1, _2, _3, _4, _5 );

    // Test Halley object.
    HalleyRootFinder halleyRootFinder( terminationConditionFunction );

    // Let the Halley method search for the root.
    const double root = halleyRootFinder.execute( testFunction, testFunction->getInitialGuess( ) );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_SMALL( root, 1.0e-100 );
    BOOST_CHECK_SMALL( testFunction->evaluate( root ), 1.0e-200 );
}

BOOST_AUTO_TEST_SUITE_END( ) // testsuite_rootfinders

} // namespace unit_tests
} // namespace tudat
