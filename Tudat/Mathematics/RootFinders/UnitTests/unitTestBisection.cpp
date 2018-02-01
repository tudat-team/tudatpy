/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      This unit test does not include testing with the function testFunctionWithZeroRoot, because
 *      this test function is not suited for the Bisection root finder. The Bisection method
 *      requires that the interval brackets the solution, and the functions values on both sides of
 *      the root have an opposite sign. Although the interval of the test function [-3,3] brackets
 *      the root, the function values of the lower bound and the upper bound are both positive,
 *      hence the root finder throws a runtime error.
 *
 */

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/RootFinders/bisection.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunction1.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunction2.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunction3.h"
#include "Tudat/Mathematics/RootFinders/UnitTests/testFunctionWithLargeRootDifference.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( testsuite_rootfinders )

using namespace tudat;
using namespace root_finders;
using namespace root_finders::termination_conditions;

//! Check if Bisection method converges on test function #1 (TestFunction1).
BOOST_AUTO_TEST_CASE( test_bisection_testFunction1 )
{
    // Create object containing the test function.
    boost::shared_ptr< TestFunction1 > testFunction = boost::make_shared< TestFunction1 >( 1 );

    // The termination condition.
    Bisection::TerminationFunction terminationConditionFunction =
            boost::bind( &RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                         boost::make_shared< RootAbsoluteToleranceTerminationCondition< double > >(
                             testFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );

    // Test Bisection object.
    Bisection bisection( terminationConditionFunction,
                         testFunction->getLowerBound( ), testFunction->getUpperBound( ) );

    // Let the Bisection method search for the root.
    const double root = bisection.execute( testFunction );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Bisection method converges on test function #2 (TestFunction2).
BOOST_AUTO_TEST_CASE( test_bisection_testFunction2 )
{
    // Create object containing the test function.
    boost::shared_ptr< TestFunction2 > testFunction = boost::make_shared< TestFunction2 >( 1 );

    // The termination condition.
    Bisection::TerminationFunction terminationConditionFunction =
            boost::bind( &RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                         boost::make_shared< RootAbsoluteToleranceTerminationCondition< double > >(
                             testFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );

    // Test Bisection object.
    Bisection bisection( terminationConditionFunction,
                         testFunction->getLowerBound( ), testFunction->getUpperBound( ) );

    // Let the Bisection method search for the root.
    const double root = bisection.execute( testFunction );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Bisection method converges on test function #3 (TestFunction3).
BOOST_AUTO_TEST_CASE( test_bisection_testFunction3 )
{
    // Create object containing the test function.
    boost::shared_ptr< TestFunction3 > testFunction = boost::make_shared< TestFunction3 >( 1 );

    // The termination condition.
    Bisection::TerminationFunction terminationConditionFunction =
            boost::bind( &RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                         boost::make_shared< RootAbsoluteToleranceTerminationCondition< double > >(
                             testFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );

    // Test Bisection object.
    Bisection bisection( terminationConditionFunction,
                         testFunction->getLowerBound( ), testFunction->getUpperBound( ) );

    // Let the Bisection method search for the root.
    const double root = bisection.execute( testFunction );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( root, testFunction->getTrueRootLocation( ), 1.0e-15 );
    BOOST_CHECK_LT( testFunction->evaluate( root ), testFunction->getTrueRootAccuracy( ) );
}

//! Check if Bisection method converges on function with large root difference
//! (testFunctionWithLargeRootDifference).
// Not the best test case. Inheritance from old code.
BOOST_AUTO_TEST_CASE( test_bisection_testFunctionWithLargeRootDifference )
{
    // Declare tolerance.
    const double tolerance = 1.0e-10;

    // Declare expected roots.
    const double expectedRootLowCase = 1.00000000793634;
    const double expectedRootHighCase = 7937.3386333591;

    // Declare interval for the bisection method.
    const double lowerBound = 1.0;
    const double upperBound = 10000.0;

    // Create objects containing the test functions. Values were obtained during a limit case
    // gravity assist calculation (while evaluating Cassini-1 trajectory).
    boost::shared_ptr< TestFunctionWithLargeRootDifference > testFunctionLowCase =
            boost::make_shared< TestFunctionWithLargeRootDifference >
            ( 1, -3.24859999867635e18, -3248600.0, 1.5707963267949 );
    boost::shared_ptr< TestFunctionWithLargeRootDifference > testFunctionHighCase =
            boost::make_shared< TestFunctionWithLargeRootDifference >
            ( 1, -3248600.0, -3.24859999867635e18, 1.5707963267949 );

    // The termination condition.
    Bisection::TerminationFunction terminationConditionFunction
            = boost::bind( &RootRelativeToleranceTerminationCondition< >::checkTerminationCondition,
                           boost::make_shared< RootRelativeToleranceTerminationCondition< > >(
                               1.0e-10 ), _1, _2, _3, _4, _5 );

    // Test Bisection object, per case.
    Bisection bisectionLowCase( terminationConditionFunction, lowerBound, upperBound );
    Bisection bisectionHighCase( terminationConditionFunction, lowerBound, upperBound );

    // Let the Bisection method search for the root.
    const double rootLowCase = bisectionLowCase.execute( testFunctionLowCase );
    const double rootHighCase = bisectionHighCase.execute( testFunctionHighCase );

    // Check if the result is within the requested accuracy.
    BOOST_CHECK_CLOSE_FRACTION( rootLowCase, expectedRootLowCase, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( rootHighCase, expectedRootHighCase, tolerance );
}

//! Check expected functionality. Test interval that does not bracket the root.
BOOST_AUTO_TEST_CASE( test_bisection_wrongBracket )
{
    // Create object containing the test functions.
    boost::shared_ptr< TestFunction1 > testFunction = boost::make_shared< TestFunction1 >( 1 );

    // The termination condition.
    Bisection::TerminationFunction terminationConditionFunction =
            boost::bind( &RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                         boost::make_shared< RootAbsoluteToleranceTerminationCondition< double > >(
                             testFunction->getTrueRootAccuracy( ) ), _1, _2, _3, _4, _5 );

    // Test Bisection object. The input interval does not bracket the solution.
    Bisection bisection( terminationConditionFunction,
                         testFunction->getLowerBound( ), 0.0 );

    // Check if a runtime error is thrown if the root finder is executed.
    BOOST_CHECK_THROW( bisection.execute( testFunction ),
                       std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( ) // testsuite_rootfinders

} // namespace unit_tests
} // namespace tudat
