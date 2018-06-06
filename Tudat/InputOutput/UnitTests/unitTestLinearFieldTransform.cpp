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
 *      This unit test does not test the parseStream functionality of TextParser. Such a test must
 *      be implemented at a later stage.
 *
 */

#define BOOST_TEST_MAIN

#include <string>

#include <boost/make_shared.hpp>
#include <memory>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/InputOutput/linearFieldTransform.h"

namespace tudat
{
namespace unit_tests
{

// Define Boost test suite.
BOOST_AUTO_TEST_SUITE( test_Linear_Field_Transform )

//! Test the linear transformation functionality.
BOOST_AUTO_TEST_CASE( testLinearFieldTransform1 )
{
    // Using declaration.
    using namespace input_output;

    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create test string.
    std::string testString( "1.5362" );

    // Create a linear field transform object: y = 2 * x + 36.
    LinearFieldTransform testLinearFieldTransform1( 2.0, 36.0 );

    // Define expected value: 2 * 1.5362 + 36 = 39.0724.
    const double expectedValue = 39.0724;

    // Transform test string.
    std::shared_ptr< std::string > returnedValue =
            testLinearFieldTransform1.transform( testString );

    // Check that returned and expected strings are identical.
    BOOST_CHECK_CLOSE_FRACTION( std::stod( *returnedValue ), expectedValue,
                                tolerance );
}

//! Test the linear transformation functionality with no slope (a=0).
BOOST_AUTO_TEST_CASE( testLinearFieldTransform2 )
{
    // Using declaration.
    using namespace input_output;

    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create test string.
    std::string testString( "1.5362" );

    // Create a linear field transform object: y = 0 * x + 36.
    LinearFieldTransform testLinearFieldTransform2( 0.0, 36.0 );

    // Define expected value: 0 * 1.5362 + 36 = 36.
    const double expectedValue = 36.0;

    // Transform test string.
    std::shared_ptr< std::string > returnedValue =
            testLinearFieldTransform2.transform( testString );

    // Check that returned and expected strings are identical.
    BOOST_CHECK_CLOSE_FRACTION( std::stod( *returnedValue ), expectedValue,
                                tolerance );
}

//! Test the linear transformation functionality with no intercept (b=0).
BOOST_AUTO_TEST_CASE( testLinearFieldTransform3 )
{
    // Using declaration.
    using namespace input_output;

    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create test string.
    std::string testString( "1.5362" );

    // Create a linear field transform object: y = -56 * x + 0.
    LinearFieldTransform testLinearFieldTransform3( -56.0, 0.0 );

    // Define expected value: -56 * 1.5362 + 0 = -86.0272.
    const double expectedValue = -86.0272;

    // Transform test string.
    std::shared_ptr< std::string > returnedValue =
            testLinearFieldTransform3.transform( testString );

    // Check that returned and expected strings are identical.
    BOOST_CHECK_CLOSE_FRACTION( std::stod( *returnedValue ), expectedValue,
                                tolerance );
}

//! Test the linear transformation functionality with no slope nor intercept (a=b=0).
BOOST_AUTO_TEST_CASE( testLinearFieldTransform4 )
{
    // Using declaration.
    using namespace input_output;

    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create test string.
    std::string testString( "1.5362" );

    // Create a linear field transform object: y = 0 * x + 0.
    LinearFieldTransform testLinearFieldTransform4( 0.0, 0.0 );

    // Transform test string.
    std::shared_ptr< std::string > returnedValue =
            testLinearFieldTransform4.transform( testString );

    // Check that returned and expected strings are identical.
    BOOST_CHECK_SMALL( std::stod( *returnedValue ), tolerance );
}

// Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
