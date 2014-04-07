/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      120606    T. Secretin       File created.
 *      120608    E. Heeren         Placed tests in unit_tests namespace.
 *      130301    S. Billemont      Updated tests to new FieldValue definition.
 *
 *    References
 *
 *    Notes
 *      This unit test does not test the parseStream functionality of TextParser. Such a test must
 *      be implemented at a later stage.
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
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
    using namespace tudat::input_output;

    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create test string.
    std::string testString( "1.5362" );

    // Create a linear field transform object: y = 2 * x + 36.
    LinearFieldTransform testLinearFieldTransform1( 2.0, 36.0 );

    // Define expected value: 2 * 1.5362 + 36 = 39.0724.
    const double expectedValue = 39.0724;

    // Transform test string.
    boost::shared_ptr< std::string > returnedValue =
            testLinearFieldTransform1.transform( testString );

    // Check that returned and expected strings are identical.
    BOOST_CHECK_CLOSE_FRACTION( boost::lexical_cast< double >( *returnedValue ), expectedValue,
                                tolerance );
}

//! Test the linear transformation functionality with no slope (a=0).
BOOST_AUTO_TEST_CASE( testLinearFieldTransform2 )
{
    // Using declaration.
    using namespace tudat::input_output;

    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create test string.
    std::string testString( "1.5362" );

    // Create a linear field transform object: y = 0 * x + 36.
    LinearFieldTransform testLinearFieldTransform2( 0.0, 36.0 );

    // Define expected value: 0 * 1.5362 + 36 = 36.
    const double expectedValue = 36.0;

    // Transform test string.
    boost::shared_ptr< std::string > returnedValue =
            testLinearFieldTransform2.transform( testString );

    // Check that returned and expected strings are identical.
    BOOST_CHECK_CLOSE_FRACTION( boost::lexical_cast< double >( *returnedValue ), expectedValue,
                                tolerance );
}

//! Test the linear transformation functionality with no intercept (b=0).
BOOST_AUTO_TEST_CASE( testLinearFieldTransform3 )
{
    // Using declaration.
    using namespace tudat::input_output;

    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create test string.
    std::string testString( "1.5362" );

    // Create a linear field transform object: y = -56 * x + 0.
    LinearFieldTransform testLinearFieldTransform3( -56.0, 0.0 );

    // Define expected value: -56 * 1.5362 + 0 = -86.0272.
    const double expectedValue = -86.0272;

    // Transform test string.
    boost::shared_ptr< std::string > returnedValue =
            testLinearFieldTransform3.transform( testString );

    // Check that returned and expected strings are identical.
    BOOST_CHECK_CLOSE_FRACTION( boost::lexical_cast< double >( *returnedValue ), expectedValue,
                                tolerance );
}

//! Test the linear transformation functionality with no slope nor intercept (a=b=0).
BOOST_AUTO_TEST_CASE( testLinearFieldTransform4 )
{
    // Using declaration.
    using namespace tudat::input_output;

    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create test string.
    std::string testString( "1.5362" );

    // Create a linear field transform object: y = 0 * x + 0.
    LinearFieldTransform testLinearFieldTransform4( 0.0, 0.0 );

    // Transform test string.
    boost::shared_ptr< std::string > returnedValue =
            testLinearFieldTransform4.transform( testString );

    // Check that returned and expected strings are identical.
    BOOST_CHECK_SMALL( boost::lexical_cast< double >( *returnedValue ), tolerance );
}

// Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
