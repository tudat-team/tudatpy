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
 *      130301    S. Billemont      Update to fieldValue definition.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/InputOutput/fieldValue.h"
#include "Tudat/InputOutput/fieldTransform.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::input_output;

//! Test transform struct used to test FieldTransform.
struct TestTransform : public tudat::input_output::FieldTransform
{
    static const boost::shared_ptr< std::string > result;
    boost::shared_ptr< std::string > transform( const std::string& input ) { return result; }
};

//! Set result of transform in TestTransform to a string.
const boost::shared_ptr< std::string > TestTransform::result( new std::string( "Transformed!") );

// Define Boost test suite.
BOOST_AUTO_TEST_SUITE( test_field_value )

//! Test that the raw data is returned.
BOOST_AUTO_TEST_CASE( testFieldValueGetRawFunction )
{
    // Create raw string data.
    const std::string testString = "raw";

    // Create a test transformer.
    boost::shared_ptr< TestTransform > transform( new TestTransform( ) );

    // Create a field value object.
    FieldValue testFieldValue( field_types::state::inclination, testString, transform );

    // Retrieve raw data.
    const std::string returnedData = testFieldValue.getRaw( );

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( returnedData, testString );
}

//! Test that the raw transformed data is returned.
BOOST_AUTO_TEST_CASE( testFieldValueGetTransformedFunction )
{
    // Create raw string data.
    const std::string testString = "raw";

    // Create a test transformer.
    boost::shared_ptr< TestTransform > transform( new TestTransform( ) );

    // Create a field value object.
    FieldValue testFieldValue( field_types::state::inclination, testString, transform );

    // Retrieve raw data.
    const std::string returnedData = testFieldValue.getTransformed( );

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( returnedData, *TestTransform::result );
}

//! Test that the transformed data is returned from get.
BOOST_AUTO_TEST_CASE( testFieldValueGetFunction )
{
    // Create raw string data.
    const std::string testString = "raw";

    // Create a test transformer.
    boost::shared_ptr< TestTransform > transform( new TestTransform( ) );

    // Test 1: Test the get function with a transform.
    // Create a field value object with transform.
    FieldValue testFieldValue( field_types::state::inclination, testString, transform );

    // Retrieve (transformed) data.
    const std::string returnedData = testFieldValue.get< std::string >( );

    // Expected data.
    const std::string expectedData = *TestTransform::result;

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( returnedData, expectedData );

    // Test 2: Test the get function without transform.
    // Create a field value object.
    FieldValue testFieldValue2( field_types::state::inclination, testString );

    // Retrieve (raw) data.
    const std::string returnedRawData = testFieldValue2.get< std::string >( );

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( returnedRawData, testString );
}

//! Test that the transformed data is returned from getPointer.
BOOST_AUTO_TEST_CASE( testFieldValueGetPointerFunction )
{
    // Create raw string data.
    const std::string testString = "raw";

    // Create a test transformer.
    boost::shared_ptr< TestTransform > transform( new TestTransform( ) );

    // Test 1: Test the get function with a transform.
    // Create a field value object with transform.
    FieldValue testFieldValue( field_types::state::inclination, testString, transform );

    // Retrieve (transformed) data.
    boost::shared_ptr< std::string > returnedData = testFieldValue.getPointer< std::string >( );

    // Expected data.
    const std::string expectedData = *TestTransform::result;

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( *returnedData, expectedData );

    // Test 2: Test the get function without transform.
    // Create a field value object.
    FieldValue testFieldValue2( field_types::state::inclination, testString );

    // Retrieve (raw) data.
    boost::shared_ptr< std::string > returnedRawData
            = testFieldValue2.getPointer< std::string >( );

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( *returnedRawData, testString );
}

//! Test that the field type is stored correctly.
BOOST_AUTO_TEST_CASE( testFieldValueType )
{
    // Inserted and expected resulting field type.
    FieldType testType = field_types::general::name;

    // Test constructor without transform.
    FieldValue testFieldValue1( testType, "foo" );
    BOOST_CHECK_EQUAL( testType, testFieldValue1.type );

    // Test constructor with transform.
    FieldValue testFieldValue2( testType, "foo", boost::make_shared< TestTransform >(  ) );
    BOOST_CHECK_EQUAL( testType, testFieldValue2.type );
}

// Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
