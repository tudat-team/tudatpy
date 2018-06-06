/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <string>

#include <memory>
#include <boost/test/unit_test.hpp>

#include "Tudat/InputOutput/fieldValue.h"
#include "Tudat/InputOutput/fieldTransform.h"

namespace tudat
{
namespace unit_tests
{

using namespace input_output;

//! Test transform struct used to test FieldTransform.
struct TestTransform : public input_output::FieldTransform
{
    static const std::shared_ptr< std::string > result;
    std::shared_ptr< std::string > transform( const std::string& input ) { return result; }
};

//! Set result of transform in TestTransform to a string.
const std::shared_ptr< std::string > TestTransform::result( new std::string( "Transformed!") );

// Define Boost test suite.
BOOST_AUTO_TEST_SUITE( test_field_value )

//! Test that the raw data is returned.
BOOST_AUTO_TEST_CASE( testFieldValueGetRawFunction )
{
    // Create raw string data.
    const std::string testString = "raw";

    // Create a test transformer.
    std::shared_ptr< TestTransform > transform( new TestTransform( ) );

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
    std::shared_ptr< TestTransform > transform( new TestTransform( ) );

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
    std::shared_ptr< TestTransform > transform( new TestTransform( ) );

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
    std::shared_ptr< TestTransform > transform( new TestTransform( ) );

    // Test 1: Test the get function with a transform.
    // Create a field value object with transform.
    FieldValue testFieldValue( field_types::state::inclination, testString, transform );

    // Retrieve (transformed) data.
    std::shared_ptr< std::string > returnedData = testFieldValue.getPointer< std::string >( );

    // Expected data.
    const std::string expectedData = *TestTransform::result;

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( *returnedData, expectedData );

    // Test 2: Test the get function without transform.
    // Create a field value object.
    FieldValue testFieldValue2( field_types::state::inclination, testString );

    // Retrieve (raw) data.
    std::shared_ptr< std::string > returnedRawData
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
    FieldValue testFieldValue2( testType, "foo", std::make_shared< TestTransform >(  ) );
    BOOST_CHECK_EQUAL( testType, testFieldValue2.type );
}

// Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
