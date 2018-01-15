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
 *      parsedDataVectorUtilities contains two (overloaded) functions called DUMP. These are used to
 *      print the contents of a parsed data vector to an ostream (e.g. std::cout). This
 *      functionality is not tested in this unit tests because one cannot compare streams (see
 *      http://stackoverflow.com/questions/7350392/comparing-streams). Note however that the
 *      DUMP functions were tested by visually inspecting the output produced.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace unit_tests
{

// Define Boost test suite.
BOOST_AUTO_TEST_SUITE( test_parsed_data_vector_utilities )

//! Test that the correct field is returned.
BOOST_AUTO_TEST_CASE( testParsedDataVectorUtilitiesGetFieldFunction )
{
    // Using declaration.
    using namespace input_output;

    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< FieldType, parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings.
    std::string testStringName = "TestName", testStringEpoch = "2456067";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueName(
                new FieldValue( field_types::general::name, testStringName ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEpoch(
                new FieldValue( field_types::time::epoch, testStringEpoch ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< FieldType, parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::general::name, testFieldValueName ) );
    testDataMap->insert( FieldDataPair( field_types::time::epoch, testFieldValueEpoch ) );

    // Retrieve name field value.
    std::string returnedStringName =
            parsed_data_vector_utilities::getField< std::string >(
                testDataMap, field_types::general::name );

    // Retrieve the epoch field value, as string, integer and double.
    std::string returnedStringEpoch = parsed_data_vector_utilities::getField< std::string >(
                testDataMap, field_types::time::epoch );
    int returnedIntegerEpoch = parsed_data_vector_utilities::getField< int >(
                testDataMap, field_types::time::epoch );
    double returnedDoubleEpoch = parsed_data_vector_utilities::getField< double >(
                testDataMap, field_types::time::epoch );

    // Verify that returned values are correct.
    BOOST_CHECK_EQUAL( returnedStringName, testStringName );

    BOOST_CHECK_EQUAL( returnedStringEpoch, testStringEpoch );

    BOOST_CHECK_EQUAL( returnedIntegerEpoch, 2456067 );

    BOOST_CHECK_EQUAL( returnedDoubleEpoch, 2456067.0 );
}

//! Test that the data vector is correctly filtered based on the given field types.
BOOST_AUTO_TEST_CASE( testParsedDataVectorUtilitiesFilterMapKeyFunction )
{
    // Using declaration.
    using namespace input_output;

    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< FieldType, parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings.
    std::string testStringName1 = "TestName1", testStringEpoch1 = "2456067",
            testStringName2 = "TestName2", testStringID1 = "5683",
            testStringID2 = "5869", testStringEpoch2 = "2456068";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueName1(
                new FieldValue( field_types::general::name, testStringName1 ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEpoch1(
                new FieldValue( field_types::time::epoch, testStringEpoch1 ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueName2(
                new FieldValue( field_types::general::name, testStringName2 ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueID1(
                new FieldValue( field_types::general::id, testStringID1 ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueID2(
                new FieldValue( field_types::general::id, testStringID2 ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEpoch2(
                new FieldValue( field_types::time::epoch, testStringEpoch2 ) );

    // Create data maps.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMapNameEpoch =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< FieldType, parsed_data_vector_utilities::FieldValuePtr >( ) );
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMapNameID =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< FieldType, parsed_data_vector_utilities::FieldValuePtr >( ) );
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMapIDEpoch =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< FieldType, parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data maps.
    testDataMapNameEpoch->insert( FieldDataPair(
                                      field_types::general::name, testFieldValueName1 ) );
    testDataMapNameEpoch->insert( FieldDataPair(
                                      field_types::time::epoch, testFieldValueEpoch1 ) );
    testDataMapNameID->insert( FieldDataPair( field_types::general::name, testFieldValueName2 ) );
    testDataMapNameID->insert( FieldDataPair( field_types::general::id, testFieldValueID1 ) );
    testDataMapIDEpoch->insert( FieldDataPair (field_types::general::id, testFieldValueID2 ) );
    testDataMapIDEpoch->insert( FieldDataPair (field_types::time::epoch, testFieldValueEpoch2 ) );

    // Create data map vector.
    parsed_data_vector_utilities::ParsedDataVectorPtr testDataMapVector(
                new parsed_data_vector_utilities::ParsedDataVector );

    // Populate data map vector.
    testDataMapVector->push_back( testDataMapNameEpoch );
    testDataMapVector->push_back( testDataMapNameID );
    testDataMapVector->push_back( testDataMapIDEpoch );

    // Test 1: Filter data map vector for data maps containing epoch entries.
    parsed_data_vector_utilities::ParsedDataVectorPtr testFilteredDataMapVectorEpoch
            = parsed_data_vector_utilities::filterMapKey( testDataMapVector, 1,
                                                          field_types::time::epoch );

    // Verify that the size of the vector is equal to 2.
    BOOST_CHECK_EQUAL( testFilteredDataMapVectorEpoch->size( ), 2 );

    // Verify that the data map vector is correctly filtered.
    BOOST_CHECK_EQUAL( testFilteredDataMapVectorEpoch->at( 0 )->find(
                           field_types::time::epoch )->second->getRaw( ), "2456067" );

    BOOST_CHECK_EQUAL( testFilteredDataMapVectorEpoch->at( 0 )->find(
                           field_types::general::name )->second->getRaw( ), "TestName1" );

    BOOST_CHECK_EQUAL( testFilteredDataMapVectorEpoch->at( 1 )->find(
                           field_types::general::id )->second->getRaw( ), "5869" );

    BOOST_CHECK_EQUAL( testFilteredDataMapVectorEpoch->at( 1 )->find(
                           field_types::time::epoch )->second->getRaw( ), "2456068" );

    // Test 2: Filter data map vector for data maps containing name and epoch entries.
    parsed_data_vector_utilities::ParsedDataVectorPtr testFilteredDataMapVectorNameEpoch
            = parsed_data_vector_utilities::filterMapKey( testDataMapVector, 2,
                                                          field_types::general::name,
                                                          field_types::time::epoch );

    // Verify that the size of the vector is equal to 1.
    BOOST_CHECK_EQUAL( testFilteredDataMapVectorNameEpoch->size( ), 1 );

    // Verify that the data map vector is correctly filtered.
    BOOST_CHECK_EQUAL( testFilteredDataMapVectorNameEpoch->at( 0 )->find(
                           field_types::time::epoch )->second->getRaw( ), "2456067" );

    BOOST_CHECK_EQUAL( testFilteredDataMapVectorNameEpoch->at( 0 )->find(
                           field_types::general::name )->second->getRaw( ), "TestName1" );
}

//! Test that the data vector is correctly filtered based on the given field types and
//! corresponding field values.
BOOST_AUTO_TEST_CASE( testParsedDataVectorUtilitiesFilterMapKeyValueFunction )
{
    // Using declaration.
    using namespace input_output;

    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< FieldType, parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings.
    std::string testStringName1 = "TestName1", testStringEpoch1 = "2456067",
            testStringName2 = "TestName2", testStringID1 = "5683",
            testStringID2 = "5869", testStringEpoch2 = "2456068";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueName1(
                new FieldValue( field_types::general::name, testStringName1 ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEpoch1(
                new FieldValue( field_types::time::epoch, testStringEpoch1 ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueName2(
                new FieldValue( field_types::general::name, testStringName2 ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueID1(
                new FieldValue( field_types::general::id, testStringID1 ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueID2(
                new FieldValue( field_types::general::id, testStringID2 ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEpoch2(
                new FieldValue( field_types::time::epoch, testStringEpoch2 ) );

    // Create data maps.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMapNameEpoch =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< FieldType, parsed_data_vector_utilities::FieldValuePtr >( ) );
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMapNameID =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< FieldType, parsed_data_vector_utilities::FieldValuePtr >( ) );
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMapIDEpoch =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< FieldType, parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data maps.
    testDataMapNameEpoch->insert( FieldDataPair(
                                      field_types::general::name, testFieldValueName1 ) );
    testDataMapNameEpoch->insert( FieldDataPair(
                                      field_types::time::epoch, testFieldValueEpoch1 ) );
    testDataMapNameID->insert( FieldDataPair( field_types::general::name, testFieldValueName2 ) );
    testDataMapNameID->insert( FieldDataPair( field_types::general::id, testFieldValueID1 ) );
    testDataMapIDEpoch->insert( FieldDataPair (field_types::general::id, testFieldValueID2 ) );
    testDataMapIDEpoch->insert( FieldDataPair (field_types::time::epoch, testFieldValueEpoch2 ) );

    // Create data map vector.
    parsed_data_vector_utilities::ParsedDataVectorPtr testDataMapVector(
                new parsed_data_vector_utilities::ParsedDataVector );

    // Populate data map vector.
    testDataMapVector->push_back( testDataMapNameEpoch );
    testDataMapVector->push_back( testDataMapNameID );
    testDataMapVector->push_back( testDataMapIDEpoch );

    // Test 1: Filter data map vector for data maps containing an epoch entry equal to 2456067.
    parsed_data_vector_utilities::ParsedDataVectorPtr testFilteredDataMapVectorEpoch
            = parsed_data_vector_utilities::filterMapKeyValue(
                testDataMapVector, 1, field_types::time::epoch, "2456067" );

    // Verify that the size of the vector is equal to 1.
    BOOST_CHECK_EQUAL( testFilteredDataMapVectorEpoch->size( ), 1);

    // Verify that the data map vector is correctly filtered.
    BOOST_CHECK_EQUAL( testFilteredDataMapVectorEpoch->at( 0 )->find(
                           field_types::time::epoch )->second->getRaw( ), "2456067" );

    BOOST_CHECK_EQUAL( testFilteredDataMapVectorEpoch->at( 0 )->find(
                           field_types::general::name )->second->getRaw( ), "TestName1" );

    // Test 2: Filter data map vector for data maps containing a name entry with field value within
    // a given range.
    parsed_data_vector_utilities::ParsedDataVectorPtr testFilteredDataMapVectorNameEpoch
            = parsed_data_vector_utilities::filterMapKeyValue(
                testDataMapVector, 1, field_types::general::name, "TestName[1-9]*");

    // Verify that the size of the vector is equal to 1.
    BOOST_CHECK_EQUAL( testFilteredDataMapVectorNameEpoch->size( ), 2 );

    // Verify that the data map vector is correctly filtered.
    BOOST_CHECK_EQUAL( testFilteredDataMapVectorNameEpoch->at( 0 )->find(
                           field_types::time::epoch )->second->getRaw( ), "2456067" );

    BOOST_CHECK_EQUAL( testFilteredDataMapVectorNameEpoch->at( 0 )->find(
                           field_types::general::name )->second->getRaw( ), "TestName1" );

    BOOST_CHECK_EQUAL( testFilteredDataMapVectorNameEpoch->at( 1 )->find(
                           field_types::general::id )->second->getRaw( ), "5683" );

    BOOST_CHECK_EQUAL( testFilteredDataMapVectorNameEpoch->at( 1 )->find(
                           field_types::general::name )->second->getRaw( ), "TestName2" );
}

// Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
