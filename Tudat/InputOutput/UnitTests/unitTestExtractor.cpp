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

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/InputOutput/extractor.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace input_output
{

// A dummy derived class to test the functionality of the extractor.
class DummyExtractor : public Extractor< double >
{
public:

    // Implement extract function, which simply returns a pointer to a double.
    boost::shared_ptr< double > extract( parsed_data_vector_utilities::ParsedDataLineMapPtr data )
    {
        // Create pointer to a double with value 2.0.
        boost::shared_ptr< double > dummyValue = boost::make_shared< double >( 2.0 );

        // Return pointer.
        return dummyValue;
    }

    // Function to access protected checkOptionalFieldType function.
    bool dummyCheckOptionalFieldType( ParsedDataLineMapPtr dataLineMap,
                                      int numberOfFields,
                                      input_output::FieldType firstFieldType,
                                      input_output::FieldType secondFieldType)
    {
        return checkOptionalFieldType( dataLineMap, numberOfFields,
                                       firstFieldType, secondFieldType );
    }

    // Function to access protected checkRequiredFieldType function.
    void dummyCheckRequiredFieldType( ParsedDataLineMapPtr dataLineMap,
                                      int numberOfFields,
                                      input_output::FieldType firstFieldType,
                                      input_output::FieldType secondFieldType)
   {
       checkRequiredFieldType( dataLineMap, numberOfFields, firstFieldType, secondFieldType );
   }

protected:

private:

};

} // namespace input_output

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_extractor )

//! Test the extract function.
BOOST_AUTO_TEST_CASE( extractor_Extract )
{
    // Using declaration.
    using namespace input_output;
    namespace field_types = input_output::field_types;

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Create extractor.
    DummyExtractor testExtractor;

    // Extract test data map.
    boost::shared_ptr<double> returnedValue = testExtractor.extract( testDataMap );

    // Verify that the returned value corresponds to the expected value.
    BOOST_CHECK_EQUAL( *returnedValue, 2.0 );
}

//! Test if the check for optional field types is done correctly.
BOOST_AUTO_TEST_CASE( extractor_CheckOptionalFieldTypes )
{
    // Using declaration.
    using namespace input_output;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings.
    std::string testStringName = "TestName", testStringEpoch = "2456067";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueName(
                new input_output::FieldValue( field_types::general::name,
                                                     testStringName ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEpoch(
                new input_output::FieldValue( field_types::time::epoch, testStringEpoch ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
                parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::general::name, testFieldValueName ) );
    testDataMap->insert( FieldDataPair( field_types::time::epoch, testFieldValueEpoch ) );

    // Create extractor.
    DummyExtractor testExtractor;

    // Check that the extractor finds the optional fields.
    bool areNameAndEpochFound = testExtractor.dummyCheckOptionalFieldType(
                testDataMap, 2, field_types::general::name, field_types::time::epoch );

    BOOST_CHECK_EQUAL( areNameAndEpochFound, true );

    // Check that the extractor does not find absent fields.
    bool areIdAndInclinationFound = testExtractor.dummyCheckOptionalFieldType(
                testDataMap, 2, field_types::general::id, field_types::state::inclination );

    BOOST_CHECK_EQUAL( areIdAndInclinationFound, false );

    // Check that the function returns false if one of the fields is absent.
    bool areNameAndIdFound = testExtractor.dummyCheckOptionalFieldType(
                testDataMap, 2, field_types::general::name, field_types::general::id );

    BOOST_CHECK_EQUAL( areNameAndIdFound, false );
}

//! Test if the check for required field types is done correctly.
BOOST_AUTO_TEST_CASE( extractor_CheckRequiredFieldTypes )
{
    // Using declaration.
    using namespace input_output;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings.
    std::string testStringName = "TestName", testStringEpoch = "2456067";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueName(
                new input_output::FieldValue( field_types::general::name,
                                                     testStringName ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEpoch(
                new input_output::FieldValue( field_types::time::epoch, testStringEpoch ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::general::name, testFieldValueName ) );
    testDataMap->insert( FieldDataPair( field_types::time::epoch, testFieldValueEpoch ) );

    // Create extractor.
    DummyExtractor testExtractor;

    // Set flags.
    bool areNameAndEpochFound = true, areNameAndIdFound = true, areIdAndInclinationFound = true;

    // Check for name and epoch, which should NOT result in a runtime error.
    try
    {
        // Check for required field types, name and epoch.
        testExtractor.dummyCheckRequiredFieldType( testDataMap, 2,
                                                   field_types::general::name,
                                                   field_types::time::epoch );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        areNameAndEpochFound = false;
    }

    // Check that name and epoch field types were found.
    BOOST_CHECK( areNameAndEpochFound );

    // Check for id and inclination, which should result in a runtime error.
    try
    {
        // Check for required field types, name and epoch.
        testExtractor.dummyCheckRequiredFieldType( testDataMap, 2,
                                                   field_types::general::id,
                                                   field_types::state::inclination );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        areIdAndInclinationFound = false;
    }

    // Check that name and epoch field types were found.
    BOOST_CHECK( !areIdAndInclinationFound );

    // Check for id and name, which should result in a runtime error.
    try
    {
        // Check for required field types, name and epoch.
        testExtractor.dummyCheckRequiredFieldType( testDataMap, 2,
                                                   field_types::general::id,
                                                   field_types::general::name );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        areNameAndIdFound = false;
    }

    // Check that name and epoch field types were found.
    BOOST_CHECK( !areNameAndIdFound );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
