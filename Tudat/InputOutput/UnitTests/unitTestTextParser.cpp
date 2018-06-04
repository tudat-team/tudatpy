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

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/InputOutput/textParser.h"

namespace tudat
{
namespace input_output
{

// A dummy derived class to test the functionality of the text parser.
class dummyTextParser : public TextParser
{
public:

    // Default constructor with parseAsStream == false.
    dummyTextParser( bool parseAsStream = false ) : TextParser( parseAsStream ){ }

    // Default destructor.
    ~dummyTextParser( ){ }

protected:

    // Parse line functionality.
    void parseLine( std::string& line )
    {
        // Define a new type: pair of field type and pointer to value.
        typedef std::pair< FieldType, parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

        // Create a new pointer to map of line data.
        parsed_data_vector_utilities::ParsedDataLineMapPtr currentLineData
                = std::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                    std::map< FieldType, parsed_data_vector_utilities::FieldValuePtr >( ) );

        // Register the data line with the global current parsed data vector.
        parsedData->push_back( currentLineData );

        // Store the field-value string.
        parsed_data_vector_utilities::FieldValuePtr value(
                    new FieldValue( field_types::general::name, line ) );

        // Store the type and value in the current line data.
        currentLineData->insert( FieldDataPair( field_types::general::name, value ) );
    }

private:

};

} // namespace input_output

namespace unit_tests
{

// Define Boost test suite.
BOOST_AUTO_TEST_SUITE( test_text_parser )

//! Test the line parsing functionality.
BOOST_AUTO_TEST_CASE( textParser_parseLine )
{
    // Using declaration.
    using namespace input_output::parsed_data_vector_utilities;

    // Create a dummy text parser.
    input_output::dummyTextParser testTextParser;

    // Create a string.
    std::string testString( " Test string " );

    // Parse the data.
    ParsedDataVectorPtr testResult = testTextParser.parse( testString );

    // Check that only one line is parsed.
    BOOST_CHECK_EQUAL( testResult->size( ), 1 );

    // Retrieve the single line.
    ParsedDataLineMapPtr testLineData = testResult->at( 0 );

    // Check that it parsed 1 field.
    BOOST_CHECK_EQUAL( testLineData->size( ), 1 );

    // Check if the data was correctly separated.
    BOOST_CHECK_EQUAL( testLineData->find(
                           input_output::field_types::general::name )->second->getRaw( ),
                       testString );

    // Check if it fails to find a field that is not passed in the constructor.
    BOOST_CHECK( testLineData->find( input_output::field_types::state::trueAnomaly )
                 == testLineData->end( ) );
}

//! Test the parsing of a stream on a line-per-line basis.
BOOST_AUTO_TEST_CASE( textParser_parseStreamAsLine )
{
    // Using declaration.
    using namespace input_output::parsed_data_vector_utilities;

    // Create a dummy text parser.
    input_output::dummyTextParser testTextParser;

    // Create a string.
    std::string testString( " Test string " );

    // Create stringstream object.
    std::stringstream testStringStream( testString );

    // Parse the data.
    ParsedDataVectorPtr testResult = testTextParser.parse( testStringStream );

    // Check that only one line is parsed.
    BOOST_CHECK_EQUAL( testResult->size( ), 1 );

    // Retrieve the single line.
    ParsedDataLineMapPtr testLineData = testResult->at( 0 );

    // Check that it parsed 1 field.
    BOOST_CHECK_EQUAL( testLineData->size( ), 1 );

    // Check if the data was correctly separated.
    BOOST_CHECK_EQUAL( testLineData->find(
                           input_output::field_types::general::name )->second->getRaw( ),
                       testString );

    // Check if it fails to find a field that is not passed in the constructor.
    BOOST_CHECK( testLineData->find( input_output::field_types::state::trueAnomaly )
                 == testLineData->end( ) );
}

// Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
