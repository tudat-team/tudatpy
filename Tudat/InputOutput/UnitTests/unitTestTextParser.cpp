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
                = boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
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
