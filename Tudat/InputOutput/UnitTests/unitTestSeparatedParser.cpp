/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      111103    S. Billemont      First creation of code.
 *      111130    S. Billemont      Rewrote tests in boost.test form.
 *      120718    A. Ronse          Added unit test from "whitespaceparser".
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <sstream>
#include <string>

#include <boost/assign.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Tudat/InputOutput/fieldType.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/InputOutput/separatedParser.h"
#include "Tudat/InputOutput/linearFieldTransform.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( testsuite_separated_parser )

using namespace tudat::input_output::parsed_data_vector_utilities;

//! Create a parser that will be used in all tests (see fixtures).
struct test_separatedparser_fixture
{
public:

    //! Default constructor.
    test_separatedparser_fixture( )
        : parser( std::string( ", " ), 4,
                  tudat::input_output::field_types::general::id,
                  tudat::input_output::field_types::state::cartesianXCoordinate,
                  tudat::input_output::field_types::state::cartesianYCoordinate,
                  tudat::input_output::field_types::state::cartesianZCoordinate )
    { }

    //! Default destructor.
    ~test_separatedparser_fixture( ) { }

    //! A separated parser.
    tudat::input_output::SeparatedParser parser;

protected:

private:
};

// Create a new test suite for all the seperated parser tests, with the parser fixture.
BOOST_FIXTURE_TEST_SUITE( testsuite_separatedParser, test_separatedparser_fixture )

BOOST_AUTO_TEST_CASE( separatedParser_singleLine )
{
    using namespace tudat::input_output::field_types;

    // Create test data.
    std::string csvText( " 12 , 1, 2, 3,2" );

    // Parse the data.
    ParsedDataVectorPtr result = parser.parse(csvText);
    
    // Check that only one line is parsed.
    BOOST_CHECK_EQUAL( result->size( ), 1 ); // result->size() == 1
    
    // Retrieve the single line.
    ParsedDataLineMapPtr myLineData = result->at( 0 );

    // Check that it parsed 4 fields.
    BOOST_CHECK_EQUAL( myLineData->size( ), 4 );

    // Check if the data was correcly separated.
    BOOST_CHECK_EQUAL( *( myLineData->find( general::id )->second->getRaw( ) ),
                       ( "12" ) ); // Trim additional spaces is enabled, should no be padded.
    BOOST_CHECK_EQUAL( *( myLineData->find( state::cartesianXCoordinate )->second->getRaw( ) ),
                       ( "1" ) );
    BOOST_CHECK_EQUAL( *( myLineData->find( state::cartesianYCoordinate )->second->getRaw( ) ),
                       ( "2" ) );
    BOOST_CHECK_EQUAL( *( myLineData->find( state::cartesianZCoordinate)->second->getRaw( ) ),
                       ( "3,2" ) ); // No space behind the comma, so it should not be split.

    // Check if it fails to find a field that is not passed in the constructor
    BOOST_CHECK( myLineData->find( time::epoch ) == myLineData->end( ) );
}

BOOST_AUTO_TEST_CASE( separatedParser_multiLine )
{
    using namespace tudat::input_output::parsed_data_vector_utilities;
    using namespace tudat::input_output::field_types;

    // Create test data.
    std::string csvText( " 12 , 1, 2, 3\n24, 5, 6, 7, 8\na, b c, d, e" );

    // Declare stream to capture std::cerr warning that the number of fields defined (5) is
    // greater than the number of fields provided in the parser (4).
    std::stringstream warningStream;

    // Redirect std::cerr output.
    std::streambuf* standardStreamBuffer = std::cerr.rdbuf( warningStream.rdbuf( ) );

    // Parse the data.
    parser.setTrim( false );
    ParsedDataVectorPtr result = parser.parse( csvText );

    // Hand std::cerr output back to standard logger.
    std::cerr.rdbuf( standardStreamBuffer );

    // Set expected warning message.
    std::stringstream expectedWarning;
    expectedWarning << "Number of elements in the line (5) "
                    << "does not match the specified number of fields (4)" << std::endl;

    // Check that warning outputted matches expected warning.
    BOOST_CHECK( warningStream.str( ).compare( expectedWarning.str( ) ) == 0 );
    
    // Check that three lines are parsed.
    BOOST_CHECK_EQUAL( result->size( ), 3 );
    
    // Retrieve the first line.
    ParsedDataLineMapPtr myLineData = result->at( 0 );

    // Check that it parsed 4 fields.
    BOOST_CHECK_EQUAL( myLineData->size( ), 4 );

    // Check if the data was correcly separated.
    BOOST_CHECK_EQUAL( *( myLineData->find( general::id )->second->getRaw( ) ),
                       ( " 12 " ) ); // Trim additional spaces is disabled, so it should give
                                     // padding spaces.

    BOOST_CHECK_EQUAL( *( myLineData->find( state::cartesianXCoordinate )->second->getRaw( ) ),
                       ( "1" ) );
    BOOST_CHECK_EQUAL( *( myLineData->find( state::cartesianYCoordinate )->second->getRaw( ) ),
                       ( "2" ) );
    BOOST_CHECK_EQUAL( *(myLineData->find( state::cartesianZCoordinate )->second->getRaw( ) ),
                       ( "3" ) );

    // Check if it fails to find a field that is not passed in the constructor.
    BOOST_CHECK( myLineData->find( time::epoch ) == myLineData->end( ) );

    // Similar for line 2.
    myLineData = result->at( 1 );

    // There should only be four fields even though 5 are provided. It should ignore the last
    // field.
    BOOST_CHECK_EQUAL( myLineData->size( ), 4 );
    BOOST_CHECK_EQUAL( *( myLineData->find( general::id )->second->getRaw( ) ), ( "24" ) );
    BOOST_CHECK_EQUAL( *( myLineData->find( state::cartesianXCoordinate )->second->getRaw( ) ),
                       ( "5" ) );
    BOOST_CHECK_EQUAL( *( myLineData->find( state::cartesianYCoordinate )->second->getRaw( ) ),
                       ( "6" ) );
    BOOST_CHECK_EQUAL( *( myLineData->find( state::cartesianZCoordinate )->second->getRaw( ) ),
                       ( "7" ) );

    // Test last line.
    myLineData = result->at( 2 );
    BOOST_CHECK_EQUAL( myLineData->size( ), 4 );
    BOOST_CHECK_EQUAL( *(myLineData->find( general::id )->second->getRaw( ) ), ( "a" ) );
    BOOST_CHECK_EQUAL( *( myLineData->find( state::cartesianXCoordinate )->second->getRaw( ) ),
                       ( "b c" ) );
    BOOST_CHECK_EQUAL( *( myLineData->find( state::cartesianYCoordinate )->second->getRaw( ) ),
                       ( "d" ) );
    BOOST_CHECK_EQUAL( *( myLineData->find( state::cartesianZCoordinate )->second->getRaw( ) ),
                       ( "e" ) );
}

BOOST_AUTO_TEST_CASE( separatedParser_whitespace )
{
    using namespace tudat::input_output::parsed_data_vector_utilities;
    using namespace tudat::input_output::field_types;

    // Create GTOC2 white space parser.
    tudat::input_output::SeparatedParser
            testWhiteSpaceParser( " ",
                                  8,
                                  general::id,
                                  state::semiMajorAxis,
                                  state::eccentricity,
                                  state::inclination,
                                  state::longitudeOfAscendingNode,
                                  state::argumentOfPeriapsis,
                                  state::meanAnomaly,
                                  time::epoch );

    // Create test data (First line in GTOC2 problem data file with additional whitespaces).
    std::string testString(
                "  2011542 3.9501468 0.2391642   6.87574   16.88982  48.9603 229.49648 54000   " );

    // Parse the data.
    ParsedDataVectorPtr testResult = testWhiteSpaceParser.parse( testString );

    // Check that only one line is parsed.
    BOOST_CHECK_EQUAL( testResult->size( ), 1 );

    // Retrieve the single line.
    ParsedDataLineMapPtr testLineData = testResult->at( 0 );

    // Check that it parsed 9 fields.
    BOOST_CHECK_EQUAL( testLineData->size( ), 8 );

    // Check if the data was correcly separated.
    BOOST_CHECK_EQUAL( *( testLineData->find( general::id )->second->getRaw( ) ), ( "2011542" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find( state::semiMajorAxis )->second->getRaw( ) ),
                       ( "3.9501468" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find( state::eccentricity )->second->getRaw( ) ),
                       ( "0.2391642" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find( state::inclination )->second->getRaw( ) ),
                       ( "6.87574" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find( state::longitudeOfAscendingNode
                                              )->second->getRaw( ) ), ( "16.88982" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find( state::argumentOfPeriapsis )->second->getRaw( ) ),
                       ( "48.9603" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find( state::meanAnomaly )->second->getRaw( ) ),
                       ( "229.49648" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find( time::epoch )->second->getRaw( ) ), ( "54000" ) );

    // Check if it fails to find a field that is not passed in the constructor.
    BOOST_CHECK( testLineData->find( state::trueAnomaly ) == testLineData->end( ) );
}

BOOST_AUTO_TEST_CASE( separatedParser_fieldTransform )
{
    using namespace tudat::input_output;
    using namespace tudat::input_output;

    // Create GTOC2 white space parser.
    tudat::input_output::SeparatedParser
            testWhiteSpaceParser( " ",
                                  8,
                                  field_types::general::id,
                                  field_types::state::semiMajorAxis,
                                  field_types::state::eccentricity,
                                  field_types::state::inclination,
                                  field_types::state::longitudeOfAscendingNode,
                                  field_types::state::argumentOfPeriapsis,
                                  field_types::state::meanAnomaly,
                                  field_types::time::epoch );

    // Create test data (First line in GTOC2 problem data file with additional whitespaces).
    std::string testString(
                "  2011542 3.9501468 0.2391642   6.87574   16.88982  48.9603 229.49648 54000   " );

    // Create unit transformation map
    std::map< FieldType, boost::shared_ptr< FieldTransform > > unitTransformationMap =
            boost::assign::map_list_of
            ( field_types::state::semiMajorAxis, boost::shared_ptr< FieldTransform >(
                  new LinearFieldTransform(
                      tudat::unit_conversions::convertAstronomicalUnitsToMeters< double >( 1.0 ),
                      0.0 ) ) )
            ( field_types::state::inclination, boost::shared_ptr< FieldTransform >(
                  new LinearFieldTransform(
                      tudat::unit_conversions::convertDegreesToRadians< double >( 1.0 ),
                      0.0 ) ) )
            ( field_types::state::longitudeOfAscendingNode, boost::shared_ptr< FieldTransform >(
                  new LinearFieldTransform(
                      tudat::unit_conversions::convertDegreesToRadians< double >( 1.0 ),
                      0.0 ) ) )
            ( field_types::state::argumentOfPeriapsis, boost::shared_ptr< FieldTransform >(
                  new LinearFieldTransform(
                      tudat::unit_conversions::convertDegreesToRadians< double >( 1.0 ),
                      0.0 ) ) )
            ( field_types::state::meanAnomaly, boost::shared_ptr< FieldTransform >(
                  new LinearFieldTransform(
                      tudat::unit_conversions::convertDegreesToRadians< double >( 1.0 ),
                      0.0 ) ) )
            ( field_types::time::epoch, boost::shared_ptr< FieldTransform >(
                  new LinearFieldTransform( 1.0, 2400000.5 ) ) );

    // Pass unit transformation map to parser.
    testWhiteSpaceParser.setUnitTransformationMap( unitTransformationMap );

    // Parse the data.
    ParsedDataVectorPtr testResult = testWhiteSpaceParser.parse( testString );

    // Check that only one line is parsed.
    BOOST_CHECK_EQUAL( testResult->size( ), 1 );

    // Retrieve the single line.
    ParsedDataLineMapPtr testLineData = testResult->at( 0 );

    // Check that it parsed 9 fields.
    BOOST_CHECK_EQUAL( testLineData->size( ), 8 );

    // Check if the data was correcly separated.
    BOOST_CHECK_EQUAL( *( testLineData->find( field_types::general::id )->second->getRaw( ) ),
                       ( "2011542" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find(
                              field_types::state::semiMajorAxis )->second->get( ) ),
                       ( "590933550196.86743" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find( field_types::state::eccentricity )->second->get( ) ),
                       ( "0.2391642" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find( field_types::state::inclination )->second->get( ) ),
                       ( "0.12000430151107493" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find(
                              field_types::state::longitudeOfAscendingNode )->second->get( ) ),
                       ( "0.29478296906918866" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find(
                              field_types::state::argumentOfPeriapsis )->second->get( ) ),
                       ( "0.85451843776417968" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find( field_types::state::meanAnomaly )->second->get( ) ),
                       ( "4.005469197737316" ) );
    BOOST_CHECK_EQUAL( *( testLineData->find( field_types::time::epoch )->second->get( ) ),
                       ( "2454000.5" ) );

    // Check if it fails to find a field that is not passed in the constructor.
    BOOST_CHECK( testLineData->find( field_types::state::trueAnomaly ) == testLineData->end( ) );
}

BOOST_AUTO_TEST_SUITE_END( ) // testsuite_separatedParser

BOOST_AUTO_TEST_SUITE_END( ) // testsuite_ephemeris

} // namespace unit_tests
} // namespace tudat
