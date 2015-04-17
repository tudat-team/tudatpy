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
 *      111209    D.J. Gondelach    File created.
 *      120615    D.J. Gondelach    Adapted code to Boost test suite.
 *      120718    A. Ronse          Added multiline and non-trim unit test.
 *      130301    S. Billemont      Updated tests to new FieldValue definition.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <string>

#include <boost/test/unit_test.hpp>

#include "Tudat/InputOutput/fieldType.h"
#include "Tudat/InputOutput/fixedWidthParser.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_fixedwidth_parser )

//! Test if a single line is correctly parsed.
BOOST_AUTO_TEST_CASE( testFixedWidthParserSingleLine )
{
    using namespace input_output::parsed_data_vector_utilities;
    using namespace input_output::field_types;

    // Create GTOC2 fixed width parser.
    input_output::FixedWidthParser
            testFixedWidthParser( 9,
                                  general::name,
                                  time::epoch,
                                  state::periapsisDistance,
                                  state::eccentricity,
                                  state::inclination,
                                  state::argumentOfPeriapsis,
                                  state::longitudeOfAscendingNode,
                                  state::timeOfPeriapsisPassage,
                                  general::id,
                                  43, 8, 12, 11, 10, 10, 10, 15, 13 );

    // Create test data.
    std::string testString( "  9P/Tempel 1                                 54951  1.50908111 "
                            "0.51694646  10.52513 178.92837  68.92947 20110112.22958 JPL 154     "
                            );

    // Parse the data.
    ParsedDataVectorPtr testResult = testFixedWidthParser.parse( testString );

    // Check that only one line is parsed.
    BOOST_CHECK_EQUAL( testResult->size( ), 1 );

    // Retrieve the single line.
    ParsedDataLineMapPtr testLineData = testResult->at( 0 );

    // Check that it parsed 9 fields.
    BOOST_CHECK_EQUAL( testLineData->size( ), 9 );

    // Check if the data was correctly separated.
    BOOST_CHECK_EQUAL( testLineData->find( general::name )->second->getRaw( ), "9P/Tempel 1" );
    BOOST_CHECK_EQUAL( testLineData->find( time::epoch )->second->getRaw( ), "54951" );
    BOOST_CHECK_EQUAL( testLineData->find( state::periapsisDistance )->second->getRaw( ),
                       "1.50908111" );
    BOOST_CHECK_EQUAL( testLineData->find( state::eccentricity )->second->getRaw( ),
                       "0.51694646" );
    BOOST_CHECK_EQUAL( testLineData->find( state::inclination )->second->getRaw( ),
                       "10.52513" );
    BOOST_CHECK_EQUAL( testLineData->find( state::argumentOfPeriapsis )->second->getRaw( ),
                       "178.92837" );
    BOOST_CHECK_EQUAL( testLineData->find(state::longitudeOfAscendingNode )->second->getRaw( ),
                       "68.92947" );
    BOOST_CHECK_EQUAL( testLineData->find( state::timeOfPeriapsisPassage )->second->getRaw( ),
                       "20110112.22958" );
    BOOST_CHECK_EQUAL( testLineData->find( general::id )->second->getRaw( ), "JPL 154" );

    // Check if it fails to find a field that is not passed in the constructor.
    BOOST_CHECK( testLineData->find( state::trueAnomaly ) == testLineData->end( ) );
}

BOOST_AUTO_TEST_CASE( testFixedWidthParserMultiLineWithoutTrim )
{
    using namespace input_output::parsed_data_vector_utilities;
    using namespace input_output::field_types;

    // Create GTOC2AsteroidEphemerides fixed width parser.
    input_output::FixedWidthParser
        testFixedWidthParser( 8,
                              general::id,
                              state::semiMajorAxis,
                              state::eccentricity,
                              state::inclination,
                              state::longitudeOfAscendingNode,
                              state::argumentOfPeriapsis,
                              state::meanAnomaly,
                              time::epoch,
                              11, 18, 18, 15, 16, 15, 17, 8 );

    // Disable trim.
    testFixedWidthParser.setTrim( false );

    // Create test data.
    std::string testString( "2011542          3.9501468       0.2391642         6.87574        "
                            "16.88982         48.9603       229.49648       54000\n2001038       "
                            "   3.9619932        0.227483         9.22988        58.20488       "
                            "307.19412       201.38868       54000" );

    // Parse the data.
    ParsedDataVectorPtr testResult = testFixedWidthParser.parse( testString );

    // Check that two lines are parsed.
    BOOST_CHECK_EQUAL( testResult->size( ), 2 );

    // Retrieve the first line.
    ParsedDataLineMapPtr Line1Data = testResult->at( 0 );

    // Check that it parsed 9 fields.
    BOOST_CHECK_EQUAL( Line1Data->size( ), 8 );

    // Check if the data was correctly separated.
    BOOST_CHECK_EQUAL( Line1Data->find( general::id )->second->getRaw( ), "2011542    " );
    BOOST_CHECK_EQUAL( Line1Data->find( state::semiMajorAxis )->second->getRaw( ),
                       "      3.9501468   " );
    BOOST_CHECK_EQUAL( Line1Data->find( state::eccentricity )->second->getRaw( ),
                       "    0.2391642     " );
    BOOST_CHECK_EQUAL( Line1Data->find( state::inclination )->second->getRaw( ),
                       "    6.87574    " );
    BOOST_CHECK_EQUAL( Line1Data->find( state::longitudeOfAscendingNode )->second->getRaw( ),
                       "    16.88982    " );
    BOOST_CHECK_EQUAL( Line1Data->find( state::argumentOfPeriapsis )->second->getRaw( ),
                       "     48.9603   " );
    BOOST_CHECK_EQUAL( Line1Data->find( state::meanAnomaly )->second->getRaw( ),
                       "    229.49648    " );
    BOOST_CHECK_EQUAL( Line1Data->find( time::epoch )->second->getRaw( ), "   54000" );

    // Check if it fails to find a field that is not passed in the constructor.
    BOOST_CHECK( Line1Data->find( state::trueAnomaly ) == Line1Data->end( ) );

    // Similar for line 2.
    ParsedDataLineMapPtr Line2Data = testResult->at( 1 );

    // Check that it parsed 9 fields.
    BOOST_CHECK_EQUAL( Line2Data->size( ), 8 );

    // Check if the data was correcly separated.
    BOOST_CHECK_EQUAL( Line2Data->find( general::id )->second->getRaw( ), "2001038    " );
    BOOST_CHECK_EQUAL( Line2Data->find( state::semiMajorAxis )->second->getRaw( ),
                       "      3.9619932   " );
    BOOST_CHECK_EQUAL( Line2Data->find( state::eccentricity )->second->getRaw( ),
                       "     0.227483     " );
    BOOST_CHECK_EQUAL( Line2Data->find( state::inclination )->second->getRaw( ),
                       "    9.22988    " );
    BOOST_CHECK_EQUAL( Line2Data->find( state::longitudeOfAscendingNode )->second->getRaw( ),
                       "    58.20488    " );
    BOOST_CHECK_EQUAL( Line2Data->find( state::argumentOfPeriapsis )->second->getRaw( ),
                       "   307.19412   " );
    BOOST_CHECK_EQUAL( Line2Data->find( state::meanAnomaly )->second->getRaw( ),
                       "    201.38868    " );
    BOOST_CHECK_EQUAL( Line2Data->find( time::epoch )->second->getRaw( ), ( "   54000" ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
