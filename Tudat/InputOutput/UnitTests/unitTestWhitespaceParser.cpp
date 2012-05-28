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
 *      111103    S. Billemont      Creation of code.
 *      111209    T. Secretin       Adapted code to Boost test suite.
 *      120424    T. Secretin       Code-check, layout changes.
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "Tudat/InputOutput/fieldType.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/InputOutput/whiteSpaceParser.h"

// Add groupNumber to fieldtypes namespace.
namespace tudat
{
namespace input_output
{
namespace field_types
{
namespace general
{
    //! Group number [-].
    static const tudat::input_output::FieldType groupNumber
        = hash_constructor( "General: Group_Number" );

} // namespace general
} // namespace field_types
} // namespace input_output
} // namespace tudat

// Define Boost test suite.
BOOST_AUTO_TEST_SUITE( test_whitespace_parser )

//! Test if a single line is correctly parsed.
BOOST_AUTO_TEST_CASE( whiteSpaceParser_singleLine )
{
    // Using declaration.
    using namespace tudat::input_output::parsed_data_vector_utilities;

    // Create GTOC2 white space parser.
    tudat::input_output::WhiteSpaceParser
            testWhiteSpaceParser(
                9,
                tudat::input_output::field_types::general::id,
                tudat::input_output::field_types::state::semiMajorAxis,
                tudat::input_output::field_types::state::eccentricity,
                tudat::input_output::field_types::state::inclination,
                tudat::input_output::field_types::state::longitudeOfAscendingNode,
                tudat::input_output::field_types::state::argumentOfPeriapsis,
                tudat::input_output::field_types::state::meanAnomaly,
                tudat::input_output::field_types::time::epoch,
                tudat::input_output::field_types::general::groupNumber );

    // Create test data (First line in GTOC2 problem data file with additional whitespaces).
    std::string testString(
              "  2011542 3.9501468 0.2391642   6.87574   16.88982  48.9603 229.49648 54000 1   " );

    // Parse the data.
    ParsedDataVectorPtr testResult = testWhiteSpaceParser.parse( testString );

    // Check that only one line is parsed.
    BOOST_CHECK_EQUAL( testResult->size( ), 1);

    // Retrieve the single line.
    ParsedDataLineMapPtr testLineData = testResult->at( 0 );

    // Check that it parsed 9 fields.
    BOOST_CHECK_EQUAL( testLineData->size( ), 9 );

    // Check if the data was correcly separated.
    BOOST_CHECK_EQUAL(
        *( testLineData->find(
               tudat::input_output::field_types::general::id )->second->getRaw( ) ),
        ( "2011542" ) );
    BOOST_CHECK_EQUAL(
        *( testLineData->find(
              tudat::input_output::field_types::state::semiMajorAxis )->second->getRaw( ) ),
        ( "3.9501468" ) );
    BOOST_CHECK_EQUAL(
        *( testLineData->find(
              tudat::input_output::field_types::state::eccentricity )->second->getRaw( ) ),
        ( "0.2391642" ) );
    BOOST_CHECK_EQUAL(
        *( testLineData->find(
               tudat::input_output::field_types::state::inclination )->second->getRaw( ) ),
        ( "6.87574" ) );
    BOOST_CHECK_EQUAL(
        *(testLineData->find(
              tudat::input_output::field_types::state::
              longitudeOfAscendingNode )->second->getRaw( ) ),
        ( "16.88982" ) );
    BOOST_CHECK_EQUAL(
        *( testLineData->find(
               tudat::input_output::field_types::state::argumentOfPeriapsis )->second->getRaw( ) ),
        ( "48.9603" ) );
    BOOST_CHECK_EQUAL(
        *( testLineData->find(
               tudat::input_output::field_types::state::meanAnomaly )->second->getRaw( ) ),
        ( "229.49648" ) );
    BOOST_CHECK_EQUAL(
        *( testLineData->find(
               tudat::input_output::field_types::time::epoch )->second->getRaw( ) ),
        ( "54000" ) );
    BOOST_CHECK_EQUAL(
        *( testLineData->find(
               tudat::input_output::field_types::general::groupNumber )->second->getRaw( ) ),
        ( "1" ) );

    // Check if it fails to find a field that is not passed in the constructor.
    BOOST_CHECK(
        testLineData->find( tudat::input_output::field_types::state::trueAnomaly )
                == testLineData->end( ) );
}

// Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )
