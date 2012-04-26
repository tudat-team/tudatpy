/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      111103    S. Billemont      First creation of code.
 *      111209    T. Secretin       Adapted code to Boost test suite.
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <boost/test/unit_test.hpp>

#include "Tudat/InputOutput/fieldType.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/InputOutput/whiteSpaceParser.h"

// Add groupNumber to fieldtypes namespace
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
        = hash_constructor("General: Group_Number");

} } } }

// Define Boost test suite.
BOOST_AUTO_TEST_SUITE( test_whitespace_parser )

//! Test if a single line is correctly parsed.
BOOST_AUTO_TEST_CASE( whiteSpaceParser_singleLine )
{
    // Using declaration.
    using namespace tudat::input_output::parsed_data_vector_utilities;

    // Create GTOC2 white space parser.
    tudat::input_output::WhiteSpaceParser
            testWhiteSpaceParser(9,
                                 tudat::input_output::field_types::general::id,
                                 tudat::input_output::field_types::state::semiMajorAxis,
                                 tudat::input_output::field_types::state::eccentricity,
                                 tudat::input_output::field_types::state::inclination,
                                 tudat::input_output::field_types::state::longitudeOfAscendingNode,
                                 tudat::input_output::field_types::state::argumentOfPeriapsis,
                                 tudat::input_output::field_types::state::meanAnomaly,
                                 tudat::input_output::field_types::time::epoch,
                                 tudat::input_output::field_types::general::groupNumber);

    // Create test data (First line in GTOC2 problem data file with additional whitespaces).
    std::string testString(
                "  2011542 3.9501468 0.2391642   6.87574   16.88982  48.9603 229.49648 54000 1   ");

    // Parse the data.
    ParsedDataVectorPtr testResult = testWhiteSpaceParser.parse(testString);

    // Check that only one line is parsed.
    BOOST_CHECK_EQUAL(testResult->size(), 1);

    // Retrieve the single line.
    ParsedDataLineMapPtr testLineData = testResult->at(0);

    // Check that it parsed 9 fields.
    BOOST_CHECK_EQUAL(testLineData->size(), 9);

    // Check if the data was correcly separated.
    BOOST_CHECK_EQUAL(
        *(testLineData->find(tudat::input_output::field_types::general::id)->second->getRaw()),
        ("2011542") );
    BOOST_CHECK_EQUAL(
        *(testLineData->find(tudat::input_output::field_types::state::semiMajorAxis)->second->getRaw()),
        ("3.9501468") );
    BOOST_CHECK_EQUAL(
        *(testLineData->find(tudat::input_output::field_types::state::eccentricity)->second->getRaw()),
        ("0.2391642") );
    BOOST_CHECK_EQUAL(
        *(testLineData->find(tudat::input_output::field_types::state::inclination)->second->getRaw()),
        ("6.87574") );
    BOOST_CHECK_EQUAL(
        *(testLineData->find(tudat::input_output::field_types::state::longitudeOfAscendingNode)->second->getRaw()),
        ("16.88982") );
    BOOST_CHECK_EQUAL(
        *(testLineData->find(tudat::input_output::field_types::state::argumentOfPeriapsis)->second->getRaw()),
        ("48.9603") );
    BOOST_CHECK_EQUAL(
        *(testLineData->find(tudat::input_output::field_types::state::meanAnomaly)->second->getRaw()),
        ("229.49648") );
    BOOST_CHECK_EQUAL(
        *(testLineData->find(tudat::input_output::field_types::time::epoch)->second->getRaw()),
        ("54000") );
    BOOST_CHECK_EQUAL(
        *(testLineData->find(tudat::input_output::field_types::general::groupNumber)->second->getRaw()),
        ("1") );

    // Check if it fails to find a field that is not passed in the constructor.
    BOOST_CHECK(
        testLineData->find(tudat::input_output::field_types::state::trueAnomaly)
                == testLineData->end());
}

// Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END()
