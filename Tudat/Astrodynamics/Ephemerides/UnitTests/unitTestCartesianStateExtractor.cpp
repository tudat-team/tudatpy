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
 *      120607    T. Secretin       First creation of code.
 *
 *    References
 *      NASA/JPL HORIZONS Web-Interface, http://ssd.jpl.nasa.gov/horizons.cgi [1].
 *
 *     Notes
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Ephemerides/cartesianStateExtractor.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace unit_tests
{

// Short-hand notation.
namespace parsed_data_vector_utilities = tudat::input_output::parsed_data_vector_utilities;

BOOST_AUTO_TEST_SUITE( test_cartesian_state_extractor )

//! Test the extract function.
BOOST_AUTO_TEST_CASE( cartesianStateExtractor_Extract )
{
    // Using declaration.
    using namespace tudat::ephemerides;
    namespace field_types = tudat::input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< tudat::input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testYCoordinate = "-138702420.4967065";
    std::string testZCoordinate = "2356.083772440332";
    std::string testXVelocity = "26.73418833094654";
    std::string testYVelocity = "-12.25189877608797";
    std::string testZVelocity = "0.0005017881278709447";

    // Convert strings to doubles.
    double expectedXCoordinate = boost::lexical_cast< double >( testXCoordinate );
    double expectedYCoordinate = boost::lexical_cast< double >( testYCoordinate );
    double expectedZCoordinate = boost::lexical_cast< double >( testZCoordinate );
    double expectedXVelocity = boost::lexical_cast< double >( testXVelocity );
    double expectedYVelocity = boost::lexical_cast< double >( testYVelocity );
    double expectedZVelocity = boost::lexical_cast< double >( testZVelocity );

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                                     testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                                     testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                                     testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianXVelocity,
                                                     testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianYVelocity,
                                                     testYVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianZVelocity,
                                                     testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< tudat::input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXCoordinate,
                                        testFieldValueXCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYCoordinate,
                                        testFieldValueYCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZCoordinate,
                                        testFieldValueZCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXVelocity,
                                        testFieldValueXVelocity ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYVelocity,
                                        testFieldValueYVelocity ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZVelocity,
                                        testFieldValueZVelocity ) );

    // Create Cartesian state extractor.
    CartesianStateExtractor testCartesianStateExtractor;

    // Extract test data map.
    boost::shared_ptr< CartesianElements > returnedCartesianElements =
            testCartesianStateExtractor.extract( testDataMap );

    // Verify that the returned value corresponds to the expected value.
    BOOST_CHECK_EQUAL( returnedCartesianElements->getCartesianElementX( ), expectedXCoordinate );
    BOOST_CHECK_EQUAL( returnedCartesianElements->getCartesianElementY( ), expectedYCoordinate );
    BOOST_CHECK_EQUAL( returnedCartesianElements->getCartesianElementZ( ), expectedZCoordinate );
    BOOST_CHECK_EQUAL( returnedCartesianElements->getCartesianElementXDot( ), expectedXVelocity );
    BOOST_CHECK_EQUAL( returnedCartesianElements->getCartesianElementYDot( ), expectedYVelocity );
    BOOST_CHECK_EQUAL( returnedCartesianElements->getCartesianElementZDot( ), expectedZVelocity );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( cartesianStateExtractor_MissingX )
{
    // Using declaration.
    using namespace tudat::ephemerides;
    namespace field_types = tudat::input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< tudat::input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testYCoordinate = "-138702420.4967065";
    std::string testZCoordinate = "2356.083772440332";
    std::string testXVelocity = "26.73418833094654";
    std::string testYVelocity = "-12.25189877608797";
    std::string testZVelocity = "0.0005017881278709447";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                                     testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                                     testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianXVelocity,
                                                     testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianYVelocity,
                                                     testYVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianZVelocity,
                                                     testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< tudat::input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map, except x-coordinate.
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYCoordinate,
                                        testFieldValueYCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZCoordinate,
                                        testFieldValueZCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXVelocity,
                                        testFieldValueXVelocity ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYVelocity,
                                        testFieldValueYVelocity ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZVelocity,
                                        testFieldValueZVelocity ) );

    // Create Cartesian state extractor.
    CartesianStateExtractor testCartesianStateExtractor;

    // Set flag.
    bool isXCoordinateFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< CartesianElements > returnedCartesianElements =
                testCartesianStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isXCoordinateFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isXCoordinateFound );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( cartesianStateExtractor_MissingY )
{
    // Using declaration.
    using namespace tudat::ephemerides;
    namespace field_types = tudat::input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< tudat::input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testZCoordinate = "2356.083772440332";
    std::string testXVelocity = "26.73418833094654";
    std::string testYVelocity = "-12.25189877608797";
    std::string testZVelocity = "0.0005017881278709447";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                                     testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                                     testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianXVelocity,
                                                     testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianYVelocity,
                                                     testYVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianZVelocity,
                                                     testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< tudat::input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map, except x-coordinate.
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXCoordinate,
                                        testFieldValueXCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZCoordinate,
                                        testFieldValueZCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXVelocity,
                                        testFieldValueXVelocity ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYVelocity,
                                        testFieldValueYVelocity ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZVelocity,
                                        testFieldValueZVelocity ) );

    // Create Cartesian state extractor.
    CartesianStateExtractor testCartesianStateExtractor;

    // Set flag.
    bool isYCoordinateFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< CartesianElements > returnedCartesianElements =
                testCartesianStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isYCoordinateFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isYCoordinateFound );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( cartesianStateExtractor_MissingZ )
{
    // Using declaration.
    using namespace tudat::ephemerides;
    namespace field_types = tudat::input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< tudat::input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testYCoordinate = "-138702420.4967065";
    std::string testXVelocity = "26.73418833094654";
    std::string testYVelocity = "-12.25189877608797";
    std::string testZVelocity = "0.0005017881278709447";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                                     testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                                     testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianXVelocity,
                                                     testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianYVelocity,
                                                     testYVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianZVelocity,
                                                     testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< tudat::input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map, except x-coordinate.
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXCoordinate,
                                        testFieldValueXCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYCoordinate,
                                        testFieldValueYCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXVelocity,
                                        testFieldValueXVelocity ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYVelocity,
                                        testFieldValueYVelocity ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZVelocity,
                                        testFieldValueZVelocity ) );

    // Create Cartesian state extractor.
    CartesianStateExtractor testCartesianStateExtractor;

    // Set flag.
    bool isZCoordinateFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< CartesianElements > returnedCartesianElements =
                testCartesianStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isZCoordinateFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isZCoordinateFound );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( cartesianStateExtractor_MissingXDot )
{
    // Using declaration.
    using namespace tudat::ephemerides;
    namespace field_types = tudat::input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< tudat::input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testYCoordinate = "-138702420.4967065";
    std::string testZCoordinate = "2356.083772440332";
    std::string testYVelocity = "-12.25189877608797";
    std::string testZVelocity = "0.0005017881278709447";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                                     testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                                     testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                                     testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianYVelocity,
                                                     testYVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianZVelocity,
                                                     testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< tudat::input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map, except x-coordinate.
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXCoordinate,
                                        testFieldValueXCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYCoordinate,
                                        testFieldValueYCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZCoordinate,
                                        testFieldValueZCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYVelocity,
                                        testFieldValueYVelocity ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZVelocity,
                                        testFieldValueZVelocity ) );

    // Create Cartesian state extractor.
    CartesianStateExtractor testCartesianStateExtractor;

    // Set flag.
    bool isXVelocityFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< CartesianElements > returnedCartesianElements =
                testCartesianStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isXVelocityFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isXVelocityFound );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( cartesianStateExtractor_MissingYDot )
{
    // Using declaration.
    using namespace tudat::ephemerides;
    namespace field_types = tudat::input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< tudat::input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testYCoordinate = "-138702420.4967065";
    std::string testZCoordinate = "2356.083772440332";
    std::string testXVelocity = "26.73418833094654";
    std::string testZVelocity = "0.0005017881278709447";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                                     testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                                     testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                                     testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianXVelocity,
                                                     testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianZVelocity,
                                                     testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< tudat::input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map, except x-coordinate.
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXCoordinate,
                                        testFieldValueXCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYCoordinate,
                                        testFieldValueYCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZCoordinate,
                                        testFieldValueZCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXVelocity,
                                        testFieldValueXVelocity ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZVelocity,
                                        testFieldValueZVelocity ) );

    // Create Cartesian state extractor.
    CartesianStateExtractor testCartesianStateExtractor;

    // Set flag.
    bool isYVelocityFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< CartesianElements > returnedCartesianElements =
                testCartesianStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isYVelocityFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isYVelocityFound );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( cartesianStateExtractor_MissingZDot )
{
    // Using declaration.
    using namespace tudat::ephemerides;
    namespace field_types = tudat::input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< tudat::input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testYCoordinate = "-138702420.4967065";
    std::string testZCoordinate = "2356.083772440332";
    std::string testXVelocity = "26.73418833094654";
    std::string testYVelocity = "-12.25189877608797";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                                     testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                                     testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new tudat::input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                                     testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianXVelocity,
                                                     testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new tudat::input_output::FieldValue( field_types::state::cartesianYVelocity,
                                                     testYVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< tudat::input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map, except x-coordinate.
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXCoordinate,
                                        testFieldValueXCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYCoordinate,
                                        testFieldValueYCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianZCoordinate,
                                        testFieldValueZCoordinate ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianXVelocity,
                                        testFieldValueXVelocity ) );
    testDataMap->insert( FieldDataPair( field_types::state::cartesianYVelocity,
                                        testFieldValueYVelocity ) );

    // Create Cartesian state extractor.
    CartesianStateExtractor testCartesianStateExtractor;

    // Set flag.
    bool isZVelocityFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< CartesianElements > returnedCartesianElements =
                testCartesianStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isZVelocityFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isZVelocityFound );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
