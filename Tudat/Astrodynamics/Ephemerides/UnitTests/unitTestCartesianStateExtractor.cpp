/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      NASA/JPL HORIZONS Web-Interface, http://ssd.jpl.nasa.gov/horizons.cgi [1].
 *
 */

#define BOOST_TEST_MAIN

#include <stdexcept>
#include <string>
#include <utility>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Ephemerides/cartesianStateExtractor.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace unit_tests
{

// Short-hand notation.
namespace parsed_data_vector_utilities = input_output::parsed_data_vector_utilities;
using Eigen::Vector6d;

BOOST_AUTO_TEST_SUITE( test_cartesian_state_extractor )

//! Test the extract function.
BOOST_AUTO_TEST_CASE( cartesianStateExtractor_Extract )
{
    // Using declaration.
    using namespace ephemerides;
    using namespace orbital_element_conversions;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testYCoordinate = "-138702420.4967065";
    std::string testZCoordinate = "2356.083772440332";
    std::string testXVelocity = "26.73418833094654";
    std::string testYVelocity = "-12.25189877608797";
    std::string testZVelocity = "0.0005017881278709447";

    // Convert strings to doubles.
    const double expectedXCoordinate = std::stod( testXCoordinate );
    const double expectedYCoordinate = std::stod( testYCoordinate );
    const double expectedZCoordinate = std::stod( testZCoordinate );
    const double expectedXVelocity = std::stod( testXVelocity );
    const double expectedYVelocity = std::stod( testYVelocity );
    const double expectedZVelocity = std::stod( testZVelocity );

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                              testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                              testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                              testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new input_output::FieldValue( field_types::state::cartesianXVelocity,
                                              testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new input_output::FieldValue( field_types::state::cartesianYVelocity,
                                              testYVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new input_output::FieldValue( field_types::state::cartesianZVelocity,
                                              testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            std::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
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
    std::shared_ptr< Vector6d > returnedCartesianElements =
            testCartesianStateExtractor.extract( testDataMap );

    // Verify that the returned value corresponds to the expected value.
    BOOST_CHECK_EQUAL( ( *returnedCartesianElements )( xCartesianPositionIndex ),
                       expectedXCoordinate );
    BOOST_CHECK_EQUAL( ( *returnedCartesianElements )( yCartesianPositionIndex ),
                       expectedYCoordinate );
    BOOST_CHECK_EQUAL( ( *returnedCartesianElements )( zCartesianPositionIndex ),
                       expectedZCoordinate );
    BOOST_CHECK_EQUAL( ( *returnedCartesianElements )( xCartesianVelocityIndex ),
                       expectedXVelocity );
    BOOST_CHECK_EQUAL( ( *returnedCartesianElements )( yCartesianVelocityIndex ),
                       expectedYVelocity );
    BOOST_CHECK_EQUAL( ( *returnedCartesianElements )( zCartesianVelocityIndex ),
                       expectedZVelocity );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( cartesianStateExtractor_MissingX )
{
    // Using declaration.
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testYCoordinate = "-138702420.4967065";
    std::string testZCoordinate = "2356.083772440332";
    std::string testXVelocity = "26.73418833094654";
    std::string testYVelocity = "-12.25189877608797";
    std::string testZVelocity = "0.0005017881278709447";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                                     testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                                     testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new input_output::FieldValue( field_types::state::cartesianXVelocity,
                                                     testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new input_output::FieldValue( field_types::state::cartesianYVelocity,
                                                     testYVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new input_output::FieldValue( field_types::state::cartesianZVelocity,
                                                     testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            std::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
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
        std::shared_ptr< Vector6d > returnedCartesianElements =
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
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testZCoordinate = "2356.083772440332";
    std::string testXVelocity = "26.73418833094654";
    std::string testYVelocity = "-12.25189877608797";
    std::string testZVelocity = "0.0005017881278709447";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                                     testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                                     testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new input_output::FieldValue( field_types::state::cartesianXVelocity,
                                                     testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new input_output::FieldValue( field_types::state::cartesianYVelocity,
                                                     testYVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new input_output::FieldValue( field_types::state::cartesianZVelocity,
                                                     testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            std::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
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
        std::shared_ptr< Vector6d > returnedCartesianElements =
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
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testYCoordinate = "-138702420.4967065";
    std::string testXVelocity = "26.73418833094654";
    std::string testYVelocity = "-12.25189877608797";
    std::string testZVelocity = "0.0005017881278709447";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                                     testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                                     testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new input_output::FieldValue( field_types::state::cartesianXVelocity,
                                                     testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new input_output::FieldValue( field_types::state::cartesianYVelocity,
                                                     testYVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new input_output::FieldValue( field_types::state::cartesianZVelocity,
                                                     testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            std::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
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
        std::shared_ptr< Vector6d > returnedCartesianElements =
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
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testYCoordinate = "-138702420.4967065";
    std::string testZCoordinate = "2356.083772440332";
    std::string testYVelocity = "-12.25189877608797";
    std::string testZVelocity = "0.0005017881278709447";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                                     testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                                     testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                                     testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new input_output::FieldValue( field_types::state::cartesianYVelocity,
                                                     testYVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new input_output::FieldValue( field_types::state::cartesianZVelocity,
                                                     testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            std::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
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
        std::shared_ptr< Vector6d > returnedCartesianElements =
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
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testYCoordinate = "-138702420.4967065";
    std::string testZCoordinate = "2356.083772440332";
    std::string testXVelocity = "26.73418833094654";
    std::string testZVelocity = "0.0005017881278709447";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                                     testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                                     testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                                     testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new input_output::FieldValue( field_types::state::cartesianXVelocity,
                                                     testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZVelocity(
                new input_output::FieldValue( field_types::state::cartesianZVelocity,
                                                     testZVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            std::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
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
        std::shared_ptr< Vector6d > returnedCartesianElements =
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
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth position and velocity at JD = 2456074.5 in SI units [1].
    std::string testXCoordinate = "-62163980.78633109";
    std::string testYCoordinate = "-138702420.4967065";
    std::string testZCoordinate = "2356.083772440332";
    std::string testXVelocity = "26.73418833094654";
    std::string testYVelocity = "-12.25189877608797";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXCoordinate(
                new input_output::FieldValue( field_types::state::cartesianXCoordinate,
                                                     testXCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYCoordinate(
                new input_output::FieldValue( field_types::state::cartesianYCoordinate,
                                                     testYCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueZCoordinate(
                new input_output::FieldValue( field_types::state::cartesianZCoordinate,
                                                     testZCoordinate ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueXVelocity(
                new input_output::FieldValue( field_types::state::cartesianXVelocity,
                                                     testXVelocity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueYVelocity(
                new input_output::FieldValue( field_types::state::cartesianYVelocity,
                                                     testYVelocity ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            std::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
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
        std::shared_ptr< Vector6d > returnedCartesianElements =
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
