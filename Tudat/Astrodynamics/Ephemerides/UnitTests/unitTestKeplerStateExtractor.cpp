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

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerStateExtractor.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace unit_tests
{

// Short-hand notation.
namespace parsed_data_vector_utilities = input_output::parsed_data_vector_utilities;
using Eigen::Vector6d;

BOOST_AUTO_TEST_SUITE( test_keplerStateExtractor )

//! Test the extract function.
BOOST_AUTO_TEST_CASE( keplerStateExtractor_Extract )
{
    // Using declaration.
    using namespace ephemerides;
    using namespace orbital_element_conversions;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
            parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth orbital elements at JD = 2456074.5 in SI units [1].
    std::string testSemiMajorAxis = "149641767.7265875";
    std::string testEccentricity = "0.01625818315929578";
    std::string testInclination = "0.000038294684687687365297820321195993";
    std::string testLongitudeOfAscendingNode = "3.3677017965694967895481414920601";
    std::string testArgumentOfPeriapsis = "4.7066684303635934184870904755709";
    std::string testTrueAnomaly = "2.5013767559662625175504405244245";

    // Convert strings to doubles.
    const double expectedSemiMajorAxis = std::stod( testSemiMajorAxis );
    const double expectedEccentricity = std::stod( testEccentricity );
    const double expectedInclination = std::stod( testInclination );
    const double expectedLongitudeOfAscendingNode = std::stod(
                testLongitudeOfAscendingNode );
    const double expectedArgumentOfPeriapsis = std::stod(
                testArgumentOfPeriapsis );
    const double expectedTrueAnomaly = std::stod( testTrueAnomaly );

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueSemiMajorAxis(
                new input_output::FieldValue( field_types::state::semiMajorAxis,
                                              testSemiMajorAxis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEccentricity(
                new input_output::FieldValue( field_types::state::eccentricity,
                                              testEccentricity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueInclination(
                new input_output::FieldValue( field_types::state::inclination,
                                              testInclination ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueLongitudeOfAscendingNode(
                new input_output::FieldValue( field_types::state::longitudeOfAscendingNode,
                                              testLongitudeOfAscendingNode ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueArgumentOfPeriapsis(
                new input_output::FieldValue( field_types::state::argumentOfPeriapsis,
                                              testArgumentOfPeriapsis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueTrueAnomaly(
                new input_output::FieldValue( field_types::state::trueAnomaly,
                                              testTrueAnomaly ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
                parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::state::semiMajorAxis,
                                        testFieldValueSemiMajorAxis ) );
    testDataMap->insert( FieldDataPair( field_types::state::eccentricity,
                                        testFieldValueEccentricity ) );
    testDataMap->insert( FieldDataPair( field_types::state::inclination,
                                        testFieldValueInclination ) );
    testDataMap->insert( FieldDataPair( field_types::state::longitudeOfAscendingNode,
                                        testFieldValueLongitudeOfAscendingNode ) );
    testDataMap->insert( FieldDataPair( field_types::state::argumentOfPeriapsis,
                                        testFieldValueArgumentOfPeriapsis ) );
    testDataMap->insert( FieldDataPair( field_types::state::trueAnomaly,
                                        testFieldValueTrueAnomaly ) );

    // Create Cartesian state extractor.
    KeplerStateExtractor testKeplerStateExtractor;

    // Extract test data map.
    boost::shared_ptr< Vector6d > returnedKeplerianElements =
            testKeplerStateExtractor.extract( testDataMap );

    // Verify that the returned values correspond to the expected values.
    BOOST_CHECK_EQUAL( ( *returnedKeplerianElements )( semiMajorAxisIndex ),
                       expectedSemiMajorAxis );
    BOOST_CHECK_EQUAL( ( *returnedKeplerianElements )( eccentricityIndex ),
                       expectedEccentricity );
    BOOST_CHECK_EQUAL( ( *returnedKeplerianElements )( inclinationIndex ),
                       expectedInclination );
    BOOST_CHECK_EQUAL( ( *returnedKeplerianElements )( argumentOfPeriapsisIndex ),
                       expectedArgumentOfPeriapsis );
    BOOST_CHECK_EQUAL( ( *returnedKeplerianElements )( longitudeOfAscendingNodeIndex ),
                       expectedLongitudeOfAscendingNode );
    BOOST_CHECK_EQUAL( ( *returnedKeplerianElements )( trueAnomalyIndex ),
                       expectedTrueAnomaly );
}

//! Test the extract function with mean anomaly as input.
BOOST_AUTO_TEST_CASE( keplerStateExtractor_ExtractWithMeanAnomaly )
{
    // Using declaration.
    using namespace ephemerides;
    using namespace orbital_element_conversions;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
            parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth orbital elements at JD = 2456074.5 in SI units [1].
    std::string testSemiMajorAxis = "149641767.7265875";
    std::string testEccentricity = "0.01625818315929578";
    std::string testInclination = "0.000038294684687687365297820321195993";
    std::string testLongitudeOfAscendingNode = "3.3677017965694967895481414920601";
    std::string testArgumentOfPeriapsis = "4.7066684303635934184870904755709";
    std::string testMeanAnomaly = "2.4817611918827110773177514320574";

    // Convert strings to doubles.
    double expectedSemiMajorAxis = std::stod( testSemiMajorAxis );
    double expectedEccentricity = std::stod( testEccentricity );
    double expectedInclination = std::stod( testInclination );
    double expectedLongitudeOfAscendingNode = std::stod(
                testLongitudeOfAscendingNode );
    double expectedArgumentOfPeriapsis = std::stod( testArgumentOfPeriapsis );

    // Set expected true anomaly.
    double expectedTrueAnomaly = 2.5013767559662625175504405244245;

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueSemiMajorAxis(
                new input_output::FieldValue( field_types::state::semiMajorAxis,
                                              testSemiMajorAxis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEccentricity(
                new input_output::FieldValue( field_types::state::eccentricity,
                                              testEccentricity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueInclination(
                new input_output::FieldValue( field_types::state::inclination,
                                              testInclination ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueLongitudeOfAscendingNode(
                new input_output::FieldValue( field_types::state::longitudeOfAscendingNode,
                                              testLongitudeOfAscendingNode ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueArgumentOfPeriapsis(
                new input_output::FieldValue( field_types::state::argumentOfPeriapsis,
                                              testArgumentOfPeriapsis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueMeanAnomaly(
                new input_output::FieldValue( field_types::state::meanAnomaly,
                                              testMeanAnomaly ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
                parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::state::semiMajorAxis,
                                        testFieldValueSemiMajorAxis ) );
    testDataMap->insert( FieldDataPair( field_types::state::eccentricity,
                                        testFieldValueEccentricity ) );
    testDataMap->insert( FieldDataPair( field_types::state::inclination,
                                        testFieldValueInclination ) );
    testDataMap->insert( FieldDataPair( field_types::state::longitudeOfAscendingNode,
                                        testFieldValueLongitudeOfAscendingNode ) );
    testDataMap->insert( FieldDataPair( field_types::state::argumentOfPeriapsis,
                                        testFieldValueArgumentOfPeriapsis ) );
    testDataMap->insert( FieldDataPair( field_types::state::meanAnomaly,
                                        testFieldValueMeanAnomaly ) );

    // Create Cartesian state extractor.
    KeplerStateExtractor testKeplerStateExtractor;

    // Extract test data map.
    boost::shared_ptr< Vector6d > returnedKeplerianElements =
            testKeplerStateExtractor.extract( testDataMap );

    // Verify that the returned values correspond to the expected values.
    BOOST_CHECK_EQUAL( ( *returnedKeplerianElements )( semiMajorAxisIndex ),
                       expectedSemiMajorAxis );
    BOOST_CHECK_EQUAL( ( *returnedKeplerianElements )( eccentricityIndex ),
                       expectedEccentricity );
    BOOST_CHECK_EQUAL( ( *returnedKeplerianElements )( inclinationIndex ),
                       expectedInclination );
    BOOST_CHECK_EQUAL( ( *returnedKeplerianElements )( argumentOfPeriapsisIndex ),
                       expectedArgumentOfPeriapsis );
    BOOST_CHECK_EQUAL( ( *returnedKeplerianElements )( longitudeOfAscendingNodeIndex ),
                       expectedLongitudeOfAscendingNode );

    BOOST_CHECK_CLOSE_FRACTION( ( *returnedKeplerianElements )( trueAnomalyIndex ),
                                expectedTrueAnomaly,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( keplerStateExtractor_MissingSemiMajorAxis )
{
    // Using declaration.
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
            parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth orbital elements at JD = 2456074.5 in SI units [1].
    std::string testEccentricity = "0.01625818315929578";
    std::string testInclination = "0.000038294684687687365297820321195993";
    std::string testLongitudeOfAscendingNode = "3.3677017965694967895481414920601";
    std::string testArgumentOfPeriapsis = "4.7066684303635934184870904755709";
    std::string testTrueAnomaly = "2.4817611918827110773177514320574";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEccentricity(
                new input_output::FieldValue( field_types::state::eccentricity,
                                              testEccentricity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueInclination(
                new input_output::FieldValue( field_types::state::inclination,
                                              testInclination ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueLongitudeOfAscendingNode(
                new input_output::FieldValue( field_types::state::longitudeOfAscendingNode,
                                              testLongitudeOfAscendingNode ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueArgumentOfPeriapsis(
                new input_output::FieldValue( field_types::state::argumentOfPeriapsis,
                                              testArgumentOfPeriapsis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueTrueAnomaly(
                new input_output::FieldValue( field_types::state::trueAnomaly,
                                              testTrueAnomaly ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
                parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::state::eccentricity,
                                        testFieldValueEccentricity ) );
    testDataMap->insert( FieldDataPair( field_types::state::inclination,
                                        testFieldValueInclination ) );
    testDataMap->insert( FieldDataPair( field_types::state::longitudeOfAscendingNode,
                                        testFieldValueLongitudeOfAscendingNode ) );
    testDataMap->insert( FieldDataPair( field_types::state::argumentOfPeriapsis,
                                        testFieldValueArgumentOfPeriapsis ) );
    testDataMap->insert( FieldDataPair( field_types::state::trueAnomaly,
                                        testFieldValueTrueAnomaly ) );

    // Create Cartesian state extractor.
    KeplerStateExtractor testKeplerStateExtractor;

    // Set flag.
    bool isSemiMajorAxisFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< Vector6d > returnedKeplerianElements =
                testKeplerStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isSemiMajorAxisFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isSemiMajorAxisFound );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( keplerStateExtractor_MissingEccentricity )
{
    // Using declaration.
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth orbital elements at JD = 2456074.5 in SI units [1].
    std::string testSemiMajorAxis = "149641767.7265875";
    std::string testInclination = "0.000038294684687687365297820321195993";
    std::string testLongitudeOfAscendingNode = "3.3677017965694967895481414920601";
    std::string testArgumentOfPeriapsis = "4.7066684303635934184870904755709";
    std::string testTrueAnomaly = "2.4817611918827110773177514320574";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueSemiMajorAxis(
                new input_output::FieldValue( field_types::state::semiMajorAxis,
                                                     testSemiMajorAxis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueInclination(
                new input_output::FieldValue( field_types::state::inclination,
                                                     testInclination ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueLongitudeOfAscendingNode(
                new input_output::FieldValue( field_types::state::longitudeOfAscendingNode,
                                                     testLongitudeOfAscendingNode ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueArgumentOfPeriapsis(
                new input_output::FieldValue( field_types::state::argumentOfPeriapsis,
                                                     testArgumentOfPeriapsis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueTrueAnomaly(
                new input_output::FieldValue( field_types::state::trueAnomaly,
                                                     testTrueAnomaly ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::state::semiMajorAxis,
                                        testFieldValueSemiMajorAxis ) );
    testDataMap->insert( FieldDataPair( field_types::state::inclination,
                                        testFieldValueInclination ) );
    testDataMap->insert( FieldDataPair( field_types::state::longitudeOfAscendingNode,
                                        testFieldValueLongitudeOfAscendingNode ) );
    testDataMap->insert( FieldDataPair( field_types::state::argumentOfPeriapsis,
                                        testFieldValueArgumentOfPeriapsis ) );
    testDataMap->insert( FieldDataPair( field_types::state::trueAnomaly,
                                        testFieldValueTrueAnomaly ) );

    // Create Cartesian state extractor.
    KeplerStateExtractor testKeplerStateExtractor;

    // Set flag.
    bool isEccentricityFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< Vector6d > returnedKeplerianElements =
                testKeplerStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isEccentricityFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isEccentricityFound );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( keplerStateExtractor_MissingInclination )
{
    // Using declaration.
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth orbital elements at JD = 2456074.5 in SI units [1].
    std::string testSemiMajorAxis = "149641767.7265875";
    std::string testEccentricity = "0.01625818315929578";
    std::string testLongitudeOfAscendingNode = "3.3677017965694967895481414920601";
    std::string testArgumentOfPeriapsis = "4.7066684303635934184870904755709";
    std::string testTrueAnomaly = "2.4817611918827110773177514320574";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueSemiMajorAxis(
                new input_output::FieldValue( field_types::state::semiMajorAxis,
                                                     testSemiMajorAxis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEccentricity(
                new input_output::FieldValue( field_types::state::eccentricity,
                                                     testEccentricity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueLongitudeOfAscendingNode(
                new input_output::FieldValue( field_types::state::longitudeOfAscendingNode,
                                                     testLongitudeOfAscendingNode ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueArgumentOfPeriapsis(
                new input_output::FieldValue( field_types::state::argumentOfPeriapsis,
                                                     testArgumentOfPeriapsis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueTrueAnomaly(
                new input_output::FieldValue( field_types::state::trueAnomaly,
                                                     testTrueAnomaly ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::state::semiMajorAxis,
                                        testFieldValueSemiMajorAxis ) );
    testDataMap->insert( FieldDataPair( field_types::state::eccentricity,
                                        testFieldValueEccentricity ) );
    testDataMap->insert( FieldDataPair( field_types::state::longitudeOfAscendingNode,
                                        testFieldValueLongitudeOfAscendingNode ) );
    testDataMap->insert( FieldDataPair( field_types::state::argumentOfPeriapsis,
                                        testFieldValueArgumentOfPeriapsis ) );
    testDataMap->insert( FieldDataPair( field_types::state::trueAnomaly,
                                        testFieldValueTrueAnomaly ) );

    // Create Cartesian state extractor.
    KeplerStateExtractor testKeplerStateExtractor;

    // Set flag.
    bool isInclinationFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< Vector6d > returnedKeplerianElements =
                testKeplerStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isInclinationFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isInclinationFound );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( keplerStateExtractor_MissingLongitudeOfAscendingNode )
{
    // Using declaration.
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth orbital elements at JD = 2456074.5 in SI units [1].
    std::string testSemiMajorAxis = "149641767.7265875";
    std::string testEccentricity = "0.01625818315929578";
    std::string testInclination = "0.000038294684687687365297820321195993";
    std::string testArgumentOfPeriapsis = "4.7066684303635934184870904755709";
    std::string testTrueAnomaly = "2.4817611918827110773177514320574";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueSemiMajorAxis(
                new input_output::FieldValue( field_types::state::semiMajorAxis,
                                                     testSemiMajorAxis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEccentricity(
                new input_output::FieldValue( field_types::state::eccentricity,
                                                     testEccentricity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueInclination(
                new input_output::FieldValue( field_types::state::inclination,
                                                     testInclination ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueArgumentOfPeriapsis(
                new input_output::FieldValue( field_types::state::argumentOfPeriapsis,
                                                     testArgumentOfPeriapsis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueTrueAnomaly(
                new input_output::FieldValue( field_types::state::trueAnomaly,
                                                     testTrueAnomaly ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::state::semiMajorAxis,
                                        testFieldValueSemiMajorAxis ) );
    testDataMap->insert( FieldDataPair( field_types::state::eccentricity,
                                        testFieldValueEccentricity ) );
    testDataMap->insert( FieldDataPair( field_types::state::inclination,
                                        testFieldValueInclination ) );
    testDataMap->insert( FieldDataPair( field_types::state::argumentOfPeriapsis,
                                        testFieldValueArgumentOfPeriapsis ) );
    testDataMap->insert( FieldDataPair( field_types::state::trueAnomaly,
                                        testFieldValueTrueAnomaly ) );

    // Create Cartesian state extractor.
    KeplerStateExtractor testKeplerStateExtractor;

    // Set flag.
    bool isLongitudeOfAscendingNodeFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< Vector6d > returnedKeplerianElements =
                testKeplerStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isLongitudeOfAscendingNodeFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isLongitudeOfAscendingNodeFound );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( keplerStateExtractor_MissingArgumentOfPeriapsis )
{
    // Using declaration.
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth orbital elements at JD = 2456074.5 in SI units [1].
    std::string testSemiMajorAxis = "149641767.7265875";
    std::string testEccentricity = "0.01625818315929578";
    std::string testInclination = "0.000038294684687687365297820321195993";
    std::string testLongitudeOfAscendingNode = "3.3677017965694967895481414920601";
    std::string testTrueAnomaly = "2.4817611918827110773177514320574";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueSemiMajorAxis(
                new input_output::FieldValue( field_types::state::semiMajorAxis,
                                                     testSemiMajorAxis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEccentricity(
                new input_output::FieldValue( field_types::state::eccentricity,
                                                     testEccentricity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueInclination(
                new input_output::FieldValue( field_types::state::inclination,
                                                     testInclination ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueLongitudeOfAscendingNode(
                new input_output::FieldValue( field_types::state::longitudeOfAscendingNode,
                                                     testLongitudeOfAscendingNode ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueTrueAnomaly(
                new input_output::FieldValue( field_types::state::trueAnomaly,
                                                     testTrueAnomaly ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::state::semiMajorAxis,
                                        testFieldValueSemiMajorAxis ) );
    testDataMap->insert( FieldDataPair( field_types::state::eccentricity,
                                        testFieldValueEccentricity ) );
    testDataMap->insert( FieldDataPair( field_types::state::inclination,
                                        testFieldValueInclination ) );
    testDataMap->insert( FieldDataPair( field_types::state::longitudeOfAscendingNode,
                                        testFieldValueLongitudeOfAscendingNode ) );
    testDataMap->insert( FieldDataPair( field_types::state::trueAnomaly,
                                        testFieldValueTrueAnomaly ) );

    // Create Cartesian state extractor.
    KeplerStateExtractor testKeplerStateExtractor;

    // Set flag.
    bool isArgumentOfPeriapsisFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< Vector6d > returnedKeplerianElements =
                testKeplerStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isArgumentOfPeriapsisFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isArgumentOfPeriapsisFound );
}

//! Test if the extract function throws the necessary exceptions.
BOOST_AUTO_TEST_CASE( keplerStateExtractor_MissingTrueOrMeanAnomaly )
{
    // Using declaration.
    using namespace ephemerides;
    namespace field_types = input_output::field_types;

    // Create parsed data line map pointer.
    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< input_output::FieldType,
                       parsed_data_vector_utilities::FieldValuePtr > FieldDataPair;

    // Create test strings, based on Earth orbital elements at JD = 2456074.5 in SI units [1].
    std::string testSemiMajorAxis = "149641767.7265875";
    std::string testEccentricity = "0.01625818315929578";
    std::string testInclination = "0.000038294684687687365297820321195993";
    std::string testLongitudeOfAscendingNode = "3.3677017965694967895481414920601";
    std::string testArgumentOfPeriapsis = "4.7066684303635934184870904755709";

    // Store strings as field values.
    parsed_data_vector_utilities::FieldValuePtr testFieldValueSemiMajorAxis(
                new input_output::FieldValue( field_types::state::semiMajorAxis,
                                                     testSemiMajorAxis ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueEccentricity(
                new input_output::FieldValue( field_types::state::eccentricity,
                                                     testEccentricity ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueInclination(
                new input_output::FieldValue( field_types::state::inclination,
                                                     testInclination ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueLongitudeOfAscendingNode(
                new input_output::FieldValue( field_types::state::longitudeOfAscendingNode,
                                                     testLongitudeOfAscendingNode ) );
    parsed_data_vector_utilities::FieldValuePtr testFieldValueArgumentOfPeriapsis(
                new input_output::FieldValue( field_types::state::argumentOfPeriapsis,
                                                     testArgumentOfPeriapsis ) );

    // Create a new pointer to data map.
    parsed_data_vector_utilities::ParsedDataLineMapPtr testDataMap =
            boost::make_shared< parsed_data_vector_utilities::ParsedDataLineMap >(
                std::map< input_output::FieldType,
                          parsed_data_vector_utilities::FieldValuePtr >( ) );

    // Store field values in data map.
    testDataMap->insert( FieldDataPair( field_types::state::semiMajorAxis,
                                        testFieldValueSemiMajorAxis ) );
    testDataMap->insert( FieldDataPair( field_types::state::eccentricity,
                                        testFieldValueEccentricity ) );
    testDataMap->insert( FieldDataPair( field_types::state::inclination,
                                        testFieldValueInclination ) );
    testDataMap->insert( FieldDataPair( field_types::state::longitudeOfAscendingNode,
                                        testFieldValueLongitudeOfAscendingNode ) );
    testDataMap->insert( FieldDataPair( field_types::state::argumentOfPeriapsis,
                                        testFieldValueArgumentOfPeriapsis ) );

    // Create Cartesian state extractor.
    KeplerStateExtractor testKeplerStateExtractor;

    // Set flag.
    bool isTrueAnomalyFound = true;

    // Try to extract, which should result in a runtime error.
    try
    {
        // Extract test data map.
        boost::shared_ptr< Vector6d > returnedKeplerianElements =
                testKeplerStateExtractor.extract( testDataMap );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isTrueAnomalyFound = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isTrueAnomalyFound );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
