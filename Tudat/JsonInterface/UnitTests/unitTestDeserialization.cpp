/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include "Tudat/JsonInterface/UnitTests/unitTestSupport.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"
#include "Tudat/JsonInterface/Support/deserialization.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_deserialization )

// Test 1: value access
BOOST_AUTO_TEST_CASE( test_json_valueAccess )
{
    using namespace json_interface;

    const nlohmann::json dog = parseJSONFile( INPUT( "valueAccess" ) );

    // Numbers
    BOOST_CHECK_EQUAL( getValue< unsigned int >( dog, "age" ), 11 );
    BOOST_CHECK_EQUAL( getValue< int >( dog, "age" ), 11 );
    BOOST_CHECK_EQUAL( getValue< double >( dog, "mass" ), 19.5 );
    BOOST_CHECK_EQUAL( getValue< float >( dog, "mass" ), 19.5 );
    BOOST_CHECK_EQUAL( getValue< long double >( dog, "mass" ), 19.5 );

    // Strings
    BOOST_CHECK_EQUAL( getValue< std::string >( dog, "name" ), "Bumper" );

    // Arrays
    const std::vector< std::string > hobbies = { "eat", "sleep" };
    const std::string hobbiesKey = "hobbies";
    BOOST_CHECK( getValue< std::vector< std::string > >( dog, hobbiesKey ) == hobbies );
    BOOST_CHECK_EQUAL( getValue< std::string >( dog, hobbiesKey / 0 ), hobbies.at( 0 ) );
    BOOST_CHECK_EQUAL( getValue< std::string >( dog, hobbiesKey / 1 ), hobbies.at( 1 ) );

    // Context: one level
    const nlohmann::json enemy = getValue< std::vector< nlohmann::json > >( dog, "enemies" ).front( );
    BOOST_CHECK_EQUAL( getRootObject( enemy ), dog );
    BOOST_CHECK_EQUAL( getValue< double >( enemy, "mass" ), 2.6 );
    BOOST_CHECK_EQUAL( getValue< double >( enemy, SpecialKeys::up / SpecialKeys::up / "mass" ), 19.5 );
    BOOST_CHECK_EQUAL( getValue< double >( enemy, SpecialKeys::root / "mass" ), 19.5 );

    // Context: several levels
    const nlohmann::json valencia =
            getValue< nlohmann::json >( enemy, std::string( "mother" ) / "birthplace" / "city" );
    BOOST_CHECK_EQUAL( valencia.at( "name" ), "Valencia" );
    BOOST_CHECK_EQUAL( getValue< double >( valencia, SpecialKeys::root / "mass" ), 19.5 );
    BOOST_CHECK_EQUAL( getValue< double >( valencia, SpecialKeys::up / SpecialKeys::up / SpecialKeys::up /
                                           SpecialKeys::up / SpecialKeys::up / "mass" ), 19.5 );
    BOOST_CHECK_EQUAL( getValue< double >( valencia, SpecialKeys::up / "continent" / "temperatureRange" / 0 ), -15 );
    BOOST_CHECK_EQUAL( getValue< double >( valencia, SpecialKeys::up / "continent" / "temperatureRange" / 1 ), 45 );

    // Eigen
    const Eigen::Matrix3d matrix = getValue< Eigen::Matrix3d >( dog, "orientation" );
    const Eigen::Matrix3d matrix2 = ( Eigen::Matrix3d( ) << 1.0, 0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 3.0, -1.0 ).finished( );
    BOOST_CHECK_EQUAL( matrix, matrix2 );

    // Map with non-char keys
    const std::map< double, std::string > food = getValue< std::map< double, std::string > >( dog, "food" );
    const std::map< double, std::string > food2 = { { 7, "feed" }, { 12, "meat" }, { 15, "feed" }, { 19, "feed" } };
    BOOST_CHECK( food == food2 );
}

// Test 2: modular
BOOST_AUTO_TEST_CASE( test_json_modular )
{
    using namespace json_interface;

    const nlohmann::json modular = getDeserializedJSON( getPathForJSONFile( INPUT( "modular" ) ) );

    nlohmann::json simulation;
    simulation[ "bodies" ] =
            R"(
            {
            "Earth": {
            "useDefaultSettings": true
            },
            "satellite": {
            "mass": 500
            }
            }
            )"_json;

    simulation[ "propagators" ] =
            R"(
            [
            {
            "centralBodies": [
            "Earth"
            ],
            "bodiesToPropagate": [
            "satellite"
            ]
            }
            ]
            )"_json;

    simulation[ "integrator" ] =
            R"(
            {
            "type": "rungeKutta4",
            "stepSize": 30
            }
            )"_json;

    simulation[ "export" ] =
            R"(
            [
            {
            "variables": [
            {
            "type": "independent"
            }
            ]
            },
            {
            "variables": [
            {
            "body": "satellite",
            "dependentVariableType": "relativePosition",
            "relatieToBody": "Earth"
            },
            {
            "body": "satellite",
            "dependentVariableType": "relativeVelocity",
            "relatieToBody": "Earth"
            }
            ]
            }
            ]
            )"_json;

    simulation[ "export" ][ 0 ][ "file" ] =
            ( boost::filesystem::path( "export" ) / "../outputs/epochs.txt" ).string( );
    simulation[ "export" ][ 1 ][ "file" ] =
            ( boost::filesystem::path( "export" ) / "../states.txt" ).string( );

    BOOST_CHECK_EQUAL( modular, simulation );
}

// Test 3: mergeable
BOOST_AUTO_TEST_CASE( test_json_mergeable )
{
    using namespace json_interface;
    using namespace boost::filesystem;

    const nlohmann::json merged1 = getDeserializedJSON(
                getPathForJSONFile( INPUT( ( path( "mergeable" ) / "inputs" / "merge1" ).string( ) ) ) );

    const nlohmann::json manual1 =
            R"(
            {
            "type": "rungeKutta4",
            "stepSize": 20
            }
            )"_json;

    BOOST_CHECK_EQUAL( merged1, manual1 );


    const nlohmann::json merged2 = getDeserializedJSON(
                getPathForJSONFile( INPUT( ( path( "mergeable" ) / "inputs" / "merge2" ).string( ) ) ) );

    const nlohmann::json manual2 =
            R"(
            {
            "integrator": {
            "type": "rungeKutta4",
            "stepSize": 20,
            "initialTimes": [
            0,
            86400
            ]
            }
            }
            )"_json;

    BOOST_CHECK_EQUAL( merged2, manual2 );


    const nlohmann::json merged3 = getDeserializedJSON(
                getPathForJSONFile( INPUT( ( path( "mergeable" ) / "inputs" / "merge3" ).string( ) ) ) );

    const nlohmann::json manual3 =
            R"(
            {
            "integrator": {
            "type": "rungeKutta4",
            "stepSize": 20,
            "initialTimes": [
            0,
            43200,
            86400
            ]
            },
            "spice": {
            "useStandardKernels": true,
            "preloadEpehemeris": false
            }
            }
            )"_json;

    BOOST_CHECK_EQUAL( merged3, manual3 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
