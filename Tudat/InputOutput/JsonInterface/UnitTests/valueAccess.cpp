/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include "unitTestSupport.h"
#include <Tudat/InputOutput/JsonInterface/Support/valueAccess.h>
#include <Tudat/InputOutput/JsonInterface/Support/valueConversions.h>

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_valueAccess )

BOOST_AUTO_TEST_CASE( test_json_valueAccess_object )
{
    using namespace json_interface;

    const nlohmann::json dog = parseJSONFile( INPUT( "object" ) );

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
    const std::vector< nlohmann::json > enemies = getValue< std::vector< nlohmann::json > >( dog, "enemies" );
    const nlohmann::json enemy = enemies.front( );

    // Context: one level
    BOOST_CHECK_EQUAL( getRootObject( enemy ), dog );
    BOOST_CHECK_EQUAL( getRootObject( enemies.at( 0 ) ), dog );
    BOOST_CHECK_EQUAL( getParentKey( enemy ), "enemies" );
    BOOST_CHECK_EQUAL( getValue< double >( enemy, "mass" ), 2.6 );
    BOOST_CHECK_EQUAL( getValue< double >( enemy, SpecialKeys::up / "mass" ), 19.5 );
    BOOST_CHECK_EQUAL( getValue< double >( enemy, SpecialKeys::root / "mass" ), 19.5 );

    // Context: several levels
    const nlohmann::json valencia = getValue< nlohmann::json >( enemy, std::string( "mother" ) / "birthplace" / "city" );
    BOOST_CHECK_EQUAL( valencia.at( "name" ), "Valencia" );
    BOOST_CHECK_EQUAL( getValue< double >( valencia, SpecialKeys::root / "mass" ), 19.5 );
    BOOST_CHECK_EQUAL( getValue< double >( valencia, SpecialKeys::up / SpecialKeys::up /
                                           SpecialKeys::up / SpecialKeys::up / "mass" ), 19.5 );
    BOOST_CHECK_EQUAL( getValue< double >( valencia, SpecialKeys::up / "continent" / "temperatureRange" / 0 ), -15 );
    BOOST_CHECK_EQUAL( getValue< double >( valencia, SpecialKeys::up / "continent" / "temperatureRange" / 1 ), 45 );

    // Eigen
    const Eigen::Matrix3d matrix = getValue< Eigen::Matrix3d >( dog, "orientation" );
    // BOOST_CHECK_EQUAL( matrix,  );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
