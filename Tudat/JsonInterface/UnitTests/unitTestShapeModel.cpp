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
#include "Tudat/JsonInterface/Environment/shapeModel.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_shapeModel )

// Test 1: radiation pressure types
BOOST_AUTO_TEST_CASE( test_json_shapeModel_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            simulation_setup::bodyShapeTypes,
                            simulation_setup::unsupportedBodyShapeTypes );
}

// Test 2: spherical shape model
BOOST_AUTO_TEST_CASE( test_json_shapeModel_spherical )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create BodyShapeSettings from JSON file
    const boost::shared_ptr< BodyShapeSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< BodyShapeSettings > >( INPUT( "spherical" ) );

    // Create BodyShapeSettings manually
    const double radius = 6.4e6;
    const boost::shared_ptr< BodyShapeSettings > manualSettings =
            boost::make_shared< SphericalBodyShapeSettings >( radius );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 3: spherical Spice shape model
BOOST_AUTO_TEST_CASE( test_json_shapeModel_sphericalSpice )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create BodyShapeSettings from JSON file
    const boost::shared_ptr< BodyShapeSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< BodyShapeSettings > >( INPUT( "sphericalSpice" ) );

    // Create BodyShapeSettings manually
    const boost::shared_ptr< BodyShapeSettings > manualSettings =
            boost::make_shared< BodyShapeSettings >( spherical_spice );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 4: oblate spherical shape model
BOOST_AUTO_TEST_CASE( test_json_shapeModel_oblateSpherical )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create BodyShapeSettings from JSON file
    const boost::shared_ptr< BodyShapeSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< BodyShapeSettings > >( INPUT( "oblateSpherical" ) );

    // Create BodyShapeSettings manually
    const double equatorialRadius = 6.378e6;
    const double flattening = 0.0034;
    const boost::shared_ptr< BodyShapeSettings > manualSettings =
            boost::make_shared< OblateSphericalBodyShapeSettings >( equatorialRadius,
                                                                    flattening );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
