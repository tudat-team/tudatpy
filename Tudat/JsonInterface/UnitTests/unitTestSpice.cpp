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
#include "Tudat/JsonInterface/Environment/spice.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_spice )

// Test 1: standard kernels
BOOST_AUTO_TEST_CASE( test_json_spice_standard )
{
    using namespace json_interface;

    // Create SpiceSettings from JSON file
    const boost::shared_ptr< SpiceSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< SpiceSettings > >( INPUT( "standard" ) );

    // Create SpiceSettings manually
    boost::shared_ptr< SpiceSettings > manualSettings = boost::make_shared< SpiceSettings >( );
    manualSettings->useStandardKernels_ = true;
    manualSettings->preloadEphemeris_ = true;
    manualSettings->interpolationOffsets_ = { 10.0, 400.0 };

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 2: alternative kernels
BOOST_AUTO_TEST_CASE( test_json_spice_alternative )
{
    using namespace json_interface;

    // Create SpiceSettings from JSON file
    const boost::shared_ptr< SpiceSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< SpiceSettings > >( INPUT( "alternative" ) );

    // Create SpiceSettings manually
    boost::shared_ptr< SpiceSettings > manualSettings = boost::make_shared< SpiceSettings >( );
    manualSettings->useStandardKernels_ = true;
    manualSettings->alternativeKernels_ = { "foo.txt", "oof.txt" };
    manualSettings->preloadEphemeris_ = false;

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 3: custom kernels
BOOST_AUTO_TEST_CASE( test_json_spice_custom )
{
    using namespace json_interface;

    // Create SpiceSettings from JSON file
    const boost::shared_ptr< SpiceSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< SpiceSettings > >( INPUT( "custom" ) );

    // Create SpiceSettings manually
    boost::shared_ptr< SpiceSettings > manualSettings = boost::make_shared< SpiceSettings >( );
    manualSettings->useStandardKernels_ = false;
    manualSettings->kernels_ = { "foo.txt" };
    manualSettings->preloadEphemeris_ = true;
    manualSettings->interpolationStep_ = 9.0;

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
