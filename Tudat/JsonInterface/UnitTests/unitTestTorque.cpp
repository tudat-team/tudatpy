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
#include "Tudat/JsonInterface/Propagation/torque.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_torque )

// Test 1: torque types
BOOST_AUTO_TEST_CASE( test_json_torque_ )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            basic_astrodynamics::torqueTypes,
                            basic_astrodynamics::unsupportedTorqueTypes );
}

// Test 2: second order gravitational torque
BOOST_AUTO_TEST_CASE( test_json_gravityField_secondOrderGravitational )
{
    using namespace basic_astrodynamics;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create TorqueSettings from JSON file
    const boost::shared_ptr< TorqueSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< TorqueSettings > >( INPUT( "secondOrderGravitational" ) );

    // Create TorqueSettings manually
    const boost::shared_ptr< TorqueSettings > manualSettings =
            boost::make_shared< TorqueSettings >( second_order_gravitational_torque );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 3: aerodynamic torque
BOOST_AUTO_TEST_CASE( test_json_gravityField_aerodynamic )
{
    using namespace basic_astrodynamics;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create TorqueSettings from JSON file
    const boost::shared_ptr< TorqueSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< TorqueSettings > >( INPUT( "aerodynamic" ) );

    // Create TorqueSettings manually
    const boost::shared_ptr< TorqueSettings > manualSettings =
            boost::make_shared< TorqueSettings >( aerodynamic_torque );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
