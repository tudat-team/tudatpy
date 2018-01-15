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
#include "Tudat/JsonInterface/Propagation/variable.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_variable )

// Test 1: variable types
BOOST_AUTO_TEST_CASE( test_json_variable_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            propagators::variableTypes,
                            propagators::unsupportedVariableTypes );
}

// Test 2: dependent variable types
BOOST_AUTO_TEST_CASE( test_json_variable_dependentTypes )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "dependentTypes" ),
                            propagators::dependentVariableTypes,
                            propagators::unsupportedDependentVariableTypes );
}

// Test 3: dependent variable
BOOST_AUTO_TEST_CASE( test_json_variable_dependent )
{
    using namespace propagators;
    using namespace json_interface;

    // Create SingleDependentVariableSaveSettings from JSON file
    const boost::shared_ptr< SingleDependentVariableSaveSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< SingleDependentVariableSaveSettings > >( INPUT( "dependent" ) );

    // Create SingleDependentVariableSaveSettings manually
    const boost::shared_ptr< SingleDependentVariableSaveSettings > manualSettings =
            boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "body", "Earth" );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 4: acceleration variable
BOOST_AUTO_TEST_CASE( test_json_variable_acceleration )
{
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace json_interface;

    // Create SingleDependentVariableSaveSettings from JSON file
    const boost::shared_ptr< SingleDependentVariableSaveSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< SingleDependentVariableSaveSettings > >( INPUT( "acceleration" ) );

    // Create SingleDependentVariableSaveSettings manually
    const boost::shared_ptr< SingleDependentVariableSaveSettings > manualSettings =
            boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                thrust_acceleration, "body", "body", true );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 5: torque variable
BOOST_AUTO_TEST_CASE( test_json_variable_torque )
{
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace json_interface;

    // Create SingleDependentVariableSaveSettings from JSON file
    const boost::shared_ptr< SingleDependentVariableSaveSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< SingleDependentVariableSaveSettings > >( INPUT( "torque" ) );

    // Create SingleDependentVariableSaveSettings manually
    const boost::shared_ptr< SingleDependentVariableSaveSettings > manualSettings =
            boost::make_shared< SingleTorqueDependentVariableSaveSettings >(
                second_order_gravitational_torque, "Moon", "Earth", false, 2 );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 6: intermediate aerodynamic rotation matrix variable
BOOST_AUTO_TEST_CASE( test_json_variable_intermediateAerodynamicRotationMatrix )
{
    using namespace basic_astrodynamics;
    using namespace reference_frames;
    using namespace propagators;
    using namespace json_interface;

    // Create SingleDependentVariableSaveSettings from JSON file
    const boost::shared_ptr< SingleDependentVariableSaveSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< SingleDependentVariableSaveSettings > >(
                INPUT( "intermediateAerodynamicRotationMatrix" ) );

    // Create SingleDependentVariableSaveSettings manually
    const boost::shared_ptr< SingleDependentVariableSaveSettings > manualSettings =
            boost::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                "vehicle", body_frame, inertial_frame, 8 );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 7: relative body aerodynamic orientation angle variable
BOOST_AUTO_TEST_CASE( test_json_variable_relativeBodyAerodynamicOrientationAngle )
{
    using namespace basic_astrodynamics;
    using namespace reference_frames;
    using namespace propagators;
    using namespace json_interface;

    // Create SingleDependentVariableSaveSettings from JSON file
    const boost::shared_ptr< SingleDependentVariableSaveSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< SingleDependentVariableSaveSettings > >(
                INPUT( "relativeBodyAerodynamicOrientationAngle" ) );

    // Create SingleDependentVariableSaveSettings manually
    const boost::shared_ptr< SingleDependentVariableSaveSettings > manualSettings =
            boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >( "vehicle", bank_angle );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
