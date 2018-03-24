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
#include "Tudat/JsonInterface/Propagation/termination.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_termination )

// Test 1: single condition
BOOST_AUTO_TEST_CASE( test_json_termination_single )
{
    using namespace propagators;
    using namespace json_interface;

    // Create PropagationTerminationSettings from JSON file
    const boost::shared_ptr< PropagationTerminationSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< PropagationTerminationSettings > >( INPUT( "single" ) );

    // Create PropagationTerminationSettings manually
    const boost::shared_ptr< SingleDependentVariableSaveSettings > variable1 =
            boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "body", "Earth" );
    const boost::shared_ptr< PropagationTerminationSettings > condition1 =
            boost::make_shared< PropagationDependentVariableTerminationSettings >( variable1, 105.0E+3, true );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, condition1 );
}

// Test 2: multiple conditions (any of)
BOOST_AUTO_TEST_CASE( test_json_termination_anyOf )
{
    using namespace propagators;
    using namespace json_interface;

    // Create PropagationTerminationSettings from JSON file
    const boost::shared_ptr< PropagationTerminationSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< PropagationTerminationSettings > >( INPUT( "anyOf" ) );

    // Create PropagationTerminationSettings manually
    const boost::shared_ptr< SingleDependentVariableSaveSettings > variable1 =
            boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "body", "Earth" );
    const boost::shared_ptr< PropagationTerminationSettings > condition1 =
            boost::make_shared< PropagationDependentVariableTerminationSettings >( variable1, 105.0E+3, true );

    const boost::shared_ptr< PropagationTerminationSettings > condition2 =
            boost::make_shared< PropagationCPUTimeTerminationSettings >( 60.0 );

    const boost::shared_ptr< PropagationTerminationSettings > manualSettings =
            boost::make_shared< PropagationHybridTerminationSettings >(
                std::vector< boost::shared_ptr< PropagationTerminationSettings > >( { condition1, condition2 } ),
                true );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 3: multiple conditions (all of)
BOOST_AUTO_TEST_CASE( test_json_termination_allOf )
{
    using namespace propagators;
    using namespace json_interface;

    // Create PropagationTerminationSettings from JSON file
    const boost::shared_ptr< PropagationTerminationSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< PropagationTerminationSettings > >( INPUT( "allOf" ) );

    // Create PropagationTerminationSettings manually
    const boost::shared_ptr< SingleDependentVariableSaveSettings > variable1 =
            boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "body", "Earth" );
    const boost::shared_ptr< PropagationTerminationSettings > condition1 =
            boost::make_shared< PropagationDependentVariableTerminationSettings >( variable1, 105.0E+3, true );

    const boost::shared_ptr< SingleDependentVariableSaveSettings > variable3 =
            boost::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "body" );
    const boost::shared_ptr< PropagationTerminationSettings > condition3 =
            boost::make_shared< PropagationDependentVariableTerminationSettings >( variable3, 1.0, false );

    const boost::shared_ptr< PropagationTerminationSettings > manualSettings =
            boost::make_shared< PropagationHybridTerminationSettings >(
                std::vector< boost::shared_ptr< PropagationTerminationSettings > >( { condition1, condition3 } ),
                false );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 4: multiple conditions (combined)
BOOST_AUTO_TEST_CASE( test_json_termination_combined )
{
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace json_interface;

    // Create PropagationTerminationSettings from JSON file
    const boost::shared_ptr< PropagationTerminationSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< PropagationTerminationSettings > >( INPUT( "combined" ) );

    // Create PropagationTerminationSettings manually
    const boost::shared_ptr< SingleDependentVariableSaveSettings > variable1 =
            boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "body", "Earth" );
    const boost::shared_ptr< PropagationTerminationSettings > condition1 =
            boost::make_shared< PropagationDependentVariableTerminationSettings >( variable1, 105.0E+3, true );

    const boost::shared_ptr< PropagationTerminationSettings > condition2 =
            boost::make_shared< PropagationCPUTimeTerminationSettings >( 60.0 );

    const boost::shared_ptr< SingleDependentVariableSaveSettings > variable3 =
            boost::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "body" );
    const boost::shared_ptr< PropagationTerminationSettings > condition3 =
            boost::make_shared< PropagationDependentVariableTerminationSettings >( variable3, 1.0, false );

    const boost::shared_ptr< PropagationTerminationSettings > condition4 =
            boost::make_shared< PropagationTimeTerminationSettings >( 86400.0 );

    const boost::shared_ptr< SingleAccelerationDependentVariableSaveSettings > variable5 =
            boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                aerodynamic, "body", "Earth", false, 0 );
    const boost::shared_ptr< PropagationTerminationSettings > condition5 =
            boost::make_shared< PropagationDependentVariableTerminationSettings >( variable5, 0.0, true );

    const boost::shared_ptr< SingleAccelerationDependentVariableSaveSettings > variable6 =
            boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                aerodynamic, "body", "Earth", false, 1 );
    const boost::shared_ptr< PropagationTerminationSettings > condition6 =
            boost::make_shared< PropagationDependentVariableTerminationSettings >( variable6, 0.0, true );

    const boost::shared_ptr< PropagationTerminationSettings > c1or5or6 =
            boost::make_shared< PropagationHybridTerminationSettings >(
                std::vector< boost::shared_ptr< PropagationTerminationSettings > >(
    { condition1, condition5, condition6 } ), true );

    const boost::shared_ptr< PropagationTerminationSettings > c1or5or6and3 =
            boost::make_shared< PropagationHybridTerminationSettings >(
                std::vector< boost::shared_ptr< PropagationTerminationSettings > >( { c1or5or6, condition3 } ), false );

    const boost::shared_ptr< PropagationTerminationSettings > manualSettings =
            boost::make_shared< PropagationHybridTerminationSettings >(
                std::vector< boost::shared_ptr< PropagationTerminationSettings > >(
    { condition4, condition2, c1or5or6and3 } ), true );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
