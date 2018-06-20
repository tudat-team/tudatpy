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
#include "Tudat/JsonInterface/Propagation/export.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_export )

// Test 1: full result
BOOST_AUTO_TEST_CASE( test_json_export_full_result )
{
    using namespace tudat::propagators;
    using namespace tudat::json_interface;

    // Create ExportSettings from JSON file
    const std::shared_ptr< ExportSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< ExportSettings > >( INPUT( "fullResult" ) );

    // Create ExportSettings manually
    const std::string outputFile = "full.txt";
    const std::vector< std::shared_ptr< VariableSettings > > variables =
    {
        std::make_shared< VariableSettings >( independentVariable ),
        std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "body", "Earth" ),
    };
    std::shared_ptr< ExportSettings > manualSettings =
            std::make_shared< ExportSettings >( outputFile, variables );
    manualSettings->header_ = "Foo\n";
    manualSettings->epochsInFirstColumn_ = false;

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 2: partial result
BOOST_AUTO_TEST_CASE( test_json_export_partial_result )
{
    using namespace tudat::propagators;
    using namespace tudat::json_interface;

    // Create ExportSettings from JSON file
    const std::shared_ptr< ExportSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< ExportSettings > >( INPUT( "partialResult" ) );

    // Create ExportSettings manually
    const std::string outputFile = "partial.txt";
    const std::vector< std::shared_ptr< VariableSettings > > variables =
    {
        std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "body", "Earth" ),
    };
    std::shared_ptr< ExportSettings > manualSettings =
            std::make_shared< ExportSettings >( outputFile, variables );
    manualSettings->epochsInFirstColumn_ = true;
    manualSettings->onlyInitialStep_ = true;
    manualSettings->onlyFinalStep_ = true;
    manualSettings->numericalPrecision_ = 6;

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
