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
#include "Tudat/JsonInterface/Propagation/massRateModel.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_massRateModel )

// Test 1: mass rate types
BOOST_AUTO_TEST_CASE( test_json_massRateModel_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            basic_astrodynamics::massRateTypes,
                            basic_astrodynamics::unsupportedMassRateType );
}

// Test 2: from thrust mass rate model
BOOST_AUTO_TEST_CASE( test_json_massRateModel_fromThrust )
{
    using namespace basic_astrodynamics;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create MassRateModelSettings from JSON file
    const std::shared_ptr< MassRateModelSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< MassRateModelSettings > >( INPUT( "fromThrust" ) );

    // Create MassRateModelSettings manually
    const bool useAllThrustModels = false;
    const std::string associatedThrustSource = "booster2";
    const std::shared_ptr< MassRateModelSettings > manualSettings =
            std::make_shared< FromThrustMassModelSettings >( useAllThrustModels, associatedThrustSource );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
