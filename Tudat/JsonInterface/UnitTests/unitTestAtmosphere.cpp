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
#include "Tudat/JsonInterface/Environment/atmosphere.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_atmosphere )

// Test 1: atmosphere models
BOOST_AUTO_TEST_CASE( test_json_atmosphere_models )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "models" ),
                            simulation_setup::atmosphereTypes,
                            simulation_setup::unsupportedAtmosphereTypes );
}

// Test 2: exponential atmosphere
BOOST_AUTO_TEST_CASE( test_json_atmosphere_exponential )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create AtmosphereSettings from JSON file
    const std::shared_ptr< AtmosphereSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< AtmosphereSettings > >( INPUT( "exponential" ) );

    // Create AtmosphereSettings manually
    const double densityScaleHeight = 5.0;
    const double constantTemperature = 290.0;
    const double densityAtZeroAltitude = 1.0;
    const double specificGasConstant = 3.0;
    const double ratioOfSpecificHeats = -10.0;
    const std::shared_ptr< AtmosphereSettings > manualSettings =
            std::make_shared< ExponentialAtmosphereSettings >( densityScaleHeight,
                                                               constantTemperature,
                                                               densityAtZeroAltitude,
                                                               specificGasConstant,
                                                               ratioOfSpecificHeats );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 3: tabulated atmosphere
BOOST_AUTO_TEST_CASE( test_json_atmosphere_tabulated )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create AtmosphereSettings from JSON file
    const std::shared_ptr< AtmosphereSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< AtmosphereSettings > >( INPUT( "tabulated" ) );

    // Create AtmosphereSettings manually
    std::map< int, std::string > atmosphereTableFile;
    atmosphereTableFile[ 0 ] = "atmosphereTable.foo";
    std::vector< AtmosphereIndependentVariables > independentVariablesNames = { altitude_dependent_atmosphere, time_dependent_atmosphere };
    std::vector< AtmosphereDependentVariables > dependentVariablesNames = { temperature_dependent_atmosphere,
                                                                            gas_constant_dependent_atmosphere,
                                                                            density_dependent_atmosphere,
                                                                            pressure_dependent_atmosphere };
    const double specificGasConstant = 10.0;
    const double ratioOfSpecificHeats = 3.0;
    std::vector< interpolators::BoundaryInterpolationType > boundaryHandling = { interpolators::use_boundary_value,
                                                                                 interpolators::throw_exception_at_boundary };
    const std::shared_ptr< AtmosphereSettings > manualSettings = std::make_shared< TabulatedAtmosphereSettings >(
                atmosphereTableFile, independentVariablesNames, dependentVariablesNames,
                specificGasConstant, ratioOfSpecificHeats, boundaryHandling );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 4: NRLMSISE00 atmosphere
BOOST_AUTO_TEST_CASE( test_json_atmosphere_nrlmsise00 )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create AtmosphereSettings from JSON file
    const std::shared_ptr< AtmosphereSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< AtmosphereSettings > >( INPUT( "nrlmsise00" ) );

    // Create AtmosphereSettings manually
    const std::shared_ptr< AtmosphereSettings > manualSettings =
            std::make_shared< AtmosphereSettings >( nrlmsise00 );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 5: NRLMSISE00 atmosphere (custom space weather file)
BOOST_AUTO_TEST_CASE( test_json_atmosphere_nrlmsise00_custom )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create AtmosphereSettings from JSON file
    const std::shared_ptr< AtmosphereSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< AtmosphereSettings > >( INPUT( "nrlmsise00_custom" ) );

    // Create AtmosphereSettings manually
    const std::shared_ptr< AtmosphereSettings > manualSettings =
            std::make_shared< NRLMSISE00AtmosphereSettings >( "spaceWeatherFile.foo" );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
