/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      To run this unit tests, a number of spice kernels need to be placed in the
 *      Spice kernel folder, by default External/SpiceInterface/Kernels or the
 *      SPICE_KERNEL_CUSTOM_FOLDER folder set as an argument to CMake or in UserSetings.txt.
 *      The required kernels are:
 *           de421.bsp
 *           pck00009.tpc
 *           naif0009.tls
 *           de-403-masses.tpc
 *      They can be found in a single zip file on the wiki at
 *      http://tudat.tudelft.nl/projects/tudat/wiki/SpiceInterface/ on the Tudat website or,
 *      alternatively, on the NAIF server at ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/.
 *
 */

#define BOOST_TEST_MAIN

#include "unitTestSupport.h"
#include <Tudat/External/JsonInterface/Environment/radiationPressure.h>

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_radiationPressure )

// Test 1: radiation pressure types
BOOST_AUTO_TEST_CASE( test_json_radiationPressure_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            simulation_setup::radiationPressureTypes,
                            simulation_setup::unsupportedRadiationPressureTypes );
}

// Test 2: cannon ball radiation pressure
BOOST_AUTO_TEST_CASE( test_json_radiationPressure_cannonBall )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create RadiationPressureInterfaceSettings from JSON file
    const boost::shared_ptr< RadiationPressureInterfaceSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< RadiationPressureInterfaceSettings > >( INPUT( "cannonBall" ) );

    // Create RadiationPressureInterfaceSettings manually
    const std::string sourceBody = "Sun";
    const double referenceArea = 2.0;
    const double radiationPressureCoefficient = 1.5;
    const std::vector< std::string > occultingBodies = { "Earth", "Moon" };
    const boost::shared_ptr< RadiationPressureInterfaceSettings > manualSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >( sourceBody,
                                                                                referenceArea,
                                                                                radiationPressureCoefficient,
                                                                                occultingBodies );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
