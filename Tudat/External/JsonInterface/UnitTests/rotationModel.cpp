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
#include <Tudat/External/JsonInterface/Environment/rotationModel.h>

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_rotationModel )

// Test 1: rotation model types
BOOST_AUTO_TEST_CASE( test_json_rotationModel_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            simulation_setup::rotationModelTypes,
                            simulation_setup::unsupportedRotationModelTypes );
}

// Test 2: simple rotation model
BOOST_AUTO_TEST_CASE( test_json_rotationModel_simple )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create RotationModelSettings from JSON file
    const boost::shared_ptr< RotationModelSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< RotationModelSettings > >( INPUT( "simple" ) );

    // Create RotationModelSettings manually
    const std::string originalFrame = "A";
    const std::string targetFrame = "B";
    const Eigen::Quaterniond initialOrientation( ( Eigen::Matrix3d( ) <<
                                                   1.0,  0.0,  0.0,
                                                   0.0,  1.0, -1.0,
                                                   0.0, -1.0,  1.0 ).finished( ) );
    const double initialTime = 42.0;
    const double rotationRate = 2.0e-5;
    const boost::shared_ptr< RotationModelSettings > manualSettings =
            boost::make_shared< SimpleRotationModelSettings >( originalFrame,
                                                               targetFrame,
                                                               initialOrientation,
                                                               initialTime,
                                                               rotationRate );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 3: Spice rotation model
BOOST_AUTO_TEST_CASE( test_json_rotationModel_spice )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create RotationModelSettings from JSON file
    const boost::shared_ptr< RotationModelSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< RotationModelSettings > >( INPUT( "spice" ) );

    // Create RotationModelSettings manually
    const std::string originalFrame = "foo";
    const std::string targetFrame = "oof";
    const boost::shared_ptr< RotationModelSettings > manualSettings =
            boost::make_shared< RotationModelSettings >( spice_rotation_model,
                                                         originalFrame,
                                                         targetFrame );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
