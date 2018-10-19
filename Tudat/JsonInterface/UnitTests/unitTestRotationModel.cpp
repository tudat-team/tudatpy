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
#include "Tudat/JsonInterface/Environment/rotationModel.h"

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
    using namespace spice_interface;
    using namespace json_interface;

    // Load spice kernels (needed for computeRotationQuaternionBetweenFrames)
    loadStandardSpiceKernels( );

    // Create RotationModelSettings from JSON file
    const std::shared_ptr< RotationModelSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< RotationModelSettings > >( INPUT( "simple" ) );

    // Create RotationModelSettings manually
    const std::string originalFrame = "ECLIPJ2000";
    const std::string targetFrame = "IAU_Earth";
    const double initialTime = 42.0;
    const Eigen::Quaterniond initialOrientation =
            spice_interface::computeRotationQuaternionBetweenFrames( originalFrame, targetFrame, initialTime );
    const double rotationRate = 2.0e-5;
    const std::shared_ptr< RotationModelSettings > manualSettings =
            std::make_shared< SimpleRotationModelSettings >( originalFrame,
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
    const std::shared_ptr< RotationModelSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< RotationModelSettings > >( INPUT( "spice" ) );

    // Create RotationModelSettings manually
    const std::string originalFrame = "foo";
    const std::string targetFrame = "oof";
    const std::shared_ptr< RotationModelSettings > manualSettings =
            std::make_shared< RotationModelSettings >( spice_rotation_model,
                                                         originalFrame,
                                                         targetFrame );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 3: Spice GCRS<->ITRS rotation model
BOOST_AUTO_TEST_CASE( test_json_rotationModel_gcrs_itrs )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create RotationModelSettings from JSON file
    const std::shared_ptr< RotationModelSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< RotationModelSettings > >( INPUT( "gcrsItrs" ) );

    // Create RotationModelSettings manually
    const std::string originalFrame = "J2000";
    const basic_astrodynamics::IAUConventions iauConvention = basic_astrodynamics::iau_2000_b;

    const std::shared_ptr< RotationModelSettings > manualSettings =
            std::make_shared< GcrsToItrsRotationModelSettings >(
                iauConvention, originalFrame );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
