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
#include "Tudat/JsonInterface/Environment/ephemeris.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_ephemeris )

// Test 1: ephemeris types
BOOST_AUTO_TEST_CASE( test_json_ephemeris_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            simulation_setup::ephemerisTypes,
                            simulation_setup::unsupportedEphemerisTypes );
}

// Test 2: bodies with ephemeris data
BOOST_AUTO_TEST_CASE( test_json_ephemeris_bodiesWithEphemerisData )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "bodiesWithEphemerisData" ),
                            ephemerides::bodiesWithEphemerisData,
                            ephemerides::unsupportedBodiesWithEphemerisData );
}

// Test 3: approximate planet position ephemeris
BOOST_AUTO_TEST_CASE( test_json_ephemeris_approximatePlanetPositions )
{
    using namespace ephemerides;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create EphemerisSettings from JSON file
    const std::shared_ptr< EphemerisSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< EphemerisSettings > >( INPUT( "approximatePlanetPositions" ) );

    // Create EphemerisSettings manually
    const ApproximatePlanetPositionsBase::BodiesWithEphemerisData bodyIdentifier =
            ApproximatePlanetPositionsBase::earthMoonBarycenter;
    const bool useCircularCoplanarApproximation = false;
    const std::shared_ptr< EphemerisSettings > manualSettings =
            std::make_shared< ApproximatePlanetPositionSettings >( bodyIdentifier,
                                                                     useCircularCoplanarApproximation );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 4: direct Spice ephemeris
BOOST_AUTO_TEST_CASE( test_json_ephemeris_directSpice )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create EphemerisSettings from JSON file
    const std::shared_ptr< EphemerisSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< EphemerisSettings > >( INPUT( "directSpice" ) );

    // Create EphemerisSettings manually
    const std::string frameOrigin = "Foo";
    const std::string frameOrientation = "FOO";
    const bool correctForStellarAberration = true;
    const bool correctForLightTimeAberration = false;
    const bool convergeLighTimeAberration = true;
    const std::shared_ptr< EphemerisSettings > manualSettings =
            std::make_shared< DirectSpiceEphemerisSettings >( frameOrigin,
                                                                frameOrientation,
                                                                correctForStellarAberration,
                                                                correctForLightTimeAberration,
                                                                convergeLighTimeAberration );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 5: tabulated ephemeris
BOOST_AUTO_TEST_CASE( test_json_ephemeris_tabulated )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create EphemerisSettings from JSON file
    const std::shared_ptr< EphemerisSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< EphemerisSettings > >( INPUT( "tabulated" ) );

    // Create EphemerisSettings manually
    std::map< double, Eigen::Vector6d > bodyStateHistory;
    bodyStateHistory[ 0.0 ] = ( Eigen::Vector6d( ) << 1.0, 0.0, 0.0, 0.0, -0.4, 0.0 ).finished( );
    bodyStateHistory[ 1.0 ] = ( Eigen::Vector6d( ) << 3.0, 0.0, 0.0, 0.0, -0.2, 0.0 ).finished( );
    bodyStateHistory[ 2.0 ] = ( Eigen::Vector6d( ) << 4.0, 0.0, 0.0, 0.0, -0.1, 0.0 ).finished( );


    const std::shared_ptr< EphemerisSettings > manualSettings =
            std::make_shared< TabulatedEphemerisSettings >( bodyStateHistory );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 6: interpolated Spice ephemeris
BOOST_AUTO_TEST_CASE( test_json_ephemeris_interpolatedSpice )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create EphemerisSettings from JSON file
    const std::shared_ptr< EphemerisSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< EphemerisSettings > >( INPUT( "interpolatedSpice" ) );

    // Create EphemerisSettings manually
    const double initialTime = 2.0;
    const double finalTime = 100.0;
    const double timeStep = 10.0;
    const std::string frameOrigin = "Foo";
    const std::string frameOrientation = "FOO";
    const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 4 );
    std::shared_ptr< EphemerisSettings > manualSettings =
            std::make_shared< InterpolatedSpiceEphemerisSettings >( initialTime,
                                                                      finalTime,
                                                                      timeStep,
                                                                      frameOrigin,
                                                                      frameOrientation,
                                                                      interpolatorSettings );
    manualSettings->resetMakeMultiArcEphemeris( true );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 7: constant ephemeris
BOOST_AUTO_TEST_CASE( test_json_ephemeris_constant )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create EphemerisSettings from JSON file
    const std::shared_ptr< EphemerisSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< EphemerisSettings > >( INPUT( "constant" ) );

    // Create EphemerisSettings manually
    const Eigen::Vector6d constantState = ( Eigen::Vector6d( ) << 0.0, 1.0, 0.0, -0.1, 0.0, 0.0 ).finished( );
    const std::shared_ptr< EphemerisSettings > manualSettings =
            std::make_shared< ConstantEphemerisSettings >( constantState );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 8: kepler ephemeris
BOOST_AUTO_TEST_CASE( test_json_ephemeris_kepler )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create EphemerisSettings from JSON file
    const std::shared_ptr< EphemerisSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< EphemerisSettings > >( INPUT( "kepler" ) );

    // Create EphemerisSettings manually
    const Eigen::Vector6d initialStateInKeplerianElements =
            ( Eigen::Vector6d( ) << 7.0e6, 0.1, 0.0, 0.0, 0.0, 0.0 ).finished( );
    const double epochOfInitialState= -4.0e4;
    const double centralBodyGravitationalParameter = 4.0e14;
    const std::string referenceFrameOrigin = "Foo";
    const std::string referenceFrameOrientation = "FOO";
    const double rootFinderAbsoluteTolerance = 1.0e-9;
    const double rootFinderMaximumNumberOfIterations = 100.0;
    const std::shared_ptr< EphemerisSettings > manualSettings =
            std::make_shared< KeplerEphemerisSettings >( initialStateInKeplerianElements,
                                                           epochOfInitialState,
                                                           centralBodyGravitationalParameter,
                                                           referenceFrameOrigin,
                                                           referenceFrameOrientation,
                                                           rootFinderAbsoluteTolerance,
                                                           rootFinderMaximumNumberOfIterations );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
