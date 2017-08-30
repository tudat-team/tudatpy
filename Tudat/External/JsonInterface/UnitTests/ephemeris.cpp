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
#include <Tudat/External/JsonInterface/Environment/ephemeris.h>

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_json_ephemeris )

// Test 1: ephemeris types
BOOST_AUTO_TEST_CASE( test_json_ephemeris_types )
{
    BOOST_CHECK_EQUAL_ENUM( "ephemeris_types",
                            simulation_setup::ephemerisTypes,
                            simulation_setup::unsupportedEphemerisTypes );
}

// Test 2: bodies with ephemeris data
BOOST_AUTO_TEST_CASE( test_json_ephemeris_bodiesWithEphemerisData )
{
    BOOST_CHECK_EQUAL_ENUM( "ephemeris_bodiesWithEphemerisData",
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
    const boost::shared_ptr< EphemerisSettings > fromFileSettings =
            readInputFile< boost::shared_ptr< EphemerisSettings > >( "ephemeris_approximatePlanetPositions" );

    // Create EphemerisSettings manually
    const ApproximatePlanetPositionsBase::BodiesWithEphemerisData bodyIdentifier =
            ApproximatePlanetPositionsBase::earthMoonBarycenter;
    const bool useCircularCoplanarApproximation = false;
    const boost::shared_ptr< EphemerisSettings > manualSettings =
            boost::make_shared< ApproximatePlanetPositionSettings >( bodyIdentifier,
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
    const boost::shared_ptr< EphemerisSettings > fromFileSettings =
            readInputFile< boost::shared_ptr< EphemerisSettings > >( "ephemeris_directSpice" );

    // Create EphemerisSettings manually
    const std::string frameOrigin = "Foo";
    const std::string frameOrientation = "FOO";
    const bool correctForStellarAbberation = true;
    const bool correctForLightTimeAbberation = false;
    const bool convergeLighTimeAbberation = true;
    const boost::shared_ptr< EphemerisSettings > manualSettings =
            boost::make_shared< DirectSpiceEphemerisSettings >( frameOrigin,
                                                                frameOrientation,
                                                                correctForStellarAbberation,
                                                                correctForLightTimeAbberation,
                                                                convergeLighTimeAbberation );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 5: tabulated ephemeris
BOOST_AUTO_TEST_CASE( test_json_ephemeris_tabulated )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create EphemerisSettings from JSON file
    const boost::shared_ptr< EphemerisSettings > fromFileSettings =
            readInputFile< boost::shared_ptr< EphemerisSettings > >( "ephemeris_tabulated" );

    // Create EphemerisSettings manually
    const std::map< double, Eigen::Vector6d > bodyStateHistory =
    {
        { 0.0, ( Eigen::Vector6d( ) << 1.0, 0.0, 0.0, 0.0, -0.4, 0.0 ).finished( ) },
        { 1.0, ( Eigen::Vector6d( ) << 3.0, 0.0, 0.0, 0.0, -0.2, 0.0 ).finished( ) },
        { 2.0, ( Eigen::Vector6d( ) << 4.0, 0.0, 0.0, 0.0, -0.1, 0.0 ).finished( ) }
    };
    const boost::shared_ptr< EphemerisSettings > manualSettings =
            boost::make_shared< TabulatedEphemerisSettings >( bodyStateHistory );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 6: interpolated Spice ephemeris
BOOST_AUTO_TEST_CASE( test_json_ephemeris_interpolatedSpice )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create EphemerisSettings from JSON file
    const boost::shared_ptr< EphemerisSettings > fromFileSettings =
            readInputFile< boost::shared_ptr< EphemerisSettings > >( "ephemeris_interpolatedSpice" );

    // Create EphemerisSettings manually
    const double initialTime = 2.0;
    const double finalTime = 100.0;
    const double timeStep = 10.0;
    const std::string frameOrigin = "Foo";
    const std::string frameOrientation = "FOO";
    const boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            boost::make_shared< interpolators::LagrangeInterpolatorSettings >( 4 );
    boost::shared_ptr< EphemerisSettings > manualSettings =
            boost::make_shared< InterpolatedSpiceEphemerisSettings >( initialTime,
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
    const boost::shared_ptr< EphemerisSettings > fromFileSettings =
            readInputFile< boost::shared_ptr< EphemerisSettings > >( "ephemeris_constant" );

    // Create EphemerisSettings manually
    const Eigen::Vector6d constantState = ( Eigen::Vector6d( ) << 0.0, 1.0, 0.0, -0.1, 0.0, 0.0 ).finished( );
    const boost::shared_ptr< EphemerisSettings > manualSettings =
            boost::make_shared< ConstantEphemerisSettings >( constantState );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 8: kepler ephemeris
BOOST_AUTO_TEST_CASE( test_json_ephemeris_kepler )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create EphemerisSettings from JSON file
    const boost::shared_ptr< EphemerisSettings > fromFileSettings =
            readInputFile< boost::shared_ptr< EphemerisSettings > >( "ephemeris_kepler" );

    // Create EphemerisSettings manually
    const Eigen::Vector6d initialStateInKeplerianElements =
            ( Eigen::Vector6d( ) << 7.0e6, 0.1, 0.0, 0.0, 0.0, 0.0 ).finished( );
    const double epochOfInitialState= -4.0e4;
    const double centralBodyGravitationalParameter = 4.0e14;
    const std::string referenceFrameOrigin = "Foo";
    const std::string referenceFrameOrientation = "FOO";
    const double rootFinderAbsoluteTolerance = 1.0e-9;
    const double rootFinderMaximumNumberOfIterations = 100.0;
    const boost::shared_ptr< EphemerisSettings > manualSettings =
            boost::make_shared< KeplerEphemerisSettings >( initialStateInKeplerianElements,
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
