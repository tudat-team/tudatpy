/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/testLightTimeCorrections.h"

#include <limits>
#include <string>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/basics/testMacros.h"

#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;
using namespace observation_models;
using namespace spice_interface;

BOOST_AUTO_TEST_SUITE( test_light_time )

//! Test light-time calculator.
BOOST_AUTO_TEST_CASE( testLightWithSpice )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define names of bodies and frames.
    const std::string earth = "Earth";
    const std::string moon = "Moon";
    const std::string frame = "ECLIPJ2000";

    // Create ephemerides of Earth and Moon, with data from Spice.
    std::shared_ptr< SpiceEphemeris > earthEphemeris = std::make_shared< SpiceEphemeris >(
                earth, "SSB", false, false, false, frame );
    std::shared_ptr< SpiceEphemeris > moonEphemeris = std::make_shared< SpiceEphemeris >(
                moon, "SSB", false, false, false, frame );

    // Create light-time calculator, Earth center transmitter, Moon center receiver.
    std::shared_ptr< LightTimeCalculator< > > lightTimeEarthToMoon =
            std::make_shared< LightTimeCalculator< > >
            ( std::bind( &Ephemeris::getCartesianState, earthEphemeris, std::placeholders::_1 ),
              std::bind( &Ephemeris::getCartesianState, moonEphemeris, std::placeholders::_1 ) );

    // Define input time for tests.
    const double testTime = 1.0E6;

    // Define Spice output variables.
    double spiceOutputState[ 6 ] = { };
    double spiceMoonLightTime = 0.0;

    // Calculate observed (i.e. relative) position of Earth, and 'light time' at 'testTime' on
    // Moon, using spice. (Reception case with converged Newtonian light-time correction.)
    spkezr_c( earth.c_str( ), testTime, frame.c_str( ), std::string( "CN" ).c_str( ),
              moon.c_str( ), spiceOutputState, &spiceMoonLightTime );
    Eigen::Vector3d spiceMoonToEarthVector = Eigen::Vector3d::Zero( );
    for( int i = 0; i < 3; i++ )
    {
        // Convert from kilometers to meters.
        spiceMoonToEarthVector( i ) = -1000.0 * spiceOutputState[ i ];
    }

    // Calculate light time, with as input time the reception time, using light-time calculator.
    double testMoonLightTime = lightTimeEarthToMoon->calculateLightTime( testTime, true );
    BOOST_CHECK_CLOSE_FRACTION( testMoonLightTime, spiceMoonLightTime, 1.0E-9 );

    // Calculate relativeRange vector, with as input time the reception time, using light-time
    // calculator.
    const Eigen::Vector3d testMoonToEarthVector =
            lightTimeEarthToMoon->calculateRelativeRangeVector( testTime, true );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMoonToEarthVector, spiceMoonToEarthVector, 1.0E-12 );

    // Calculate observed (i.e. relative) position of Earth, and 'light time' at
    // 'testTime+light time' on Moon, using spice. (Transmission case with converged Newtonian
    // light-time correction.)
    spkezr_c( moon.c_str( ), testTime, frame.c_str( ), std::string( "XCN" ).c_str( ),
              earth.c_str( ), spiceOutputState, &spiceMoonLightTime );
    Eigen::Vector3d spiceEarthToMoonVector = Eigen::Vector3d::Zero( );
    for( int i = 0; i < 3; i++ )
    {
        // Convert from kilometers to meters.
        spiceEarthToMoonVector( i ) = 1000.0 * spiceOutputState[ i ];
    }

    // Calculate light time, with as input time the transmission time, using light-time calculator.
    testMoonLightTime = lightTimeEarthToMoon->calculateLightTime( testTime, false );
    BOOST_CHECK_CLOSE_FRACTION( testMoonLightTime, spiceMoonLightTime, 1.0E-9 );

    // Calculate relativeRange vector, with as input time the transmission time, using light-time
    // calculator.
    const Eigen::Vector3d testEarthToMoonVector =
            lightTimeEarthToMoon->calculateRelativeRangeVector( testTime, false );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testEarthToMoonVector, spiceEarthToMoonVector, 1.0E-10 );

    // Test light time and link end state functions.
    double testOutputTime = 0.0;
    Eigen::Vector6d testEarthState = Eigen::Vector6d::Zero( );
    Eigen::Vector6d testMoonState = Eigen::Vector6d::Zero( );
    Eigen::Vector6d spiceEarthState = Eigen::Vector6d::Zero( );
    Eigen::Vector6d spiceMoonState = Eigen::Vector6d::Zero( );

    // Get link end states, assuming input time is transmission time.
    // SSB = Solar system barycenter.
    testOutputTime = lightTimeEarthToMoon->calculateLightTimeWithLinkEndsStates(
                testMoonState, testEarthState, testTime, false );
    spiceEarthState = getBodyCartesianStateAtEpoch( earth, "SSB", "ECLIPJ2000", "NONE", testTime );
    spiceMoonState = getBodyCartesianStateAtEpoch(
                moon, "SSB", "ECLIPJ2000", "NONE", testTime + testOutputTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                spiceEarthState, testEarthState, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                spiceMoonState, testMoonState, std::numeric_limits< double >::epsilon( ) );

    // Get link end states, assuming input time is reception time.
    testOutputTime = lightTimeEarthToMoon->calculateLightTimeWithLinkEndsStates(
                testMoonState, testEarthState, testTime, true );
    spiceEarthState = getBodyCartesianStateAtEpoch( earth, "SSB", "ECLIPJ2000", "NONE",
                                                    testTime - testOutputTime );
    spiceMoonState = getBodyCartesianStateAtEpoch( moon, "SSB", "ECLIPJ2000", "NONE", testTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                spiceEarthState, testEarthState, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                spiceMoonState, testMoonState, std::numeric_limits< double >::epsilon( ) );

    // Test light time with corrections.

    // Set single light-time correction function.
    std::vector< LightTimeCorrectionFunction > lightTimeCorrections;
    lightTimeCorrections.push_back( &getTimeDifferenceLightTimeCorrection );

    // Create light-time object with correction.
    std::shared_ptr< LightTimeCalculator< > > lightTimeEarthToMoonWithCorrection =
            std::make_shared< LightTimeCalculator< > >
            ( std::bind( &Ephemeris::getCartesianState, earthEphemeris, std::placeholders::_1 ),
              std::bind( &Ephemeris::getCartesianState, moonEphemeris, std::placeholders::_1 ),
              lightTimeCorrections, true );

    // Calculate newtonian light time.
    double newtonianLightTime = lightTimeEarthToMoonWithCorrection->calculateRelativeRangeVector(
                testTime, true ).norm( ) / physical_constants::SPEED_OF_LIGHT;

    // Calculate light time (including correction), at reception.
    testMoonLightTime = lightTimeEarthToMoonWithCorrection->calculateLightTime( testTime, true );

    // Calculate expected correction.
    double expectedCorrection = getTimeDifferenceLightTimeCorrection(
                Eigen::Vector6d::Zero( ), Eigen::Vector6d::Zero( ),
                testTime - testMoonLightTime, testTime );

    // Test whether results are approximately equal.
    BOOST_CHECK_CLOSE_FRACTION( newtonianLightTime + expectedCorrection,
                                testMoonLightTime,
                                1E-14 );

    // Create light-time object with correction, without iterating light-time corrections.
    lightTimeEarthToMoonWithCorrection =
            std::make_shared< LightTimeCalculator< > >
            ( std::bind( &Ephemeris::getCartesianState, earthEphemeris, std::placeholders::_1 ),
              std::bind( &Ephemeris::getCartesianState, moonEphemeris, std::placeholders::_1 ),
              lightTimeCorrections, false );

    // Calculate newtonian light time.
    newtonianLightTime = lightTimeEarthToMoonWithCorrection->calculateRelativeRangeVector(
                testTime, true ).norm( ) / physical_constants::SPEED_OF_LIGHT;

    // Calculate light time (including correction), at reception.
    testMoonLightTime = lightTimeEarthToMoonWithCorrection->calculateLightTime( testTime, true );

    // Calculate expected correction.
    expectedCorrection = getTimeDifferenceLightTimeCorrection(
                Eigen::Vector6d::Zero( ), Eigen::Vector6d::Zero( ),
                testTime - testMoonLightTime, testTime );

    // Test whether results are approximately equal.
    BOOST_CHECK_CLOSE_FRACTION( newtonianLightTime + expectedCorrection,
                                testMoonLightTime,
                                1E-14 );

    // Add two more light-time corrections.
    lightTimeCorrections.push_back( &getVelocityDifferenceLightTimeCorrection );
    lightTimeCorrections.push_back( &getPositionDifferenceLightTimeCorrection );

    // Create light-time object with multiple corrections.
    lightTimeEarthToMoonWithCorrection =
            std::make_shared< LightTimeCalculator< > >
            ( std::bind( &Ephemeris::getCartesianState, earthEphemeris, std::placeholders::_1 ),
              std::bind( &Ephemeris::getCartesianState, moonEphemeris, std::placeholders::_1 ),
              lightTimeCorrections, true );

    // Calculate newtonian light time.
    newtonianLightTime = lightTimeEarthToMoonWithCorrection->calculateRelativeRangeVector(
                testTime, true ).norm( ) / physical_constants::SPEED_OF_LIGHT;

    // Calculate light time (including corrections), at reception.
    testMoonLightTime = lightTimeEarthToMoonWithCorrection->calculateLightTimeWithLinkEndsStates(
                testMoonState, testEarthState, testTime, true );

    // Calculate and sum expected correction.
    expectedCorrection = getTimeDifferenceLightTimeCorrection(
                testEarthState, testMoonState, testTime - testMoonLightTime, testTime );
    expectedCorrection += getPositionDifferenceLightTimeCorrection(
                testEarthState, testMoonState, testTime - testMoonLightTime, testTime );
    expectedCorrection += getVelocityDifferenceLightTimeCorrection(
                testEarthState, testMoonState, testTime - testMoonLightTime, testTime );

    // Test whether results are approximately equal.
    BOOST_CHECK_CLOSE_FRACTION( newtonianLightTime + expectedCorrection,
                                testMoonLightTime,
                                1E-14 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
