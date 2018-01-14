/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      US Standard Atmosphere 1976,
 *          http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539_1977009539.pdf.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Aerodynamics/exponentialAtmosphere.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_exponential_atmosphere )

// Summary of tests.
// Test 1: Test set- and get-functions of constants.
// Test 2: Test exponential atmosphere at sea level.
// Test 3: Test exponential atmosphere at 10 km altitude.
// Test 4: Test if the position-independent functions work.

//! Test set- and get-functions of constants.
BOOST_AUTO_TEST_CASE( testExponentialAtmosphereGetSet )
{
    // Initialize constants that need to be set.
    const double constantTemperature = 288.16;
    const double densityAtZeroAltitude = 1.225;
    const double scaleHeight = 7.050e3;

    // Create an exponential atmosphere object.
    aerodynamics::ExponentialAtmosphere exponentialAtmosphere(
                scaleHeight, constantTemperature, densityAtZeroAltitude );

    BOOST_CHECK_EQUAL( constantTemperature, exponentialAtmosphere.getConstantTemperature( ) );
    BOOST_CHECK_EQUAL( densityAtZeroAltitude, exponentialAtmosphere.getDensityAtZeroAltitude( ) );
    BOOST_CHECK_EQUAL( scaleHeight, exponentialAtmosphere.getScaleHeight( ) );
}

//! Check whether the atmosphere is calculated correctly at sea level.
// Values from (US Standard Atmosphere, 1976).
BOOST_AUTO_TEST_CASE( testExponentialAtmosphereSeaLevel )
{
    // Initialize constants that need to be set.
    const double constantTemperature = 288.16;
    const double densityAtZeroAltitude = 1.225;
    const double scaleHeight = 7.050e3;
    const double altitude = 0.0;

    // Initialize the pressure at zero altitude.
    const double pressureAtZeroAltitude = 101325.0;

    // Create an exponential atmosphere object.
    aerodynamics::ExponentialAtmosphere exponentialAtmosphere(
                scaleHeight, constantTemperature, densityAtZeroAltitude );

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    BOOST_CHECK_CLOSE_FRACTION( constantTemperature, exponentialAtmosphere.
                                getTemperature( altitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( densityAtZeroAltitude, exponentialAtmosphere.
                                getDensity( altitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( pressureAtZeroAltitude, exponentialAtmosphere.
                                getPressure( altitude ), 0.002 * pressureAtZeroAltitude );
}

//! Test exponential atmosphere at 10 km altitude.
// The given value for pressure was calculated off-line with a calculator.
BOOST_AUTO_TEST_CASE( testExponentialAtmosphereAt10km )
{
    // Initialize constants that need to be set.
    const double constantTemperature = 288.16;
    const double densityAtZeroAltitude = 1.225;
    const double scaleHeight = 7.050e3;
    const double altitude = 10.0e3;

    // Longitude, latitude and time to check overloading.
    double longitude = 0.0;
    double latitude = 0.0;
    double time = 0.0;

    // Create an exponential atmosphere object.
    aerodynamics::ExponentialAtmosphere exponentialAtmosphere(
                scaleHeight, constantTemperature, densityAtZeroAltitude );

    // Declare and set expected density.
    const double expectedDensity  = densityAtZeroAltitude * std::exp ( -altitude / scaleHeight );
    const double expectedPressure = 24526.24934607106;

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );


    BOOST_CHECK_CLOSE_FRACTION( constantTemperature, exponentialAtmosphere.
                                getTemperature( altitude, longitude, latitude, time ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedDensity, exponentialAtmosphere.
                                getDensity( altitude, longitude, latitude, time ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedPressure, exponentialAtmosphere.
                                getPressure( altitude, longitude, latitude, time ), tolerance );
}

//! Test if the position-independent functions work.
BOOST_AUTO_TEST_CASE( testExponentialAtmospherePositionIndependentFunctions )
{
    // Initialize constants that need to be set.
    const double constantTemperature = 288.16;
    const double densityAtZeroAltitude = 1.225;
    const double scaleHeight = 7.050e3;
    const double altitude = 10.0e3;

    // Longitude, latitude and time to check overloading.
    double longitude = 0.0;
    double latitude = 0.0;
    double time = 0.0;

    // Create an exponential atmosphere object.
    aerodynamics::ExponentialAtmosphere exponentialAtmosphere(
                scaleHeight, constantTemperature, densityAtZeroAltitude );

    const double density1 = exponentialAtmosphere.getDensity( altitude );
    const double density2 = exponentialAtmosphere.getDensity( altitude, longitude, latitude,
                                                              time );

    const double pressure1 = exponentialAtmosphere.getPressure( altitude );
    const double pressure2 = exponentialAtmosphere.getPressure( altitude, longitude, latitude,
                                                                time );

    const double temperature1 = exponentialAtmosphere.getTemperature( altitude );
    const double temperature2 = exponentialAtmosphere.getTemperature( altitude, longitude,
                                                                      latitude, time );

    BOOST_CHECK_EQUAL( density1, density2 );
    BOOST_CHECK_EQUAL( pressure1, pressure2 );
    BOOST_CHECK_EQUAL( temperature1, temperature2 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
