/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110316    F.M. Engelen      File created.
 *      110324    J. Melman         Time taken out of the equation.
 *      110427    F.M. Engelen      Changed input arguments for functions.
 *      110629    F.M. Engelen      Simplified unit test, removed dependency on other classes.
 *      110705    F.M. Engelen      Update to relative numerical errors.
 *      111128    B. Tong Minh      Added location-independent function test.
 *      111211    K. Kumar          Minor corrections to location-independent function test.
 *      120618    A. Ronse          Boostified unit test.
 *      120627    P. Musegaas       Changed scope of some variable + minor corrections.
 *
 *    References
 *      US Standard Atmosphere 1976,
 *          http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539_1977009539.pdf.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>

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
    tudat::aerodynamics::ExponentialAtmosphere exponentialAtmosphere;

    // Initialize the exponential atmosphere.
    exponentialAtmosphere.setConstantTemperature( constantTemperature );
    exponentialAtmosphere.setDensityAtZeroAltitude( densityAtZeroAltitude );
    exponentialAtmosphere.setScaleHeight( scaleHeight );
    exponentialAtmosphere.setSpecificGasConstant(
                    tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR );

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
    tudat::aerodynamics::ExponentialAtmosphere exponentialAtmosphere;

    // Initialize the exponential atmosphere.
    exponentialAtmosphere.setConstantTemperature( constantTemperature );
    exponentialAtmosphere.setDensityAtZeroAltitude( densityAtZeroAltitude );
    exponentialAtmosphere.setScaleHeight( scaleHeight );
    exponentialAtmosphere.setSpecificGasConstant(
                    tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR );

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
    tudat::aerodynamics::ExponentialAtmosphere exponentialAtmosphere;

    // Initialize the exponential atmosphere.
    exponentialAtmosphere.setConstantTemperature( constantTemperature );
    exponentialAtmosphere.setDensityAtZeroAltitude( densityAtZeroAltitude );
    exponentialAtmosphere.setScaleHeight( scaleHeight );
    exponentialAtmosphere.setSpecificGasConstant(
                    tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR );

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
    tudat::aerodynamics::ExponentialAtmosphere exponentialAtmosphere;

    // Initialize the exponential atmosphere.
    exponentialAtmosphere.setConstantTemperature( constantTemperature );
    exponentialAtmosphere.setDensityAtZeroAltitude( densityAtZeroAltitude );
    exponentialAtmosphere.setScaleHeight( scaleHeight );
    exponentialAtmosphere.setSpecificGasConstant(
                    tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR );

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
