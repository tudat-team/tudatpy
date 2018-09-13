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
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/Astrodynamics/Aerodynamics/customConstantTemperatureAtmosphere.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_custom_constant_temperature_atmosphere )

//! Test set- and get-functions of constants.
BOOST_AUTO_TEST_CASE( testCustomConstantTemperatureAtmosphereGetSet )
{
    // Initialize constants that need to be set.
    const double constantTemperature = 288.16;
    std::vector< double > modelSpecificParameters;
    modelSpecificParameters.push_back( 0.0 ); // reference altitude
    modelSpecificParameters.push_back( 1.225 ); // density at reference altitude
    modelSpecificParameters.push_back( 7.050e3 ); // scale height

    // Create a custom atmosphere object.
    aerodynamics::CustomConstantTemperatureAtmosphere customAtmosphere(
                aerodynamics::exponential_atmosphere_model, constantTemperature,
                physical_constants::SPECIFIC_GAS_CONSTANT_AIR, 1.4,
                modelSpecificParameters );

    BOOST_CHECK_EQUAL( constantTemperature, customAtmosphere.getConstantTemperature( ) );
}

//! Check whether the atmosphere is calculated correctly at sea level.
// Values from (US Standard Atmosphere, 1976).
BOOST_AUTO_TEST_CASE( testCustomConstantTemperatureAtmosphereSeaLevel )
{
    // Initialize constants that need to be set.
    const double constantTemperature = 288.16;
    std::vector< double > modelSpecificParameters;
    modelSpecificParameters.push_back( 0.0 ); // reference altitude
    modelSpecificParameters.push_back( 1.225 ); // density at reference altitude
    modelSpecificParameters.push_back( 7.050e3 ); // scale height
    const double altitude = 0.0;

    // Initialize the pressure at zero altitude.
    const double pressureAtZeroAltitude = 101325.0;

    // Create a custom atmosphere object.
    aerodynamics::CustomConstantTemperatureAtmosphere customAtmosphere(
                aerodynamics::exponential_atmosphere_model, constantTemperature,
                physical_constants::SPECIFIC_GAS_CONSTANT_AIR, 1.4,
                modelSpecificParameters );

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    BOOST_CHECK_CLOSE_FRACTION( constantTemperature, customAtmosphere.
                                getTemperature( altitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( modelSpecificParameters.at( 1 ), customAtmosphere.
                                getDensity( altitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( pressureAtZeroAltitude, customAtmosphere.
                                getPressure( altitude ), 0.002 * pressureAtZeroAltitude );
}

//! Test custom constant temperature atmosphere at 10 km altitude.
// The given value for pressure was calculated off-line with a calculator.
BOOST_AUTO_TEST_CASE( testCustomConstantTemperatureAtmosphereAt10km )
{
    {
        // Initialize constants that need to be set.
        const double constantTemperature = 288.16;
        std::vector< double > modelSpecificParameters;
        modelSpecificParameters.push_back( 0.0 ); // reference altitude
        modelSpecificParameters.push_back( 1.225 ); // density at reference altitude
        modelSpecificParameters.push_back( 7.050e3 ); // scale height
        const double altitude = 10.0e3;

        // Longitude, latitude and time to check overloading.
        double longitude = 0.0;
        double latitude = 0.0;
        double time = 0.0;

        // Create a custom atmosphere object.
        aerodynamics::CustomConstantTemperatureAtmosphere customAtmosphere(
                    aerodynamics::exponential_atmosphere_model, constantTemperature,
                    physical_constants::SPECIFIC_GAS_CONSTANT_AIR, 1.4,
                    modelSpecificParameters );

        // Declare and set expected density.
        const double expectedDensity  = modelSpecificParameters.at( 1 ) *
                std::exp( - altitude / modelSpecificParameters.at( 2 ) );
        const double expectedPressure = 24526.24934607106;

        // Declare tolerance used for Boost tests.
        const double tolerance = std::numeric_limits< double >::epsilon( );

        BOOST_CHECK_CLOSE_FRACTION( constantTemperature, customAtmosphere.
                                    getTemperature( altitude, longitude, latitude, time ), tolerance );
        BOOST_CHECK_CLOSE_FRACTION( expectedDensity, customAtmosphere.
                                    getDensity( altitude, longitude, latitude, time ), tolerance );
        BOOST_CHECK_CLOSE_FRACTION( expectedPressure, customAtmosphere.
                                    getPressure( altitude, longitude, latitude, time ), tolerance );
    }

    // Test three longitudinal waves atmosphere
    {
        // Initialize constants that need to be set.
        const double constantTemperature = 288.16;
        std::vector< double > modelSpecificParameters;
        modelSpecificParameters.push_back( 0.0 ); // reference altitude
        modelSpecificParameters.push_back( 1.225 ); // density at reference altitude
        modelSpecificParameters.push_back( 7.050e3 ); // scale height
        modelSpecificParameters.push_back( 1.0 ); // uncertainty factor
        modelSpecificParameters.push_back( 0.0 ); // dust storm factor
        const double altitude = 10.0e3;

        // Longitude, latitude and time to check overloading.
        double longitude = 0.0;
        double latitude = 0.0;
        double time = 0.0;

        // Create a custom atmosphere object.
        aerodynamics::CustomConstantTemperatureAtmosphere customAtmosphere(
                    aerodynamics::three_wave_atmosphere_model, constantTemperature,
                    physical_constants::SPECIFIC_GAS_CONSTANT_AIR, 1.4,
                    modelSpecificParameters );

        // Declare and set expected density.
        const double expectedDensity = modelSpecificParameters.at( 0 ) *
                std::exp( - altitude / modelSpecificParameters.at( 2 ) ) *
                ( modelSpecificParameters.at( 3 ) + modelSpecificParameters.at( 4 ) + // atmosphere parameters
                  0.1 * std::sin( 1.0 * longitude ) + // first longitudinal wave
                  0.2 * std::sin( 2.0 * ( longitude -
                                          unit_conversions::convertDegreesToRadians( 50.0 ) ) ) + // second longitudinal wave
                  0.1 * std::sin( 3.0 * ( longitude -
                                          unit_conversions::convertDegreesToRadians( 55.0 ) ) ) ); // third longitudinal wave
        const double expectedPressure = expectedDensity *
                physical_constants::SPECIFIC_GAS_CONSTANT_AIR * constantTemperature;

        // Declare tolerance used for Boost tests.
        const double tolerance = std::numeric_limits< double >::epsilon( );

        BOOST_CHECK_CLOSE_FRACTION( constantTemperature, customAtmosphere.
                                    getTemperature( altitude, longitude, latitude, time ), tolerance );
        BOOST_CHECK_CLOSE_FRACTION( expectedDensity, customAtmosphere.
                                    getDensity( altitude, longitude, latitude, time ), tolerance );
        BOOST_CHECK_CLOSE_FRACTION( expectedPressure, customAtmosphere.
                                    getPressure( altitude, longitude, latitude, time ), tolerance );
    }

    // Test three term atmosphere, where cosine and sine terms are zero (i.e., the model is the
    // conventional exponential atmosphere)
    {
        // Initialize constants that need to be set.
        const double constantTemperature = 288.16;
        std::vector< double > modelSpecificParameters;
        modelSpecificParameters.push_back( 0.0 ); // reference altitude
        modelSpecificParameters.push_back( 1.225 ); // density at reference altitude
        modelSpecificParameters.push_back( 7.050e3 ); // scale height
        modelSpecificParameters.push_back( -1.0 ); // weight for exponential term
        modelSpecificParameters.push_back( 0.0 ); // weight for cosine term
        modelSpecificParameters.push_back( 0.0 ); // weight for cosine term
        const double altitude = 10.0e3;

        // Longitude, latitude and time to check overloading.
        double longitude = 0.0;
        double latitude = 0.0;
        double time = 0.0;

        // Create a custom atmosphere object.
        aerodynamics::CustomConstantTemperatureAtmosphere customAtmosphere(
                    aerodynamics::three_term_atmosphere_model, constantTemperature,
                    physical_constants::SPECIFIC_GAS_CONSTANT_AIR, 1.4,
                    modelSpecificParameters );

        // Declare and set expected density.
        const double expectedDensity = modelSpecificParameters.at( 1 ) *
                std::exp( - altitude / modelSpecificParameters.at( 2 ) );
        const double expectedPressure = 24526.24934607106;

        // Declare tolerance used for Boost tests.
        const double tolerance = std::numeric_limits< double >::epsilon( );

        BOOST_CHECK_CLOSE_FRACTION( constantTemperature, customAtmosphere.
                                    getTemperature( altitude, longitude, latitude, time ), tolerance );
        BOOST_CHECK_CLOSE_FRACTION( expectedDensity, customAtmosphere.
                                    getDensity( altitude, longitude, latitude, time ), tolerance );
        BOOST_CHECK_CLOSE_FRACTION( expectedPressure, customAtmosphere.
                                    getPressure( altitude, longitude, latitude, time ), tolerance );
    }
}

//! Test if the position-independent functions work.
BOOST_AUTO_TEST_CASE( testCustomConstantTemperatureAtmospherePositionIndependentFunctions )
{
    // Initialize constants that need to be set.
    const double constantTemperature = 288.16;
    std::vector< double > modelSpecificParameters;
    modelSpecificParameters.push_back( 0.0 ); // reference altitude
    modelSpecificParameters.push_back( 1.225 ); // density at reference altitude
    modelSpecificParameters.push_back( 7.050e3 ); // scale height
    const double altitude = 10.0e3;

    // Longitude, latitude and time to check overloading.
    double longitude = 0.0;
    double latitude = 0.0;
    double time = 0.0;

    // Create a custom atmosphere object.
    aerodynamics::CustomConstantTemperatureAtmosphere customAtmosphere(
                aerodynamics::exponential_atmosphere_model, constantTemperature,
                physical_constants::SPECIFIC_GAS_CONSTANT_AIR, 1.4,
                modelSpecificParameters );

    const double density1 = customAtmosphere.getDensity( altitude );
    const double density2 = customAtmosphere.getDensity( altitude, longitude, latitude, time );

    const double pressure1 = customAtmosphere.getPressure( altitude );
    const double pressure2 = customAtmosphere.getPressure( altitude, longitude, latitude, time );

    const double temperature1 = customAtmosphere.getTemperature( altitude );
    const double temperature2 = customAtmosphere.getTemperature( altitude, longitude, latitude, time );

    BOOST_CHECK_EQUAL( density1, density2 );
    BOOST_CHECK_EQUAL( pressure1, pressure2 );
    BOOST_CHECK_EQUAL( temperature1, temperature2 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
