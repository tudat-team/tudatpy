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
 *      Introduction to Flight, Fifth edition, Appendix A, John D. Anderson Jr., McGraw Hill, 2005.
 *      US Standard Atmosphere 1976,
 *          http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539_1977009539.pdf.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_tabulated_atmosphere )

// Summary of tests.
// Test 1: Test tabulated atmosphere at sea level.
// Test 2: Test tabulated atmosphere at 10.0 km altitude including passing arbitrary longitude
//         and latitude.
// Test 3: Test tabulated atmosphere at 10.05 km altitude when just passing the altitude..
// Test 4: Test tabulated atmosphere at 1000 km altitude with table.
// Test 5: Test if the atmosphere file can be read multiple times.
// Test 6: Test if the position-independent functions work.

//! Check if the atmosphere is calculated correctly at sea level.
// Values from (US Standard Atmosphere, 1976).
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereAtSeaLevel )
{
    // Create a tabulated atmosphere object.
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere(
                input_output::getAtmosphereTablesPath( ) + "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    const double altitude = 0.0;
    BOOST_CHECK_CLOSE_FRACTION( 288.15, tabulatedAtmosphere.getTemperature( altitude ),
                                tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 1.225, tabulatedAtmosphere.getDensity( altitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 101325.0, tabulatedAtmosphere.getPressure( altitude ), 1.0e-4 );
}

//! Test tabulated atmosphere at 10km altitude including passing arbitrary longitude and latitude.
// The given value for pressure was obtained from table in book.
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereAt10km )
{
    // Create a tabulated atmosphere object.
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere(
                input_output::getAtmosphereTablesPath( ) + "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );
    const double altitude = 10.0e3;
    const double longitude = 0.0;
    const double latitude = 0.0;
    const double time = 0.0;

    BOOST_CHECK_SMALL( 223.26 - tabulatedAtmosphere.getTemperature( altitude,  longitude,
                                                                    latitude, time ), 1.0e-2 );
    BOOST_CHECK_SMALL( 0.41351 - tabulatedAtmosphere.getDensity( altitude,  longitude,
                                                                 latitude, time ), 1.0e-4 );
    BOOST_CHECK_SMALL( 26500.0 - tabulatedAtmosphere.getPressure( altitude,  longitude,
                                                                  latitude, time ), 1.0 );
}

//! Test tabulated atmosphere at 10.05 km altitude when just passing the altitude.
// Check whether the atmosphere is calculated correctly at 10.05 km.
// The values are linear interpolated values based on book values.
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereAt10p5km )
{
    // Create a tabulated atmosphere object.
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere(
                input_output::getAtmosphereTablesPath( ) + "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );
    const double altitude = 10.05e3;

    BOOST_CHECK_SMALL( 222.9350 - tabulatedAtmosphere.getTemperature( altitude ), 2.0e-2 );
    BOOST_CHECK_SMALL( 0.4110 - tabulatedAtmosphere.getDensity( altitude ), 1.0e-3 );
    BOOST_CHECK_SMALL( 26299.0 - tabulatedAtmosphere.getPressure( altitude ), 1.0 );
}

//! Test tabulated atmosphere at 1000 km altitude.
// Check whether the atmosphere is calculated correctly at 1000 km. Compared with the input table.
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereAt1000kmtab )
{
    // Create a tabulated atmosphere object.
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere(
                input_output::getAtmosphereTablesPath( ) + "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );
    const double altitude = 1.0e6 ;

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    BOOST_CHECK_CLOSE_FRACTION( 1000.0, tabulatedAtmosphere.getTemperature( altitude ),
                                tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 3.5618e-15, tabulatedAtmosphere.getDensity( altitude ),
                                tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 7.5158e-9, tabulatedAtmosphere.getPressure( altitude ),
                                tolerance );
}

//! Test if the position-independent functions work.
BOOST_AUTO_TEST_CASE( testTabulatedAtmospherePositionIndependentFunctions)
{
    // Create a tabulated atmosphere object.
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere(
                input_output::getAtmosphereTablesPath( ) + "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );

    const double altitude  = 10.0e3;

    // Longitude, latitude and time to check overloading.
    const double longitude = 0.0;
    const double latitude = 0.0;
    const double time = 0.0;

    const double density1 = tabulatedAtmosphere.getDensity( altitude );
    const double density2 = tabulatedAtmosphere.getDensity( altitude, longitude, latitude, time );

    const double pressure1 = tabulatedAtmosphere.getPressure( altitude );
    const double pressure2 = tabulatedAtmosphere.getPressure(
                altitude, longitude, latitude, time );

    const double temperature1 = tabulatedAtmosphere.getTemperature( altitude );
    const double temperature2 = tabulatedAtmosphere.getTemperature(
                altitude, longitude, latitude, time );

    BOOST_CHECK_EQUAL( density1, density2 );
    BOOST_CHECK_EQUAL( pressure1, pressure2 );
    BOOST_CHECK_EQUAL( temperature1, temperature2 );
}

//! Check if the atmosphere is calculated correctly when the variables are shuffled
// Values from (US Standard Atmosphere, 1976).
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereDependentVariables )
{

    std::vector< aerodynamics::TabulatedAtmosphere::AtmosphereDependentVariables > dependentVariables;
    dependentVariables.push_back( aerodynamics::TabulatedAtmosphere::pressure_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::TabulatedAtmosphere::density_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::TabulatedAtmosphere::temperature_dependent_atmosphere );

    // Create a tabulated atmosphere object.
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere(
                input_output::getAtmosphereTablesPath( ) + "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat",
                dependentVariables);

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    const double altitude = 0.0;
    BOOST_CHECK_CLOSE_FRACTION( 288.15, tabulatedAtmosphere.getTemperature( altitude ),
                                tolerance );
    // Pressure and Density are switched
    BOOST_CHECK_CLOSE_FRACTION( 101325.0, tabulatedAtmosphere.getDensity( altitude ), 1.0e-4 );
    BOOST_CHECK_CLOSE_FRACTION( 1.225, tabulatedAtmosphere.getPressure( altitude ), tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
