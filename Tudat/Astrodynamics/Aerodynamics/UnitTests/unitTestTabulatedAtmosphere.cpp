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
 *      Forget, François, et al. "Improved general circulation models of the Martian atmosphere from the surface to above 80 km."
 *          Journal of Geophysical Research: Planets 104.E10 (1999): 24155-24175.
 *      Millour, E., et al. "The Mars Climate Database (MCD version 5.2)." European Planetary Science Congress 2015. Vol. 10. EPSC, 2015.
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
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
// Test 3: Test tabulated atmosphere at 10.05 km altitude when just passing the altitude.
// Test 4: Test tabulated atmosphere at 1000 km altitude with table.
// Test 5: Test if the atmosphere file can be read multiple times.
// Test 6: Test if the position-independent functions work.

//! Check if the atmosphere is calculated correctly at sea level.
// Values from (US Standard Atmosphere, 1976).
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereAtSeaLevel )
{
    // Create a tabulated atmosphere object.
    std::string tabulatedAtmosphereFile = input_output::getAtmosphereTablesPath( ) +
            "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFile );

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    const double altitude = 0.0;
    BOOST_CHECK_CLOSE_FRACTION( 1.225, tabulatedAtmosphere.getDensity( altitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 101325.0, tabulatedAtmosphere.getPressure( altitude ), 1.0e-4 );
    BOOST_CHECK_CLOSE_FRACTION( 288.15, tabulatedAtmosphere.getTemperature( altitude ), tolerance );
}

//! Test tabulated atmosphere at 10km altitude including passing arbitrary longitude and latitude.
// The given value for pressure was obtained from table in book.
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereAt10km )
{
    // Create a tabulated atmosphere object.
    std::string tabulatedAtmosphereFile = input_output::getAtmosphereTablesPath( ) +
            "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFile );
    const double altitude = 10.0e3;
    const double longitude = 0.0;
    const double latitude = 0.0;
    const double time = 0.0;

    BOOST_CHECK_SMALL( 0.41351 - tabulatedAtmosphere.getDensity( altitude,  longitude, latitude, time ), 1.0e-4 );
    BOOST_CHECK_SMALL( 26500.0 - tabulatedAtmosphere.getPressure( altitude,  longitude, latitude, time ), 1.0 );
    BOOST_CHECK_SMALL( 223.26 - tabulatedAtmosphere.getTemperature( altitude,  longitude, latitude, time ), 1.0e-2 );
}

//! Test tabulated atmosphere at 10.05 km altitude when just passing the altitude.
// Check whether the atmosphere is calculated correctly at 10.05 km.
// The values are linear interpolated values based on book values.
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereAt10p5km )
{
    // Create a tabulated atmosphere object.
    std::string tabulatedAtmosphereFile = input_output::getAtmosphereTablesPath( ) +
            "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFile );
    const double altitude = 10.05e3;

    BOOST_CHECK_SMALL( 0.4110 - tabulatedAtmosphere.getDensity( altitude ), 1.0e-3 );
    BOOST_CHECK_SMALL( 26299.0 - tabulatedAtmosphere.getPressure( altitude ), 1.0 );
    BOOST_CHECK_SMALL( 222.9350 - tabulatedAtmosphere.getTemperature( altitude ), 2.0e-2 );
}

//! Test tabulated atmosphere at 1000 km altitude.
// Check whether the atmosphere is calculated correctly at 1000 km. Compared with the input table.
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereAt1000kmtab )
{
    // Create a tabulated atmosphere object.
    std::string tabulatedAtmosphereFiles = input_output::getAtmosphereTablesPath( ) +
            "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles );
    const double altitude = 1.0e6;

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    BOOST_CHECK_CLOSE_FRACTION( 3.5618e-15, tabulatedAtmosphere.getDensity( altitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 7.5158e-9, tabulatedAtmosphere.getPressure( altitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 1000.0, tabulatedAtmosphere.getTemperature( altitude ), tolerance );
}

//! Check if the atmosphere is calculated correctly when the variables are shuffled
// Values from (US Standard Atmosphere, 1976).
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereDependentVariables )
{
    // Create vector of dependent variables
    std::vector< aerodynamics::AtmosphereDependentVariables > dependentVariables;
    dependentVariables.push_back( aerodynamics::pressure_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::density_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::temperature_dependent_atmosphere );

    // Create a tabulated atmosphere object.
    std::string tabulatedAtmosphereFile = input_output::getAtmosphereTablesPath( ) +
            "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFile, dependentVariables );

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    const double altitude = 0.0;
    // Pressure and Density are switched
    BOOST_CHECK_CLOSE_FRACTION( 101325.0, tabulatedAtmosphere.getDensity( altitude ), 1.0e-4 );
    BOOST_CHECK_CLOSE_FRACTION( 1.225, tabulatedAtmosphere.getPressure( altitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 288.15, tabulatedAtmosphere.getTemperature( altitude ), tolerance );
}

//! Check if the atmosphere is calculated correctly when all the possible dependent variables are used and shuffled
// Values from Mars Climate Database Web Interface (http://www-mars.lmd.jussieu.fr/mcd_python/), avaraged over time, latitude and longitude.
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereAllDependentVariables )
{
    // Create vector of dependent and independent variables
    std::vector< aerodynamics::AtmosphereDependentVariables > dependentVariables;
    dependentVariables.push_back( aerodynamics::pressure_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::molar_mass_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::specific_heat_ratio_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::gas_constant_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::temperature_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::density_dependent_atmosphere );
    std::vector< aerodynamics::AtmosphereIndependentVariables > independentVariables = { aerodynamics::latitude_dependent_atmosphere };

    // Create a tabulated atmosphere object.
    std::map< int, std::string > tabulatedAtmosphereFiles;
    tabulatedAtmosphereFiles[ 0 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphere.dat";
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, independentVariables, dependentVariables );

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    const double latitude = 5e4;
    BOOST_CHECK_CLOSE_FRACTION( 4.3604106892e-02, tabulatedAtmosphere.getDensity( 0.0, 0.0, latitude ), 1.0e-4 );
    BOOST_CHECK_CLOSE_FRACTION( 7.6178752157e-05, tabulatedAtmosphere.getPressure( 0.0, 0.0, latitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 1.3794210386e+00, tabulatedAtmosphere.getTemperature( 0.0, 0.0, latitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 1.9068065814e+02, tabulatedAtmosphere.getSpecificGasConstant( 0.0, 0.0, latitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 1.6104994141e+02, tabulatedAtmosphere.getRatioOfSpecificHeats( 0.0, 0.0, latitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 2.3477457225e+00, tabulatedAtmosphere.getMolarMass( 0.0, 0.0, latitude ), 1.0e-4 );
    BOOST_CHECK_CLOSE_FRACTION( 205.817372408135, tabulatedAtmosphere.getSpeedOfSound( 0.0, 0.0, latitude ), 10.0 * tolerance );
}

//! Check if the atmosphere is calculated correctly when more than one independent variable is used, and when
//! dependent variables are shuffled.
// Values from Mars Climate Database Web Interface (http://www-mars.lmd.jussieu.fr/mcd_python/), avaraged over time.
BOOST_AUTO_TEST_CASE( testMultiDimensionalTabulatedAtmosphere )
{
    // Create vector of dependent and independent variables
    std::vector< aerodynamics::AtmosphereDependentVariables > dependentVariables = {
        aerodynamics::density_dependent_atmosphere, aerodynamics::pressure_dependent_atmosphere,
        aerodynamics::temperature_dependent_atmosphere};
    std::vector< aerodynamics::AtmosphereIndependentVariables > independentVariables = {
        aerodynamics::longitude_dependent_atmosphere, aerodynamics::latitude_dependent_atmosphere,
        aerodynamics::altitude_dependent_atmosphere };

    // Create a tabulated atmosphere object.
    std::map< int, std::string > tabulatedAtmosphereFiles;
    tabulatedAtmosphereFiles[ 0 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/density.dat";
    tabulatedAtmosphereFiles[ 1 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/pressure.dat";
    tabulatedAtmosphereFiles[ 2 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/temperature.dat";
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, independentVariables, dependentVariables );

    // Declare tolerance used for Boost tests.
    const double tolerance = 1e-7;

    // Define independent variable conditions
    const double altitude = 5.0e4;
    const double longitude = unit_conversions::convertDegreesToRadians( -180.0 );
    const double latitude = unit_conversions::convertDegreesToRadians( -90.0 );

    // Check that values matches with MATLAB interpolation (note that dependent variables are shuffled)
    BOOST_CHECK_CLOSE_FRACTION( 5.2805275e-05, tabulatedAtmosphere.getDensity( altitude, longitude, latitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 1.6627685, tabulatedAtmosphere.getPressure( altitude, longitude, latitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 151.544, tabulatedAtmosphere.getTemperature( altitude, longitude, latitude ), 1e2 * tolerance );
}

//! Check if the atmosphere is calculated correctly when more than one independent variable is used, and when both independent
//! and dependent variables are shuffled. Also tests for linear interpolation.
// Values from Mars Climate Database Web Interface (http://www-mars.lmd.jussieu.fr/mcd_python/), avaraged over time.
BOOST_AUTO_TEST_CASE( testMultiDimensionalTabulatedAtmosphereWithInterpolationAndShuffledVariables )
{
    // Create vector of dependent and independent variables
    std::vector< aerodynamics::AtmosphereDependentVariables > dependentVariables = {
        aerodynamics::specific_heat_ratio_dependent_atmosphere, aerodynamics::temperature_dependent_atmosphere,
        aerodynamics::density_dependent_atmosphere, aerodynamics::pressure_dependent_atmosphere,
        aerodynamics::gas_constant_dependent_atmosphere };
    std::vector< aerodynamics::AtmosphereIndependentVariables > independentVariables = {
        aerodynamics::longitude_dependent_atmosphere, aerodynamics::latitude_dependent_atmosphere,
        aerodynamics::altitude_dependent_atmosphere };

    // Create a tabulated atmosphere object.
    std::map< int, std::string > tabulatedAtmosphereFiles;
    tabulatedAtmosphereFiles[ 0 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/specificHeatRatio.dat";
    tabulatedAtmosphereFiles[ 1 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/temperature.dat";
    tabulatedAtmosphereFiles[ 2 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/density.dat";
    tabulatedAtmosphereFiles[ 3 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/pressure.dat";
    tabulatedAtmosphereFiles[ 4 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/gasConstant.dat";
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, independentVariables, dependentVariables );

    // Declare tolerance used for Boost tests.
    const double tolerance = 1e-5;

    // Define independent variable conditions
    const double altitude = 236.9862e3;
    const double longitude = unit_conversions::convertDegreesToRadians( 72.98632 );
    const double latitude = unit_conversions::convertDegreesToRadians( -65.9762 );

    // Check that values matches with MATLAB interpolation (note that dependent variables are shuffled)
    BOOST_CHECK_CLOSE_FRACTION( 5.50924276619674e-13, tabulatedAtmosphere.getDensity( altitude, longitude, latitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 5.05433422347663e-08, tabulatedAtmosphere.getPressure( altitude, longitude, latitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 174.841221295332, tabulatedAtmosphere.getTemperature( altitude, longitude, latitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 536.101455389577, tabulatedAtmosphere.getSpecificGasConstant( altitude, longitude, latitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 1.66646490380313, tabulatedAtmosphere.getRatioOfSpecificHeats( altitude, longitude, latitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 395.224168752852, tabulatedAtmosphere.getSpeedOfSound( altitude, longitude, latitude ), tolerance );
}

//! Check if the atmosphere is calculated correctly when default extrapolation values are used.
// Values from Mars Climate Database Web Interface (http://www-mars.lmd.jussieu.fr/mcd_python/), avaraged over time.
BOOST_AUTO_TEST_CASE( testMultiDimensionalTabulatedAtmosphereDefaultExtrapolation )
{
    // Create vector of dependent and independent variables
    std::vector< aerodynamics::AtmosphereDependentVariables > dependentVariables = {
        aerodynamics::pressure_dependent_atmosphere, aerodynamics::temperature_dependent_atmosphere,
        aerodynamics::specific_heat_ratio_dependent_atmosphere, aerodynamics::density_dependent_atmosphere };
    std::vector< aerodynamics::AtmosphereIndependentVariables > independentVariables = {
        aerodynamics::longitude_dependent_atmosphere, aerodynamics::latitude_dependent_atmosphere,
        aerodynamics::altitude_dependent_atmosphere };

    // Create a tabulated atmosphere object.
    std::map< int, std::string > tabulatedAtmosphereFiles;
    tabulatedAtmosphereFiles[ 0 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/pressure.dat";
    tabulatedAtmosphereFiles[ 1 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/temperature.dat";
    tabulatedAtmosphereFiles[ 2 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/specificHeatRatio.dat";
    tabulatedAtmosphereFiles[ 3 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/density.dat";

    // Test 1: check with single vector of default values
    {
        std::vector< interpolators::BoundaryInterpolationType > boundaryHandling = { interpolators::throw_exception_at_boundary,
                                                                                     interpolators::use_boundary_value_with_warning,
                                                                                     interpolators::use_default_value_with_warning };
        std::vector< double > defualtExtrapolationValues = { 100, 250, 1.3, 0 };
        aerodynamics::TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, independentVariables, dependentVariables,
                                                               boundaryHandling, defualtExtrapolationValues );

        // Declare tolerance used for Boost tests.
        const double tolerance = 1e-7;

        // Define independent variable conditions
        const std::vector< double > longitude = { unit_conversions::convertDegreesToRadians( 0.0 ),
                                                  unit_conversions::convertDegreesToRadians( 0.0 ) + 10.0,
                                                  unit_conversions::convertDegreesToRadians( 0.0 ) };
        const std::vector< double > latitude = { unit_conversions::convertDegreesToRadians( 0.0 ),
                                                 unit_conversions::convertDegreesToRadians( 0.0 ),
                                                 unit_conversions::convertDegreesToRadians( 0.0 ) + 10.0 };
        const std::vector< double > altitude = { 5.0e7, 5.0e4, 5.0e4 };

        // Expected results
        std::vector< double > boundaryValues = { 0.877653865782941, 163.354129386099, 1.37774305607899, 2.82975546661156e-05 };
        const std::vector< bool > expectedException = { false, true, false };
        std::vector< double > expectedResult;
        bool exception;

        // Check that values matches with MATLAB interpolation (note that dependent variables are shuffled)
        for ( unsigned int i = 0; i < 3; i++ )
        {
            exception = false;
            expectedResult = ( i == 0 ) ? defualtExtrapolationValues : boundaryValues;
            try
            {
                BOOST_CHECK_CLOSE_FRACTION( expectedResult.at( 0 ),
                                            tabulatedAtmosphere.getPressure( altitude.at( i ), longitude.at( i ), latitude.at( i ) ),
                                            tolerance );
                BOOST_CHECK_CLOSE_FRACTION( expectedResult.at( 1 ),
                                            tabulatedAtmosphere.getTemperature( altitude.at( i ), longitude.at( i ), latitude.at( i ) ),
                                            tolerance );
                BOOST_CHECK_CLOSE_FRACTION( expectedResult.at( 2 ),
                                            tabulatedAtmosphere.getRatioOfSpecificHeats( altitude.at( i ), longitude.at( i ), latitude.at( i ) ),
                                            tolerance );
                BOOST_CHECK_CLOSE_FRACTION( expectedResult.at( 3 ),
                                            tabulatedAtmosphere.getDensity( altitude.at( i ), longitude.at( i ), latitude.at( i ) ),
                                            tolerance );
            }
            catch ( std::runtime_error )
            {
                exception = true;
            }
            BOOST_CHECK_EQUAL( exception, expectedException.at( i ) );
        }
    }

    // Test 2: check with full set of default values
    {
        std::vector< interpolators::BoundaryInterpolationType > boundaryHandling =
        std::vector< interpolators::BoundaryInterpolationType >( independentVariables.size( ), interpolators::use_default_value );
        std::vector< std::vector< std::pair< double, double > > > defualtExtrapolationValues =
        {
            std::vector< std::pair< double, double > >( independentVariables.size( ), std::make_pair( -100.0, 100.0 ) ),
            std::vector< std::pair< double, double > >( independentVariables.size( ), std::make_pair( 250.0, 0.08652 ) ),
            std::vector< std::pair< double, double > >( independentVariables.size( ), std::make_pair( 1.4, 1.3 ) ),
            std::vector< std::pair< double, double > >( independentVariables.size( ), std::make_pair( 1e-3, 0.0 ) )
        };
        aerodynamics::TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, independentVariables, dependentVariables,
                                                               boundaryHandling, defualtExtrapolationValues );

        // Declare tolerance used for Boost tests.
        const double tolerance = 1e-7;

        // Define independent variable conditions
        const std::vector< double > longitude = { unit_conversions::convertDegreesToRadians( 0.0 ) - 10.0,
                                                  unit_conversions::convertDegreesToRadians( 0.0 ) + 10.0,
                                                  unit_conversions::convertDegreesToRadians( 0.0 ),
                                                  unit_conversions::convertDegreesToRadians( 0.0 ),
                                                  unit_conversions::convertDegreesToRadians( 0.0 ),
                                                  unit_conversions::convertDegreesToRadians( 0.0 ) };
        const std::vector< double > latitude = { unit_conversions::convertDegreesToRadians( 0.0 ),
                                                 unit_conversions::convertDegreesToRadians( 0.0 ),
                                                 unit_conversions::convertDegreesToRadians( 0.0 ) - 10.0,
                                                 unit_conversions::convertDegreesToRadians( 0.0 ) + 10.0,
                                                 unit_conversions::convertDegreesToRadians( 0.0 ),
                                                 unit_conversions::convertDegreesToRadians( 0.0 ) };
        const std::vector< double > altitude = { 5.0e5, 5.0e5, 5.0e5, 5.0e5, 0.0, 5.0e10 };

        // Check that values matches default values
        for ( unsigned int i = 0; i < altitude.size( ); i++ )
        {
            std::vector< double > expectedValue;
            for ( unsigned int j = 0; j < dependentVariables.size( ); j++ )
            {
                if ( ( i % 2 ) == 0.0 )
                {
                    expectedValue.push_back( defualtExtrapolationValues.at( j ).at( std::floor( i / 2.0 ) ).first );
                }
                else
                {
                    expectedValue.push_back( defualtExtrapolationValues.at( j ).at( std::floor( i / 2.0 ) ).second );
                }
            }

            BOOST_CHECK_CLOSE_FRACTION( expectedValue.at( 0 ),
                                        tabulatedAtmosphere.getPressure( altitude.at( i ), longitude.at( i ), latitude.at( i ) ),
                                        tolerance );
            BOOST_CHECK_CLOSE_FRACTION( expectedValue.at( 1 ),
                                        tabulatedAtmosphere.getTemperature( altitude.at( i ), longitude.at( i ), latitude.at( i ) ),
                                        tolerance );
            BOOST_CHECK_CLOSE_FRACTION( expectedValue.at( 2 ),
                                        tabulatedAtmosphere.getRatioOfSpecificHeats( altitude.at( i ), longitude.at( i ), latitude.at( i ) ),
                                        tolerance );
            BOOST_CHECK_CLOSE_FRACTION( expectedValue.at( 3 ),
                                        tabulatedAtmosphere.getDensity( altitude.at( i ), longitude.at( i ), latitude.at( i ) ),
                                        tolerance );
        }
    }
}

//! Check if the atmosphere is calculated correctly when heat ratio and gas constants are added.
// Values from (US Standard Atmosphere, 1976).
BOOST_AUTO_TEST_CASE( testTabulatedAtmosphereExtraVariables )
{
    // Create vector of dependent variables
    std::vector< aerodynamics::AtmosphereDependentVariables > dependentVariables;
    dependentVariables.push_back( aerodynamics::density_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::pressure_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::temperature_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::specific_heat_ratio_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::gas_constant_dependent_atmosphere );

    // Create a tabulated atmosphere object.
    aerodynamics::TabulatedAtmosphere tabulatedAtmosphere(
                input_output::getAtmosphereTablesPath( ) + "USSA1976Until100kmPer100mUntil1000kmPer1000m_wHR_GC.dat",
                dependentVariables );

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    const double altitude = 0.0;
    BOOST_CHECK_CLOSE_FRACTION( 8.0, tabulatedAtmosphere.getSpecificGasConstant( altitude ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( 1.7, tabulatedAtmosphere.getRatioOfSpecificHeats( altitude ), 1.0e-4 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
