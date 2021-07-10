/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_physical_constants )

//! Test if the physical constants have the correct relations (ratios/offsets).
BOOST_AUTO_TEST_CASE( testRelationsBetweenPhysicalConstant )
{
    using namespace physical_constants;

    // Test for the number of seconds in a year.
    BOOST_CHECK_CLOSE_FRACTION( JULIAN_YEAR, JULIAN_DAY * JULIAN_YEAR_IN_DAYS,
                                std::numeric_limits< double >::epsilon( ) );

    // Test for the number of seconds in a year.
    BOOST_CHECK_CLOSE_FRACTION( SIDEREAL_YEAR, JULIAN_DAY * SIDEREAL_YEAR_IN_DAYS,
                                std::numeric_limits< double >::epsilon( ) );

    // Test pre-computed powers of speed of light.
    BOOST_CHECK_CLOSE_FRACTION( std::pow( SPEED_OF_LIGHT, -2.0 ), INVERSE_SQUARE_SPEED_OF_LIGHT,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( std::pow( SPEED_OF_LIGHT, -3.0 ), INVERSE_CUBIC_SPEED_OF_LIGHT,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( std::pow( SPEED_OF_LIGHT, -4.0 ), INVERSE_QUARTIC_SPEED_OF_LIGHT,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( std::pow( SPEED_OF_LIGHT, -5.0 ), INVERSE_QUINTIC_SPEED_OF_LIGHT,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if physical constants have the expected value.
BOOST_AUTO_TEST_CASE( testOtherConstants )
{
    using namespace physical_constants;

    // Test for gravitational constant.
    BOOST_CHECK_CLOSE_FRACTION( GRAVITATIONAL_CONSTANT, 6.67259e-11,
                                std::numeric_limits< double >::epsilon( ) );

    // Test for speed of light.
    BOOST_CHECK_CLOSE_FRACTION( SPEED_OF_LIGHT, 299792458.0,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( SPEED_OF_LIGHT_LONG, 299792458.0L,
                                std::numeric_limits< long double >::epsilon( ) );

    // Test for astronomical unit.
    BOOST_CHECK_CLOSE_FRACTION( ASTRONOMICAL_UNIT, 1.49597870691e11,
                                std::numeric_limits< double >::epsilon( ) );

    // Test for molar gas constant.
    BOOST_CHECK_CLOSE_FRACTION( MOLAR_GAS_CONSTANT, 8.3144598,
                                std::numeric_limits< double >::epsilon( ) );

    // Test for Avogadro's number.
    BOOST_CHECK_CLOSE_FRACTION( AVOGADRO_CONSTANT, 6.022140857e23,
                                std::numeric_limits< double >::epsilon( ) );

    // Test for Planck constant.
    BOOST_CHECK_CLOSE_FRACTION( PLANCK_CONSTANT, 6.62606957E-34,
                                std::numeric_limits< double >::epsilon( ) );

    // Test for Boltzmann constant.
    BOOST_CHECK_CLOSE_FRACTION( BOLTZMANN_CONSTANT, 1.3806488E-23,
                                std::numeric_limits< double >::epsilon( ) );

    // Test permittivity/permeability fo vacuum
    BOOST_CHECK_CLOSE_FRACTION( VACUUM_PERMEABILITY, 4.0 * mathematical_constants::PI * 1.0E-7,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( VACUUM_PERMITTIVITY, 1.0 /
                                ( 4.0 * mathematical_constants::PI * 1.0E-7 * std::pow( SPEED_OF_LIGHT, 2.0 ) ),
                                std::numeric_limits< double >::epsilon( ) );

    // Test for Stefan-Boltzmann constant relation (derived from Planck and Boltzmann constants).
    BOOST_CHECK_CLOSE_FRACTION( STEFAN_BOLTZMANN_CONSTANT, 2.0 * physical_constants::compile_time_pow(
                                    mathematical_constants::PI, 5 ) *
                                physical_constants::compile_time_pow( BOLTZMANN_CONSTANT, 4 ) /
                                ( 15.0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT *
                                  PLANCK_CONSTANT * PLANCK_CONSTANT * PLANCK_CONSTANT ),
                                std::numeric_limits< double >::epsilon( ) );

    // Test for Stefan-boltzmann constant value (NIST, 2013)
    BOOST_CHECK_CLOSE_FRACTION( STEFAN_BOLTZMANN_CONSTANT, 5.670373E-8, 1.0E-7 );

    // Test time scale rate difference factors
    BOOST_CHECK_CLOSE_FRACTION( LG_TIME_RATE_TERM, 6.969290134E-10, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( LB_TIME_RATE_TERM, 1.550519768E-8, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( LG_TIME_RATE_TERM_LONG, 6.969290134E-10L, std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( LB_TIME_RATE_TERM_LONG, 1.550519768E-8L, std::numeric_limits< long  double >::epsilon( ) );


}

//! Check if the time constants have the expected values.
BOOST_AUTO_TEST_CASE( testTimeConstants )
{
    using namespace physical_constants;

    // Test for the number of Julian days in a year.
    BOOST_CHECK_CLOSE_FRACTION( JULIAN_YEAR_IN_DAYS, 365.25,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( JULIAN_YEAR_IN_DAYS_LONG, 365.25L,
                                std::numeric_limits< long double >::epsilon( ) );

    // Test for the number of sidereal days in a year.
    BOOST_CHECK_CLOSE_FRACTION( SIDEREAL_YEAR_IN_DAYS, 365.25636,
                                std::numeric_limits< double >::epsilon( ) );

    // Test for the Julian day length.
    BOOST_CHECK_CLOSE_FRACTION( SIDEREAL_DAY, 86164.09054,
                                std::numeric_limits< double >::epsilon( ) );

    // Test for the sidereal day length.
    BOOST_CHECK_CLOSE_FRACTION( JULIAN_DAY, 86400.0,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( JULIAN_DAY_LONG, 86400.0L,
                                std::numeric_limits< long double >::epsilon( ) );
}

//! Test templated get physical constant functions
BOOST_AUTO_TEST_CASE( testTemplatedConstantFunctions )
{
    using namespace physical_constants;

    // Test for the number of Julian days in a year.
    BOOST_CHECK_CLOSE_FRACTION( JULIAN_YEAR_IN_DAYS, getJulianYearInDays< double >( ),
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( JULIAN_YEAR_IN_DAYS_LONG, getJulianYearInDays< long double >( ),
                                std::numeric_limits< long double >::epsilon( ) );

    // Test for the sidereal day length.
    BOOST_CHECK_CLOSE_FRACTION( JULIAN_DAY, getJulianDay< double >( ),
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( JULIAN_DAY_LONG, getJulianDay< long double >( ),
                                std::numeric_limits< long double >::epsilon( ) );

    // Test for speed of light.
    BOOST_CHECK_CLOSE_FRACTION( SPEED_OF_LIGHT, getSpeedOfLight< double >( ),
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( SPEED_OF_LIGHT_LONG, getSpeedOfLight< long double >( ),
                                std::numeric_limits< long double >::epsilon( ) );

    // Test time scale rate difference factors
    BOOST_CHECK_CLOSE_FRACTION( LG_TIME_RATE_TERM, getLgTimeRateTerm< double >( ),
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( LB_TIME_RATE_TERM, getLbTimeRateTerm< double >( ),
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( LG_TIME_RATE_TERM_LONG, getLgTimeRateTerm< long double >( ),
                                std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( LB_TIME_RATE_TERM_LONG, getLbTimeRateTerm< long double >( ),
                                std::numeric_limits< long  double >::epsilon( ) );


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
