/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      100910    J. Melman         First creation of code.
 *      110111    J. Melman         Adapted to the offical Tudat standards.
 *      110124    J. Melman         Further adapted to the offical Tudat standards.
 *      110201    J. Melman         Made the tests for obliquity and astronomical unit more
 *                                  accurate.
 *      120127    D. Dirkx          Moved to Tudat core.
 *      120127    K. Kumar          Transferred unit tests over to Boost unit test framework.
 *      120128    K. Kumar          Changed BOOST_CHECK to BOOST_CHECK_CLOSE_FRACTION for unit test
 *                                  comparisons.
 *      130111    D. Dirkx          Added unit tests for Planck, Boltzmann, and Stefan-Boltzmann
 *                                  constants.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

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
            BOOST_CHECK_CLOSE_FRACTION( STEFAN_BOLTZMANN_CONSTANT, 2.0 * std::pow(
                                            mathematical_constants::PI, 5.0 ) *
                                        std::pow( BOLTZMANN_CONSTANT, 4.0 ) /
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
