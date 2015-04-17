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
 *      121011    M. Ganeff         File created.
 *      130319    K. Kumar          Renamed filed and updated tests for celestial body constants.
 *      130321    D. Dirkx          Added unit tests for Earth's J2 and flattening factor.
 *      130322    K. Kumar          Moved J2 and flattening factor unit tests into their own Boost
 *                                  test definitions.
 *
 *    References
 *      International Earth Rotation and Reference System Service, Conventions 2010 (IERS2010),
 *        http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html.
 *      JPL, NASA. Astrodynamic Constants, http://ssd.jpl.nasa.gov/?constants,
 *        last updated: 13 Dec, 2012, last accessed: 19th March, 2013.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_celestial_body_constants )

//! Test if the correct gravitational parameters are returned.
BOOST_AUTO_TEST_CASE( test_GravitationalParameters )
{
    // Test if Sun's gravitational parameter is defined correctly.
    {
        // Defined expected gravitational parameter, as stated in (JPL, 2012).
        const double expectedGravitationalParameter = 1.32712440018e20;

        // Test that gravitational parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                    expectedGravitationalParameter, 1.0e-15 );
    }

    // Test if Mercury's gravitational parameter is defined correctly.
    {
        // Defined expected gravitational parameter, as derived from mass ratio given by
        // (JPL, 2012).
        const double expectedGravitationalParameter =
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER
                / 6023600.0;

        // Test that gravitational parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::MERCURY_GRAVITATIONAL_PARAMETER,
                    expectedGravitationalParameter, 1.0e-15 );
    }

    // Test if Venus's gravitational parameter is defined correctly.
    {
        // Defined expected gravitational parameter, as derived from mass ratio given by
        // (JPL, 2012).
        const double expectedGravitationalParameter =
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER
                / 408523.71;

        // Test that gravitational parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::VENUS_GRAVITATIONAL_PARAMETER,
                    expectedGravitationalParameter, 1.0e-15 );
    }

    // Test if Earth's gravitational parameter is defined correctly.
    {
        // Defined expected gravitational parameter, as stated in (IERS, 2010).
        const double exptectedGravitationalParameter = 3.986004418e14;

        // Test that gravitational parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER,
                    exptectedGravitationalParameter, 1.0e-15 );
    }

    // Test if Moon's gravitational parameter is defined correctly.
    {
        // Defined expected gravitational parameter, as derived from mass ratio given by
        // (JPL, 2012).
        const double expectedGravitationalParameter =
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER
                / ( 328900.56 * ( 1.0 + 81.30059 ) );

        // Test that gravitational parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::MOON_GRAVITATIONAL_PARAMETER,
                    expectedGravitationalParameter, 1.0e-15 );
    }

    // Test if Mars's gravitational parameter is defined correctly.
    {
        // Defined expected gravitational parameter, as derived from mass ratio given by
        // (JPL, 2012).
        const double expectedGravitationalParameter =
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER
                / 3098708.0;

        // Test that gravitational parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER,
                    expectedGravitationalParameter, 1.0e-15 );
    }

    // Test if Jupiter's gravitational parameter is defined correctly.
    {
        // Defined expected gravitational parameter, as derived from mass ratio given by
        // (JPL, 2012).
        const double expectedGravitationalParameter =
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER
                / 1047.3486;

        // Test that gravitational parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::JUPITER_GRAVITATIONAL_PARAMETER,
                    expectedGravitationalParameter, 1.0e-15 );
    }

    // Test if Saturn's gravitational parameter is defined correctly.
    {
        // Defined expected gravitational parameter, as derived from mass ratio given by
        // (JPL, 2012).
        const double expectedGravitationalParameter =
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER
                / 3497.898;

        // Test that gravitational parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::SATURN_GRAVITATIONAL_PARAMETER,
                    expectedGravitationalParameter, 1.0e-15 );
    }

    // Test if Uranus's gravitational parameter is defined correctly.
    {
        // Defined expected gravitational parameter, as derived from mass ratio given by
        // (JPL, 2012).
        const double expectedGravitationalParameter =
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER
                / 22902.98;

        // Test that gravitational parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::URANUS_GRAVITATIONAL_PARAMETER,
                    expectedGravitationalParameter, 1.0e-15 );
    }

    // Test if Neptune's gravitational parameter is defined correctly.
    {
        // Defined expected gravitational parameter, as derived from mass ratio given by
        // (JPL, 2012).
        const double expectedGravitationalParameter =
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER
                / 19412.24;

        // Test that gravitational parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::NEPTUNE_GRAVITATIONAL_PARAMETER,
                    expectedGravitationalParameter, 1.0e-15 );
    }

    // Test if Pluto's gravitational parameter is defined correctly.
    {
        // Defined expected gravitational parameter, as derived from mass ratio given by
        // (JPL, 2012).
        const double expectedGravitationalParameter =
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER
                / 1.35e8;

        // Test that gravitational parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::PLUTO_GRAVITATIONAL_PARAMETER,
                    expectedGravitationalParameter, 1.0e-15 );
    }
}

//! Test if the correct spherical harmonics coefficients are returned.
BOOST_AUTO_TEST_CASE( test_SphericalHarmonicsCoefficients )
{
    // Test if Earth's J2 parameter is defined correctly.
    {
        // Defined expected J2 parameter, as derived from mass ratio given by (JPL, 2012).
        const double expectedJ2 = -0.484165143790815E-03;

        // Test that J2 parameter is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::EARTH_GEODESY_NORMALIZED_J2,
                    expectedJ2, 1.0e-15 );

    }
}

//! Test if the correct equatorial radii are returned.
BOOST_AUTO_TEST_CASE( test_Radii )
{
    // Test if Earth's equatorial radius is defined correctly.
    {
        // Defined expected radius, as stated in (IERS, 2010).
        const double exptectedEquatorialRadius = 6378136.6;

        // Test that radius is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::EARTH_EQUATORIAL_RADIUS,
                    exptectedEquatorialRadius, 1.0e-15 );
    }
}

//! Test if the correct flattening parameters are returned.
BOOST_AUTO_TEST_CASE( test_FlatteningFactors )
{
    // Test if Earth's flattening factor is defined correctly.
    {
        // Defined expected flattening factor, as stated in (IERS, 2010).
        const double exptectedFlatteningFactor = 298.25642;

        // Test that flattening factor is defined correctly.
        BOOST_CHECK_CLOSE_FRACTION(
                    celestial_body_constants::EARTH_FLATTENING_FACTOR,
                    exptectedFlatteningFactor, 1.0e-15 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
