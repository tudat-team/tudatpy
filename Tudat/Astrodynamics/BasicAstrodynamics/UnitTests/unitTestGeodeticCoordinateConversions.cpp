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
 *      130301    D. Dirkx          Migrated from personal code.
 *      130308    E.D. Brandon      Minor changes.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_geodetic_coordinate_conversions )

BOOST_AUTO_TEST_CASE( testGeodeticCoordinateConversions )
{
    using namespace tudat::basic_astrodynamics::coordinate_conversions;
    using namespace tudat::basic_astrodynamics::unit_conversions;

    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testCartesianPosition( 1917032.190, 6029782.349, -801376.113 );

    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testGeodeticPosition( -63.667,
                                                convertDegreesToRadians( -7.26654999 ),
                                                convertDegreesToRadians( 72.36312094 ) );

    // Central body characteristics (WGS84 Earth ellipsoid).
    const double flattening = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;

    // Test conversion to geodetic coordinates.
    {
        // Calculate geodetic position.
        const Eigen::Vector3d calculatedGeodeticPosition =
                convertCartesianToGeodeticCoordinates(
                    testCartesianPosition, equatorialRadius, flattening, 1.0E-4 );

        // Compare per coefficients (different tolerances).
        BOOST_CHECK_SMALL( calculatedGeodeticPosition.x( ) - testGeodeticPosition.x( ), 1.0E-4 );
        BOOST_CHECK_SMALL( calculatedGeodeticPosition.y( ) - testGeodeticPosition.y( ), 1.0E-10 );
        BOOST_CHECK_SMALL( calculatedGeodeticPosition.z( ) - testGeodeticPosition.z( ), 1.0E-10 );
    }

    // Test separate functions for altitude and geodetic latitude.
    {
        // Calculate altitude and geodetic latitude using dedicated functions.
        const double directAltitude = calculateAltitudeOverOblateSpheroid(
                    testCartesianPosition, equatorialRadius, flattening, 1.0E-4 );
        const double directGeodeticLatitude = calculateGeodeticLatitude(
                    testCartesianPosition, equatorialRadius, flattening, 1.0E-4 );

        // Compare values.
        BOOST_CHECK_SMALL( directAltitude - testGeodeticPosition.x( ), 1.0E-4 );
        BOOST_CHECK_SMALL( directGeodeticLatitude - testGeodeticPosition.y( ), 1.0E-10 );
    }

    // Test conversions from geodetic coordinates to cartesian position.
    {
        const Eigen::Vector3d calculateCartesianPosition =
                convertGeodeticToCartesianCoordinates(
                    testGeodeticPosition, equatorialRadius, flattening  );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    calculateCartesianPosition, testCartesianPosition, 1.0E-9 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
