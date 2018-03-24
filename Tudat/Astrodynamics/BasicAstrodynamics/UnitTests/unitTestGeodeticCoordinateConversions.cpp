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
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_geodetic_coordinate_conversions )

BOOST_AUTO_TEST_CASE( testGeodeticCoordinateConversions )
{
    using namespace coordinate_conversions;
    using namespace unit_conversions;

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
