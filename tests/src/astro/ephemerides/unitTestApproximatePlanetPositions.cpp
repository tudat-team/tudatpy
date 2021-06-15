/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      HORIZONS Web-Interface, http://ssd.jpl.nasa.gov/horizons.cgi, last accessed: 5 April, 2011.
 *
 *    Notes
 *      It is noted during the 120513 check that this is not a very extensive and/or precise unit
 *      test. Given that the ephemeris class will soon be updated, this is not deemed a big issue.
 *      However the unit test will have to be improved a lot in the next update. Some code was
 *      outcommented when boostifying the unit test. It is attached in commented version at the
 *      bottom of this file.
 *
 *      Also the test of the approximate planet positions (3D) was changed. It used to check the
 *      only the spherical position coordinates. It was changed to check the cartesian elements in
 *      total. The accuracy with which this is possible is very low though.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/approximatePlanetPositionsCircularCoplanar.h"

namespace tudat
{
namespace unit_tests
{

//! Test the functionality of the approximate planet position functions
BOOST_AUTO_TEST_SUITE( test_approximate_planet_positions )

//! Test the orbital elements function against the orbital elements of Mars at JD 2455626.5.
BOOST_AUTO_TEST_CASE( testOrbitalElements )
{
    using unit_conversions::convertDegreesToRadians;
    using namespace ephemerides;

    // Set tolerance.
    const double tolerance = 2.0e-2;

    // Expected result.
    Eigen::Matrix< double, 6, 1 > expectedKeplerianElements;
    expectedKeplerianElements[ 0 ] = 2.279361944126564e11;
    expectedKeplerianElements[ 1 ] = 9.338126166083623e-2;
    expectedKeplerianElements[ 2 ] = convertDegreesToRadians( 1.848907897011101 );
    expectedKeplerianElements[ 3 ] = convertDegreesToRadians( 2.866464026954701e2 );
    expectedKeplerianElements[ 4 ] = convertDegreesToRadians( 4.952419052428279e1 );
    expectedKeplerianElements[ 5 ] = convertDegreesToRadians( 3.577219707986779e2 );

    // Create Mars ephemeris.
    ApproximateJplEphemeris marsEphemeris( "Mars" );

    // Convert the expected Keplerian elements to Cartesian elements.
    Eigen::Vector6d expectedEphemeris;
    expectedEphemeris = orbital_element_conversions::
            convertKeplerianToCartesianElements(
            expectedKeplerianElements,
            marsEphemeris.getSunGravitationalParameter( ) );

    // Retrieve state of Mars in Cartesian elements at Julian date 2455626.5.
    Eigen::Vector6d marsState = marsEphemeris.getCartesianState(
                ( 2455626.5 - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY );

    // Test if the computed ephemeris matches the expected ephemeris within the tolerance set.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedEphemeris, marsState, tolerance );

    // Check that the reference frame properties are as expected.
    BOOST_CHECK_EQUAL( marsEphemeris.getReferenceFrameOrientation( ), "ECLIPJ2000" );
    BOOST_CHECK_EQUAL( marsEphemeris.getReferenceFrameOrigin( ), "Sun" );
}

//! Test the cicular coplanar function against orbital elements of Mars at JD 2455626.5.
BOOST_AUTO_TEST_CASE( testCircularCoplannar )
{
    using namespace ephemerides;

    ApproximateJplCircularCoplanarEphemeris marsEphemeris(
                "Mars" );

    Eigen::Vector6d marsStateCircularCoplanar
            = marsEphemeris.getCartesianState(
                ( 2455626.5 - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY );

    // Compute the Keplerian elements from this ephemeris.
    Eigen::Vector6d keplerianElementsCircularCoplanar;
    keplerianElementsCircularCoplanar = orbital_element_conversions::
            convertCartesianToKeplerianElements( marsStateCircularCoplanar,
                    marsEphemeris.getSunGravitationalParameter( ) + marsEphemeris.getPlanetGravitationalParameter( ) );

    // Check the eccentricity, inclination and z-component of velocity and position are 0.
    BOOST_CHECK_SMALL( keplerianElementsCircularCoplanar( 1 ), 1e-15 );
    BOOST_CHECK_SMALL( keplerianElementsCircularCoplanar( 2 ),
                       std::numeric_limits< double >::min( ) );
    BOOST_CHECK_SMALL( marsStateCircularCoplanar( 2 ), 2.0e-5 );
    BOOST_CHECK_SMALL( marsStateCircularCoplanar( 5 ),
                       std::numeric_limits< double >::min( ) );

    // Check that the reference frame properties are as expected.
    BOOST_CHECK_EQUAL( marsEphemeris.getReferenceFrameOrientation( ), "ECLIPJ2000" );
    BOOST_CHECK_EQUAL( marsEphemeris.getReferenceFrameOrigin( ), "Sun" );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

//    // Compute the difference in semi-major axis between Test 2 and
//    // the external ephemeris_data "p_elem_t2.txt".
//    double errorSemiMajorAxis = fabs( positionMars.norm( )
//                                      - convertAstronomicalUnitsToMeters( 1.52371243 ) )
//            / convertAstronomicalUnitsToMeters( 1.52371243 );

//    if ( errorSemiMajorAxis > errorTolerance_ )
//    {
//        isApproximateSolarSystemEphemerisErroneous = true;

//        // Generate error statements.
//        cerr << "The computed relative error in position of the  " << endl;
//        cerr << "coplanar circular position of Mars ( " << errorSemiMajorAxis << " )" << endl;
//        cerr << "using the ApproximateJplCircularCoplanarEphemeris class, exceeds "
//             << "the maximum expected error " << endl;
//        cerr << "( " << errorTolerance_ << " )." << endl;
//    }

//    // Check orientation of position vector by comparison of separate components.
//    // Error in position should be smaller than maximum expected offset with respect to
//    // elliptical and inclined orbits.
//    double maximumErrorPosition =   keplerianElementsTest3D.getSemiMajorAxis( ) * (
//                keplerianElementsTest3D.getEccentricity( ) + 1.0
//                - cos( keplerianElementsTest3D.getInclination( ) ) );
//    Vector3d errorPositionVector = positionMars - marsEphemeris.getPosition( );

//    if ( fabs( errorPositionVector( 0 ) ) > maximumErrorPosition
//         || fabs( errorPositionVector( 1 ) ) > maximumErrorPosition )
//    {
//        isApproximateSolarSystemEphemerisErroneous = true;

//        // Generate error statements.
//        cerr << "The computed error in position vector of the  " << endl;
//        cerr << "coplanar circular position of Mars ( "
//             << errorPositionVector << " meters )" << endl;
//        cerr << "using the ApproximateJplCircularCoplanarEphemeris class, exceeds "
//             << "the expected error (" << endl;
//        cerr << "( " << maximumErrorPosition << " meters )." << endl;
//    }

    /* FIX THIS TEST!!!
    // Check size of velocity.
    Eigen::Vector3d errorVelocity = velocityMars - marsEphemeris.segment( 3, 3 );

    // Error in scalar velocity should be smaller than maximum expected offset with respect to
    // ellipitical and inclined orbits.
    double expectedErrorVelocity = fabs(
                sqrt( predefinedSun.getGravitationalParameter( )
                      / marsEphemeris.segment( 0, 3 ).norm( ) ) *
                ( ( 1.0 - cos( keplerianElementsTest3D.getInclination( ) )
                    + sqrt( ( 1.0 - keplerianElementsTest3D.getEccentricity( ) ) /
                            ( 1.0 + keplerianElementsTest3D.getEccentricity( ) ) ) - 1.0 ) ) );


    if ( errorVelocity.norm( ) > expectedErrorVelocity )
    {
        isApproximateSolarSystemEphemerisErroneous = true;

        // Generate error statements.
        cerr << "The computed error in velocity of the " << endl;
        cerr << "coplanar circular position of Mars "
             << "( " << errorVelocity.norm( ) << " meters per second )" << endl;
        cerr << "using the ApproximateJplCircularCoplanarEphemeris class, exceeds "
             << "the expected error " << endl;
        cerr << "( " << expectedErrorVelocity << " meters per second )." << endl;
        cerr << "The computed error exceeds the expected error by: "
             << fabs( errorVelocity.norm( ) - expectedErrorVelocity )
             << " meters per second." << endl;
    }
    */
