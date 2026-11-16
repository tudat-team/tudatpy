/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      NASA GMAT, src/UnitTests/TestParam/TestBplaneParam.cpp, expected values credited there to
 *          S. Hughes (2005).
 */

#define BOOST_TEST_MAIN

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <vector>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/close_encounters/bPlane.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_b_plane )

//! Check B-plane coordinates against the GMAT unit test for B-plane parameters.
BOOST_AUTO_TEST_CASE( testBPlaneCoordinatesMatchGmat )
{
    // GMAT default gravitational parameter of the Earth [m^3/s^2]
    const double earthGravitationalParameter = 398600.4415e9;

    // Hyperbolic geocentric state, EarthMJ2000Eq [m, m/s]
    Eigen::Vector6d relativeState;
    relativeState << 233410.6846140172e3, 83651.0868276347170e3, -168884.42195943173e3, -0.4038280708568842e3, 2.0665425988121107e3,
            0.4654706868112324e3;

    // GMAT expected values [m]
    const double expectedBDotT = 365738.686341826e3;
    const double expectedBDotR = 276107.260600374e3;

    const close_encounters::BPlaneData bPlane = close_encounters::makeBPlane( relativeState, earthGravitationalParameter );

    BOOST_CHECK_CLOSE_FRACTION( bPlane.bT, expectedBDotT, 1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( bPlane.bR, expectedBDotR, 1.0E-12 );
}

//! Check the B-plane of hyperbolic Kepler orbits against closed-form conic relations, at several points along
//! the orbit and far along both asymptotes.
BOOST_AUTO_TEST_CASE( testBPlaneOfHyperbolicKeplerOrbits )
{
    using namespace close_encounters;
    using namespace orbital_element_conversions;

    const double gravitationalParameter = 398600.4415e9;
    const double semiMajorAxis = -2.0E7;
    const double pi = mathematical_constants::PI;

    // Inclined, equatorial prograde, equatorial retrograde and polar orbit
    const std::vector< double > inclinations = { 0.7, 0.0, pi, pi / 2.0 };
    const std::vector< double > eccentricities = { 1.5, 2.0, 2.5, 3.0 };
    for( unsigned int i = 0; i < inclinations.size( ); i++ )
    {
        const double eccentricity = eccentricities.at( i );
        Eigen::Vector6d keplerianElements;
        keplerianElements << semiMajorAxis, eccentricity, inclinations.at( i ), 1.1, 2.3, 0.0;

        // At periapsis, the position and velocity are along the perifocal P and Q axes
        const Eigen::Vector6d periapsisState = convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter );
        const Eigen::Vector3d pUnitVector = periapsisState.segment< 3 >( 0 ).normalized( );
        const Eigen::Vector3d qUnitVector = periapsisState.segment< 3 >( 3 ).normalized( );
        const BPlaneData bPlane = makeBPlane( periapsisState, gravitationalParameter );
        const Eigen::Vector3d outgoingAsymptote = computeOutgoingAsymptoteDirection( periapsisState, gravitationalParameter );

        // Closed-form asymptotic velocity, impact parameter (semi-minor axis), asymptotes and turn angle
        const double asymptoteQComponent = std::sqrt( 1.0 - 1.0 / ( eccentricity * eccentricity ) );
        BOOST_CHECK_CLOSE_FRACTION( bPlane.asymptoticVelocity, std::sqrt( -gravitationalParameter / semiMajorAxis ), 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( bPlane.impactParameter, -semiMajorAxis * std::sqrt( eccentricity * eccentricity - 1.0 ), 1.0E-14 );
        BOOST_CHECK_SMALL( ( bPlane.sUnitVector - ( pUnitVector / eccentricity + asymptoteQComponent * qUnitVector ) ).norm( ), 1.0E-14 );
        BOOST_CHECK_SMALL( ( outgoingAsymptote - ( -pUnitVector / eccentricity + asymptoteQComponent * qUnitVector ) ).norm( ), 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION(
                std::acos( bPlane.sUnitVector.dot( outgoingAsymptote ) ), 2.0 * std::asin( 1.0 / eccentricity ), 1.0E-14 );

        // With the default reference direction -z, B is along +T (prograde) or -T (retrograde) for an equatorial
        // orbit, and along R for a polar orbit
        if( i == 1 || i == 2 )
        {
            BOOST_CHECK_CLOSE_FRACTION( bPlane.bT, ( i == 1 ? 1.0 : -1.0 ) * bPlane.impactParameter, 1.0E-14 );
            BOOST_CHECK_SMALL( bPlane.bR, 1.0E-14 * bPlane.impactParameter );
        }
        else if( i == 3 )
        {
            BOOST_CHECK_SMALL( bPlane.bT, 1.0E-14 * bPlane.impactParameter );
        }

        // Same B-plane at other points along the hyperbola
        const double asymptoteTrueAnomaly = std::acos( -1.0 / eccentricity );
        for( const double trueAnomalyFraction : { -0.9, -0.5, 0.5, 0.9 } )
        {
            keplerianElements( trueAnomalyIndex ) = trueAnomalyFraction * asymptoteTrueAnomaly;
            const BPlaneData otherBPlane =
                    makeBPlane( convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter ), gravitationalParameter );
            BOOST_CHECK_CLOSE_FRACTION( otherBPlane.asymptoticVelocity, bPlane.asymptoticVelocity, 1.0E-14 );
            BOOST_CHECK_SMALL( ( otherBPlane.sUnitVector - bPlane.sUnitVector ).norm( ), 1.0E-14 );
            BOOST_CHECK_SMALL( ( otherBPlane.bVector - bPlane.bVector ).norm( ), 1.0E-14 * bPlane.impactParameter );
        }

        // Far along the asymptotes, the velocity is along the asymptote, and the incoming position, projected on the
        // B-plane, is the B vector (to the rounding error of the far position, about 0.5 m at 2E15 m)
        const double farDistance = -1.0E8 * semiMajorAxis;
        const double semiLatusRectum = semiMajorAxis * ( 1.0 - eccentricity * eccentricity );
        const double farTrueAnomaly = std::acos( ( semiLatusRectum / farDistance - 1.0 ) / eccentricity );

        keplerianElements( trueAnomalyIndex ) = -farTrueAnomaly;
        const Eigen::Vector6d incomingState = convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter );
        const Eigen::Vector3d incomingPosition = incomingState.segment< 3 >( 0 );
        BOOST_CHECK_SMALL( ( incomingState.segment< 3 >( 3 ).normalized( ) - bPlane.sUnitVector ).norm( ), 1.0E-14 );
        BOOST_CHECK_SMALL( ( incomingPosition - incomingPosition.dot( bPlane.sUnitVector ) * bPlane.sUnitVector - bPlane.bVector ).norm( ),
                           1.0E-7 * bPlane.impactParameter );

        keplerianElements( trueAnomalyIndex ) = farTrueAnomaly;
        const Eigen::Vector6d outgoingState = convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter );
        BOOST_CHECK_SMALL( ( outgoingState.segment< 3 >( 3 ).normalized( ) - outgoingAsymptote ).norm( ), 1.0E-14 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
