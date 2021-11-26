/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

namespace tudat
{
namespace unit_tests
{

using namespace orbital_element_conversions;
using namespace unit_conversions;
using namespace reference_frames;

BOOST_AUTO_TEST_SUITE( test_spherical_state_conversions )

//! Test Cartesian to spherical orbital state transformations, using Matlab script of Erwin Mooij to generate reference data.
//! See test_aerodynamic_angle_calculator unit test
BOOST_AUTO_TEST_CASE( testSphericalStateConversions )
{
    // Test case 1: arbitrary rotation
    {
        // Define Cartesian state
        Eigen::Vector6d cartesianState;
        cartesianState << -1656517.23153109, -5790058.28764025, -2440584.88186829,
                6526.30784888051, -2661.34558272018, 2377.09572383163;

        // Define associated spherical orbital state
        double testHeadingAngle = 1.229357188236127;
        double testFlightPathAngle = -0.024894033070522;
        double testLatitude = -0.385027359562548;
        double testLongitude = -1.849449608688977;
        double radius = cartesianState.segment( 0, 3 ).norm( );
        double speed = cartesianState.segment( 3, 3 ).norm( );

        // Convert pack and forth
        Eigen::Vector6d sphericalOrbitState  = convertCartesianToSphericalOrbitalState(
                    cartesianState );
        Eigen::Vector6d reconvertedCartesianState  = convertSphericalOrbitalToCartesianState(
                    sphericalOrbitState );

        // Check computed spherical orbital state against reference data
        BOOST_CHECK_SMALL(
                    std::fabs( radius - sphericalOrbitState( radiusIndex ) ),
                    radius * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testLatitude - sphericalOrbitState( latitudeIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testLongitude - sphericalOrbitState( longitudeIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( speed - sphericalOrbitState( speedIndex ) ),
                    speed * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testFlightPathAngle - sphericalOrbitState( flightPathIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testHeadingAngle - sphericalOrbitState( headingAngleIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );

        // Check consistency of back-and-forth conversion
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i ) - cartesianState( i ) ),
                        radius * 2.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i  + 3 ) - cartesianState( i + 3 ) ),
                        speed * 2.0 * std::numeric_limits< double >::epsilon( ) );
        }
    }

    // Test case 2: rotation with zero (heading, flight path, latitude) and half pi (longitude) angles.
    {
        // Define Cartesian state
        Eigen::Vector6d cartesianState;
        cartesianState << 0.0, 6498098.09700000, 0.0, 0.0, 0.0, 7.438147520000000e+03;

        // Define associated spherical orbital state
        double testHeadingAngle = 0.0;
        double testFlightPathAngle = 0.0;
        double testLatitude = 0.0;
        double testLongitude = mathematical_constants::PI / 2.0;
        double radius = cartesianState.segment( 0, 3 ).norm( );
        double speed = cartesianState.segment( 3, 3 ).norm( );

        // Convert pack and forth
        Eigen::Vector6d sphericalOrbitState  = convertCartesianToSphericalOrbitalState(
                    cartesianState );
        Eigen::Vector6d reconvertedCartesianState  = convertSphericalOrbitalToCartesianState(
                    sphericalOrbitState );

        // Check computed spherical orbital state against reference data
        BOOST_CHECK_SMALL(
                    std::fabs( radius - sphericalOrbitState( radiusIndex ) ),
                    radius * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testLatitude - sphericalOrbitState( latitudeIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testLongitude - sphericalOrbitState( longitudeIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( speed - sphericalOrbitState( speedIndex ) ),
                    4.0 * speed * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testFlightPathAngle - sphericalOrbitState( flightPathIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testHeadingAngle - sphericalOrbitState( headingAngleIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );

        // Check consistency of back-and-forth conversion
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i ) - cartesianState( i ) ),
                        radius * 2.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i  + 3 ) - cartesianState( i + 3 ) ),
                        speed * 4.0 * std::numeric_limits< double >::epsilon( ) );
        }
    }

    // Test case 3: rotation with zero (heading, flight path, longitude) and half pi (latitude) angles.
    {
        // Define Cartesian state
        Eigen::Vector6d cartesianState;
        cartesianState << 0.0, 0.0, 6.498098097000000e3, -7.438147520000000e3, 0.0, 0.0;

        // Define associated spherical orbital state
        double testHeadingAngle = 0.0;
        double testFlightPathAngle = 0.0;
        double testLatitude = mathematical_constants::PI / 2.0;
        double testLongitude = 0.0;
        double radius = cartesianState.segment( 0, 3 ).norm( );
        double speed = cartesianState.segment( 3, 3 ).norm( );

        // Convert pack and forth
        Eigen::Vector6d sphericalOrbitState  = convertCartesianToSphericalOrbitalState(
                    cartesianState );
        Eigen::Vector6d reconvertedCartesianState  = convertSphericalOrbitalToCartesianState(
                    sphericalOrbitState );

        // Check computed spherical orbital state against reference data
        BOOST_CHECK_SMALL(
                    std::fabs( radius - sphericalOrbitState( radiusIndex ) ),
                    radius * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testLatitude - sphericalOrbitState( latitudeIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testLongitude - sphericalOrbitState( longitudeIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( speed - sphericalOrbitState( speedIndex ) ),
                    speed * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testFlightPathAngle - sphericalOrbitState( flightPathIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testHeadingAngle - sphericalOrbitState( headingAngleIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );

        // Check consistency of back-and-forth conversion
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i ) - cartesianState( i ) ),
                        radius * 2.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i  + 3 ) - cartesianState( i + 3 ) ),
                        speed * 2.0 * std::numeric_limits< double >::epsilon( ) );
        }
    }

    // Test case 4: rotation with zero (heading, longitude) and half pi (latitude, flight path) angles.
    {
        Eigen::Vector6d cartesianState;
        cartesianState << 0.0, 0.0, 6.498098097000000e3, 0.0, 0.0, -7.438147520000000e3;

        // Define associated spherical orbital state
        double testFlightPathAngle = -mathematical_constants::PI / 2.0;
        double testLatitude = mathematical_constants::PI / 2.0;
        double testLongitude = 0.0;
        double radius = cartesianState.segment( 0, 3 ).norm( );
        double speed = cartesianState.segment( 3, 3 ).norm( );

        // Convert pack and forth
        Eigen::Vector6d sphericalOrbitState  = convertCartesianToSphericalOrbitalState(
                    cartesianState );
        Eigen::Vector6d reconvertedCartesianState  = convertSphericalOrbitalToCartesianState(
                    sphericalOrbitState );

        // Check computed spherical orbital state against reference data (heading angle not tested, because close to
        // undefined; value is numerically irrelevant for physical state).
        BOOST_CHECK_SMALL(
                    std::fabs( radius - sphericalOrbitState( radiusIndex ) ),
                    radius * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testLatitude - sphericalOrbitState( latitudeIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testLongitude - sphericalOrbitState( longitudeIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( speed - sphericalOrbitState( speedIndex ) ),
                    speed * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL(
                    std::fabs( testFlightPathAngle - sphericalOrbitState( flightPathIndex ) ),
                    2.0 * mathematical_constants::PI * std::numeric_limits< double >::epsilon( ) );


        // Check consistency of back-and-forth conversion
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i ) - cartesianState( i ) ),
                        radius * 2.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i  + 3 ) - cartesianState( i + 3 ) ),
                        speed * 2.0 * std::numeric_limits< double >::epsilon( ) );
        }
    }

    // Test case 5: rotation with undefined heading angle.
    {
        Eigen::Vector6d cartesianState;
        cartesianState << 0.0, 0.0, 6.498098097000000e3, 0.0, 0.0, -7.438147520000000e3;

        // Define associated spherical orbital state
        Eigen::Vector6d sphericalOrbitalState;
        sphericalOrbitalState( headingAngleIndex )  = TUDAT_NAN;
        sphericalOrbitalState( flightPathIndex ) = -mathematical_constants::PI / 2.0;
        sphericalOrbitalState( latitudeIndex ) = mathematical_constants::PI / 2.0;
        sphericalOrbitalState( longitudeIndex ) = 0.0;
        sphericalOrbitalState( radiusIndex )  = cartesianState.segment( 0, 3 ).norm( );
        sphericalOrbitalState( speedIndex )  = cartesianState.segment( 3, 3 ).norm( );

        Eigen::Vector6d reconvertedCartesianState  = convertSphericalOrbitalToCartesianState(
                    sphericalOrbitalState );

        // Check consistency of back-and-forth conversion
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i ) - cartesianState( i ) ),
                        sphericalOrbitalState( radiusIndex ) * 2.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i  + 3 ) - cartesianState( i + 3 ) ),
                        sphericalOrbitalState( speedIndex ) * 2.0 * std::numeric_limits< double >::epsilon( ) );
        }
    }

    // Test case 5: rotation with undefined flight path angle.
    {
        Eigen::Vector6d cartesianState;
        cartesianState << 0.0, 0.0, 6.498098097000000e3, 0.0, 0.0, 0.0;

        // Define associated spherical orbital state
        Eigen::Vector6d sphericalOrbitalState;
        sphericalOrbitalState( headingAngleIndex )  = TUDAT_NAN;
        sphericalOrbitalState( flightPathIndex ) = TUDAT_NAN;
        sphericalOrbitalState( latitudeIndex ) = mathematical_constants::PI / 2.0;
        sphericalOrbitalState( longitudeIndex ) = 0.0;
        sphericalOrbitalState( radiusIndex )  = cartesianState.segment( 0, 3 ).norm( );
        sphericalOrbitalState( speedIndex )  = cartesianState.segment( 3, 3 ).norm( );

        Eigen::Vector6d reconvertedCartesianState  = convertSphericalOrbitalToCartesianState(
                    sphericalOrbitalState );

        // Check consistency of back-and-forth conversion
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i ) - cartesianState( i ) ),
                        sphericalOrbitalState( radiusIndex ) * 2.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( reconvertedCartesianState( i  + 3 ) - cartesianState( i + 3 ) ),
                        2.0 * std::numeric_limits< double >::epsilon( ) );
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat


