/*    Copyright (c) 2010-2023, Delft University of Technology
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

#include "tudat/basics/testMacros.h"
#include "tudat/astro/gravitation/ringGravityField.h"
#include "tudat/astro/basic_astro/physicalConstants.h"

#include "tudat/astro/gravitation/centralGravityModel.h"

namespace tudat
{
namespace unit_tests
{

//! Test the functionality of the polyhedron gravity field class.
BOOST_AUTO_TEST_SUITE( test_ring_gravity_field )

//! Test getters
BOOST_AUTO_TEST_CASE( testSettingAndGettingParameters )
{
    const double tolerance = 1.0E-14;

    // Define ring gravity parameters
    const double gravitationalParameter = 2.39e21 * physical_constants::GRAVITATIONAL_CONSTANT;
    const double ringRadius = 2.7 * physical_constants::ASTRONOMICAL_UNIT;
    const bool ellipticIntegralSFromDAndB = false;

    gravitation::RingGravityField gravityField = gravitation::RingGravityField(
            gravitationalParameter, ringRadius, ellipticIntegralSFromDAndB);

    // Gravitational parameter
    double retrievedGravitationalParameter = gravityField.getGravitationalParameter();
    BOOST_CHECK_CLOSE_FRACTION( gravitationalParameter, retrievedGravitationalParameter, tolerance );

    // Ring radius
    double retrievedRadius = gravityField.getRingRadius();
    BOOST_CHECK_CLOSE_FRACTION( ringRadius, retrievedRadius, tolerance );

    // Elliptic integrals flag
    bool retrievedEllipticIntegralSFromDAndB = gravityField.getEllipticIntegralSFromDAndB();
    BOOST_CHECK_EQUAL( ellipticIntegralSFromDAndB, retrievedEllipticIntegralSFromDAndB );

}

// Test computation of elliptic integrals
BOOST_AUTO_TEST_CASE( testEllipticIntegralsComputation )
{
    const double tolerance = 1.0E-14;

    // Check values of integrals for m=0
    gravitation::RingGravityCache gravityCache1 = gravitation::RingGravityCache( 1.0, true);

    Eigen::Vector3d bodyFixedPosition = Eigen::Vector3d::Zero();
    gravityCache1.update( bodyFixedPosition );

    BOOST_CHECK_CLOSE_FRACTION(
            mathematical_constants::PI / 2.0, gravityCache1.getEllipticIntegralK(), tolerance );
    BOOST_CHECK_CLOSE_FRACTION(
            mathematical_constants::PI / 4.0, gravityCache1.getEllipticIntegralB(), tolerance );
    BOOST_CHECK_CLOSE_FRACTION(
            mathematical_constants::PI / 2.0, gravityCache1.getEllipticIntegralE(), tolerance );
    BOOST_CHECK_CLOSE_FRACTION(
            mathematical_constants::PI / 16.0, gravityCache1.getEllipticIntegralS(), tolerance );


    // Check consistency of computation of S integral between different methods for m in ]0,0.1]
    // Ring radius and position selected such that: m = 0.039211841976276834
    bodyFixedPosition << 100.0, 0.0, 0.0;
    double ringRadius = 1.0;
    double m = 0.039211841976276834;

    // Compute expected S(m) via taylor series (sec. A.1 of Fukushima, 2010; valid for [0,0.1])
    double expectedEllipticIntegralS = 0;
    std::vector< double > taylorCoefficients = {
                    0.204012532440038310, 0.159513582234205843, 0.130422818255893004, 0.111687838140976463,
                    0.098925188226691425, 0.089815348807960028, 0.083084759300136632, 0.077987984857306626,
                    0.074062924745595950, 0.071009059783923539, 0.068623059119746445, 0.066762755430661757,
                    0.065325983044110253 };
    double m0 = 0.05;
    for ( unsigned int taylorOrder = 0; taylorOrder <= 12; ++taylorOrder )
    {
        expectedEllipticIntegralS += taylorCoefficients.at( taylorOrder ) * std::pow( m - m0, taylorOrder );
    }

    for ( bool ellipticIntegralSFromDAndB : { true, false } )
    {
        gravitation::RingGravityCache gravityCache2 = gravitation::RingGravityCache(
                ringRadius, ellipticIntegralSFromDAndB );
        gravityCache2.update( bodyFixedPosition );
        
        if ( ellipticIntegralSFromDAndB )
        {
            BOOST_CHECK_CLOSE_FRACTION( expectedEllipticIntegralS, gravityCache2.getEllipticIntegralS( ), 2e-12 );
        }
        // Computation of S(m) via K(m) and E(m) is more sensitive to numerical cancellation, hence the larger tolerance.
        else
        {
            BOOST_CHECK_CLOSE_FRACTION( expectedEllipticIntegralS, gravityCache2.getEllipticIntegralS( ), 2e-12 );
        }

    }


}

//! Test computation of potential, gradient of potential
BOOST_AUTO_TEST_CASE( testGravityComputation )
{
    const double tolerance = 1.0E-14;

    // Define ring gravity parameters
    const double gravitationalParameter = 2.39e21 * physical_constants::GRAVITATIONAL_CONSTANT;
    const double ringRadius = 2.7 * physical_constants::ASTRONOMICAL_UNIT;

    for ( bool ellipticIntegralSFromDAndB : { true, false } )
    {

        gravitation::RingGravityField gravityField = gravitation::RingGravityField(
            gravitationalParameter, ringRadius, ellipticIntegralSFromDAndB );

        Eigen::Vector3d bodyFixedPosition;
        double expectedPotential;
        Eigen::Vector3d expectedGradient;

        // Point at the origin
        (bodyFixedPosition << 0.0, 0.0, 0.0).finished();
        expectedPotential = gravitationalParameter / ringRadius;
        (expectedGradient << 0.0, 0.0, 0.0).finished();

        BOOST_CHECK_CLOSE_FRACTION(
                expectedPotential, gravityField.getGravitationalPotential( bodyFixedPosition ), tolerance );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                expectedGradient, gravityField.getGradientOfPotential( bodyFixedPosition ), tolerance );

        // Point in x=y=0 line. Check potential and direction of acceleration
        (bodyFixedPosition << 0.0, 0.0, 500.0e3).finished();
        expectedPotential = gravitationalParameter / std::sqrt( std::pow(ringRadius, 2) + std::pow( bodyFixedPosition(2), 2 ) );
        (expectedGradient << 0.0, 0.0, TUDAT_NAN).finished();

        BOOST_CHECK_CLOSE_FRACTION(
                expectedPotential, gravityField.getGravitationalPotential( bodyFixedPosition ), tolerance );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                expectedGradient.segment(0, 2), gravityField.getGradientOfPotential( bodyFixedPosition ).segment(0, 2), tolerance );
        // Check if z component of acceleration is negative
        BOOST_CHECK_EQUAL(
                true, std::signbit( gravityField.getGradientOfPotential( bodyFixedPosition )(2) ) );

        // Point in z=0 plane inside ring. Check direction of acceleration
        (bodyFixedPosition << 100.0e3, 500.0e3, 0.0).finished();
        Eigen::Vector3d computedGradient = gravityField.getGradientOfPotential( bodyFixedPosition );
        BOOST_CHECK_CLOSE_FRACTION(
                bodyFixedPosition(0) / bodyFixedPosition(1),
                computedGradient(0) / computedGradient(1), tolerance );

        // Point in z=0 plane outside ring. Check direction of acceleration
        (bodyFixedPosition << - 2 * ringRadius, 1.3 * ringRadius, 0.0).finished();
        computedGradient = gravityField.getGradientOfPotential( bodyFixedPosition );
        BOOST_CHECK_CLOSE_FRACTION(
                bodyFixedPosition(0) / bodyFixedPosition(1),
                computedGradient(0) / computedGradient(1), tolerance );
        BOOST_CHECK_EQUAL(
                std::signbit( - bodyFixedPosition(0) ), std::signbit( computedGradient(0) ) );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace tudat
} // namespace unit_tests
