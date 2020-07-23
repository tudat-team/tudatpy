/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <map>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/keplerPropagatorTestData.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_keplerEphemeris )

//! Test 1: Comparison of KeplerEphemeris output with benchmark data from (Melman, 2010).
//! (see testPropagateKeplerOrbit_Eccentric_Melman).
BOOST_AUTO_TEST_CASE( testKeplerEphemerisElliptical )
{
    // Load the expected propagation history.
    // Create expected propagation history.
    PropagationHistory expectedPropagationHistory = getODTBXBenchmarkData( );

    // Set Earth gravitational parameter [m^3 s^-2].
    const double earthGravitationalParameter = 398600.4415e9;

    // Compute propagation history.
    PropagationHistory computedPropagationHistory;
    computedPropagationHistory[ 0.0 ] = expectedPropagationHistory[ 0.0 ];

    ephemerides::KeplerEphemeris keplerEphemeris(
                expectedPropagationHistory[ 0.0 ],
                0.0, earthGravitationalParameter );

    for( PropagationHistory::iterator stateIterator = expectedPropagationHistory.begin( );
         stateIterator != expectedPropagationHistory.end( ); stateIterator++ )
    {
        // Compute next entry.
        computedPropagationHistory[ stateIterator->first ] =
                orbital_element_conversions::convertCartesianToKeplerianElements(
                    keplerEphemeris.getCartesianState( stateIterator->first ),
                    earthGravitationalParameter );

        // Check that computed results match expected results.
        BOOST_CHECK_CLOSE_FRACTION(
                    computedPropagationHistory[ stateIterator->first ]( 5 ),
                    expectedPropagationHistory[ stateIterator->first ]( 5 ),
                    2.5e-14 );
    }
}

//! Test 2: Comparison of KeplerEphemeris with that of GTOP (hyperbolic).
//! (see testPropagateKeplerOrbit_hyperbolic_GTOP).
BOOST_AUTO_TEST_CASE( testKeplerEphemerisHyperbolic )
{
    // Load the expected propagation history.
    PropagationHistory expectedPropagationHistory = getGTOPBenchmarkData( );

    // Compute propagation history.
    PropagationHistory computedPropagationHistory;
    computedPropagationHistory[ 0.0 ] = expectedPropagationHistory[ 0.0 ];

    ephemerides::KeplerEphemeris keplerEphemeris(
                expectedPropagationHistory[ 0.0 ],
                0.0, getGTOPGravitationalParameter( ) );

    for( PropagationHistory::iterator stateIterator = expectedPropagationHistory.begin( );
         stateIterator != expectedPropagationHistory.end( ); stateIterator++ )
    {
        // Compute next entry.
        computedPropagationHistory[ stateIterator->first ] =
                orbital_element_conversions::convertCartesianToKeplerianElements(
                    keplerEphemeris.getCartesianState( stateIterator->first ),
                    getGTOPGravitationalParameter( ) );

        // Check that computed results match expected results.
        BOOST_CHECK_CLOSE_FRACTION(
                    computedPropagationHistory[ stateIterator->first ]( 5 ),
                    expectedPropagationHistory[ stateIterator->first ]( 5 ),
                    1.0e-15 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat

