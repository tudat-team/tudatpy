/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      150411    D. Dirkx          Migrated and updated from personal code.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <fstream>
#include <limits>
#include <map>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Basics/testMacros.h>
#include <Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/UnitTests/keplerPropagatorTestData.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/InputOutput/basicInputOutput.h"

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
                    keplerEphemeris.getCartesianStateFromEphemeris( stateIterator->first ),
                    earthGravitationalParameter );

        // Check that computed results match expected results.
        BOOST_CHECK_CLOSE_FRACTION(
                    computedPropagationHistory[ stateIterator->first ]( 5 ),
                    expectedPropagationHistory[ stateIterator->first ]( 5 ),
                    2.0e-14 );
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
                    keplerEphemeris.getCartesianStateFromEphemeris( stateIterator->first ),
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

