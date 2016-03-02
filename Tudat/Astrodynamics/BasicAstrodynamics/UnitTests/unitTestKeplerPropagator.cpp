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
 *      110216    K. Kumar          File created.
 *      110217    E. Iorfida        Minor changes made.
 *      110221    K. Kumar          Updated variable-naming to comply with protocol.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 *      120215    K. Kumar          Updated to use new propagateKeplerOrbit() function;
 *                                  Boostified unit test; added new unit tests based on ODTBX.
 *      120217    K. Kumar          Updated computeModuloForSignedValues() to computeModulo() from
 *                                  Tudat Core.
 *      120225    K. Kumar          Updated unit test using data from (Melman, 2010) to fix bug in
 *                                  Linux; data is no longer imported from input text file.
 *      120508    K. Kumar          Corrected bug in backwards propagation loop end-condition.
 *      120607    P. Musegaas       Updated unit test to new interface.
 *      120813    P. Musegaas       Updated unit test to new root finding structure.
 *      120823    P. Musegaas       Separated existing unit tests. Added unit test for hyperbolic
 *                                  kepler propagation, with accompanying data from GTOP.
 *      120904    P. Musegaas       Added unit test of a case that failed on the old modulo
 *                                  function in the Kepler propagator. Removed dedicated modulo
 *                                  unit test.
 *      121205    P. Musegaas       Updated code to final version of rootfinders.
 *      150417    D. Dirkx          Made modifications for templated element conversions.
 *
 *    References
 *      Melman, J. Propagate software, J.C.P.Melman@tudelft.nl, 2010.
 *      NASA, Goddard Spaceflight Center. Orbit Determination Toolbox (ODTBX), NASA - GSFC Open
 *          Source Software, http://opensource.gsfc.nasa.gov/projects/ODTBX/, last accessed:
 *          31st January, 2012.
 *      ESA, GTOP Toolbox, http://www.esa.int/gsp/ACT/doc/INF/Code/globopt/GTOPtoolbox.rar.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"

#include <Eigen/Core>

#include <map>
#include <limits>


#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{
namespace unit_tests
{

using namespace orbital_element_conversions;

//! Typedef for propagation history.
typedef std::map < double, basic_mathematics::Vector6d > PropagationHistory;

//! Get Earth gravitational parameter for benchmark data from (Melman, 2010).
double getMelmanEarthGravitationalParameter( )
{
    // Return Earth gravitational parameter [m^3 s^-2].
    return 3.986004415e14;
}

//! Get benchmark data from (Melman, 2010).
PropagationHistory getMelmanBenchmarkData( )
{
    // Declare benchmark pragation history.
    PropagationHistory benchmarkPropagationHistory;

    // Populate benchmark propagation history.
    basic_mathematics::Vector6d stateInCartesianElements;
    basic_mathematics::Vector6d stateInKeplerianElements;

    stateInCartesianElements << 6.75e6, 0.0, 0.0, 0.0, 8.0595973215e3, 0.0;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getMelmanEarthGravitationalParameter( ) );
    benchmarkPropagationHistory[ 0.0 ] = stateInKeplerianElements;

    stateInCartesianElements << -6.1318272067e6, 5.1974105627e6, 0.0,
            -4.7375063953e3, -4.8565484865e3, 0.0;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getMelmanEarthGravitationalParameter( ) );
    benchmarkPropagationHistory[ 86400.0 ] = stateInKeplerianElements;

    // Return benchmark data.
    return benchmarkPropagationHistory;
}

//! Get ODTBX benchmark data.
PropagationHistory getODTBXBenchmarkData( )
{
    // Declare benchmark pragation history.
    PropagationHistory benchmarkPropagationHistory;

    // Populate benchmark propagation history.
    basic_mathematics::Vector6d stateInKeplerianElements;

    stateInKeplerianElements << 42165.3431351313e3, 0.26248354351331, 0.30281462522101,
            4.71463172847351, 4.85569272927819, 2.37248926702153;

    // Set time step.
    double timeStep = 8640.0;

    for ( unsigned int i = 0; i < 11; i++ )
    {
        benchmarkPropagationHistory[ static_cast< double >( i ) * timeStep ]
                = stateInKeplerianElements;
    }

    benchmarkPropagationHistory[ 1.0 * timeStep ]( 5 ) = 2.79722436211144;
    benchmarkPropagationHistory[ 2.0 * timeStep ]( 5 ) = 3.18337407409023;
    benchmarkPropagationHistory[ 3.0 * timeStep ]( 5 ) = 3.57400974200765;
    benchmarkPropagationHistory[ 4.0 * timeStep ]( 5 ) = 4.01425565759545;
    benchmarkPropagationHistory[ 5.0 * timeStep ]( 5 ) = 4.57232665706546;
    benchmarkPropagationHistory[ 6.0 * timeStep ]( 5 ) = 5.35956850972672;
    benchmarkPropagationHistory[ 7.0 * timeStep ]( 5 ) = 0.137251905665217;
    benchmarkPropagationHistory[ 8.0 * timeStep ]( 5 ) = 1.14521863765007;
    benchmarkPropagationHistory[ 9.0 * timeStep ]( 5 ) = 1.86433634881636;
    benchmarkPropagationHistory[ 10.0 * timeStep ]( 5 ) = 2.38486787064101;

    return benchmarkPropagationHistory;
}

//! Get GTOP gravitational parameter for benchmark data from GTOP.
double getGTOPGravitationalParameter( )
{
    // Return inaccurate gravitational parameter of the Sun [m^3 s^-2].
    return 1.327e20;
}

//! Get benchmark data from GTOP.
PropagationHistory getGTOPBenchmarkData( )
{
    // Declare benchmark pragation history.
    PropagationHistory benchmarkPropagationHistory;

    // Populate benchmark propagation history. Obtained by propagating an orbit starting at
    // x = 1.5e11, V_y = 6.0e4 for a period of 100 days four times consecutively. Since GTOP does
    // not work with similar Keplerian orbital elements, the cartesian elements resulting from that
    // are converted to Keplerian elements first. These initial starting coordinates correspond to
    // a semi major axis of -7.24873e+010 meters and an eccentricity of 3.06933.
    basic_mathematics::Vector6d stateInCartesianElements, stateInKeplerianElements;

    stateInCartesianElements << 1.5e11, 0.0, 0.0, 0.0, 6.0e4, 0.0;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getGTOPGravitationalParameter( ) );
    benchmarkPropagationHistory[ 0.0 ] = stateInKeplerianElements;

    stateInCartesianElements << 50369576778.98602, 453006898372.5074, 2.156946592732799e-005,
            -14654.13750690802, 46884.94068619227, 4.665334803219454e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getGTOPGravitationalParameter( ) );
    benchmarkPropagationHistory[ 8640000.0 ] = stateInKeplerianElements;

    stateInCartesianElements << -76810236076.38216, 842661848023.4473, 6.100268297443444e-005,
            -14683.57015580225, 43917.12010513522, 4.48721854707566e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getGTOPGravitationalParameter( ) );
    benchmarkPropagationHistory[ 8640000.0 * 2.0 ] = stateInKeplerianElements;

    stateInCartesianElements << -203052258817.9893, 1216808019495.603, 9.937145023346651e-005,
            -14543.34378775917, 42828.66589049961, 4.403399939593385e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getGTOPGravitationalParameter( ) );
    benchmarkPropagationHistory[ 8640000.0 * 3.0 ] = stateInKeplerianElements;

    stateInCartesianElements << -328225472457.8796, 1584186440047.591, 0.0001371949389119038,
            -14437.813524927732, 42264.20425643964, 4.355914471377053e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getGTOPGravitationalParameter( ) );
    benchmarkPropagationHistory[ 8640000.0 * 4.0 ] = stateInKeplerianElements;

    // Return benchmark data.
    return benchmarkPropagationHistory;
}

//! Test 1: Comparison of propagateKeplerOrbit() output with benchmark data from (Melman, 2010).
BOOST_AUTO_TEST_CASE( testPropagateKeplerOrbit_Eccentric_Melman )
{
    // Load benchmark data.
    // This data originates from J. Melman and is generated by the software package Propagate.

    // Create propagation history map for benchmark data to be stored in.
    PropagationHistory benchmarkKeplerPropagationHistory = getMelmanBenchmarkData( );

    // Propagate to final state in Keplerian elements.
    basic_mathematics::Vector6d computedFinalStateInKeplerianElements
            = propagateKeplerOrbit(
                benchmarkKeplerPropagationHistory.begin( )->second,
                benchmarkKeplerPropagationHistory.rbegin( )->first -
                benchmarkKeplerPropagationHistory.begin( )->first,
                getMelmanEarthGravitationalParameter( ) );

    // Check that computed results match expected results.
    BOOST_CHECK_CLOSE_FRACTION(
                benchmarkKeplerPropagationHistory.rbegin( )->second( 5 ),
                basic_mathematics::computeModulo( computedFinalStateInKeplerianElements( 5 ),
                                                  2.0 * mathematical_constants::PI ), 1.0e-8 );
}

//! Test 2: Comparison of kepprop2b() test output from (GSFC, 2012) using modulo option.
BOOST_AUTO_TEST_CASE( testPropagateKeplerOrbit_Eccentric_kepprop2b_modulo )
{
    // Create expected propagation history.
    PropagationHistory expectedPropagationHistory = getODTBXBenchmarkData( );

    // Set Earth gravitational parameter [m^3 s^-2].
    const double earthGravitationalParameter = 398600.4415e9;

    // Set time step for ODTBX benchmark data.
    const double timeStep = 8640.0;

    // Compute propagation history.
    PropagationHistory computedPropagationHistory;
    computedPropagationHistory[ 0.0 ] = expectedPropagationHistory[ 0.0 ];

    for ( unsigned int i = 1; i < expectedPropagationHistory.size( ); i++ )
    {
        computedPropagationHistory[ static_cast< double >( i ) * timeStep ]
                = propagateKeplerOrbit(
                    computedPropagationHistory[ static_cast< double >( i - 1 ) * timeStep ],
                    timeStep, earthGravitationalParameter  );

        computedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 )
                = basic_mathematics::computeModulo(
                    computedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 ),
                    2.0 * mathematical_constants::PI );

        // Check that computed results match expected results.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    computedPropagationHistory[ static_cast< double >( i ) * timeStep ],
                    expectedPropagationHistory[ static_cast< double >( i ) * timeStep ],
                    1.0e-13 );
    }
}

//! Test 3: Comparison of kepprop2b() test output from (GSFC, 2012), propagating backwards.
BOOST_AUTO_TEST_CASE( testPropagateKeplerOrbit_Eccentric_kepprop2b_backwards )
{
    // Create expected propagation history.
    PropagationHistory expectedPropagationHistory = getODTBXBenchmarkData( );

    // Set Earth gravitational parameter [m^3 s^-2].
    const double earthGravitationalParameter = 398600.4415e9;

    // Set time step for ODTBX benchmark data.
    const double timeStep = 8640.0;

    // Compute propagation history.
    PropagationHistory computedPropagationHistory;
    computedPropagationHistory[ 10.0 * 8640.0 ] = expectedPropagationHistory[ 10.0 * 8640.0 ];

    for ( int i = expectedPropagationHistory.size( ) - 2; i >= 0; i-- )
    {
        computedPropagationHistory[ static_cast< double >( i ) * timeStep ]
                = propagateKeplerOrbit(
                    computedPropagationHistory[ static_cast< double >( i + 1 ) * timeStep ],
                    -timeStep, earthGravitationalParameter );

        computedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 )
                = basic_mathematics::computeModulo(
                    computedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 ),
                    2.0 * mathematical_constants::PI );

        // Check that computed results match expected results.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    computedPropagationHistory[ static_cast< double >( i ) * timeStep ],
                    expectedPropagationHistory[ static_cast< double >( i ) * timeStep ],
                    1.0e-13 );
    }
}

//! Test 4: Comparison of hyperbolic kepler propagation with that of GTOP.
BOOST_AUTO_TEST_CASE( testPropagateKeplerOrbit_hyperbolic_GTOP )
{
    // Load the expected propagation history.
    PropagationHistory expectedPropagationHistory = getGTOPBenchmarkData( );

    // Set the time step for the GTOP benchmark data to 100 days.
    const double timeStep = 86400.0 * 100.0;

    // Compute propagation history.
    PropagationHistory computedPropagationHistory;
    computedPropagationHistory[ 0.0 ] = expectedPropagationHistory[ 0.0 ];

    for ( unsigned int i = 1; i < expectedPropagationHistory.size( ); i++ )
    {
        // Compute next entry.
        computedPropagationHistory[ static_cast< double >( i ) * timeStep ] =
                propagateKeplerOrbit(
                    computedPropagationHistory[ static_cast< double >( i - 1 ) * timeStep ],
                    timeStep, getGTOPGravitationalParameter( ) );

        // Check that computed results match expected results.
        BOOST_CHECK_CLOSE_FRACTION(
                    computedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 ),
                    expectedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 ),
                    1.0e-15 );
    }
}

//! Test 5: Unit test that failed on versions that caused the old modulo function to crash.
BOOST_AUTO_TEST_CASE( testPropagateKeplerOrbit_FunctionFailingOnOldModuloFunction )
{
    // Set expected true anomaly.
    const double expectedTrueAnomaly = -3.1245538487052089;

    // Set the propagation time.
    const double propagationTime = 8651869.8944624383;

    // Set the gravitational parameter (of the Sun).
    const double gravitationalParameter = 1.32712428e20;

    // Set initial Keplerian elements.
    basic_mathematics::Vector6d keplerElements;
    keplerElements << 56618890355.593132, 0.99961601437304082, 1.0238269559089248,
            3.1526292818328812, 1.5807574453453865, 3.1478950321924795;

    // Propagate Keplerian elements.
    keplerElements = propagateKeplerOrbit(
                keplerElements, propagationTime, gravitationalParameter );

    // Check that computed results match expected results.
    BOOST_CHECK_CLOSE_FRACTION( keplerElements( trueAnomalyIndex ),
                                expectedTrueAnomaly,
                                1.0e-15 );
}

//! Test 6. Propagation test using ODTBX test Kepler elements.
BOOST_AUTO_TEST_CASE( testMeanAnomalyAgainstMeanMotion )
{
    std::vector< double > doubleErrors;
    // Test using double parameters.
    {
        double gravitationalParameter = 398600.4415e9;
        basic_mathematics::Vector6d initialStateInKeplerianElements;

        initialStateInKeplerianElements << 42165.3431351313e3, 0.26248354351331, 0.30281462522101,
                4.71463172847351, 4.85569272927819, 2.37248926702153;
        double timeStep = 600.0;
        double meanMotion = std::sqrt( gravitationalParameter /
                                       std::pow( initialStateInKeplerianElements( 0 ), 3.0 ) );

        double initialMeanAnomaly = convertEccentricAnomalyToMeanAnomaly(
                    convertTrueAnomalyToEccentricAnomaly(
                        initialStateInKeplerianElements( 5 ), initialStateInKeplerianElements( 1 ) ),
                    initialStateInKeplerianElements( 1 ) );

        double propagationTime, propagatedMeanAnomaly;

        basic_mathematics::Vector6d propagatedKeplerElements;


        for( int i = -25; i < 26; i++ )
        {
            propagationTime = static_cast< double >( i ) * timeStep;
            propagatedKeplerElements = propagateKeplerOrbit(
                        initialStateInKeplerianElements, propagationTime, gravitationalParameter );
            propagatedMeanAnomaly = convertEccentricAnomalyToMeanAnomaly(
                        convertTrueAnomalyToEccentricAnomaly(
                            propagatedKeplerElements( 5 ), initialStateInKeplerianElements( 1 ) ),
                        initialStateInKeplerianElements( 1 ) );
            doubleErrors.push_back( meanMotion * propagationTime - ( propagatedMeanAnomaly - initialMeanAnomaly ) );
        }
    }

    std::vector< double > longDoubleErrors;
    // Test using long double parameters.
    {
        long double gravitationalParameter = 398600.4415e9L;
        Eigen::Matrix< long double, 6, 1 > initialStateInKeplerianElements;

        initialStateInKeplerianElements << 42165.3431351313e3L, 0.26248354351331L, 0.30281462522101L,
                4.71463172847351L, 4.85569272927819L, 2.37248926702153L;
        long double timeStep = 600.0L;
        long double meanMotion = std::sqrt( gravitationalParameter /
                                            ( initialStateInKeplerianElements( 0 ) *
                                              initialStateInKeplerianElements( 0 ) *
                                              initialStateInKeplerianElements( 0 ) ) );

        long double initialMeanAnomaly = convertEccentricAnomalyToMeanAnomaly< long double >(
                    convertTrueAnomalyToEccentricAnomaly< long double >(
                        initialStateInKeplerianElements( 5 ), initialStateInKeplerianElements( 1 ) ),
                    initialStateInKeplerianElements( 1 ) );

        long double propagationTime, propagatedMeanAnomaly;

        Eigen::Matrix< long double, 6, 1 > propagatedKeplerElements;


        for( int i = -25; i < 26; i++ )
        {
            propagationTime = static_cast< long double >( i ) * timeStep;

            propagatedKeplerElements = propagateKeplerOrbit< long double >(
                        initialStateInKeplerianElements, propagationTime, gravitationalParameter );
            propagatedMeanAnomaly = convertEccentricAnomalyToMeanAnomaly< long double >(
                        convertTrueAnomalyToEccentricAnomaly< long double >(
                            propagatedKeplerElements( 5 ), initialStateInKeplerianElements( 1 ) ),
                        initialStateInKeplerianElements( 1 ) );

            longDoubleErrors.push_back( static_cast< double >(
                                            meanMotion * propagationTime -
                                            ( propagatedMeanAnomaly - initialMeanAnomaly ) ) );
        }
    }

    for( unsigned int i = 0; i < doubleErrors.size( ); i++ )
    {
        if( std::fabs( doubleErrors.at( i ) ) > 0.0 )
        {
            BOOST_CHECK_SMALL(
                        static_cast< double >( std::fabs( longDoubleErrors.at( i ) / doubleErrors.at( i ) ) ),
                        static_cast< double >( 5.0 * std::numeric_limits< long double >::epsilon( ) /
                        std::numeric_limits< double >::epsilon( ) ) );
        }
        else
        {
            BOOST_CHECK_SMALL( static_cast< double >( longDoubleErrors.at( i ) ),
                               static_cast< double >( 5.0 * std::numeric_limits< long double >::epsilon( ) ) );
        }
    }
}
} // namespace unit_tests
} // namespace tudat
