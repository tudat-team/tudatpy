/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      110216    K. Kumar          Creation of code.
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
 *
 *    References
 *      Melman, J. Propagate software, J.C.P.Melman@tudelft.nl, 2010.
 *      NASA, Goddard Spaceflight Center. Orbit Determination Toolbox (ODTBX), NASA - GSFC Open
 *          Source Software, http://opensource.gsfc.nasa.gov/projects/ODTBX/, last accessed:
 *          31st January, 2012.
 *
 */

#define BOOST_TEST_MAIN

#include <fstream>
#include <limits>
#include <map>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Basics/testMacros.h>
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

//! Typedef for propagation history.
typedef std::map < double, Eigen::VectorXd > PropagationHistory;

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
    Eigen::VectorXd stateInCartesianElements( 6 );
    Eigen::VectorXd stateInKeplerianElements( 6 );

    stateInCartesianElements << 6.75e6, 0.0, 0.0, 0.0, 8.0595973215e3, 0.0;
    stateInKeplerianElements = orbital_element_conversions::convertCartesianToKeplerianElements(
                stateInCartesianElements, getMelmanEarthGravitationalParameter( ) );
    benchmarkPropagationHistory[ 0.0 ] = stateInKeplerianElements;

    stateInCartesianElements << -6.1318272067e6, 5.1974105627e6, 0.0,
            -4.7375063953e3, -4.8565484865e3, 0.0;
    stateInKeplerianElements = orbital_element_conversions::convertCartesianToKeplerianElements(
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
    Eigen::VectorXd stateInKeplerianElements( 6 );

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

//! Test if the orbital period of a Kepler orbit is computed correctly.
BOOST_AUTO_TEST_CASE( testPropagateKeplerOrbit )
{
    // Case 1: Comparison of propagateKeplerOrbit() output with benchmark data from (Melman, 2010).
    {
        // Load benchmark data.
        // This data originates from J. Melman and is generated by the software package Propagate.

        // Create propagation history map for benchmark data to be stored in.
        PropagationHistory benchmarkKeplerPropagationHistory = getMelmanBenchmarkData( );

        // Propagate to final state in Keplerian elements.
        Eigen::VectorXd computedFinalStateInKeplerianElements
                = basic_astrodynamics::orbital_element_conversions::propagateKeplerOrbit(
                    benchmarkKeplerPropagationHistory.begin( )->second,
                    benchmarkKeplerPropagationHistory.rbegin( )->first -
                    benchmarkKeplerPropagationHistory.begin( )->first,
                    getMelmanEarthGravitationalParameter( ), 1.0e-10 );

        // Check that computed results match expected results.
        BOOST_CHECK_CLOSE_FRACTION(
                    benchmarkKeplerPropagationHistory.rbegin( )->second( 5 ),
                    mathematics::computeModulo( computedFinalStateInKeplerianElements( 5 ),
                                                2.0 * mathematics::PI ), 1.0e-8 );
    }

    // Case 2: Comparison of kepprop2b() test output from (GSFC, 2012) using modulo option.
    {
        // Create expected propagation history.
        PropagationHistory expectedPropagationHistory = getODTBXBenchmarkData( );

        // Set Earth gravitational parameter [m^3 s^-2].
        double earthGravitationalParameter = 398600.4415e9;

        // Set time step for ODTBX benchmark data.
        double timeStep = 8640.0;

        // Compute propagation history.
        PropagationHistory computedPropagationHistory;
        computedPropagationHistory[ 0.0 ] = expectedPropagationHistory[ 0.0 ];

        for ( unsigned int i = 1; i < expectedPropagationHistory.size( ); i++ )
        {
            computedPropagationHistory[ static_cast< double >( i ) * timeStep ]
                    = basic_astrodynamics::orbital_element_conversions::propagateKeplerOrbit(
                        computedPropagationHistory[ static_cast< double >( i - 1 ) * timeStep ],
                        timeStep, earthGravitationalParameter,
                        1.0e-10 );

            computedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 )
                    = mathematics::computeModulo(
                        computedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 ),
                        2.0 * mathematics::PI );

            // Check that computed results match expected results.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        computedPropagationHistory[ static_cast< double >( i ) * timeStep ],
                        expectedPropagationHistory[ static_cast< double >( i ) * timeStep ],
                        1.0e-13 );
        }
    }

    // Case 3: Comparison of kepprop2b() test output from (GSFC, 2012) without using modulo option.
    {
        // Create expected propagation history.
        PropagationHistory expectedPropagationHistory = getODTBXBenchmarkData( );

        // Set Earth gravitational parameter [m^3 s^-2].
        double earthGravitationalParameter = 398600.4415e9;

        // Set time step for ODTBX benchmark data.
        double timeStep = 8640.0;

        // Compute propagation history.
        PropagationHistory computedPropagationHistory;
        computedPropagationHistory[ 0.0 ] = expectedPropagationHistory[ 0.0 ];

        for ( unsigned int i = 1; i < expectedPropagationHistory.size( ); i++ )
        {
            computedPropagationHistory[ static_cast< double >( i ) * timeStep ]
                    = basic_astrodynamics::orbital_element_conversions::propagateKeplerOrbit(
                        computedPropagationHistory[ static_cast< double >( i - 1 ) * timeStep ],
                        timeStep, earthGravitationalParameter, false );

            computedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 )
                    = mathematics::computeModulo(
                        computedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 ),
                        2.0 * mathematics::PI );

            // Check that computed results match expected results.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        computedPropagationHistory[ static_cast< double >( i ) * timeStep ],
                        expectedPropagationHistory[ static_cast< double >( i ) * timeStep ],
                        1.0e-13 );
        }
    }

    // Case 4: Comparison of kepprop2b() test output from (GSFC, 2012), propagating backwards.
    {
        // Create expected propagation history.
        PropagationHistory expectedPropagationHistory = getODTBXBenchmarkData( );

        // Set Earth gravitational parameter [m^3 s^-2].
        double earthGravitationalParameter = 398600.4415e9;

        // Set time step for ODTBX benchmark data.
        double timeStep = 8640.0;

        // Compute propagation history.
        PropagationHistory computedPropagationHistory;
        computedPropagationHistory[ 10.0 * 8640.0 ] = expectedPropagationHistory[ 10.0 * 8640.0 ];

        for ( int i = expectedPropagationHistory.size( ) - 2; i >= 0; i-- )
        {
            computedPropagationHistory[ static_cast< double >( i ) * timeStep ]
                    = basic_astrodynamics::orbital_element_conversions::propagateKeplerOrbit(
                        computedPropagationHistory[ static_cast< double >( i + 1 ) * timeStep ],
                        -timeStep, earthGravitationalParameter, false );

            computedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 )
                    = mathematics::computeModulo(
                        computedPropagationHistory[ static_cast< double >( i ) * timeStep ]( 5 ),
                        2.0 * mathematics::PI );

            // Check that computed results match expected results.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        computedPropagationHistory[ static_cast< double >( i ) * timeStep ],
                        expectedPropagationHistory[ static_cast< double >( i ) * timeStep ],
                        1.0e-13 );
        }
    }
}

} // namespace unit_tests
} // namespace tudat
