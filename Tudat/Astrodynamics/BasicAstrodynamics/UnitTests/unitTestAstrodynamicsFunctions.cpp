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
 *      100910    K. Kumar          First creation of code.
 *      111115    K. Kumar          Altered tests to override default mass values; corrected
 *                                  bool-variable name; corrected cerr-statement.
 *      111122    K. Kumar          Altered tests to use reference data.
 *      111206    K. Kumar          Updated synodic period test and comments.
 *      120127    D. Dirkx          Moved to Tudat core.
 *      120127    K. Kumar          Transferred unit tests over to Boost unit test framework.
 *      120128    K. Kumar          Changed BOOST_CHECK to BOOST_CHECK_CLOSE_FRACTION for unit test
 *                                  comparisons.
 *
 *    References
 *      Wikipedia. Geostationary orbit, http://en.wikipedia.org/wiki/Geostationary_orbit, last
 *          accessed: 22nd November, 2011.
 *      Keefe, T.J. Synodic Period Calculator, http://www.ccri.edu/physics/keefe/synodic_calc.htm,
 *          last accessed: 6th December, 2011, last modified: 18th November, 2011.
 *
 *    Notes
 *      The tests need to be updated to check benchmark values from literature.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{
namespace unit_tests
{

//! Test suite for astrodynamics functions.
BOOST_AUTO_TEST_SUITE( test_astrodynamics_functions )

//! Test if the orbital period of a Kepler orbit is computed correctly.
BOOST_AUTO_TEST_CASE( testKeplerOrbitalPeriod )
{
    // Declare and set satellite mass [kg].
    double satelliteMass = 1.0e3;

    // Declare and set gravitational parameter of Earth [m^3 s^-2].
    double earthGravitationalParameter
            = physical_constants::GRAVITATIONAL_CONSTANT * 5.9736e24;

    // Declare and set distance between Earth center and satellite.
    double distanceBetweenSatelliteAndEarth = 4.2164e7;

    // Compute orbital period of satellite.
    double orbitalPeriod = basic_astrodynamics::computeKeplerOrbitalPeriod(
                distanceBetweenSatelliteAndEarth, earthGravitationalParameter, satelliteMass );

    // Declare and set expected orbital period [s].
    double expectedOrbitalPeriod = 86164.09054;

    // Check if computed orbital period matches expected orbital period.
    BOOST_CHECK_CLOSE_FRACTION( orbitalPeriod, expectedOrbitalPeriod, 1.0e-5 );
}

//! Test if the orbital angular momentum of a kepler orbit is computed correctly.
BOOST_AUTO_TEST_CASE( testKeplerAngularMomentum )
{
    // Reference: http://en.wikipedia.org/wiki/Geostationary_orbit.
    // Declare and set satellite mass [kg].
    double satelliteMass = 1.0e3;

    // Declare and set gravitational parameter of Earth [m^3 s^-2].
    double earthGravitationalParameter
            = physical_constants::GRAVITATIONAL_CONSTANT * 5.9736e24;

    // Declare and set distance between Earth center and satellite.
    double distanceBetweenSatelliteAndEarth = 4.2164e7;

    // Declare and set eccentricity of satellite orbit.
    double eccentricityOfSatelliteOrbit = 0.0;

    // Compute Kepler angular momentum.
    double angularMomentum = basic_astrodynamics::computeKeplerAngularMomentum(
                distanceBetweenSatelliteAndEarth, eccentricityOfSatelliteOrbit,
                earthGravitationalParameter, satelliteMass );

    // Declare and set expected angular momentum.
    // The expected angular momentum is computed using the fact that for a circular orbit,
    // H = mRV. This is an independent check of the code, which computes angular
    // momentum differently.
    double expectedAngularMomentum = satelliteMass * distanceBetweenSatelliteAndEarth
            * std::sqrt( earthGravitationalParameter / distanceBetweenSatelliteAndEarth );

    // Check if computed angular momentum matches expected angular momentum.
    BOOST_CHECK_CLOSE_FRACTION( angularMomentum, expectedAngularMomentum,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if the mean motion of a kepler orbit is computed correctly.
BOOST_AUTO_TEST_CASE( testMeanMotion )
{
    // Declare and set satellite mass [kg].
    double satelliteMass = 1.0e3;

    // Declare and set gravitational parameter of Earth [m^3 s^-2].
    double earthGravitationalParameter
            = physical_constants::GRAVITATIONAL_CONSTANT * 5.9736e24;

    // Declare and set distance between Earth center and satellite.
    double distanceBetweenSatelliteAndEarth = 4.2164e7;

    // Reference: http://en.wikipedia.org/wiki/Geostationary_orbit.
    double meanMotion = basic_astrodynamics::computeKeplerMeanMotion(
                distanceBetweenSatelliteAndEarth, earthGravitationalParameter, satelliteMass );

    // Declare and set expected mean motion [rad/s].
    double expectedMeanMotion = 7.2921e-5;

    // Check if computed mean motion matches expected mean motion.
    BOOST_CHECK_CLOSE_FRACTION( meanMotion, expectedMeanMotion, 1.0e-7 );
}

//! Test if the orbital energy of a Kepler orbit is computed correctly.
BOOST_AUTO_TEST_CASE( testKeplerEnergy )
{
    // Declare and set satellite mass [kg].
    double satelliteMass = 1.0e3;

    // Declare and set gravitational parameter of Earth [m^3 s^-2].
    double earthGravitationalParameter
            = physical_constants::GRAVITATIONAL_CONSTANT * 5.9736e24;

    // Declare and set distance between Earth center and satellite.
    double distanceBetweenSatelliteAndEarth = 4.2164e7;

    // Compute Kepler energy.
    double orbitalEnergy = basic_astrodynamics::computeKeplerEnergy(
                distanceBetweenSatelliteAndEarth, earthGravitationalParameter, satelliteMass );

    // Declare and set expected orbital energy.
    // The expected orbital energy is computed using the fact that for a circular orbit,
    // E = m ( V^2/2 - mu/R ). This is an independent check of the code, which computes orbital
    // energy differently.
    double expectedOrbitalEnergy = satelliteMass * (
                0.5 * earthGravitationalParameter / distanceBetweenSatelliteAndEarth
                -  earthGravitationalParameter / distanceBetweenSatelliteAndEarth );

    // Check if computed orbital energy matches expected orbital energy.
    BOOST_CHECK_CLOSE_FRACTION( orbitalEnergy, expectedOrbitalEnergy,
                                std::numeric_limits< double >::epsilon( ) );

}

//! Test if the synodic period between two orbits is computed correctly.
BOOST_AUTO_TEST_CASE( testSynodicPeriod )
{
    // Compute synodic period between Earth and Mars.
    double synodicPeriod = basic_astrodynamics::computeSynodicPeriod( 365.256378, 686.95 );

    // Declare and set expected synodic period.
    double expectedSynodicPeriod = 779.9746457736733;

    // Check if computed synodic period matches expected synodic period.
    BOOST_CHECK_CLOSE_FRACTION( synodicPeriod, expectedSynodicPeriod,
                                std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
