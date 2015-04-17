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
 *      111117    L. van der Ham    File created.
 *      111122    K. Kumar          Minor comment changes; include-statements updated, updated
 *                                  tests.
 *      120306    L. van der Ham    Correction of some small typos.
 *      120307    K. Kumar          Boostified unit tests; moved file.
 *
 *    References
 *      Koon, W.S., et al. Dynamical Systems, the Three-Body Problem and Space Mission Design,
 *          2006.
 *      Jet Propulsion Laboratory, NASA. Planets and Pluto: Physical Characteristics.
 *          http://ssd.jpl.nasa.gov/?planet_phys_par, last updated: 5th November, 2008, last
 *          accessed: 22nd November, 2011.
 *      Jet Propulsion Laboratory, NASA. Astrodynamics Constants.
 *          http://ssd.jpl.nasa.gov/?constants, last updated: 6th September, 2011, last accessed:
 *          22nd November, 2011.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"

namespace tudat
{
namespace unit_tests
{

//! Test if Cartesian state conversion for CRTBP is computed correctly.
BOOST_AUTO_TEST_CASE( testCartesianStateConversionCircularRestrictedThreeBodyProblem )
{
    // Test 1: conversion of dimensionless to dimensional Cartesian state vector for test particle
    // in Sun-Jupiter system [m, m/s] (Koon, 2006).
    {
        namespace crtbp = gravitation::circular_restricted_three_body_problem;

        // Set distance between Sun and Jupiter [m] [NASA, 2010].
        double distanceSunJupiter = 7.784e11;

        // Set magnitude of velocity of Sun-Jupiter system [m/s].
        double velocitySunJupiter = 13.102e3;

        // Set gravitational parameter of Jupiter [m^3 s^-2].
        double gravitationalParameterJupiter = 6.67259e-11 * 1898.13e24;

        // Set gravitational parameter of Sun [m^3 s^-2].
        double gravitationalParameterSun = 1.32712440018e20;

        // Set mass parameter of Sun-Jupiter system.
        double massParameter = gravitationalParameterJupiter
                / ( gravitationalParameterJupiter + gravitationalParameterSun );

        // Set dimensionless state of the Sun.
        using orbital_element_conversions::xCartesianPositionIndex;
        Eigen::VectorXd dimensionlessStateOfSun = Eigen::VectorXd::Zero( 6 );
        dimensionlessStateOfSun( xCartesianPositionIndex ) = -massParameter;

        // Declare and set dimensionless state of Jupiter.
        Eigen::VectorXd dimensionlessStateOfJupiter = Eigen::VectorXd::Zero( 6 );
        dimensionlessStateOfJupiter( xCartesianPositionIndex ) = 1.0 - massParameter;

        // Set dimensionless state of test particle in Sun-Jupiter system.
        using orbital_element_conversions::yCartesianVelocityIndex;
        Eigen::VectorXd dimensionlessStateOfTestParticle = Eigen::VectorXd::Zero( 6 );
        dimensionlessStateOfTestParticle( xCartesianPositionIndex ) = 1.0;
        dimensionlessStateOfTestParticle( yCartesianVelocityIndex ) = 1.0;

        // Convert dimensionless to dimensional Cartesian state vector for test particle [m, m/s].
        Eigen::VectorXd computedDimensionalStateOfTestParticle( 6 );
        computedDimensionalStateOfTestParticle =
                crtbp::convertDimensionlessCartesianStateToDimensionalUnits(
                    dimensionlessStateOfTestParticle, gravitationalParameterSun,
                    gravitationalParameterJupiter, distanceSunJupiter );

        // Check if computed position and velocity matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( distanceSunJupiter,
                                    computedDimensionalStateOfTestParticle( xCartesianPositionIndex ),
                                    std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( velocitySunJupiter,
                                    computedDimensionalStateOfTestParticle( yCartesianVelocityIndex ),
                                    1.0e-2 );
    }
}

//! Test if time conversion for CRTBP is computed correctly.
BOOST_AUTO_TEST_CASE( testTimeConversionCircularRestrictedThreeBodyProblem )
{
    namespace crtbp = gravitation::circular_restricted_three_body_problem;
    using mathematical_constants::PI;

    // Test 1: conversion of dimensionless to dimensional time for test particle in Sun-Jupiter
    // system [s] (Koon, 2006).
    {
        // Set distance between Sun and Jupiter [m] (NASA, 2010).
        double distanceSunJupiter = 7.784e11;

        // Set gravitational parameter of Jupiter [m^3 s^-2].
        double gravitationalParameterJupiter = 6.67259e-11 * 1898.13e24;

        // Set gravitational parameter of Sun [m^3 s^-2].
        double gravitationalParameterSun = 1.32712440018e20;

        // Compute dimensional time for one complete orbit.
        double computedDimensionalTime =
                crtbp::convertDimensionlessTimeToDimensionalTime(
                    2.0 * PI, gravitationalParameterSun, gravitationalParameterJupiter,
                    distanceSunJupiter );

        // Set expected orbital period of system [s].
        double expectedOrbitalPeriodOfSunJupiterSystem = 3.733e08;

        // Check if computed orbital period matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedOrbitalPeriodOfSunJupiterSystem,
                                    computedDimensionalTime, 1.0e-2 );
    }
}

} // namespace unit_tests
} // namespace tudat
