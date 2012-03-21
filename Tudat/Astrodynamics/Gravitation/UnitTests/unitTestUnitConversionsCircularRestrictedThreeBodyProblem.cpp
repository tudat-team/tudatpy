/*    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
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
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

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
        namespace crtbp = astrodynamics::gravitation::circular_restricted_three_body_problem;

        // Set distance between Sun and Jupiter [m] [NASA, 2010].
        double distanceSunJupiter = 7.784e11;

        // Set magnitude of velocity of Sun-Jupiter system [m/s].
        double velocitySunJupiter = 13.102e3;

        // Set gravitational parameter of Jupiter [kg m^3 s^-2].
        double gravitationalParameterJupiter = 6.67259e-11 * 1898.13e24;

        // Set gravitational parameter of Sun [kg m^3 s^-2].
        double gravitationalParameterSun = 1.32712440018e20;

        // Set mass parameter of Sun-Jupiter system.
        double massParameter = gravitationalParameterJupiter
                / ( gravitationalParameterJupiter + gravitationalParameterSun );

        // Set dimensionless state of the Sun.
        using tudat::orbital_element_conversions::xPositionIndex;
        Eigen::VectorXd dimensionlessStateOfSun = Eigen::VectorXd::Zero( 6 );
        dimensionlessStateOfSun( xPositionIndex ) = -massParameter;

        // Declare and set dimensionless state of Jupiter.
        Eigen::VectorXd dimensionlessStateOfJupiter = Eigen::VectorXd::Zero( 6 );
        dimensionlessStateOfJupiter( xPositionIndex ) = 1.0 - massParameter;

        // Set dimensionless state of test particle in Sun-Jupiter system.
        using tudat::orbital_element_conversions::yVelocityIndex;
        Eigen::VectorXd dimensionlessStateOfTestParticle = Eigen::VectorXd::Zero( 6 );
        dimensionlessStateOfTestParticle( xPositionIndex ) = 1.0;
        dimensionlessStateOfTestParticle( yVelocityIndex ) = 1.0;

        // Convert dimensionless to dimensional Cartesian state vector for test particle [m, m/s].
        Eigen::VectorXd computedDimensionalStateOfTestParticle( 6 );
        computedDimensionalStateOfTestParticle =
                crtbp::convertDimensionlessCartesianStateToDimensionalUnits(
                    dimensionlessStateOfTestParticle, gravitationalParameterSun,
                    gravitationalParameterJupiter, distanceSunJupiter );

        // Check if computed position and velocity matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( distanceSunJupiter,
                                    computedDimensionalStateOfTestParticle( xPositionIndex ),
                                    std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( velocitySunJupiter,
                                    computedDimensionalStateOfTestParticle( yVelocityIndex ),
                                    1.0e-2 );
    }
}

//! Test if time conversion for CRTBP is computed correctly.
BOOST_AUTO_TEST_CASE( testTimeConversionCircularRestrictedThreeBodyProblem )
{
    namespace crtbp = astrodynamics::gravitation::circular_restricted_three_body_problem;
    using mathematics::PI;

    // Test 1: conversion of dimensionless to dimensional time for test particle in Sun-Jupiter
    // system [s] (Koon, 2006).
    {
        // Set distance between Sun and Jupiter [m] (NASA, 2010).
        double distanceSunJupiter = 7.784e11;

        // Set gravitational parameter of Jupiter [kg m^3 s^-2].
        double gravitationalParameterJupiter = 6.67259e-11 * 1898.13e24;

        // Set gravitational parameter of Sun [kg m^3 s^-2].
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
