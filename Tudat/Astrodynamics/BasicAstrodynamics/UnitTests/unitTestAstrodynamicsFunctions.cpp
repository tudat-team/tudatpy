/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

//! Test if the orbital distance of a Kepler orbit is computed correctly.
BOOST_AUTO_TEST_CASE( testKeplerRadialDistance )
{
    // Declare and set Keplerian elements.
    Eigen::Vector6d keplerianElements = Eigen::Vector6d::Constant( TUDAT_NAN );
    keplerianElements[ 0 ] = 25999.683025291e3;
    keplerianElements[ 1 ] = 0.864564003552322;
    keplerianElements[ 5 ] = 0.757654217738482;

    // Compute radial distance of the satellite.
    double radialDistance1 = basic_astrodynamics::computeKeplerRadialDistance(
                keplerianElements[ 0 ], keplerianElements[ 1 ], keplerianElements[ 5 ] );
    double radialDistance2 = basic_astrodynamics::computeKeplerRadialDistance( keplerianElements );

    // Declare and set expected radial distance [m].
    double expectedRadialDistance = 4032815.56442827;

    // Check if computed distance matches expected distance.
    BOOST_CHECK_CLOSE_FRACTION( radialDistance1, expectedRadialDistance, 1.0e-5 );
    BOOST_CHECK_CLOSE_FRACTION( radialDistance2, expectedRadialDistance, 1.0e-5 );
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

//! Test if the orbital velocity of a Kepler orbit is computed correctly.
BOOST_AUTO_TEST_CASE( testKeplerOrbitalVelocity )
{
    // Declare and set Keplerian elements.
    Eigen::Vector6d keplerianElements = Eigen::Vector6d::Constant( TUDAT_NAN );
    keplerianElements[ 0 ] = 25999.683025291e3;
    keplerianElements[ 1 ] = 0.864564003552322;
    keplerianElements[ 5 ] = 0.757654217738482;

    // Declare and set gravitational parameter of Earth [m^3 s^-2].
    double earthGravitationalParameter = physical_constants::GRAVITATIONAL_CONSTANT * 5.9736e24;

    // Compute radial distance of the satellite.
    double orbitalVelocity1 = basic_astrodynamics::computeKeplerOrbitalVelocity(
                keplerianElements[ 0 ], keplerianElements[ 1 ], keplerianElements[ 5 ], earthGravitationalParameter );
    double orbitalVelocity2 = basic_astrodynamics::computeKeplerOrbitalVelocity(
                keplerianElements, earthGravitationalParameter );

    // Declare and set expected orbital velocity [m/s].
    double expectedOrbitalVelocity = 13503.4992923871;

    // Check if computed distance matches expected distance.
    BOOST_CHECK_CLOSE_FRACTION( orbitalVelocity1, expectedOrbitalVelocity, 1.0e-5 );
    BOOST_CHECK_CLOSE_FRACTION( orbitalVelocity2, expectedOrbitalVelocity, 1.0e-5 );
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

//! Test if the periapsis altitude is computed correctly.
BOOST_AUTO_TEST_CASE( testPeriapsisAltitude )
{
    // Keplerian state
    Eigen::Vector6d keplerianState;
    keplerianState << 10000.0, 0.4, 0.0, 0.0, 0.0, 3.0;

    // Cartesian state, equivalent to keplerianState
    Eigen::Vector6d cartesianState;
    cartesianState << -1.376803915331821e4, 1.962586386216818e3, 0.0, -3.074098913804636e4, -1.285214845072070e5, 0.0;

    // Body radius
    const double centralBodyRadius = 1000.0;

    // Body gravitational parameter
    const double centralBodyGravitationalParameter = 3.9860044189e14;

    // Declare and set expected periapsis altitude.
    const double expectedPeriapsisAltitude = 5000.0;

    // Compute periapsis altitude from Keplerian state.
    const double periapsisAltitudeFromKeplerian = basic_astrodynamics::computePeriapsisAltitudeFromKeplerianState(
                keplerianState, centralBodyRadius );

    // Compute periapsis altitude from Cartesian state.
    const double periapsisAltitudeFromCartesian = basic_astrodynamics::computePeriapsisAltitudeFromCartesianState(
                cartesianState, centralBodyGravitationalParameter, centralBodyRadius );

    // Check if computed periapsis altitude from Keplerian is right.
    BOOST_CHECK_CLOSE_FRACTION( periapsisAltitudeFromKeplerian, expectedPeriapsisAltitude,
                                std::numeric_limits< double >::epsilon( ) );

    // Check if computed periapsis altitude from Cartesian is right.
    BOOST_CHECK_CLOSE_FRACTION( periapsisAltitudeFromCartesian, expectedPeriapsisAltitude,
                                std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
