/*! \file unitTestAstrodynamicsFunctions.cpp
 *    This unit test will test general astrodynamics functions defined in astrodynamicsFunctions.h.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 11 November, 2011
 *    Last modified     : 6 December, 2011
 *
 *    References
 *      Wikipedia. Geostationary orbit, http://en.wikipedia.org/wiki/Geostationary_orbit, last
 *          accessed: 22nd November, 2011.
 *     Keefe, T.J. Synodic Period Calculator, http://www.ccri.edu/physics/keefe/synodic_calc.htm,
 *          last accessed: 6th December, 2011, last modified: 18th November, 2011.
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *      The tests need to be updated to check benchmark values from literature.
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      100910    K. Kumar          First creation of code.
 *      111115    K. Kumar          Altered tests to override default mass values; corrected
 *                                  bool-variable name; corrected cerr-statement.
 *      111122    K. Kumar          Altered tests to use reference data.
 *      111206    K. Kumar          Updated synodic period test and comments.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include <limits>
#include "Astrodynamics/astrodynamicsFunctions.h"
#include "Astrodynamics/physicalConstants.h"

//! Test implementation of general astrodynamics functions.
int main( )
{
    // Summary of tests.
    // Test 1: Compute orbital period for geostationary satellite test case.
    // Test 2: Compute two-body angular momentum for geostationary satellite test case.
    // Test 3: Compute two-body mean motion for geostationary satellite test case.
    // Test 4: Compute two-body orbital energy for geostationary satellite test case.
    // Test 5: Compute synodic period of two bodies in different orbits.

    // Declare and initialize test result to false.
    bool isAstrodynamicsFunctionsErroneous = false;

    // Test 1: Compute orbital period for geostationary satellite test case.
    // Reference: http://en.wikipedia.org/wiki/Geostationary_orbit.

    // Declare and set satellite mass [kg].
    double satelliteMass = 1.0e3;

    // Declare and set gravitational parameter of Earth [m^3 s^-2].
    double earthGravitationalParameter
            = tudat::PhysicalConstants::GRAVITATIONAL_CONSTANT * 5.9736e24;

    // Declare and set distance between Earth center and satellite.
    double distanceBetweenSatelliteAndEarth = 4.2164e7;

    // Compute orbital period of satellite.
    double orbitalPeriod = tudat::astrodynamics::computeTwoBodyOrbitalPeriod(
                distanceBetweenSatelliteAndEarth, earthGravitationalParameter, satelliteMass );

    // Declare and set expected orbital period [s].
    double expectedOrbitalPeriod = 86164.09054;

    // Check if computed orbital period matches expected orbital period.
    if ( std::fabs( ( orbitalPeriod - expectedOrbitalPeriod ) / expectedOrbitalPeriod ) > 1.0e-5 )
    {
        std::cerr << "The two-body orbital period function is not functioning well, as the "
                  << "computed orbital period (" << orbitalPeriod << ") varies from the expected "
                  << "orbital period (" << expectedOrbitalPeriod << ") by: "
                  << orbitalPeriod - expectedOrbitalPeriod << std::endl;
        isAstrodynamicsFunctionsErroneous = true;
    }

    // Test 2: Compute two-body angular momentum for geostationary satellite test case.

    // Declare and set eccentricity of satellite orbit.
    double eccentricityOfSatelliteOrbit = 0.0;

    double angularMomentum = tudat::astrodynamics::computeTwoBodyAngularMomentum(
                distanceBetweenSatelliteAndEarth, eccentricityOfSatelliteOrbit,
                earthGravitationalParameter, satelliteMass );

    // Declare and set expected angular momentum.
    // The expected angular momentum is computed using the fact that for a circular orbit,
    // H = mRV. This is an independent check of the code, which computes angular
    // momentum differently.
    double expectedAngularMomentum = satelliteMass * distanceBetweenSatelliteAndEarth
            * std::sqrt( earthGravitationalParameter / distanceBetweenSatelliteAndEarth );

    // Check if computed angular momentum matches expected angular momentum.
    if ( std::fabs( ( angularMomentum - expectedAngularMomentum ) / expectedAngularMomentum )
         > std::numeric_limits< double >::epsilon( ) )
    {
        std::cerr << "The two-body angular momentum function is not functioning well, as the "
                  << "computed angular momentum (" << angularMomentum << ") varies from the "
                  << "expected angular momentum (" << expectedAngularMomentum << ") by: "
                  << angularMomentum - expectedAngularMomentum << std::endl;
        isAstrodynamicsFunctionsErroneous = true;
    }

    // Test 3: Compute two-body mean motion for geostationary satellite test case.
    // Reference: http://en.wikipedia.org/wiki/Geostationary_orbit.
    double meanMotion = tudat::astrodynamics::computeTwoBodyMeanMotion(
                distanceBetweenSatelliteAndEarth, earthGravitationalParameter, satelliteMass );

    // Declare and set expected mean motion [rad/s].
    double expectedMeanMotion = 7.2921e-5;

    // Check if computed mean motion matches expected mean motion.
    if ( std::fabs( ( meanMotion - expectedMeanMotion ) / expectedMeanMotion ) > 1.0e-7 )
    {
        std::cerr << "The two-body mean motion function is not functioning well, as the "
                  << "computed mean motion (" << meanMotion << ") varies from the "
                  << "expected mean motion (" << expectedMeanMotion << ") by: "
                  << meanMotion - expectedMeanMotion << std::endl;
        isAstrodynamicsFunctionsErroneous = true;
    }

    // Test 4: Compute two-body orbital energy for geostationary satellite test case.
    double orbitalEnergy = tudat::astrodynamics::computeTwoBodyOrbitalEnergy(
                distanceBetweenSatelliteAndEarth, earthGravitationalParameter, satelliteMass );

    // Declare and set expected orbital energy.
    // The expected orbital energy is computed using the fact that for a circular orbit,
    // E = m ( V^2/2 - mu/R ). This is an independent check of the code, which computes orbital
    // energy differently.
    double expectedOrbitalEnergy = satelliteMass * (
                0.5 * earthGravitationalParameter / distanceBetweenSatelliteAndEarth
                -  earthGravitationalParameter / distanceBetweenSatelliteAndEarth );

    // Check if computed orbital energy matches expected orbital energy.
    if ( std::fabs( ( orbitalEnergy - expectedOrbitalEnergy ) / expectedOrbitalEnergy )
         > std::numeric_limits< double >::epsilon( ) )
    {
        std::cerr << "The two-body orbital energy function is not functioning well, as the "
                  << "computed orbital energy (" << orbitalEnergy << ") varies from the "
                  << "expected orbital energy (" << expectedOrbitalEnergy << ") by: "
                  << orbitalEnergy - expectedOrbitalEnergy << std::endl;
        isAstrodynamicsFunctionsErroneous = true;
    }

    // Test 5: Compute synodic period of Mars with respect to Earth.
    double synodicPeriod = tudat::astrodynamics::computeSynodicPeriod( 365.256378, 686.95 );

    // Declare and set expected synodic period.
    double expectedSynodicPeriod = 779.9746457736733;

    // Check if computed synodic period matches expected synodic period.
    if ( std::fabs( ( synodicPeriod - expectedSynodicPeriod ) / expectedSynodicPeriod )
         > std::numeric_limits< double >::epsilon( ) )
    {
        std::cerr << "The synodic period function is not functioning well, as the "
                  << "computed synodic period (" << synodicPeriod << ") varies from the "
                  << "expected synodic period (" << expectedSynodicPeriod << ") by: "
                  << synodicPeriod - expectedSynodicPeriod << std::endl;
        isAstrodynamicsFunctionsErroneous = true;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isAstrodynamicsFunctionsErroneous )
    {
        std::cerr << "testAstrodynamicsFunctions failed!" << std::endl;
    }

    return isAstrodynamicsFunctionsErroneous;
}

// End of file.
