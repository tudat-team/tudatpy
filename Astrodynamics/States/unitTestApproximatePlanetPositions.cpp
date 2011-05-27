/*! \file unitTestApproximatePlanetPositions.cpp
 *    Source file for a unit test that tests the implementation of the
 *    ApproximatePlanetPositions class in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 5 April, 2011
 *    Last modified     : 11 April, 2011
 *
 *    References
 *      HORIZONS Web-Interface, http://ssd.jpl.nasa.gov/horizons.cgi,
 *          last accessed: 5 April, 2011.
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      110405    K. Kumar          File created.
 *      110406    K. Kumar          Added cerr statements.
 *      110411    K. Kumar          Changed test to check errors in position,
 *                                  right ascension, and declination.
 */

// Include statements.
#include "unitTestApproximatePlanetPositions.h"

// Using declarations.
using mathematics::computeAbsoluteValue;
using unit_conversions::convertRadiansToDegrees;
using unit_conversions::convertDegreesToArcminutes;
using unit_conversions::convertArcminutesToArcseconds;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test ApproximatePlanetPositions class.
bool testApproximatePlanetPositions( )
{
    // Test result initialised to true.
    bool isApproximatePlanetPositionsErroneous_ = false;

    // Test of ApproximatePlanetPositions class.
    // Test 1: Get orbital elements of Mars at Julian date 2455626.5.

    // Expected result.
    KeplerianElements* marsOrbitalElementsForTest1_ = new KeplerianElements;
    marsOrbitalElementsForTest1_->setSemiMajorAxis( 2.279361944126564e11 );
    marsOrbitalElementsForTest1_->setEccentricity( 9.338126166083623e-02 );
    marsOrbitalElementsForTest1_->setInclination(
            unit_conversions::convertDegreesToRadians( 1.848907897011101 ) );
    marsOrbitalElementsForTest1_->setArgumentOfPeriapsis(
            unit_conversions::convertDegreesToRadians(
                    2.866464026954701e2 ) );
    marsOrbitalElementsForTest1_->setLongitudeOfAscendingNode(
            unit_conversions::convertDegreesToRadians( 4.952419052428279e1 ) );
    marsOrbitalElementsForTest1_->setTrueAnomaly(
            unit_conversions::convertDegreesToRadians( 3.577219707986779e2 ) );

    // Create predefined Mars.
    CelestialBody* pointerToPredefinedMars_
            = predefined_celestial_bodies::createPredefinedCelestialBody(
                    predefined_celestial_bodies::mars );

    // Create predefined Sun.
    CelestialBody* pointerToPredefinedSun_
            = predefined_celestial_bodies::createPredefinedCelestialBody(
                    predefined_celestial_bodies::sun );

    // Convert expected result to Cartesian elements.
    CartesianElements* pointerToExpectedMarsEphemeris_
            = new CartesianElements;
    pointerToExpectedMarsEphemeris_
            = orbital_element_conversions::
              convertKeplerianToCartesianElements( marsOrbitalElementsForTest1_,
                                                   pointerToPredefinedSun_ );

    // Retrieve state of Mars in Cartesian elements at Julian date 2455626.5.
    CartesianElements* pointerToMarsEphemeris_ = new CartesianElements;
    pointerToMarsEphemeris_ = pointerToPredefinedMars_
                              ->getStateFromEphemeris( 2455626.5 );

    // Compute differences in position in spherical coordinates.
    VectorXd positionInSphericalCoordinates_( 3 );
    mathematics::convertCartesianToSpherical(
            pointerToMarsEphemeris_->state.segment( 0, 3 ),
            positionInSphericalCoordinates_ );

    VectorXd expectedPositionInSphericalCoordinates_( 3 );
    mathematics::convertCartesianToSpherical(
            pointerToExpectedMarsEphemeris_->state.segment( 0, 3 ),
            expectedPositionInSphericalCoordinates_ );

    VectorXd differencesInSphericalCoordinates_( 3 );
    differencesInSphericalCoordinates_
            = positionInSphericalCoordinates_
              - expectedPositionInSphericalCoordinates_;

    // Set test result to true if the test does not match the expected result.
    // Check if computed error in radius in metres matches expected value.
    // The expected errors are taken from: http://ssd.jpl.nasa.gov/?planet_pos.
    if ( computeAbsoluteValue( differencesInSphericalCoordinates_( 0 ) ) > 3e8 )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The computed error in radius of ephemeris of Mars " << endl;
        cerr << "( " << computeAbsoluteValue(
                differencesInSphericalCoordinates_( 0 ) )  << " metres )"
             << endl;
        cerr << "using the ApproximatePlanetPositions class exceeds "
             << "the expected error "
             << endl;
        cerr << "( " << 3e8 << " metres )." << endl;
        cerr << "The computed error exceeds the expected error by: "
             << computeAbsoluteValue(
                     differencesInSphericalCoordinates_( 0 ) ) - 3e8
             << " metres." << endl;
    }

    // Check if computed error in right ascension in radians matches expected
    // value.
    double rightAscensionErrorInArcsec_
            = convertArcminutesToArcseconds( convertDegreesToArcminutes(
                    convertRadiansToDegrees( computeAbsoluteValue(
                            differencesInSphericalCoordinates_( 1 ) ) ) ) );

    if ( rightAscensionErrorInArcsec_ > 100.0 )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The computed error in right ascension of ephemeris of Mars "
             << endl;
        cerr << "( " << rightAscensionErrorInArcsec_ << " arcsec )" << endl;
        cerr << "using the ApproximatePlanetPositions class exceeds "
                << "the expected error " << endl;
        cerr << "( " << 100.0 << " arcsec )." << endl;
        cerr << "The computed error exceeds the expected error by: "
                << rightAscensionErrorInArcsec_ - 100.0 << " arcsec." << endl;
    }

    // Check if computed error in declination inclination in radians matches
    // expected value.
    double declinationErrorInArcsec_
            = convertArcminutesToArcseconds( convertDegreesToArcminutes(
                    convertRadiansToDegrees( computeAbsoluteValue(
                            differencesInSphericalCoordinates_( 2 ) ) ) ) );

    if ( declinationErrorInArcsec_ > 40.0 )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The computed error in declination of ephemeris of Mars "
                << endl;
        cerr << "( " << declinationErrorInArcsec_ << " arcsec )" << endl;
        cerr << "using the ApproximatePlanetPositions class exceeds "
                << "the expected error " << endl;
        cerr << "( " << 40.0 << " arcsec )." << endl;
        cerr << "The computed error exceeds the expected error by: "
                << rightAscensionErrorInArcsec_ - 40.0 << " arcsec." << endl;
    }

    std::cout << rightAscensionErrorInArcsec_ << std::endl;
    std::cout << declinationErrorInArcsec_ << std::endl;

    // Deallocate dynamically allocated variables.
    delete marsOrbitalElementsForTest1_;
    delete pointerToMarsEphemeris_;

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isApproximatePlanetPositionsErroneous_;
}

}

// End of file.
