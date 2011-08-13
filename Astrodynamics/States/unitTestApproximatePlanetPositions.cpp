/*! \file unitTestApproximatePlanetPositions.cpp
 *    Source file for a unit test that tests the implementation of the
 *    ApproximatePlanetPositions class in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : L. van der Ham
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.vanderHam@student.tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 5 April, 2011
 *    Last modified     : 10 August, 2011
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
 *      110701    K. Kumar          Updated to use new predefined planets
 *                                  architecture; removed dynamic memory
 *                                  allocation.
 *      110714    L. van der Ham    Added circular coplanar case.
 */

// Include statements.
#include "unitTestApproximatePlanetPositions.h"

// Using declarations.
using mathematics::computeAbsoluteValue;
using mathematics::convertCartesianToSpherical;
using mathematics::MACHINE_PRECISION_DOUBLES;
using unit_conversions::convertDegreesToRadians;
using unit_conversions::convertRadiansToDegrees;
using unit_conversions::convertDegreesToArcminutes;
using unit_conversions::convertArcminutesToArcseconds;
using unit_conversions::convertAstronomicalUnitsToMeters;
using orbital_element_conversions::convertCartesianToKeplerianElements;

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
    KeplerianElements marsOrbitalElementsForTest1;
    marsOrbitalElementsForTest1.setSemiMajorAxis( 2.279361944126564e11 );
    marsOrbitalElementsForTest1.setEccentricity( 9.338126166083623e-02 );
    marsOrbitalElementsForTest1.setInclination(
                convertDegreesToRadians( 1.848907897011101 ) );
    marsOrbitalElementsForTest1.setArgumentOfPeriapsis(
                convertDegreesToRadians( 2.866464026954701e2 ) );
    marsOrbitalElementsForTest1.setLongitudeOfAscendingNode(
            convertDegreesToRadians( 4.952419052428279e1 ) );
    marsOrbitalElementsForTest1.setTrueAnomaly(
            convertDegreesToRadians( 3.577219707986779e2 ) );

    // Create predefined Mars.
    Planet predefinedMars;
    predefinedMars.setPredefinedPlanetSettings( Planet::mars );

    // Create predefined Sun.
    Planet predefinedSun;
    predefinedSun.setPredefinedPlanetSettings( Planet::sun );

    // Convert expected result to Cartesian elements.
    CartesianElements expectedMarsEphemeris;
    expectedMarsEphemeris
            = orbital_element_conversions::
              convertKeplerianToCartesianElements( &marsOrbitalElementsForTest1,
                                                   &predefinedSun );

    // Retrieve state of Mars in Cartesian elements at Julian date 2455626.5.
    CartesianElements marsEphemeris;
    marsEphemeris = *predefinedMars.getStateFromEphemeris( 2455626.5 );

    // Compute absolute differences in position in spherical coordinates.
    VectorXd positionInSphericalCoordinates( 3 );
    convertCartesianToSpherical( marsEphemeris.state.segment( 0, 3 ),
                                 positionInSphericalCoordinates );

    VectorXd expectedPositionInSphericalCoordinates( 3 );
    convertCartesianToSpherical( expectedMarsEphemeris.state.segment( 0, 3 ),
                                 expectedPositionInSphericalCoordinates );

    VectorXd differenceInSphericalCoordinates( 3 );
    differenceInSphericalCoordinates
            = positionInSphericalCoordinates
                   - expectedPositionInSphericalCoordinates;

    // Set test result to true if the test does not match the expected result.
    // Check if computed relative error in radius in metres matches expected value.
    // The error bounds are taken from: http://ssd.jpl.nasa.gov/?planet_pos.
    double errorBoundRadius_ = 3.0e8;
    if ( abs( differenceInSphericalCoordinates( 0 ) ) > errorBoundRadius_ )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The computed error in radius of ephemeris of Mars " << endl;
        cerr << "( " << abs( differenceInSphericalCoordinates( 0 ) )
             << " metres )" << endl;
        cerr << "using the ApproximatePlanetPositions class exceeds "
             << "the error bound "
             << endl;
        cerr << "( " << errorBoundRadius_ << " metres )." << endl;
        cerr << "The computed error exceeds the error bound by: "
             << abs( differenceInSphericalCoordinates( 0 ) ) - errorBoundRadius_
             << " metres." << endl;
    }

    // Check if computed error in right ascension in radians matches expected
    // value.
    double rightAscensionErrorInArcsec_ = convertArcminutesToArcseconds(
                convertDegreesToArcminutes( convertRadiansToDegrees(
                    abs( differenceInSphericalCoordinates( 1 ) ) ) ) );
    double errorBoundRightAscensionInArcsec_ = 100.0;

    if ( rightAscensionErrorInArcsec_ > errorBoundRightAscensionInArcsec_ )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The computed error in right ascension of ephemeris of Mars "
             << endl;
        cerr << "( " << rightAscensionErrorInArcsec_ << " arcsec )" << endl;
        cerr << "using the ApproximatePlanetPositions class exceeds "
                << "the expected error " << endl;
        cerr << "( " << errorBoundRightAscensionInArcsec_ << " arcsec )." << endl;
        cerr << "The computed error exceeds the expected error by: "
                << rightAscensionErrorInArcsec_ - errorBoundRightAscensionInArcsec_
                << " arcsec." << endl;
    }

    // Check if computed error in declination inclination in radians matches
    // expected value.
    double declinationErrorInArcsec_ = convertArcminutesToArcseconds(
                convertDegreesToArcminutes( convertRadiansToDegrees(
                    abs( differenceInSphericalCoordinates( 2 ) ) ) ) );
    double errorBoundDeclinationInArcsec_ = 40.0;

    if ( errorBoundDeclinationInArcsec_ > 40.0 )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The computed error in declination of ephemeris of Mars "
                << endl;
        cerr << "( " << declinationErrorInArcsec_ << " arcsec )" << endl;
        cerr << "using the ApproximatePlanetPositions class exceeds "
                << "the expected error " << endl;
        cerr << "( " << errorBoundDeclinationInArcsec_ << " arcsec )." << endl;
        cerr << "The computed error exceeds the expected error by: "
                << declinationErrorInArcsec_ - errorBoundDeclinationInArcsec_
                << " arcsec." << endl;
    }

    // Test 2: Get circular coplanar position of Mars at Julian date 2455626.5.

    // Define tolerance.
    double errorTolerance_ = MACHINE_PRECISION_DOUBLES;

    // Set Julian Date of test data.
    double julianDate = 2455626.5;

    // Create predefined Mars.
    CelestialBody predefinedMarsCelestialBody;

    // Get Mars state in Cartesian elements for circular coplanar case at given Julian date.
    ApproximatePlanetPositionsCircularCoplanar approximatePlanetPositionsCircularCoplanar_;
    predefinedMarsCelestialBody.setEphemeris( &approximatePlanetPositionsCircularCoplanar_ );
    CartesianElements marsEphemerisCircularCoplanar;
    marsEphemerisCircularCoplanar = *predefinedMarsCelestialBody.
                                    getStateFromEphemeris( julianDate );

    // Get Cartesian position from state.
    Vector3d positionMars;
    positionMars = marsEphemerisCircularCoplanar.getPosition( );

    // Get Cartesian velocity from state.
    Vector3d velocityMars;
    velocityMars = marsEphemerisCircularCoplanar.getVelocity( );

    // Compute Keplerian elements of Test 2.
    KeplerianElements keplerianElementsTest;
    keplerianElementsTest = convertCartesianToKeplerianElements(
                &marsEphemerisCircularCoplanar, &predefinedSun );

    // Convert Cartesian to Keplerian elements of Test 1.
    // Use these elements as reference values in Test 2.
    KeplerianElements keplerianElementsTest3D;
    keplerianElementsTest3D = convertCartesianToKeplerianElements(
                 &marsEphemeris, &predefinedSun );

    // Compute the difference in semi-major axis between Test 2 and
    // the external EphemerisData "p_elem_t2.txt".
    double errorSemiMajorAxis = abs( positionMars.norm( )
                                     - convertAstronomicalUnitsToMeters( 1.52371243 ) )
                                / convertAstronomicalUnitsToMeters( 1.52371243 ) ;

    if ( errorSemiMajorAxis > errorTolerance_ )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The computed relative error in position of the  "<< endl;
        cerr << "coplanar circular position of Mars ( "
             << errorSemiMajorAxis << " )"<< endl;
        cerr << "using the ApproximatePlanetPositionsCircularCoplanar class, exceeds "
             << "the maximum expected error " << endl;
        cerr << "( " << errorTolerance_ << " )." << endl;
    }

    // Check orientation of position vector by comparison of separate components.
    // Error in position should be smaller than maximum expected offset with respect to
    // elliptical and inclined orbits.
    double maximumErrorPosition =   keplerianElementsTest3D.getSemiMajorAxis( ) * (
                                                keplerianElementsTest3D.getEccentricity( ) + 1.0
                                                - cos( keplerianElementsTest3D.getInclination( ) ) );
    Vector3d errorPositionVector = positionMars - marsEphemeris.getPosition( );

    if ( abs( errorPositionVector( 0 ) ) > maximumErrorPosition
         || abs( errorPositionVector( 1 ) ) > maximumErrorPosition )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The computed error in position vector of the  "<< endl;
        cerr << "coplanar circular position of Mars ( "
             << errorPositionVector << " meters )" << endl;
        cerr << "using the ApproximatePlanetPositionsCircularCoplanar class, exceeds "
             << "the expected error (" << endl;
        cerr << "( " << maximumErrorPosition << " meters )." << endl;
    }

    // Check size of velocity.
    Vector3d errorVelocity = velocityMars - marsEphemeris.getVelocity( );

    // Error in scalar velocity should be smaller than maximum expected offset with respect to
    // ellipitical and inclined orbits.
    double expectedErrorVelocity =  abs( sqrt( predefinedSun.getGravitationalParameter( )
                                               / marsEphemeris.getPosition( ).norm( ) ) *
            ( ( 1.0 - cos( keplerianElementsTest3D.getInclination( ) )
                + sqrt( ( 1.0 - keplerianElementsTest3D.getEccentricity( ) ) /
                        ( 1.0 + keplerianElementsTest3D.getEccentricity( ) ) ) - 1.0 ) ) );

    if ( errorVelocity.norm( ) > expectedErrorVelocity )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The computed error in velocity of the "<< endl;
        cerr << "coplanar circular position of Mars "
             << "( " << errorVelocity.norm( ) << " meters per second )" << endl;
        cerr << "using the ApproximatePlanetPositionsCircularCoplanar class, exceeds "
             << "the expected error " << endl;
        cerr << "( " << expectedErrorVelocity << " meters per second )." << endl;
        cerr << "The computed error exceeds the expected error by: "
             << abs( errorVelocity.norm( ) - expectedErrorVelocity )
             << " meters per second." << endl;
    }

    // Check if the computed eccentricity equals zero.
    double errorEccentricity = keplerianElementsTest.getEccentricity( );

    if ( errorEccentricity > errorTolerance_ )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The computed error in eccentricity of the  "<< endl;
        cerr << " coplanar circular position of Mars is " <<endl;
        cerr << errorEccentricity << ", which exceeds the tolerated error "
             << errorTolerance_ << "." << endl;
    }

    // Check if the computed inclination equals zero.
    double relativeErrorInclination_ = abs( keplerianElementsTest.getInclination( ) - 0.0 )
                                       / ( 2.0 * M_PI );
    if ( relativeErrorInclination_ > errorTolerance_ )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The relative error in the inclination of the  "<< endl;
        cerr << "coplanar circular position of Mars is " <<endl;
        cerr << unit_conversions::convertRadiansToDegrees( relativeErrorInclination_ )
             << " degrees." << endl;
    }

    // Check if the computed position component in Z direction equals zero.
    double relativeErrorPositionComponentZ_ = positionMars( 2 ) / positionMars.norm( );

    if ( relativeErrorPositionComponentZ_ > errorTolerance_ )
    {
        isApproximatePlanetPositionsErroneous_ = true;

        // Generate error statements.
        cerr << "The relative error in vertical position of the  "<< endl;
        cerr << "coplanar circular position of Mars is " <<endl;
        cerr << relativeErrorPositionComponentZ_ << " meters." << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isApproximatePlanetPositionsErroneous_;
}

}

// End of file.
