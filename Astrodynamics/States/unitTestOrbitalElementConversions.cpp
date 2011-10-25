/*! \file unitTestOrbitalElementConversions.cpp
 *    Source file of unit test for the orbitalElementConversion, from Cartesian to Keplerian and
 *    viceversa. The first part of the code tests the code for elliptical, parabolic, hyperbolic
 *    and circular orbits. SI units are used. The second part of the code tests the code from
 *    Cartesian to Keplerian with the example 3.4 pag. 63 of the book "Fondamenti di Meccanica del
 *    Volo Spaziale" (G. Mengali, A.A. Quarta). In this part of the code, canonical units are used.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 13
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 3 December, 2010
 *    Last modified     : 10 May, 2011
 *
 *    References
 *      http://www.astro.uu.nl/~strous/AA/en/reken/kepler.html, last accessed: 16th February, 2011.
 *      Vallado, D. A., McClain, W. D. Fundamentals of astrodynamics and applications, 2nd Edition,
 *          Kluwer Academic Publishers, The Netherlands, 2004.
 *      Fortescue, P. W., et al. Spacecraft systems engineering, Third Edition,
 *          Wiley, England, 2003.
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
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
 *      101203    E. Iorfida        First creation of the code.
 *      101208    E. Iorfida        Fulfillment of the code with the elliptical case.
 *      101208    E. Iorfida        Modified punctuation.
 *      101215    E. Iorfida        Added tolerance, added parabolic, circular and hyperbolic
 *                                  cases.
 *      101217    E. Iorfida        Added computeAbsoluteValue( ) in the errors computation,
 *                                  modified punctuation.
 *      101219    J. Melman         Put gravitational parameters in one place, changed first right
 *                                  ascension to 15.0 * pi / 8.0, thereby exposing a possible
 *                                  error.
 *      110107    E. Iorfida        orbitalConversionBookExampleUnitTest.test added to this file,
 *                                  to have a unique unit test file for the conversion code. Also
 *                                  some punctuation modifications have been made.
 *      110109    J. Melman         Included test for semi-latus rectum of circular case. Gave the
 *                                  orbital angles less trivial values, and not almost exclusively
 *                                  in the first quadrant.
 *      110111    E. Iorfida        Updated to the new format of unitTest file and added hyperbolic
 *                                  equatorial case.
 *      110204    K. Kumar          Removed "vector" from naming.
 *      110216    K. Kumar          Added unit tests for new orbital element conversion functions.
 *      110310    K. Kumar          Changed right ascension of ascending node to longitude of
 *                                  ascending node.
 *      110510    K. Kumar          Updated to use new orbital element conversion functions and
 *                                  removed dynamic memory allocation.
 */

// Include statements.
#include <cmath>
#include "Astrodynamics/Bodies/celestialBody.h"
#include "Astrodynamics/Bodies/planet.h"
#include "Astrodynamics/EnvironmentModels/gravityFieldModel.h"
#include "Astrodynamics/EnvironmentModels/sphericalHarmonicsGravityField.h"
#include "Astrodynamics/States/convertMeanAnomalyToEccentricAnomaly.h"
#include "Astrodynamics/States/convertMeanAnomalyToHyperbolicEccentricAnomaly.h"
#include "Astrodynamics/States/orbitalElementConversions.h"
#include "Astrodynamics/States/unitTestOrbitalElementConversions.h"
#include "Mathematics/basicMathematicsFunctions.h"
#include "Mathematics/unitConversions.h"
#include "Mathematics/RootFindingMethods/newtonRaphson.h"

// Using declarations.
using std::cerr;
using std::endl;
using std::fabs;
using std::pow;
using mathematics::MACHINE_PRECISION_DOUBLES;
using orbital_element_conversions::convertCartesianToKeplerianElements;
using orbital_element_conversions::convertKeplerianToCartesianElements;
using orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly;
using orbital_element_conversions::ConvertMeanAnomalyToHyperbolicEccentricAnomaly;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test of orbitalElementConversion code.
bool testOrbitalElementConversions( )
{
    // Test of orbital element conversion methods imeplemented in Tudat.
    // Test 1: Test of Cartesian-to-Keplerian elements conversion and
    //         Keplerian-to-Cartesian elements conversion.
    // Test 2: Test of true anomaly to eccentric anomaly conversion.
    // Test 3: Test of eccentric anomaly to true anomaly conversion.
    // Test 4: Test of true anomaly to hyperbolic eccentric anomaly conversion.
    // Test 5: Test of hyperbolic eccentric anomaly to true anomaly conversion.
    // Test 6: Test of eccentric anomaly to mean anomaly conversion.
    // Test 7: Test of mean anomaly to eccentric anomaly conversion.
    // Test 8: Test of hyperbolic eccentric anomaly to mean anomaly conversion.
    // Test 9: Test of mean anomaly to hyperbolic eccentric anomaly conversion.
    // Test 10: Test of elapsed time to mean anomaly for elliptical orbits.
    // Test 11: Test of mean anomaly to elapsed time for elliptical orbits.
    // Test 12: Test of elapsed time to mean anomaly for hyperbolic orbits.
    // Test 13: Test of mean anomaly to elapsed time for hyperbolic orbits.

    // Test 1: Test of Cartesian-to-Keplerian elements conversion and
    //         Keplerian-to-Cartesian elements conversion.

    // Initialize unit test result to false.
    bool isOrbitalElementConversionErroneous = false;

    // Define tolerance.
    double errorTolerance_ = 1.0e2 * MACHINE_PRECISION_DOUBLES;

    // Create predefind Earth and set different gravitational parameter value.
    Planet predefinedEarth;
    predefinedEarth.setPredefinedPlanetSettings( Planet::earth );
    GravityFieldModel* pointerToEarthGravityField = predefinedEarth.getGravityFieldModel( );
    pointerToEarthGravityField->setGravitationalParameter( 398600.4418e9 );

    // Create predefined Mars.
    Planet predefinedMars;
    predefinedMars.setPredefinedPlanetSettings( Planet::mars );

    // Create custom-defined Sun with central gravity field.
    CelestialBody customDefinedSun;
    SphericalHarmonicsGravityField sunCentralGravity;
    sunCentralGravity.setGravitationalParameter( 1.32712440018e20 );
    sunCentralGravity.setDegreeOfExpansion( 0 );
    sunCentralGravity.setOrderOfExpansion( 0 );
    customDefinedSun.setGravityFieldModel( &sunCentralGravity );

    // Create custom-defined central body.
    CelestialBody customDefinedBody;
    SphericalHarmonicsGravityField customBodyCentralGravity;
    customBodyCentralGravity.setGravitationalParameter( 1.0 );
    customBodyCentralGravity.setDegreeOfExpansion( 0 );
    customBodyCentralGravity.setOrderOfExpansion( 0 );
    customDefinedBody.setGravityFieldModel( &customBodyCentralGravity );

    // *************************************************************************
    // Elliptical orbit case around the Earth.
    // *************************************************************************

    // From Keplerian to Cartesian.
    KeplerianElements keplerianEllipticalElements1;

    // Define Keplerian elements.
    keplerianEllipticalElements1.setSemiMajorAxis(
                unit_conversions::convertAstronomicalUnitsToMeters( 0.3 ) );
    keplerianEllipticalElements1.setEccentricity( 0.2 );
    keplerianEllipticalElements1.setInclination( M_PI / 4.0 );
    keplerianEllipticalElements1.setArgumentOfPeriapsis( 4.0 * M_PI / 3.0 );
    keplerianEllipticalElements1.setLongitudeOfAscendingNode( M_PI
                                                              / 8.0 );
    keplerianEllipticalElements1.setTrueAnomaly( M_PI / 3.0 );
    keplerianEllipticalElements1.setSemiLatusRectum(
                keplerianEllipticalElements1.getSemiMajorAxis( )
                * ( 1.0 - pow( keplerianEllipticalElements1.getEccentricity( ), 2.0 ) ) );

    // Compute Cartesian elements.
    CartesianElements cartesianEllipticalElements;

    cartesianEllipticalElements = convertKeplerianToCartesianElements(
                &keplerianEllipticalElements1, &predefinedEarth );

    // From Cartesian to Keplerian.
    // Compute Keplerian elements.
    KeplerianElements keplerianEllipticalElements2;

    keplerianEllipticalElements2 = convertCartesianToKeplerianElements(
                &cartesianEllipticalElements, &predefinedEarth );

    // Set test result to false if the output Keplerian elements of the
    // conversion from Cartesian to Keplerian are equal to the input Keplerian
    // elements of the conversion from Keplerian to Cartesian, within a
    // tolerance limit.
    if ( fabs( ( keplerianEllipticalElements2.getSemiMajorAxis( ) -
                 keplerianEllipticalElements1.getSemiMajorAxis( ) ) /
               keplerianEllipticalElements1.getSemiMajorAxis( ) ) >= errorTolerance_ ||
         fabs( keplerianEllipticalElements2.getEccentricity( ) -
               keplerianEllipticalElements1.getEccentricity( ) ) >= errorTolerance_ ||
         fabs( keplerianEllipticalElements2.getInclination( ) -
               keplerianEllipticalElements1.getInclination( ) ) >= errorTolerance_ ||
         fabs( keplerianEllipticalElements2.getArgumentOfPeriapsis( ) -
               keplerianEllipticalElements1.getArgumentOfPeriapsis( ) ) >= errorTolerance_ ||
         fabs( keplerianEllipticalElements2.getLongitudeOfAscendingNode( ) -
               keplerianEllipticalElements1.getLongitudeOfAscendingNode( ) ) >= errorTolerance_ ||
         fabs( keplerianEllipticalElements2.getTrueAnomaly( ) -
               keplerianEllipticalElements1.getTrueAnomaly( ) ) >= errorTolerance_ ||
         fabs( ( keplerianEllipticalElements2.getSemiLatusRectum( ) -
                 keplerianEllipticalElements1.getSemiLatusRectum( ) ) /
               keplerianEllipticalElements1.getSemiLatusRectum( ) ) >= errorTolerance_ )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The orbital element conversion for an elliptical orbit is erroneous." << endl;
    }

    // *************************************************************************
    // Parabolic orbit case around Mars.
    // *************************************************************************

    // From Keplerian to Cartesian.
    KeplerianElements keplerianParabolicElements1;

    // Define Keplerian elements.
    keplerianParabolicElements1.setSemiLatusRectum(
            unit_conversions::convertAstronomicalUnitsToMeters( 4.0 ) );
    keplerianParabolicElements1.setEccentricity( 1.0 );
    keplerianParabolicElements1.setInclination( M_PI / 6.0 );
    keplerianParabolicElements1.setArgumentOfPeriapsis( M_PI / 8.0 );
    keplerianParabolicElements1.setLongitudeOfAscendingNode( 8.0 * M_PI / 7.0 );
    keplerianParabolicElements1.setTrueAnomaly( 7.0 * M_PI / 4.0 );

    // Compute Cartesian elements.
    CartesianElements cartesianParabolicElements = convertKeplerianToCartesianElements(
                &keplerianParabolicElements1, &predefinedMars );

    // From Cartesian to Keplerian.
    // Compute Keplerian elements.
    KeplerianElements keplerianParabolicElements2 =  convertCartesianToKeplerianElements(
                &cartesianParabolicElements, &predefinedMars );

    // Set test result to false if the output Keplerian elements of the
    // conversion from Cartesian to Keplerian are equal to the input Keplerian
    // elements of the conversion from Keplerian to Cartesian, within a
    // tolerance limit.
    if ( fabs( ( keplerianParabolicElements2.getSemiLatusRectum( ) -
                 keplerianParabolicElements1.getSemiLatusRectum( ) ) /
               keplerianParabolicElements1.getSemiLatusRectum( ) ) >=  errorTolerance_ ||
         fabs( keplerianParabolicElements2.getEccentricity( ) -
               keplerianParabolicElements1.getEccentricity( ) ) >= errorTolerance_ ||
         fabs( keplerianParabolicElements2.getInclination( ) -
               keplerianParabolicElements1.getInclination( ) ) >= errorTolerance_ ||
         fabs( keplerianParabolicElements2.getArgumentOfPeriapsis( ) -
               keplerianParabolicElements1.getArgumentOfPeriapsis( ) ) >= errorTolerance_ ||
         fabs( keplerianParabolicElements2.getLongitudeOfAscendingNode( ) -
               keplerianParabolicElements1.getLongitudeOfAscendingNode( ) ) >= errorTolerance_ ||
         fabs( keplerianParabolicElements2.getTrueAnomaly( ) -
               keplerianParabolicElements1.getTrueAnomaly( ) ) >= errorTolerance_ )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The orbital element conversion for a parabolic orbit is erroneous." << endl;
    }

    // *************************************************************************
    // Circular equatorial orbit case around the Earth.
    // *************************************************************************

    // From Keplerian to Cartesian.
    KeplerianElements keplerianCircularElements1;

    // Define Keplerian elements.
    keplerianCircularElements1.setSemiMajorAxis(
            unit_conversions::convertKilometersToMeters( 7.0e3 ) );
    keplerianCircularElements1.setEccentricity( 0.0 );
    keplerianCircularElements1.setInclination( 0.0 );
    keplerianCircularElements1.setArgumentOfPeriapsis( 0.0 );
    keplerianCircularElements1.setLongitudeOfAscendingNode( 0.0 );
    keplerianCircularElements1.setTrueAnomaly( M_PI / 4.0 );

    // Compute Cartesian elements.
    CartesianElements cartesianCircularElements = convertKeplerianToCartesianElements(
                &keplerianCircularElements1, &predefinedEarth );

    // From Cartesian to Keplerian.
    // Compute Keplerian elements.
    KeplerianElements keplerianCircularElements2 = convertCartesianToKeplerianElements(
                &cartesianCircularElements, &predefinedEarth );

    // Set test result to false if the output Keplerian elements of the
    // conversion from Cartesian to Keplerian are equal to the input Keplerian
    // elements of the conversion from Keplerian to Cartesian, within a
    // tolerance limit.

    if ( fabs( ( keplerianCircularElements2.getSemiMajorAxis( ) -
                 keplerianCircularElements1.getSemiMajorAxis( ) ) /
               keplerianCircularElements1.getSemiMajorAxis( ) ) >= errorTolerance_ ||
         fabs( keplerianCircularElements2.getEccentricity( ) -
               keplerianCircularElements1.getEccentricity( ) ) >= errorTolerance_ ||
         fabs( keplerianCircularElements2.getInclination( ) -
               keplerianCircularElements1.getInclination( ) ) >=  errorTolerance_ ||
         fabs( keplerianCircularElements2.getArgumentOfPeriapsis( ) -
               keplerianCircularElements1.getArgumentOfPeriapsis( ) ) >= errorTolerance_ ||
         fabs( keplerianCircularElements2.getLongitudeOfAscendingNode( ) -
               keplerianCircularElements1.getLongitudeOfAscendingNode( ) ) >= errorTolerance_ ||
         fabs( keplerianCircularElements2.getTrueAnomaly( ) -
               keplerianCircularElements1.getTrueAnomaly( ) ) >= errorTolerance_ )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The orbital element conversion for a circular orbit is erroneous." << endl;
    }

    // *************************************************************************
    // Hyperbolic equatorial orbit case around the Sun.
    // *************************************************************************

    // From Keplerian to Cartesian.
    KeplerianElements keplerianHyperbolicElements1;

    // Define Keplerian elements.
    keplerianHyperbolicElements1.setSemiMajorAxis(
                unit_conversions::convertAstronomicalUnitsToMeters( -3.0 ) );
    keplerianHyperbolicElements1.setEccentricity( 2.0 );
    keplerianHyperbolicElements1.setInclination( 0.0 );
    keplerianHyperbolicElements1.setArgumentOfPeriapsis( 11.0 * M_PI / 8.0 );
    keplerianHyperbolicElements1.setLongitudeOfAscendingNode( 0.0 );
    keplerianHyperbolicElements1.setTrueAnomaly( 9.0 * M_PI / 16.0 );

    // Compute Cartesian elements.
    CartesianElements cartesianHyperbolicElements = convertKeplerianToCartesianElements(
                &keplerianHyperbolicElements1, &customDefinedSun );

    // From Cartesian to Keplerian.
    // Compute Keplerian elements.
    KeplerianElements keplerianHyperbolicElements2 = convertCartesianToKeplerianElements(
                &cartesianHyperbolicElements, &customDefinedSun );

    // Set test result to false if the output Keplerian elements of the
    // conversion from Cartesian to Keplerian are equal to the input Keplerian
    // elements of the conversion from Keplerian to Cartesian, within a
    // tolerance limit.

    if ( fabs( ( keplerianHyperbolicElements2.getSemiMajorAxis( ) -
                 keplerianHyperbolicElements1.getSemiMajorAxis( ) ) /
               keplerianHyperbolicElements1.getSemiMajorAxis( ) ) >= errorTolerance_ ||
         fabs( keplerianHyperbolicElements2.getEccentricity( ) -
               keplerianHyperbolicElements1.getEccentricity( ) ) >= errorTolerance_ ||
         fabs( keplerianHyperbolicElements2.getInclination( ) -
               keplerianHyperbolicElements1.getInclination( ) ) >= errorTolerance_ ||
         fabs( keplerianHyperbolicElements2.getArgumentOfPeriapsis( ) -
               keplerianHyperbolicElements1.getArgumentOfPeriapsis( ) ) >= errorTolerance_ ||
         fabs( keplerianHyperbolicElements2.getLongitudeOfAscendingNode( ) -
               keplerianHyperbolicElements1.getLongitudeOfAscendingNode( ) ) >= errorTolerance_ ||
         fabs( keplerianHyperbolicElements2.getTrueAnomaly( ) -
               keplerianHyperbolicElements1.getTrueAnomaly( ) ) >= errorTolerance_ )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The orbital element conversion for a hyperbolic orbit is erroneous." << endl;
    }

    // *************************************************************************
    // Book example.
    // *************************************************************************

    // Define tolerance, related to the precision of the values in the book.
    double errorToleranceBookExample_ = 1.0e-04;

    // From Cartesian to Keplerian.
    CartesianElements cartesianElements;

    // Define Cartesian elements.
    // Position expressed in canonical units.
    cartesianElements.setCartesianElementX( 1.0 );
    cartesianElements.setCartesianElementY( 2.0 );
    cartesianElements.setCartesianElementZ( 1.0 );

    // Velocity expressed in canonical units.
    cartesianElements.setCartesianElementXDot( -0.25 );
    cartesianElements.setCartesianElementYDot( -0.25 );
    cartesianElements.setCartesianElementZDot( 0.5 );

    // Define Keplerian elements.
    KeplerianElements keplerianElements;

    // Convert Cartesian to Keplerian elements.
    // Gravitational parameter is equal to 1 in the applied units.
    keplerianElements = convertCartesianToKeplerianElements( &cartesianElements,
                                                             &customDefinedBody );

    // Set test result to false if the output Keplerian elements of the
    // conversion from Cartesian to Keplerian are equal to the output Keplerian
    // elements of the exercise, within a tolerance limit.
    if ( fabs( keplerianElements.getSemiMajorAxis( ) - 2.265 ) >= errorToleranceBookExample_ ||
         fabs( keplerianElements.getEccentricity( ) - 0.185 ) >= errorToleranceBookExample_ ||
         fabs( keplerianElements.getInclination( ) - 1.401 ) >= errorToleranceBookExample_ ||
         fabs( keplerianElements.getArgumentOfPeriapsis( ) - 2.6143 ) >=
         errorToleranceBookExample_ ||
         fabs( keplerianElements.getLongitudeOfAscendingNode( ) - 1.0304 ) >=
         errorToleranceBookExample_ ||
         fabs( keplerianElements.getTrueAnomaly( ) - 4.0959 ) >= errorToleranceBookExample_ )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The orbital element conversion for the book example is erroneous." << endl;
    }

    // Test 2: Test of true anomaly to eccentric anomaly conversion.
    // Source: http://www.astro.uu.nl/~strous/AA/en/reken/kepler.html.

    // Set tolerance for conversion.
    double toleranceOrbitalElementConversion = 1e-8;

    // Set eccentricity.
    double eccentricity = 0.01671;

    // Set true anomaly.
    double trueAnomaly = unit_conversions::convertDegreesToRadians( 61.6755418 );

    // Compute eccentric anomaly.
    double eccentricAnomaly = orbital_element_conversions::convertTrueAnomalyToEccentricAnomaly(
                trueAnomaly, eccentricity );

    // Check if computed eccentric anomaly is equal to reference value.
    if ( fabs( eccentricAnomaly - 1.061789204 ) > toleranceOrbitalElementConversion )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of true anomaly to eccentric anomaly is "
             << "erroneous as the computed eccentric anomaly after applying the conversion ( "
             << unit_conversions::convertRadiansToDegrees( eccentricAnomaly )
             << " ) does not match the expected value of the eccentric anomaly ( "
             << unit_conversions::convertRadiansToDegrees( 1.061789204 ) << " ) " << endl;
    }

    // Test 3: Test of eccentric anomaly to true anomaly conversion.
    // Source: http://www.astro.uu.nl/~strous/AA/en/reken/kepler.html.

    // Set tolerance for conversion.
    toleranceOrbitalElementConversion = 1e-8;

    // Set eccentricity.
    eccentricity = 0.01671;

    // Set eccentric anomaly.
    eccentricAnomaly = 1.061789204;

    // Compute true anomaly.
    trueAnomaly = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                eccentricAnomaly, eccentricity );

    // Check if computed true anomaly is equal to reference value.
    if ( fabs( trueAnomaly - unit_conversions::convertDegreesToRadians( 61.6755418 ) )
         > toleranceOrbitalElementConversion )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of eccentric anomaly to true anomaly is "
             << "erroneous as the computed true anomaly after applying the conversion ( "
             << unit_conversions::convertRadiansToDegrees( trueAnomaly )
             << " ) does not match the expected value of the true anomaly "
             << "( " << 61.6755418 << " ) " << endl;
    }

    // Test 4: Test of true anomaly to hyperbolic eccentric anomaly conversion.
    // Source: ( Fortescue, 2003 ).

    // Set tolerance for orbital element conversion.
    toleranceOrbitalElementConversion = 1e-4;

    // Set eccentricity.
    eccentricity = 3.0;

    // Set true anomaly.
    trueAnomaly = 0.5291;

    // Compute hyperbolic eccentric anomaly.
    double hyperbolicEccentricAnomaly = orbital_element_conversions::
            convertTrueAnomalyToHyperbolicEccentricAnomaly( trueAnomaly, eccentricity );

    // Check if computed hyperbolic eccentric anomaly is equal to reference
    // value.
    if ( fabs( hyperbolicEccentricAnomaly - 0.3879 ) > toleranceOrbitalElementConversion )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of true anomaly to hyperbolic eccentric "
             << "anomaly is erroneous as the computed hyperbolic eccentric "
             << "anomaly after applying the conversion ( "
             << unit_conversions::convertRadiansToDegrees( hyperbolicEccentricAnomaly )
             << " ) does not match the expected value of the hyperbolic eccentric anomaly ( "
             << unit_conversions::convertRadiansToDegrees( 0.3879 ) << " ) " << endl;
    }

    // Test 5: Test of hyperbolic eccentric anomaly to true anomaly conversion.
    // Source: ( Fortescue, 2003 ).

    // Set tolerance for orbital element conversion.
    toleranceOrbitalElementConversion = 1e-4;

    // Set eccentricity.
    eccentricity = 3.0;

    // Set hyperbolic eccentric anomaly.
    hyperbolicEccentricAnomaly = 0.3879;

    // Compute true anomaly.
    trueAnomaly = orbital_element_conversions::convertHyperbolicEccentricAnomalyToTrueAnomaly(
                hyperbolicEccentricAnomaly, eccentricity );

    // Check if computed true anomaly is equal to reference value.
    if ( fabs( trueAnomaly - 0.5291 ) > toleranceOrbitalElementConversion )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of hyperbolic eccentric anomaly to true "
             << "anomaly is erroneous as the computed true anomaly after "
             << "applying the conversion ( "
             << unit_conversions::convertRadiansToDegrees( trueAnomaly )
             << " ) does not match the expected value of the true anomaly ( "
             << unit_conversions::convertRadiansToDegrees( 0.5291 ) << " ) " << endl;
    }

    // Test 6: Test of eccentric anomaly to mean anomaly conversion.
    // Source: ( Vallado, 2004 ).

    // Set tolerance for conversion.
    toleranceOrbitalElementConversion = 1e-8;

    // Set eccentricity.
    eccentricity = 0.01671;

    // Set eccentric anomaly.
    eccentricAnomaly = 1.061789204;

    // Compute mean anomaly.
    double meanAnomaly = orbital_element_conversions::convertEccentricAnomalyToMeanAnomaly(
                eccentricAnomaly, eccentricity );

    // Check if computed mean anomaly is equal to reference value.
    if ( fabs( meanAnomaly - unit_conversions::convertDegreesToRadians( 60.0 ) )
         > toleranceOrbitalElementConversion )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of eccentric anomaly to mean anomaly is "
             << "erroneous as the computed mean anomaly after applying the conversion ( "
             << unit_conversions::convertRadiansToDegrees( meanAnomaly )
             << " ) does not match the expected value of the mean anomaly "
             << "( " << 60.0 << " ) " << endl;
    }

    // Test 7: Test of mean anomaly to eccentric anomaly conversion.
    // Source: ( Vallado, 2004 ).

    // Set tolerance for conversion.
    toleranceOrbitalElementConversion = 1e-8;

    // Set eccentricity.
    eccentricity = 0.01671;

    // Set mean anomaly.
    meanAnomaly = unit_conversions::convertDegreesToRadians( 60.0 );

    // Create object for mean anomaly to eccentric anomaly conversion.
    orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly
            convertMeanAnomalyToEccentricAnomaly;

    // Create pointer to Newton-Raphson object.
    NewtonRaphson* pointerToNewtonRaphson = new NewtonRaphson;

    // Set eccentricity.
    convertMeanAnomalyToEccentricAnomaly.setEccentricity( 0.01671 );

    // Set mean anomaly.
    convertMeanAnomalyToEccentricAnomaly.setMeanAnomaly( meanAnomaly );

    // Set Newton-Raphson method.
    convertMeanAnomalyToEccentricAnomaly.setNewtonRaphson( pointerToNewtonRaphson );

    // Compute eccentric anomaly.
    eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly.convert( );

    // Check if computed eccentric anomaly is equal to reference value.
    if ( fabs( eccentricAnomaly - 1.061789204 ) > toleranceOrbitalElementConversion )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of mean anomaly to eccentric anomaly is "
             << "erroneous as the computed eccentric anomaly after applying the conversion ( "
             << unit_conversions::convertRadiansToDegrees( eccentricAnomaly )
             << " ) does not match the expected value of the eccentric anomaly ( "
             << unit_conversions::convertRadiansToDegrees( 1.061789204 ) << " ) " << endl;
    }

    // Test 8: Test of hyperbolic eccentric anomaly to mean anomaly conversion.
    // Source: ( Vallado, 2004 ).

    // Set tolerance for conversion.
    toleranceOrbitalElementConversion = 1e-8;

    // Set eccentricity.
    eccentricity = 2.4;

    // Set hyperbolic eccentric anomaly.
    hyperbolicEccentricAnomaly = 1.6013761449;

    // Compute mean anomaly.
    meanAnomaly = orbital_element_conversions::convertHyperbolicEccentricAnomalyToMeanAnomaly(
                hyperbolicEccentricAnomaly, eccentricity );

    // Check if computed mean anomaly is equal to reference value.
    if ( fabs( meanAnomaly - unit_conversions::convertDegreesToRadians( 235.4  ) )
         > toleranceOrbitalElementConversion )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of hyperbolic eccentric anomaly to mean "
             << "anomaly is erroneous as the computed mean anomaly after "
             << "applying the conversion ( "
             << unit_conversions::convertRadiansToDegrees( hyperbolicEccentricAnomaly )
             << " ) does not match the expected value of the mean anomaly "
             << "( " << 235.4 << " ) " << endl;
    }

    // Test 9: Test of mean anomaly to hyperbolic eccentric anomaly conversion.
    // Source: ( Vallado, 2004 ).

    // Set tolerance for conversion.
    toleranceOrbitalElementConversion = 1e-8;

    // Set eccentricity.
    eccentricity = 2.4;

    // Set mean anomaly.
    meanAnomaly = unit_conversions::convertDegreesToRadians( 235.4 );

    // Create object for mean anomaly to hyperbolic eccentric anomaly
    // conversion.
    orbital_element_conversions::ConvertMeanAnomalyToHyperbolicEccentricAnomaly
            convertMeanAnomalyToHyperbolicEccentricAnomaly;

    // Set eccentricity.
    convertMeanAnomalyToHyperbolicEccentricAnomaly.setEccentricity( 2.4 );

    // Set mean anomaly.
    convertMeanAnomalyToHyperbolicEccentricAnomaly.setMeanAnomaly( meanAnomaly );

    // Set Newton-Raphson method.
    convertMeanAnomalyToHyperbolicEccentricAnomaly.setNewtonRaphson( pointerToNewtonRaphson );

    // Compute hyperbolic eccentric anomaly.
    hyperbolicEccentricAnomaly = convertMeanAnomalyToHyperbolicEccentricAnomaly.convert( );

    // Check if computed hyperbolic eccentric anomaly is equal to reference
    // value.
    if ( fabs( hyperbolicEccentricAnomaly - 1.6013761449 ) > toleranceOrbitalElementConversion )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of mean anomaly to hyperbolic eccentric "
             << "anomaly is erroneous as the computed hyperbolic eccentric "
             << "anomaly after applying the conversion ( "
             << unit_conversions::convertRadiansToDegrees( hyperbolicEccentricAnomaly )
             << " ) does not match the expected value of the hyperbolic "
             << "eccentric anomaly ( "
             << unit_conversions::convertRadiansToDegrees( 1.6013761449 ) << " ) " << endl;
    }

    // Test 10: Test of elapsed time to mean anomaly for elliptical orbits.

    // Set tolerance for conversion.
    toleranceOrbitalElementConversion = 1.0e-11;

    // Expected mean anomaly value;
    double expectedMeanAnomalyForTest10 = 20.203139666972554;

    // Set elapsed time.
    double expectedElapsedTime = 4000.0;

    // Set semi-major axis.
    double semiMajorAxis = unit_conversions::convertKilometersToMeters( 2500.0 );

    // Compute mean anomaly.
    meanAnomaly = orbital_element_conversions::convertElapsedTimeToMeanAnomalyForEllipticalOrbits(
                expectedElapsedTime, &predefinedEarth, semiMajorAxis );

    // Declare and compute absolute and relative errors.
    double absoluteDifference = fabs( meanAnomaly - expectedMeanAnomalyForTest10 );

    double relativeDifference = absoluteDifference / expectedMeanAnomalyForTest10;

    // Check if relative error is too large.
    if ( relativeDifference > MACHINE_PRECISION_DOUBLES )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of elapsed time to mean anomaly is erroneous "
             << "as the computed mean anomaly after applying the conversion ( "
             << unit_conversions::convertRadiansToDegrees( meanAnomaly )
             << " ) does not match the expected value of the mean anomaly ( "
             << unit_conversions::convertRadiansToDegrees( expectedMeanAnomalyForTest10 )
             << " ) " << endl;
    }

    // Test 11: Test of mean anomaly to elapsed time for elliptical orbits.
    //          Reversal of computation for Test 10.

    // Set tolerance for conversion.
    toleranceOrbitalElementConversion = 1e-11;

    // Set mean anomaly.
    meanAnomaly = expectedMeanAnomalyForTest10;

    // Declare and compute elapsed time.
    double elapsedTime = orbital_element_conversions::
            convertMeanAnomalyToElapsedTimeForEllipticalOrbits(
                meanAnomaly, &predefinedEarth, semiMajorAxis );

    // Compute absolute and relative errors.
    absoluteDifference = fabs( elapsedTime - expectedElapsedTime );

    relativeDifference = absoluteDifference / expectedElapsedTime;

    // Check if computed elapsed time is equal to reference value.
    if ( relativeDifference > toleranceOrbitalElementConversion )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of mean anomaly to elapsed time is erroneous "
             << "as the computed elapsed time after applying the conversion ( " << elapsedTime
             << " ) does not match the expected value of the elapsed time ( "
             << expectedElapsedTime << " ) " << endl;
    }

    // Test 12: Test of elapsed time to mean anomaly for hyperbolic orbits.

    // Set tolerance for conversion.
    toleranceOrbitalElementConversion = 1e-11;

    // Set elapsed time.
    elapsedTime = 1000.0;

    // Set semi-major axis.
    semiMajorAxis = unit_conversions::convertKilometersToMeters( -40000.0 );

    // Compute mean anomaly.
    meanAnomaly = orbital_element_conversions::
            convertElapsedTimeToMeanAnomalyForHyperbolicOrbits(
                elapsedTime, &predefinedEarth, semiMajorAxis );

    // Check if computed mean anomaly is equal to reference value.
    if ( fabs( meanAnomaly - 0.078918514324112 ) > toleranceOrbitalElementConversion )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of elapsed time to mean anomaly is erroneous "
             << "as the computed mean anomaly after applying the conversion ( "
             << unit_conversions::convertRadiansToDegrees( meanAnomaly )
             << " ) does not match the expected value of the mean anomaly ( "
             << unit_conversions::convertRadiansToDegrees( 0.078918514324112 ) << " ) " << endl;
    }

    // Test 13: Test of mean anomaly to elapsed time for hyperbolic orbits.

    // Set tolerance for conversion.
    toleranceOrbitalElementConversion = 1e-11;

    // Set mean anomaly.
    meanAnomaly = 0.078918514324112;

    // Set semi-major axis.
    semiMajorAxis = unit_conversions::convertKilometersToMeters( -40000.0 );

    // Compute elapsed time.
    elapsedTime = orbital_element_conversions::
            convertMeanAnomalyToElapsedTimeForHyperbolicOrbits(
                meanAnomaly, &predefinedEarth, semiMajorAxis );

    // Check if computed elapsed time is equal to reference value.
    if ( fabs( elapsedTime - 1000.0 ) > toleranceOrbitalElementConversion )
    {
        isOrbitalElementConversionErroneous = true;

        cerr << "The conversion of mean anomaly to elapsed time is erroneous "
             << "as the computed elapsed time after applying the conversion ( " << elapsedTime
             << " ) does not match the expected value of the elapsed time ( "
             << 1000.0 << " ) " << endl;
    }

    return isOrbitalElementConversionErroneous;
}

}

// End of file.
