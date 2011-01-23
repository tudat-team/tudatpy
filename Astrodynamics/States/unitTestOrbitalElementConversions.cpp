/*! \file unitTestOrbitalElementConversion.cpp
 *    Source file of unit test for the orbitalElementConversion, from Cartesian to
 *    Keplerian and viceversa.
 *    The first part of the code tests the code for elliptical, parabolic,
 *    hyperbolic and circular orbits. SI units are used.
 *    The second part of the code tests the code from Cartesian to Keplerian with
 *    the example 3.4 pag. 63 of the book "Fondamenti di Meccanica del Volo
 *    Spaziale" (G. Mengali, A.A. Quarta). In this part of the code, canonical
 *    units are used.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 9
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 03 December, 2010
 *    Last modified     : 11 January, 2011
 *
 *    References
 *
 *    Notes
 *      Test calls function and verifies result against expected value.
 *      If the test is successful, the main function returns an integer value
 *      0; if the test fails, the main function returns an integer value 1.
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
 *      YYMMDD    author         comment
 *      101203    E. Iorfida     First creation of the code.
 *      101208    E. Iorfida     Fulfillment of the code with the elliptical case.
 *      101208    E. Iorfida     Modified punctuation.
 *      101215    E. Iorfida     Added tolerance, added parabolic, circular and
 *                               hyperbolic cases.
 *      101217    E. Iorfida     Added fabs( ) in the errors computation,
 *                               modified punctuation.
 *      101219    J. Melman      Put gravitational parameters in one place,
 *                               changed first right ascension to 15.0 * pi / 8.0,
 *                               thereby exposing a possible error.
 *      110107    E. Iorfida     orbitalConversionBookExampleUnitTest.test added to
 *                               this file, to have a unique unit test file for the
 *                               conversion code. Also some punctuation modifications
 *                               have been made.
 *      110109    J. Melman      Included test for semi-latus rectum of circular case.
 *                               Gave the orbital angles less trivial values, and not
 *                               almost exclusively in the first quadrant.
 *      110111    E. Iorfida     Updated to the new format of unitTest file and
 *                               added hyperbolic equatorial case.
 */

// Include statements.
#include "orbitalElementConversions.h"
#include "unitConversions.h"

// Using directives.
using std::cerr;
using std::endl;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test of orbitalElementConversion code.
bool testOrbitalElementConversion( )
{
    // Test result initialised to false.
    bool isOrbitalElementConversionErroneous = false;

    // Define tolerance.
    double errorTolerance_ = 1.0e2 * mathematics::MACHINE_PRECISION_DOUBLES;

    // Define gravitational parameters.
    double gravitationalParameterEarth = 398600.4415e9;
    double gravitationalParameterMars  = 42828e9;
    double gravitationalParameterSun   = 132712440018e8;

    // *************************************************************************
    // Elliptical orbit case around the Earth.
    // *************************************************************************

    // From Keplerian to Cartesian.
    KeplerianElements myKeplerianEllipticalElements1_;

    // Define Keplerian elements.
    myKeplerianEllipticalElements1_.setSemiMajorAxis(
            unit_conversions::convertAstronomicalUnitsToMeters( 0.3 ) );
    myKeplerianEllipticalElements1_.setEccentricity( 0.2 );
    myKeplerianEllipticalElements1_.setInclination( M_PI / 4.0  );
    myKeplerianEllipticalElements1_.setArgumentOfPeriapsis( 4.0 * M_PI / 3.0 );
    myKeplerianEllipticalElements1_.
            setRightAscensionOfAscendingNode( 1.0 * M_PI / 8.0 );
    myKeplerianEllipticalElements1_.setTrueAnomaly( 1.0 * M_PI / 3.0 );
    myKeplerianEllipticalElements1_.setSemiLatusRectum(
            myKeplerianEllipticalElements1_.getSemiMajorAxis( ) * ( 1.0 -
            mathematics::raiseToIntegerPower(
                    myKeplerianEllipticalElements1_.getEccentricity( ), 2 ) ) );

    // Compute Cartesian elements.
    CartesianElements myCartesianEllipticalElements_;

    myCartesianEllipticalElements_ =
            orbital_element_conversions::convertKeplerianToCartesianElements(
            myKeplerianEllipticalElements1_, gravitationalParameterEarth );

    // From Cartesian to Keplerian.
    // Compute Keplerian elements.
    KeplerianElements myKeplerianEllipticalElements2_;

    myKeplerianEllipticalElements2_ =
            orbital_element_conversions::convertCartesianToKeplerianElements(
            myCartesianEllipticalElements_, gravitationalParameterEarth );

    // Set test result to false if the output Keplerian elements of the
    // conversion from Cartesian to Keplerian are equal to the input Keplerian
    // elements of the conversion from Keplerian to Cartesian, within a
    // tolerance limit.
    if ( fabs( ( myKeplerianEllipticalElements2_.getSemiMajorAxis( ) -
                 myKeplerianEllipticalElements1_.getSemiMajorAxis( ) ) /
               myKeplerianEllipticalElements1_.getSemiMajorAxis( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianEllipticalElements2_.getEccentricity( ) -
               myKeplerianEllipticalElements1_.getEccentricity( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianEllipticalElements2_.getInclination( ) -
               myKeplerianEllipticalElements1_.getInclination( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianEllipticalElements2_.getArgumentOfPeriapsis( ) -
               myKeplerianEllipticalElements1_.getArgumentOfPeriapsis( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianEllipticalElements2_.
               getRightAscensionOfAscendingNode( ) -
               myKeplerianEllipticalElements1_.
               getRightAscensionOfAscendingNode( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianEllipticalElements2_.getTrueAnomaly( ) -
               myKeplerianEllipticalElements1_.getTrueAnomaly( ) ) >=
         errorTolerance_ ||
         fabs( ( myKeplerianEllipticalElements2_.getSemiLatusRectum( ) -
                 myKeplerianEllipticalElements1_.getSemiLatusRectum( ) ) /
               myKeplerianEllipticalElements1_.getSemiLatusRectum( ) ) >=
         errorTolerance_ )
    {
        isOrbitalElementConversionErroneous = true;
        cerr << "The orbital element conversion for an elliptical orbit is "
             << "erroneous." << endl;
    }

    // *************************************************************************
    // Parabolic orbit case around Mars.
    // *************************************************************************

    // From Keplerian to Cartesian.
    KeplerianElements myKeplerianParabolicElements1_;

    // Define Keplerian elements.
    myKeplerianParabolicElements1_.setSemiLatusRectum(
            unit_conversions::convertAstronomicalUnitsToMeters( 4.0 ) );
    myKeplerianParabolicElements1_.setEccentricity( 1.0 );
    myKeplerianParabolicElements1_.setInclination( M_PI / 6.0 );
    myKeplerianParabolicElements1_.setArgumentOfPeriapsis( M_PI / 8.0 );
    myKeplerianParabolicElements1_.setRightAscensionOfAscendingNode(
            8.0 * M_PI / 7.0 );
    myKeplerianParabolicElements1_.setTrueAnomaly( 7.0 * M_PI / 4.0 );

    // Compute Cartesian elements.
    CartesianElements myCartesianParabolicElements_;

    myCartesianParabolicElements_ =
            orbital_element_conversions::convertKeplerianToCartesianElements(
            myKeplerianParabolicElements1_, gravitationalParameterMars );

    // From Cartesian to Keplerian.
    // Compute Keplerian elements.
    KeplerianElements myKeplerianParabolicElements2_;

    myKeplerianParabolicElements2_ =
            orbital_element_conversions::convertCartesianToKeplerianElements(
            myCartesianParabolicElements_, gravitationalParameterMars );

    // Set test result to false if the output Keplerian elements of the
    // conversion from Cartesian to Keplerian are equal to the input Keplerian
    // elements of the conversion from Keplerian to Cartesian, within a
    // tolerance limit.
    if ( fabs( ( myKeplerianParabolicElements2_.getSemiLatusRectum( ) -
                 myKeplerianParabolicElements1_.getSemiLatusRectum( ) ) /
               myKeplerianParabolicElements1_.getSemiLatusRectum( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianParabolicElements2_.getEccentricity( ) -
               myKeplerianParabolicElements1_.getEccentricity( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianParabolicElements2_.getInclination( ) -
               myKeplerianParabolicElements1_.getInclination( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianParabolicElements2_.getArgumentOfPeriapsis( ) -
               myKeplerianParabolicElements1_.getArgumentOfPeriapsis( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianParabolicElements2_.getRightAscensionOfAscendingNode( ) -
               myKeplerianParabolicElements1_.getRightAscensionOfAscendingNode( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianParabolicElements2_.getTrueAnomaly( ) -
               myKeplerianParabolicElements1_.getTrueAnomaly( ) ) >=
         errorTolerance_ )
    {
        isOrbitalElementConversionErroneous = true;
        cerr << "The orbital element conversion for a parabolic orbit is "
             << "erroneous." << endl;
    }

    // *************************************************************************
    // Circular equatorial orbit case around the Earth.
    // *************************************************************************

    // From Keplerian to Cartesian.
    KeplerianElements myKeplerianCircularElements1_;

    // Define Keplerian elements.
    myKeplerianCircularElements1_.setSemiMajorAxis(
            unit_conversions::convertAstronomicalUnitsToMeters( 0.1 ) );
    myKeplerianCircularElements1_.setEccentricity( 0.0 );
    myKeplerianCircularElements1_.setInclination( 0.0 );
    myKeplerianCircularElements1_.setArgumentOfPeriapsis( 0.0 );
    myKeplerianCircularElements1_.setRightAscensionOfAscendingNode( 0.0 );
    myKeplerianCircularElements1_.setTrueAnomaly( M_PI / 4.0 );

    // Compute Cartesian elements.
    CartesianElements myCartesianCircularElements_;

    myCartesianCircularElements_ =
            orbital_element_conversions::convertKeplerianToCartesianElements(
            myKeplerianCircularElements1_, gravitationalParameterEarth );

    // From Cartesian to Keplerian.
    // Compute Keplerian elements.
    KeplerianElements myKeplerianCircularElements2_;

    myKeplerianCircularElements2_ =
            orbital_element_conversions::convertCartesianToKeplerianElements(
            myCartesianCircularElements_, gravitationalParameterEarth );

    // Set test result to false if the output Keplerian elements of the
    // conversion from Cartesian to Keplerian are equal to the input Keplerian
    // elements of the conversion from Keplerian to Cartesian, within a
    // tolerance limit.
    if ( fabs( ( myKeplerianCircularElements2_.getSemiMajorAxis( ) -
                 myKeplerianCircularElements1_.getSemiMajorAxis( ) ) /
               myKeplerianCircularElements1_.getSemiMajorAxis( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianCircularElements2_.getEccentricity( ) -
               myKeplerianCircularElements1_.getEccentricity( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianCircularElements2_.getInclination( ) -
               myKeplerianCircularElements1_.getInclination( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianCircularElements2_.getArgumentOfPeriapsis( ) -
               myKeplerianCircularElements1_.getArgumentOfPeriapsis( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianCircularElements2_.getRightAscensionOfAscendingNode( ) -
               myKeplerianCircularElements1_.getRightAscensionOfAscendingNode( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianCircularElements2_.getTrueAnomaly( ) -
               myKeplerianCircularElements1_.getTrueAnomaly( ) ) >=
         errorTolerance_ )
    {
        isOrbitalElementConversionErroneous = true;
        cerr << "The orbital element conversion for a circular orbit is "
             << "erroneous." << endl;
    }

    // *************************************************************************
    // Hyperbolic equatorial orbit case around the Sun.
    // *************************************************************************

    // From Keplerian to Cartesian.
    KeplerianElements myKeplerianHyperbolicElements1_;

    // Define Keplerian elements.
    myKeplerianHyperbolicElements1_.setSemiMajorAxis(
            unit_conversions::convertAstronomicalUnitsToMeters( -3.0 ) );
    myKeplerianHyperbolicElements1_.setEccentricity( 2.0 );
    myKeplerianHyperbolicElements1_.setInclination( 0.0 );
    myKeplerianHyperbolicElements1_.setArgumentOfPeriapsis(
            11.0 * M_PI / 8.0 );
    myKeplerianHyperbolicElements1_.setRightAscensionOfAscendingNode(
            0.0 );
    myKeplerianHyperbolicElements1_.setTrueAnomaly( 9.0 * M_PI / 16.0 );

    // Compute Cartesian elements.
    CartesianElements myCartesianHyperbolicElements_;

    myCartesianHyperbolicElements_ =
            orbital_element_conversions::convertKeplerianToCartesianElements(
            myKeplerianHyperbolicElements1_, gravitationalParameterSun );

    // From Cartesian to Keplerian.
    // Compute Keplerian elements.
    KeplerianElements myKeplerianHyperbolicElements2_;

    myKeplerianHyperbolicElements2_ =
            orbital_element_conversions::convertCartesianToKeplerianElements(
            myCartesianHyperbolicElements_, gravitationalParameterSun );

    // Set test result to false if the output Keplerian elements of the
    // conversion from Cartesian to Keplerian are equal to the input Keplerian
    // elements of the conversion from Keplerian to Cartesian, within a
    // tolerance limit.
    if ( fabs( ( myKeplerianHyperbolicElements2_.getSemiMajorAxis( ) -
                 myKeplerianHyperbolicElements1_.getSemiMajorAxis( ) ) /
               myKeplerianHyperbolicElements1_.getSemiMajorAxis( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianHyperbolicElements2_.getEccentricity( ) -
               myKeplerianHyperbolicElements1_.getEccentricity( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianHyperbolicElements2_.getInclination( ) -
               myKeplerianHyperbolicElements1_.getInclination( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianHyperbolicElements2_.getArgumentOfPeriapsis( ) -
               myKeplerianHyperbolicElements1_.getArgumentOfPeriapsis( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianHyperbolicElements2_.getRightAscensionOfAscendingNode( ) -
               myKeplerianHyperbolicElements1_.getRightAscensionOfAscendingNode( ) ) >=
         errorTolerance_ ||
         fabs( myKeplerianHyperbolicElements2_.getTrueAnomaly( ) -
               myKeplerianHyperbolicElements1_.getTrueAnomaly( ) ) >=
         errorTolerance_ )
    {
        isOrbitalElementConversionErroneous = true;
        cerr << "The orbital element conversion for a hyperbolic orbit is "
             << "erroneous." << endl;
    }

    // *************************************************************************
    // Book example.
    // *************************************************************************

    // Define tolerance, related to the precision of the values in the book.
    double errorToleranceBookExample_ = 1.0e-04;

    // From Cartesian to Keplerian.
    CartesianElements myCartesianElements_;

    // Define Cartesian elements.
    // Position vector expressed in canonical units.
    myCartesianElements_.setCartesianElementX( 1.0 );
    myCartesianElements_.setCartesianElementY( 2.0 );
    myCartesianElements_.setCartesianElementZ( 1.0 );

    // Velocity vector expressed in canonical units.
    myCartesianElements_.setCartesianElementXDot( -0.25 );
    myCartesianElements_.setCartesianElementYDot( -0.25 );
    myCartesianElements_.setCartesianElementZDot( 0.5 );

    // Define Keplerian elements.
    KeplerianElements myKeplerianElements_;

    // Convert Cartesian to Keplerian elements.
    // Gravitational parameter is equal to 1 in the applied units.
    myKeplerianElements_ =
            orbital_element_conversions::convertCartesianToKeplerianElements(
                    myCartesianElements_, 1.0 );

    // Set test result to false if the output Keplerian elements of the
    // conversion from Cartesian to Keplerian are equal to the output Keplerian
    // elements of the exercise, within a tolerance limit.
    if ( fabs( myKeplerianElements_.getSemiMajorAxis( ) - 2.265 ) >=
         errorToleranceBookExample_ ||
         fabs( myKeplerianElements_.getEccentricity( ) - 0.185 ) >=
         errorToleranceBookExample_ ||
         fabs( myKeplerianElements_.getInclination( ) - 1.401 ) >=
         errorToleranceBookExample_ ||
         fabs( myKeplerianElements_.getArgumentOfPeriapsis( ) - 2.6143 ) >=
         errorToleranceBookExample_ ||
         fabs( myKeplerianElements_.getRightAscensionOfAscendingNode( ) -
               1.0304 ) >= errorToleranceBookExample_ ||
         fabs( myKeplerianElements_.getTrueAnomaly( ) - 4.0959 ) >=
         errorToleranceBookExample_ )
    {
        isOrbitalElementConversionErroneous = true;
        cerr << "The orbital element conversion for the book example is "
             << "erroneous." << endl;
    }

    return isOrbitalElementConversionErroneous;
}

}

// End of file.
