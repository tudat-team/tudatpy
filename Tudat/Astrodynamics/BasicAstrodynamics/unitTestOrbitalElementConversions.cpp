/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *
 *    References
 *      http://www.astro.uu.nl/~strous/AA/en/reken/kepler.html, last accessed: 16th February, 2011.
 *      Vallado, D. A., McClain, W. D. Fundamentals of astrodynamics and applications, 2nd Edition,
 *          Kluwer Academic Publishers, The Netherlands, 2004.
 *      Fortescue, P. W., et al. Spacecraft systems engineering, Third Edition,
 *          Wiley, England, 2003.
 *
 */

// Temporary notes (move to class/function doxygen):
// Test runs code and verifies result against expected value.
// If the tested code is erroneous, the test function returns a boolean
// true; if the code is correct, the function returns a boolean false.
// 

#include <cmath>
#include <iostream>
#include <limits>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToHyperbolicEccentricAnomaly.h"
#include "Tudat/Astrodynamics/Bodies/celestialBody.h"
#include "Tudat/Astrodynamics/Bodies/planet.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"

//! Test orbital element conversion code.
int main( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;
    using std::fabs;
    using std::pow;
    using tudat::orbital_element_conversions::convertCartesianToKeplerianElements;
    using tudat::orbital_element_conversions::convertKeplerianToCartesianElements;
    using tudat::orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly;
    using tudat::orbital_element_conversions::ConvertMeanAnomalyToHyperbolicEccentricAnomaly;
    using namespace tudat;

    // Test of orbital element conversion methods imeplemented in Tudat.
    // Test 1: Test of mean anomaly to eccentric anomaly conversion.

    // Initialize unit test result to false.
    bool isOrbitalElementConversionErroneous = false;

    // Test 1: Test of mean anomaly to eccentric anomaly conversion.
    // Source: ( Vallado, 2004 ).

    // Set tolerance for conversion.
    double toleranceOrbitalElementConversion = 1e-8;

    // Set eccentricity.
    double eccentricity = 0.01671;

    // Set mean anomaly.
    double meanAnomaly = unit_conversions::convertDegreesToRadians( 60.0 );

    // Create Newton-Raphson object.
    NewtonRaphson newtonRaphson;

    // Create object for mean anomaly to eccentric anomaly conversion.
    orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly
            convertMeanAnomalyToEccentricAnomaly( eccentricity, meanAnomaly, &newtonRaphson );

    // Compute eccentric anomaly.
    double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly.convert( );

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

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isOrbitalElementConversionErroneous )
    {
        cerr << "testOrbitalElementConversions failed!" << endl;
    }

    return isOrbitalElementConversionErroneous;
}
