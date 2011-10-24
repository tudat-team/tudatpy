/*! \file unitTestPhysicalConstants.cpp
 *    This unit test will test the physical constants that are defined in physicalConstants.h.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Date created      : 10 September, 2010
 *    Last modified     : 1 February, 2011
 *
 *    References
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
 *      100910    J. Melman         First creation of code.
 *      110111    J. Melman         Adapted to the offical Tudat standards.
 *      110124    J. Melman         Further adapted to the offical Tudat standards.
 *      110201    J. Melman         Made the tests for obliquity and astronomical unit more
 *                                  accurate.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Astrodynamics/physicalConstants.h"
#include "Astrodynamics/unitTestPhysicalConstants.h"
#include "Mathematics/unitConversions.h"

// Using directives.
using std::cerr;
using std::endl;
using std::fabs;
using mathematics::MACHINE_PRECISION_DOUBLES;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test of physical constants header file.
bool testPhysicalConstants( )
{
    // Test result initialised to false.
    bool isPhysicalConstantsErroneous = false;

    // Test for time constants.
    // Test for the number of days in a year.
    if ( fabs( PhysicalConstants::JULIAN_YEAR_IN_DAYS - 365.25 ) > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The Julian year in days is not set correctly." << endl;
        isPhysicalConstantsErroneous = true;
    }

    // Test for the number of seconds in a year.
    if ( fabs( PhysicalConstants::JULIAN_YEAR - PhysicalConstants::JULIAN_DAY
               * PhysicalConstants::JULIAN_YEAR_IN_DAYS ) > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The Julian year in seconds does not correspond to the Julian "
             << "day in seconds multiplied with the number of days per year."
             << endl;
        isPhysicalConstantsErroneous = true;
    }

    // Test for gravitational constant.
    if ( fabs( PhysicalConstants::GRAVITATIONAL_CONSTANT - 6.67259e-11 )
         > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The gravitational constant is not set correctly." << endl;
        isPhysicalConstantsErroneous = true;
    }

    // Test for the obliquity of the ecliptic.
    // Test for its absolute size (23.5 degrees, from the top of my head).
    if ( fabs( PhysicalConstants::OBLIQUITY_ECLIPTIC_IN_DEGREES - 23.439281 )
         > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The obliquity of the ecliptic is not set correctly." << endl;
        isPhysicalConstantsErroneous = true;
    }

    // Test for its relative size.
    if ( fabs( PhysicalConstants::OBLIQUITY_ECLIPTIC_IN_DEGREES
               - PhysicalConstants::OBLIQUITY_ECLIPTIC_IN_ARCSECONDS / ( 60.0 * 60.0 ) ) > 1.0e-6 )
    {
        cerr << "The obliquity of the ecliptic in arcseconds is not "
             << "set correctly." << endl;
        isPhysicalConstantsErroneous = true;
    }

    // Test for astronomical unit.
    // As expected, indeed approximately equal to 150 million kilometers.
    if ( fabs( PhysicalConstants::ASTRONOMICAL_UNIT - 1.49597870691e11 )
         > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The astronomical unit is not set correctly." << endl;
        isPhysicalConstantsErroneous = true;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isPhysicalConstantsErroneous;
}

}

// End of file.
