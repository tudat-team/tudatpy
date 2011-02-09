/*! \file unitTestPhysicalConstants.h
 *    This unit test will test the physical constants that are
 *    defined in physicalConstants.h.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Date created      : 10 September, 2010
 *    Last modified     : 12 January, 2011
 *
 *    References
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
 *      YYMMDD    author              comment
 *      110111    J. Melman           First creation of code.
 *      110111    J. Melman           Adapted to the offical Tudat standards.
 */

// Include statements.
#include <iostream>
#include "physicalConstants.h"
#include "unitConversions.h"
#include "basicMathematicsFunctions.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test of physical constants header file.
/*!
 * Test of physical constants header file.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testPhysicalConstants( );

}

// End of file.
