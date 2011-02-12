/*! \file unitTestGravityAssist.h
 *    Header file of unit test file of gravity assist code. This unit test file
 *    will test the gravity assist code.
 *
 *    Path              : /Astrodynamics/MissionsSegments/GravityAssist/
 *    Version           : 2
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
 *    Date created      : 17 January, 2011
 *    Last modified     : 5 February, 2011
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
 *      YYMMDD    Author            Comment
 *      110117    E. Iorfida        First creation of the code.
 *      110205    J. Melman         Removed 'Powered'.
 */

#ifndef UNITTESTGRAVITYASSIST_H
#define UNITTESTGRAVITYASSIST_H

// Include statements.
#include "gravityAssist.h"
#include "unitConversions.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test of gravity assist code.
/*!
 * Test of gravity assist code.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testGravityAssist( );

}

#endif // UNITTESTGRAVITYASSIST_H

// End of file.
