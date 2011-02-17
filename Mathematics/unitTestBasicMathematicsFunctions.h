/*! \file unitTestBasicMathematicsFunctions.h
 *    Source file that defines the unitTestBasicMathematicsFunctions unit test,
 *    containing all basic mathematics functions contained in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : B. Romgens *    Affiliation       : Delft University of Technology
 *    E-mail address    : bart.romgens@gmail.com
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 7 February, 2011
 *    Last modified     : 15 February, 2011
 *
 *    References
 *
 *    Notes
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
 *      110207    B. Romgens       File created.
 *      110215    K. Kumar          Minor modifications.
 */

#ifndef UNITTESTBASICMATHEMATICSFUNCTIONS_H
#define UNITTESTBASICMATHEMATICSFUNCTIONS_H

// Include statements.
#include <iostream>
#include <cmath>
#include "basicMathematicsFunctions.h"
#include "cartesianPositionElements.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test of implementation of the basic mathematics functions.
/*!
 * Test of implementation of the basic mathematics functions in
 * basicMathematicsFunctions.cpp.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testBasicMathematicsFunctions( );

}

#endif // UNITTESTBASICMATHEMATICSFUNCTIONS_H

// End of file.
