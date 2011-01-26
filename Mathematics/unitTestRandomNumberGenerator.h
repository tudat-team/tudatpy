/*! \file unitTestRandomNumberGenerator.h
 *    Header file that defines a unit test that tests the random number
 *    generator included in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Date created      : 7 January, 2010
 *    Last modified     : 7 January, 2010
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
 *      YYMMDD    author        comment
 *      110107    K. Kumar      First creation of code.
 */

#ifndef UNITTESTRANDOMNUMBERGENERATOR_H
#define UNITTESTRANDOMNUMBERGENERATOR_H

// Include statements.
#include <ctime>
#include <cmath>
#include "randomNumberGenerator.h"
#include "basicMathematicsFunctions.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test of implementation of random number generator class.
/*!
 * Test of implementation of random number generator class.
 * \return Integer indicating success of test ( 0 = successful; 1 = failed ).
 */
int testRandomNumberGenerator( );

}

#endif // UNITTESTRANDOMNUMBERGENERATOR_H

// End of file.
