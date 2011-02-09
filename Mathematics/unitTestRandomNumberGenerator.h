/*! \file unitTestRandomNumberGenerator.h
 *    Header file that defines a unit test that tests the random number
 *    generator included in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 2
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
 *    Date created      : 7 January, 2011
 *    Last modified     : 7 February, 2011
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
 *      110107    K. Kumar          First creation of code.
 *      110207    K. Kumar          Updated code to conform to protocol.
 */

#ifndef UNITTESTRANDOMNUMBERGENERATOR_H
#define UNITTESTRANDOMNUMBERGENERATOR_H

// Include statements.
#include <ctime>
#include <cmath>
#include <iostream>
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
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testRandomNumberGenerator( );

}

#endif // UNITTESTRANDOMNUMBERGENERATOR_H

// End of file.
