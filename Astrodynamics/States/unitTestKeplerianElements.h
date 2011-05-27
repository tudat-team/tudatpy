/*! \file unitTestKeplerianElements.h
 *    Header file for a unit test that tests the implementation of the
 *    Keplerian elements state class in Tudat.
 *
 *    Path              : /Astrodynamics/States/
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
 *    Date created      : 10 January, 2010
 *    Last modified     : 10 January, 2010
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
 *      110110    K. Kumar          File created.
 */

#ifndef UNITTESTKEPLERIANELEMENTS_H
#define UNITTESTKEPLERIANELEMENTS_H

// Include statements.
#include "basicMathematicsFunctions.h"
#include "keplerianElements.h"
#include "linearAlgebra.h"
#include "unitConversions.h"

// Using statements.
using unit_conversions::convertDegreesToRadians;

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test Keplerian elements state class.
/*!
 * Tests Keplerian elements state class.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testKeplerianElements( );

}

// End of file.

#endif // UNITTESTKEPLERIANELEMENTS_H
