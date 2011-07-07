/*! \file unitTestExponentialAtmosphere.h
 *    Header file that defines the exponential atmosphere unit test included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/AtmosphereModel/
 *    Version           : 1
 *    Check status      : Unchecked
 *
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 15 March, 2011
 *    Last modified     : 15 March, 2011
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
 *      140311    F.M. Engelen      File created.
 */

#ifndef UNITTESTEXPONENTIALATMOSPHERE_H
#define UNITTESTEXPONENTIALATMOSPHERE_H

// Include statements.
#include "exponentialAtmosphere.h"
#include <cmath>
#include "basicMathematicsFunctions.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test of implementation of exponential atmosphere class.
/*!
 * Test of implementation of exponential atmosphere class.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testExponentialAtmosphere( );

}

#endif // UNITTESTEXPONENTIALATMOSPHERE_H

// End of file.
