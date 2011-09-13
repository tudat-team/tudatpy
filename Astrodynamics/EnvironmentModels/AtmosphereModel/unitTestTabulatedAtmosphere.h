/*! \file unitTestTabulatedAtmosphere.h
 *    Header file that defines the Tabulated atmosphere unit test included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/AtmosphereModel/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 13 July, 2011
 *    Last modified     : 13 July, 2011
 *
 *    References
 *
 *    Notes
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
 *      110713    F.M. Engelen      File created.
 */

#ifndef UNITTESTTABULATEDATMOSPHERE_H
#define UNITTESTTABULATEDATMOSPHERE_H

// Include statements.
#include <cmath>
#include "basicMathematicsFunctions.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test of implementation of Tabulated atmosphere class.
/*!
 * Test of implementation of Tabulated atmosphere class.
 * \return Boolean indicating success of test.
 * ( false = successful; true = failed ).
 */
bool testTabulatedAtmosphere( );

}

#endif // UNITTESTTABULATEDATMOSPHERE_H

// End of file.
