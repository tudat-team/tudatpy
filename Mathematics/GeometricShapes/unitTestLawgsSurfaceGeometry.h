/*! \file unitTestLawgsSurfaceGeometry.h
 *    This file contains the definition of the LawgsPartGeometry unit test.
 *
 *    Path              : /Mathematics/GeometricShapes/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 6 February, 2011
 *    Last modified     : 8 February, 2011
 *
 *    References
 *      Craidon, C.B. A Desription of the Langley Wireframe Geometry Standard
 *          (LaWGS) format, NASA TECHNICAL MEMORANDUM 85767.
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
 *      YYMMDD    Author           Comment
 *      110206    D. Dirkx         First version of file.
 *      110208    D. Dirkx         Finalized for code check.
 */

#ifndef UNITTESTLAWGSSURFACEGEOMETRY_H
#define UNITTESTLAWGSSURFACEGEOMETRY_H

// Include statements.
#include "lawgsPartGeometry.h"
#include "sphereSegment.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test of implementation of Lawgs surface geometry class.
/*!
 * Test of implementation of Lawgs surface geometry class.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testLawgsSurfaceGeometry( );

}

#endif // UNITTESTLAWGSSURFACEGEOMETRY_H

// End of file.
