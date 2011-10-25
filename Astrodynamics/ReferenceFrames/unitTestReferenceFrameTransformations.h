/*! \file unitTestReferenceFrameTransformations.h
 *    This file contains the definition of the reference frame transformations unit test included
 *    in Tudat.
 *
 *    Path              : /Astrodynamics/ReferenceFrames/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Checker           : F. M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 25 May, 2011
 *    Last modified     : 1 July, 2011
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
 *      110525    F.M. Engelen      First creation of code.
 *      110628    K. Kumar          Minor comment and layout modifications.
 *      110701    K. Kumar          Updated file path.
 */

#ifndef UNITTESTFRAMETRANSFORMATION_H
#define UNITTESTFRAMETRANSFORMATION_H

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test reference frame transformations.
/*!
 * Tests reference frame transformations.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testReferenceFrameTransformations( );

}

#endif // UNITTESTFRAMETRANSFORMATION_H

// End of file.
