/*! \file unitTestCubicSplineInterpolation.h
 *    Header file that defines the unit test for the cubic spline interplation
 *    included in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : e.a.g.heeren@student.tudelft.nl
 *
 *    Date created      : 21 June, 2011
 *    Last modified     : 14 July, 2011
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
 *      110621    F.M. Engelen      File created.
 *      110707    E.A.G. Heeren     Minor spelling/lay-out corrections.
 *      110714    E.A.G. Heeren     Minor spelling/lay-out corrections.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

#ifndef UNITTESTCUBICSPLINEINTERPOLATION_H
#define UNITTESTCUBICSPLINEINTERPOLATION_H

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test implementation of cubic spline class.
/*!
 * Tests implementation of cubic spline class.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testCubicSplineInterpolation( );

}

#endif // UNITTESTCUBICSPLINEINTERPOLATION_H

// End of file.
