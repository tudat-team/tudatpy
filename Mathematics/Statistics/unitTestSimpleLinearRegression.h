/*! \file unitTestSimpleLinearRegression.h
 *    Header file that defines the unitTestSimpleLinearRegression unit test,
 *    which tests the simple linear regression method implemented in Tudat.
 *
 *    Path              : /Mathematics/Statistics/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 2 July, 2011
 *    Last modified     : 5 September, 2011
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
 *      110702    K. Kumar          First creation of code.
 *      110726    K. Kumar          Changed filename and class name.
 *      110802    K .Kumar          Changed filename
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

#ifndef UNITTESTSIMPLELINEARREGRESSION_H
#define UNITTESTSIMPLELINEARREGRESSION_H

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test implementation of simple linear regression method.
/*!
 * Tests implementation of simple linear regression method.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testSimpleLinearRegression( );

}

#endif // UNITTESTSIMPLELINEARREGRESSION_H

// End of file.
