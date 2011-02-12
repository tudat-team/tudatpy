/*! \file unitTestAerodynamicsNamespace.h
 *    This file contains the definition of the aerodynamics namespace unit
 *    test.
 *
 *    Path              : /Astrodynamics/ForceModels/Aerothermodynamics/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : Dominic Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : L. Abdulkadir
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.Abdulkadir@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 8 November, 2010
 *    Last modified     : 11 February, 2011
 *
 *    References
 *        Anderson Jr., J.D. , Fundamentals of Aerodynamics, 3rd edition,
 *            McGraw Hill, 2001.
 *        Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd
 *            edition, AIAA Education Series, 2006.
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
 *      110208    D. Dirkx          First version of file.
 *      110210    L. Abdulkadir     Code check.
 *      110211    K. Kumar          Corrected file header; corrected Doxygen
 *                                  comments
 */

#ifndef UNITTESTAERODYNAMICDNAMESPACE_H
#define UNITTESTAERODYNAMICDNAMESPACE_H

// Include statements.
#include <iostream>
#include "aerodynamics.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test of aerodynamics namespace.
/*!
 * Test of aerodynamics namespace.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testAerodynamicsNameSpace( );

}

#endif // UNITTESTAERODYNAMICDNAMESPACE_H

// End of file.
