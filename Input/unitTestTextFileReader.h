/*! \file unitTestTextFileReader.h
 *    This header file contains the definition of a unit test for the text file
 *    reader class.
 *
 *    Path              : /Input/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : J. Leloux
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.Leloux@tudelft.nl, J.Leloux@student.tudelft.nl,
 *                        jleloux@gmail.com
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 24 March, 2011
 *    Last modified     : 29 March, 2011
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
 *      110324    J. Leloux         First creation of code.
 *      110329    K. Kumar          Minor corrections.
 */

#ifndef UNITTESTTEXTFILEREADER_H
#define UNITTESTTEXTFILEREADER_H

// Include statements.
#include <cmath>
#include "textFileReader.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test implementation of textFileReader class.
/*!
 * Tests implementation of textFileReader class.
 * A test input file has been written and is located in the Input folder.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testTextFileReader( );

}

#endif // UNITTESTTEXTFILEREADER_H

// End of file.
