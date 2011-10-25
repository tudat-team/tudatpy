/*! \file unitTestDeepSpaceManeuver.h
 *    Header file for unit test that tests the Deep Space Maneuver (DSM) class implemented in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/MissionSegments/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 27 May, 2011
 *    Last modified     : 27 May, 2011
 *
 *    References
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
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
 *      110527    K. Kumar        First creation of the code.
 */

#ifndef UNITTESTDEEPSPACEMANEUVER_H
#define UNITTESTDEEPSPACEMANEUVER_H

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test deep space maneuver.
/*!
 * Tests deep space maneuver.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testDeepSpaceManeuver( );

}

#endif // UNITTESTDEEPSPACEMANEUVER_H

// End of file.
