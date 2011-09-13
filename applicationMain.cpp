/*! \file applicationMain.cpp
 *    Source file for a main application executable that runs all the
 *    application.
 *
 *    Path              : /
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Date created      : 1 July, 2011
 *    Last modified     : 1 July, 2011
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
 *      110701    K. Kumar          File created.
 */

// Include statements.
#include "exampleEarthOrbitingSatellite.h"

//! Execute all unit tests.
/*!
 * Executes all unit tests.
 */
int main( )
{
    // exampleEarthOrbitingSatellite: Runs an example of an Earth-orbiting
    // satellite. Outputs the final state to console, and generates the
    // propagation history to file.
    applications::executeEarthOrbitingSatelliteExample( );

    // Return integer to indicate completion of main executable file.
    return 0;
}

// End of file.
