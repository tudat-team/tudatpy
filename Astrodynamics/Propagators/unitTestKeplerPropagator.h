/*! \file unitTestKeplerPropagator.h
 *    Header file that defines a unit test that tests the Kepler propagator
 *    included in Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 16 February, 2011
 *    Last modified     : 16 Feburary, 2011
 *
 *    References
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *      Currently, this file makes use of benchmark data provided by J. Melman.
 *      In future, it is desirable that the benchmark data is the direct output
 *      of a commercial package such as STK, where are initial conditions of
 *      the simulation are known.
 *
 *      The basic_functions::ROOT_PATH variable does not seem to work
 *      correctly. Currently, this unit test will only work if the executable
 *      is built in the Tudat root directory. A generic method to handle paths
 *      within Tudat must be added in future and must be implemented here.
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
 *      110216    K. Kumar          First creation of code.
 */

#ifndef UNITTESTKEPLERPROPAGATOR_H
#define UNITTESTKEPLERPROPAGATOR_H

// Include statements.
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include "basicFunctions.h"
#include "basicMathematicsFunctions.h"
#include "cartesianElements.h"
#include "celestialBody.h"
#include "keplerPropagator.h"
#include "newtonRaphson.h"
#include "predefinedPlanets.h"
#include "state.h"
#include "unitConversions.h"
#include "vehicle.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test of implementation of Kepler propagator class.
/*!
 * Test of implementation of Kepler propagator class.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testKeplerPropagator( );

}

#endif // UNITTESTKEPLERPROPAGATOR_H

// End of file.
