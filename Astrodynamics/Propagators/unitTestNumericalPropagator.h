/*! \file unitTestNumericalPropagator.h
 *    Header file that defines a unit test that tests the numerical propagator
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
 *    Last modified     : 16 February, 2011
 *
 *    References
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *      The unit test at present only checks that the code is internally
 *      consistent and doesn't check the result against benchmark data.
 *      Validation of the code against benchmark data should be added to
 *      ensure that the output of both simulation cases tested is correct.
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

#ifndef UNITTESTNUMERICALPROPAGATOR_H
#define UNITTESTNUMERICALPROPAGATOR_H

// Include statements.
#include <iostream>
#include <map>
#include "cartesianElements.h"
#include "cartesianPositionElements.h"
#include "cartesianStateNumericalPropagator.h"
#include "celestialBody.h"
#include "gravitationalForceModel.h"
#include "numericalPropagator.h"
#include "planet.h"
#include "rungeKutta4thOrderFixedStepsize.h"
#include "sphericalHarmonicsGravityField.h"
#include "state.h"
#include "unitConversions.h"
#include "vehicle.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test of implementation of numerical propagator class.
/*!
 * Test of implementation of numerical propagator class.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testNumericalPropagator( );

}

#endif // UNITTESTNUMERICALPROPAGATOR_H

// End of file.
