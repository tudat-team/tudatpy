/*! \file unitTestApproximatePlanetPositions.h
 *    Header file for a unit test that tests the implementation of the
 *    ApproximatePlanetPositions class in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 5 April, 2011
 *    Last modified     : 5 April, 2011
 *
 *    References
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
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
 *      110405    K. Kumar          File created.
 */

#ifndef UNITTESTAPPROXIMATEPLANETPOSITIONS_H
#define UNITTESTAPPROXIMATEPLANETPOSITIONS_H

// Include statements.
#include "basicMathematicsFunctions.h"
#include "celestialBody.h"
#include "keplerianElements.h"
#include "linearAlgebra.h"
#include "orbitalElementConversions.h"
#include "predefinedCelestialBodies.h"
#include "unitConversions.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test ApproximatePlanetPositions class.
/*!
 * Tests ApproximatePlanetPositions class.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testApproximatePlanetPositions( );

}

#endif // UNITTESTAPPROXIMATEPLANETPOSITIONS_H

// End of file.
