/*! \file unitTestLibrationPoints.h
 *    Header file of unit test file of libration point code. This unit test
 *    file will test the determintation of the locations of the libration
 *    points in the Circular Restricted Three-Body Problem (CRTBP). Computation
 *    of the mass parameter is also tested.
 *
 *    Path              : /Astrodynamics/MissionSegments/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : L. van der Ham
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.vanderHam@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 1 June, 2011
 *    Last modified     : 10 July, 2011
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
 *      110601    L. van der Ham    First creation of code.
 *      110710    K. Kumar          Restructured code; added subtests.
 */

#ifndef UNITTESTLIBRATIONPOINTS_H
#define UNITTESTLIBRATIONPOINTS_H

#include <cmath>
#include <iostream>
#include "basicMathematicsFunctions.h"
#include "celestialBody.h"
#include "librationPoint.h"
#include "planet.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test determination of libration point locations.
/*!
 * Tests determination of libration point locations. Also tests computation of
 * dimensionless mass parameter. Computations are done in the Earth-Moon
 * system.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testLibrationPointLocations( );

//! Test determination of L1 location.
/*!
 * Tests determination of L1 location. Computation is done in the Earth-Moon
 * system.
 * \param isLibrationPointComputationErroneous Flag to indicate if test is erroneous.
 * \param massParameter Mass parameter.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testL1LibrationPointLocation( bool isLibrationPointComputationErroneous,
                                   double massParameter );

//! Test determination of L2 location.
/*!
 * Tests determination of L2 location. Computation is done in the Earth-Moon
 * system.
 * \param isLibrationPointComputationErroneous Flag to indicate if test is erroneous.
 * \param massParameter Mass parameter.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testL2LibrationPointLocation( bool isLibrationPointComputationErroneous,
                                   double massParameter );


//! Test determination of L3 location.
/*!
 * Tests determination of L3 location. Computation is done in the Earth-Moon
 * system.
 * \param isLibrationPointComputationErroneous Flag to indicate if test is erroneous.
 * \param massParameter Mass parameter.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testL3LibrationPointLocation( bool isLibrationPointComputationErroneous,
                                   double massParameter );

//! Test determination of L4 location.
/*!
 * Tests determination of L4 location. Computation is done in the Earth-Moon
 * system.
 * \param isLibrationPointComputationErroneous Flag to indicate if test is erroneous.
 * \param massParameter Mass parameter.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testL4LibrationPointLocation( bool isLibrationPointComputationErroneous,
                                   double massParameter );

//! Test determination of L5 location.
/*!
 * Tests determination of L5 location. Computation is done in the Earth-Moon
 * system.
 * \param isLibrationPointComputationErroneous Flag to indicate if test is erroneous.
 * \param massParameter Mass parameter.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testL5LibrationPointLocation( bool isLibrationPointComputationErroneous,
                                   double massParameter );

}

#endif // UNITTESTLIBRATIONPOINTS_H

// End of file.
