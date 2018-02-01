/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Gravitation/jacobiEnergy.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
namespace tudat
{
namespace unit_tests
{

using namespace gravitation;
using namespace orbital_element_conversions;

//! Test if Jacobi energy is computed correctly.
BOOST_AUTO_TEST_CASE( testJacobiEnergy )
{
    // Test 1: test Jacobi energy at L1.
    {

        // Set mass parameter for Earth-moon system. Value from Table 3.1 (Wakker, 2007).
        double massParameter = 0.01215;

        // Initialize position L1, from Table 3.4 (Wakker, 2007).
        Eigen::VectorXd stateAtL1 = Eigen::VectorXd::Zero( 6 );
        stateAtL1( xCartesianPositionIndex ) = 0.836914;

        // Set expected value of Jacobi energy at L1.
        double expectedJacobiEnergy = 3.1883;

        // Compute Jacobi energy.
        double computedJacobiEnergy = computeJacobiEnergy( massParameter, stateAtL1 );

        // Check if expected Jacobi energy matches computed.
        BOOST_CHECK_CLOSE_FRACTION( expectedJacobiEnergy,  computedJacobiEnergy, 1.0e-4 );
    }

    // Test 2: test Jacobi energy at L4.
    {

        // Set mass parameter for Earth-moon system. Value from Table 3.1 (Wakker, 2007).
        double massParameter = 0.01215;

        // Initialize position L4, from Table 3.4 (Wakker, 2007).
        Eigen::VectorXd stateAtL4 = Eigen::VectorXd::Zero( 6 );
        stateAtL4( xCartesianPositionIndex ) = 0.487849;
        stateAtL4( yCartesianPositionIndex ) = 0.866025;

        // Set expected value of Jacobi energy at L4.
        double expectedJacobiEnergy = 2.9880;

        // Compute Jacobi energy.
        double computedJacobiEnergy = computeJacobiEnergy( massParameter, stateAtL4 );

        // Check if expected Jacobi energy matches computed.
        BOOST_CHECK_CLOSE_FRACTION( expectedJacobiEnergy,  computedJacobiEnergy, 1.0e-6 );
    }
}

} // namespace unit_tests
} // namespace tudat
