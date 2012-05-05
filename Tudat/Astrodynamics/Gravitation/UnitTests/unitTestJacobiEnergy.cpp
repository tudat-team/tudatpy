/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      110901    L. van der Ham    First creation of code.
 *      120221    L. van der Ham    Fixed bug and up to date with latest standards.
 *      120227    K. Kumar          Boostified unit test; renamed file.
 *      120307    K. Kumar          Moved file.
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

namespace tudat
{
namespace unit_tests
{

//! Test if Jacobi energy is computed correctly.
BOOST_AUTO_TEST_CASE( testJacobiEnergy )
{
    // Test 1: test Jacobi energy at L1.
    {
        namespace crtbp = astrodynamics::gravitation::circular_restricted_three_body_problem;

        // Set mass parameter for Earth-moon system. Value from Table 3.1 (Wakker, 2007).
        double massParameter = 0.01215;

        // Initialize position L1, from Table 3.4 (Wakker, 2007).
        Eigen::VectorXd stateAtL1 = Eigen::VectorXd::Zero( 6 );
        stateAtL1( crtbp::xPositionIndex ) = 0.836914;

        // Set expected value of Jacobi energy at L1.
        double expectedJacobiEnergy = 3.1883;

        // Compute Jacobi energy.
        double computedJacobiEnergy = crtbp::computeJacobiEnergy( massParameter, stateAtL1 );

        // Check if expected Jacobi energy matches computed.
        BOOST_CHECK_CLOSE_FRACTION( expectedJacobiEnergy,  computedJacobiEnergy, 1.0e-4 );
    }

    // Test 2: test Jacobi energy at L4.
    {
        namespace crtbp = astrodynamics::gravitation::circular_restricted_three_body_problem;

        // Set mass parameter for Earth-moon system. Value from Table 3.1 (Wakker, 2007).
        double massParameter = 0.01215;

        // Initialize position L4, from Table 3.4 (Wakker, 2007).
        Eigen::VectorXd stateAtL4 = Eigen::VectorXd::Zero( 6 );
        stateAtL4( crtbp::xPositionIndex ) = 0.487849;
        stateAtL4( crtbp::yPositionIndex ) = 0.866025;

        // Set expected value of Jacobi energy at L4.
        double expectedJacobiEnergy = 2.9880;

        // Compute Jacobi energy.
        double computedJacobiEnergy = crtbp::computeJacobiEnergy( massParameter, stateAtL4 );

        // Check if expected Jacobi energy matches computed.
        BOOST_CHECK_CLOSE_FRACTION( expectedJacobiEnergy,  computedJacobiEnergy, 1.0e-6 );
    }
}

} // namespace unit_tests
} // namespace tudat
