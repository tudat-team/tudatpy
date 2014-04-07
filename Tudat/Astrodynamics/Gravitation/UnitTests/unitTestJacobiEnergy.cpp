/*    Copyright (c) 2010-2014, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110901    L. van der Ham    File created.
 *      120221    L. van der Ham    Fixed bug and up to date with latest standards.
 *      120227    K. Kumar          Boostified unit test; renamed file.
 *      120307    K. Kumar          Moved file.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 *    Notes
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
        namespace crtbp = gravitation::circular_restricted_three_body_problem;

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
        namespace crtbp = gravitation::circular_restricted_three_body_problem;

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
