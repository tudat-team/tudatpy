/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      120319    K. Kumar          Creation of code.
 *
 *    References
 *        Wakker, K.F. "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *        Howell, K.C. Three-dimensional, periodic, 'Halo' orbits, Celestial Mechanics, 32. 53-71,
 *          1984.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/bind.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>
#include <TudatCore/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

#include "Tudat/Astrodynamics/Gravitation/stateDerivativeCircularRestrictedThreeBodyProblem.h"

namespace tudat
{
namespace unit_tests
{

//! Test if state derivative for circular restricted three-body problem is computed correctly.
BOOST_AUTO_TEST_CASE( testStateDerivativeCircularRestrictedThreeBodyProblem )
{
    namespace crtbp = gravitation::circular_restricted_three_body_problem;
    using crtbp::StateDerivativeCircularRestrictedThreeBodyProblem;
    using crtbp::xPositionIndex;
    using crtbp::zPositionIndex;
    using crtbp::yVelocityIndex;

    // Test 1: test state derivative at L1.
    {
        // Set mass parameter for Earth-moon system. Value from Table 3.1 (Wakker, 2007).
        double massParameter = 0.01215;

        // Initialize position L1, from Table 3.4 (Wakker, 2007).
        Eigen::VectorXd stateAtL1 = Eigen::VectorXd::Zero( 6 );
        stateAtL1( xPositionIndex ) = 0.836914;

        // Declare state derivative object.
        StateDerivativeCircularRestrictedThreeBodyProblem stateDerivative( massParameter );

        // Compute state derivative at L1.
        Eigen::VectorXd stateDerivativeAtL1 = stateDerivative.computeStateDerivative(
                    0.0, stateAtL1 );

        // Check if expected state derivative (all zeros) matches computed.
        TUDAT_CHECK_MATRIX_BASE( Eigen::VectorXd::Zero( 6 ), stateDerivativeAtL1 )
            BOOST_CHECK_SMALL( stateDerivativeAtL1.coeff( row, col ), 1.0e-4 );
    }

    // Test 2: test Halo orbit around L1. Note that this is actually a system test, since it
    //         requires that the Runge-Kutta 4 numerical integrator is functioning correctly.
    {
        // Set mass parameter. Value from Table I (Howell, 1984).
        double massParameter = 0.04;

        // Initialize position on Halo orbit around L1 (Howell, 1984).
        Eigen::VectorXd initialStateOnHaloOrbit = Eigen::VectorXd::Zero( 6 );
        initialStateOnHaloOrbit( xPositionIndex ) = 0.723268;
        initialStateOnHaloOrbit( zPositionIndex ) = 0.04;
        initialStateOnHaloOrbit( yVelocityIndex ) = 0.198019;

        // Declare state derivative object.
        StateDerivativeCircularRestrictedThreeBodyProblem stateDerivative( massParameter );

        // Declare Runge-Kutta 4 integrator.
        numerical_integrators::RungeKutta4IntegratorXd rungeKutta4Integrator(
                    boost::bind(
                        &StateDerivativeCircularRestrictedThreeBodyProblem::computeStateDerivative,
                        &stateDerivative, _1, _2 ),
                    0.0, initialStateOnHaloOrbit );

        // Integrate Halo orbit over one period.
        Eigen::VectorXd finalStateOnHaloOrbit = rungeKutta4Integrator.integrateTo(
                    2.0 * 1.300177, 0.1 );

        // Check if expected final state (same as initial) matches computed.
        BOOST_CHECK_CLOSE_FRACTION( initialStateOnHaloOrbit( 0 ),
                                    finalStateOnHaloOrbit( 0 ), 1.0e-3 );
        BOOST_CHECK_SMALL( finalStateOnHaloOrbit( 1 ), 1.0e-3 );
        BOOST_CHECK_CLOSE_FRACTION( initialStateOnHaloOrbit( 2 ),
                                    finalStateOnHaloOrbit( 2 ), 1.0e-3 );
        BOOST_CHECK_SMALL( finalStateOnHaloOrbit( 3 ), 1.0e-2 );
        BOOST_CHECK_CLOSE_FRACTION( initialStateOnHaloOrbit( 4 ),
                                    finalStateOnHaloOrbit( 4 ), 1.0e-2 );
        BOOST_CHECK_SMALL( finalStateOnHaloOrbit( 5 ), 1.0e-3 );
    }
}

} // namespace unit_tests
} // namespace tudat
