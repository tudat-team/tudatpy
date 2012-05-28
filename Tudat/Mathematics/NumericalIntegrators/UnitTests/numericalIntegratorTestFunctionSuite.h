/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      120328    K. Kumar          File created.
 *
 *    References
 *
 */

#ifndef TUDAT_NUMERICAL_INTEGRATOR_TEST_FUNCTION_SUITE_H
#define TUDAT_NUMERICAL_INTEGRATOR_TEST_FUNCTION_SUITE_H

#include <cmath>

#include <Eigen/Core>

namespace tudat
{
namespace unit_tests
{
namespace numerical_integrator_test_function_suite
{

//! Compute state derivative that always returns zero.
/*!
 * Computes state derivative function that always returns a zero vector with length equal to the
 * input.
 * \param time Time at which the state derivative needs to be evaluated.
 * \param state State at which the state derivative needs to be evaluated.
 * \return Zero vector with length equal to state.
 */
Eigen::VectorXd computeZeroStateDerivative( const double time, const Eigen::VectorXd& state )
{
    return Eigen::VectorXd::Zero( state.rows( ) );
}

//! Compute analytical state for zero-state derivative.
/*!
 * Computes analytical state for zero-state derivative function. This function always returns the
 * initial state with length equal to input.
 * \param time Time at which the state needs to be evaluated.
 * \param initialState Initial state for analytical solution.
 * \return Constant state equal to initial state.
 */
Eigen::VectorXd computeAnalyticalStateForZeroStateDerivative( const double time,
                                                              const Eigen::VectorXd& initialState )
{
    return initialState;
}

//! Compute van der Pol oscillator state derivative.
/*!
 * Computes the van der Pol state derivative function, returning a state of dimension 2.
 * The van der Pol oscillator model is described by the following set of ordinary differential
 * equations:
 *
 * \f{eqnarray*}{
 *      \frac{ dx }{ dt }( 0 ) &=& x( 1 ) \\
 *      \frac{ dx }{ dt }( 1 ) &=& \mu * ( 1 - x( 0 )^{ 2 } ) * x( 1 ) + x( 0 ) \\
 * \f}
 *
 * The scaling parameter \f$\mu\f$ , which sets the damping and nonlinearity of the model is set to
 * 1.0.
 *
 * \param time Time at which the state derivative needs to be evaluated.
 * \param state State at which the state derivative needs to be evaluated.
 * \return Computed state derivative.
 */
Eigen::VectorXd computeVanDerPolStateDerivative( const double time, const Eigen::VectorXd& state )
{
    // Declare van der Pol state derivative vector.
    Eigen::VectorXd vanDerPolStateDerivative( state.rows( ) );

    // Compute state derivative.
    vanDerPolStateDerivative( 0 ) = state( 1 );
    vanDerPolStateDerivative( 1 ) = ( 1.0 - std::pow( state( 0 ), 2.0 ) ) * state( 1 )
            + state( 0 );

    // Return computed state derivative.
    return vanDerPolStateDerivative;
}

////! Compute analytical state for van der Pol-state derivative.
///*!
// * Computes analytical state for zero-state derivative function. This function always returns the
// * initial state with length equal to input.
// * \param time Time at which the state needs to be evaluated.
// * \param initialState Initial state for analytical solution.
// * \return Constant state equal to initial state.
// */
//std::make_pair< double, double, Eigen::VectorXd, Eigen::VectorXd >
//getAnalyticalSolutionForZeroStateDerivative( )
//{
//    return std::make_pair( 0.0, 1.0, Eigen::VectorXd)
//}

//! Compute non-autonomous model state derivative.
/*!
 * Computes the state derivative for a non-autonomous ordinary differential equation. The state
 * derivative function defined corresponds to Example 3, pg. 278 in (Burden and Faires, 2001).
 * The initial-value problem is:
 *
 * \f[
 *      y' = y - t^{ 2 } + 1
 * \f]
 *
 * \param time Time at which the state derivative needs to be evaluated.
 * \param state State, length of which should be 1, at which the state derivative needs to
 *          be evaluated.
 * \return State derivative, length equal to state values according to above expression.
 */
Eigen::VectorXd computeNonAutonomousModelStateDerivative( const double time,
                                                          const Eigen::VectorXd& state )
{
    const double stateDerivative = state( 0 ) - std::pow( time, 2.0 ) + 1.0;
    return Eigen::VectorXd::Constant( 1, stateDerivative );
}

} // numerical_integrator_test_function_suite
} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_NUMERICAL_INTEGRATOR_TEST_FUNCTION_SUITE_H
