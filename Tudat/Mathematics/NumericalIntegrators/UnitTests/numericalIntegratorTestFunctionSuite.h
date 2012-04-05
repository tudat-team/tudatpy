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
