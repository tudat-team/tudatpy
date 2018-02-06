/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_NUMERICAL_INTEGRATOR_TEST_FUNCTIONS_H
#define TUDAT_NUMERICAL_INTEGRATOR_TEST_FUNCTIONS_H

#include <cmath>

#include <Eigen/Core>

namespace tudat
{
namespace unit_tests
{
namespace numerical_integrator_test_functions
{

//! Compute state derivative that always returns zero.
/*!
 * Computes state derivative function that always returns a zero vector with length equal to the
 * input.
 * \param time Time at which the state derivative needs to be evaluated.
 * \param state State at which the state derivative needs to be evaluated.
 * \return Zero vector with length equal to state.
 */
static inline Eigen::VectorXd computeZeroStateDerivative( const double time,
                                                          const Eigen::VectorXd& state )
{
    return Eigen::VectorXd::Zero( state.rows( ) );
}

static inline Eigen::VectorXd computeConstantStateDerivative( const double time,
                                                              const Eigen::VectorXd& state )
{
    Eigen::VectorXd stateDerivative = Eigen::VectorXd::Zero( state.rows( ) );
    stateDerivative( 0 ) = 1.0;
    //std::cout<<"State der: "<<time<<" "<<stateDerivative.transpose( )<<" "<<std::endl<<state.transpose( )<<std::endl;
    return stateDerivative;
}

//! Compute analytical state for zero-state derivative.
/*!
 * Computes analytical state for zero-state derivative function. This function always returns the
 * initial state with length equal to input.
 * \param time Time at which the state needs to be evaluated.
 * \param initialState Initial state for analytical solution.
 * \return Constant state equal to initial state.
 */
static inline Eigen::VectorXd computeAnalyticalStateForZeroStateDerivative(
        const double time, const Eigen::VectorXd& initialState )
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
static inline Eigen::VectorXd computeVanDerPolStateDerivative( const double time,
                                                               const Eigen::VectorXd& state )
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

//! Computes state derivative of a Logarithmic 2-dimensional ODE
/*!
 * Computes the state derivative of a 2-dimensional ODE that is autonomous and includes log functions.
 * This function is used by Fehlberg for testing the RK5(6) integrator and has an analytical solution.
 *
 * \f{eqnarray*}{
 *      \frac{ dx }{ dt }( 0 ) &=& -2 * t * x( 0 ) * log( x(1) ) \\
 *      \frac{ dx }{ dt }( 0 ) &=&  2 * t * x( 1 ) * log( x(0) ) \\
 * \f}
 *
 * Reference:   Fehlberg, E. (1968). Classical Fifth-, Sixth-, Seventh- and Eigth-Order Runge-Kutta
 *              Formulas with Stepsize Control (page 30)
 *
 * \param time Time at which the state derivative needs to be evaluated.
 * \param state State at which the state derivative needs to be evaluated.
 * \return Computed state derivative.
 */
static inline Eigen::VectorXd
    computeFehlbergLogirithmicTestODEStateDerivative( const double time,
                                                      const Eigen::VectorXd& state){
    // Declare state derivative vector
    Eigen::VectorXd stateDerivative( state.rows() );

    // Compute state derivative
    stateDerivative( 0 ) = -2.0 * time * state(0) * std::log( state(1) );
    stateDerivative( 1 ) =  2.0 * time * state(1) * std::log( state(0) );

    // Return computed state derivative.
    return stateDerivative;
}

//! Computes the analytical solution of the state for the Logarithmic test ODE of Fehlberg
/*!
 * Analytical solution:
 * \f{eqnarray*}{
 *      x(0) = exp( cos( t^{2} ) ) \\
 *      x(1) = exp( sin( t^{2} ) ) \\
 * \f}
 *
 * with initial conditions:
 * \f{eqnarray*}{
 *      x_{0}(0) = exp(1) \\
 *      x_{0}(1) = 1
 * \f}
 *
 * Reference: Fehlberg, E. (1968). Classical Fifth-, Sixth-, Seventh- and Eigth-Order Runge-Kutta
 *            Formulas with Stepsize Control (page 30)
 *
 * \param time Time at which the state derivative needs to be evaluated.
 * \param initialState State at which the state derivative needs to be evaluated.
 * \return Computed state.
 */

static inline Eigen::VectorXd computeAnalyticalStateFehlbergODE( const double time,
                                                                 const Eigen::VectorXd& initialState){
    // Declare state derivative vector
    Eigen::VectorXd state( initialState.rows() );

    // Compute state derivative
    state(0) = exp( cos( pow(time,2.0) ) ) ; // x(0) = exp( cos( t^{2} ) )
    state(1) = exp( sin( pow(time,2.0) ) ) ; // x(1) = exp( sin( t^{2} ) )

    // Return computed state derivative.
    return state;
}

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
static inline Eigen::VectorXd computeNonAutonomousModelStateDerivative( const double time,
                                                                        const Eigen::VectorXd& state )
{
    const double stateDerivative = state( 0 ) - std::pow( time, 2.0 ) + 1.0;
    return Eigen::VectorXd::Constant( 1, stateDerivative );
}

//! Compute another non-autonomous model state derivative.
/*!
 * Computes the state derivative for a non-autonomous ordinary differential equation, recommended
 * by R. Hofstenge for testing the integrateTo function of the variable step size integrator
 *
 * \f[
 *      y' = y * ( 2 - t ) * t + t - 1
 * \f]
 *
 * \param time Time at which the state derivative needs to be evaluated.
 * \param state State, length of which should be 1, at which the state derivative needs to
 *          be evaluated.
 * \return State derivative, length equal to state values according to above expression.
 */
static inline Eigen::VectorXd computeSecondNonAutonomousModelStateDerivative(
        const double time, const Eigen::VectorXd& state )
{
    const double stateDerivative = state( 0 ) * ( 2.0 - time ) * time + time - 1.0;
    return Eigen::VectorXd::Constant( 1, stateDerivative );
}

} // numerical_integrator_test_functions
} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_NUMERICAL_INTEGRATOR_TEST_FUNCTIONS_H
