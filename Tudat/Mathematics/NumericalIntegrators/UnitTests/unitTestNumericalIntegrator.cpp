/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      120127    B. Tong Minh      File created.
 *      120128    D. Dirkx          Minor changes during code check.
 *      120206    K. Kumar          Minor comment corrections.
 *      120207    K. Kumar          Updated to Boost unit test framework.
 *      120213    K. Kumar          Modified getCurrentInterval() to getIndependentVariable().
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <iostream>
#include <limits>

#include <boost/bind.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/utilityMacros.h"
#include "Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h"

namespace tudat
{
namespace unit_tests
{

//! Test suit for numerical integrators.
BOOST_AUTO_TEST_SUITE( test_numerical_integrator )

//! Using declaration of the NumericalIntegrator.
using numerical_integrators::NumericalIntegrator;

//! State derivative function that always returns zero.
/*!
 * State derivative function that always returns a zero vector with length equal to the input.
 * \param time Time at which the state derivative needs to be evaluated.
 * \param state State at which the state derivative needs to be evalated.
 * \return Zero vector with length equal to state.
 */
Eigen::VectorXd computeZeroStateDerivative( const double time,
                                            const Eigen::VectorXd& state )
{
    TUDAT_UNUSED_PARAMETER( time );
    return Eigen::VectorXd::Zero( state.rows( ) );
}

//! Dummy numerical integrator.
/*!
 * Dummy numerical integrator that keeps track of the amount of steps made.
 */
class DummyNumericalIntegrator :
        public NumericalIntegrator< double, Eigen::VectorXd, Eigen::VectorXd >
{
public:

    //! Default constructor.
    /*!
     * Default constructor, setting zeroDerivative as the state derivative function.
     */
    DummyNumericalIntegrator( const double intervalStart, const Eigen::VectorXd& initialState )
        : NumericalIntegrator< >( &computeZeroStateDerivative ),
          numberOfSteps( 0 ),
          currentIndependentVariable_( intervalStart ),
          currentState_( initialState )
    { }

    //! Get step size of the next step.
    /*!
     * Returns the step size of the next step.
     * \return Step size to be used for the next step.
     */
    virtual double getNextStepSize( ) const { return stepSize_; }

    //! Get current state.
    /*!
     * Returns the current state stored interally by the integrator.
     * \return Current state.
     */
    virtual Eigen::VectorXd getCurrentState( ) const { return currentState_; }

    //! Get current independent variable.
    /*!
     * Returns the current value of the independent variable stored interally by the integrator.
     * \return Current independent variable.
     */
    virtual double getCurrentIndependentVariable( ) const { return currentIndependentVariable_; }

    //! Rollback the internal state to the last state.
    /*!
     * Performs rollback of the internal state to the last state. This function can only be called
     * once after calling integrateTo( ) or performIntegrationStep( ) unless specified otherwise by
     * implementations, and can not be called before any of these functions have been called. Will
     * return true if the rollback was succesful, and false otherwise.
     * \return True if the rollback was successful.
     */
    virtual bool rollbackToPreviousState( ) { return true; }

    //! Perform a single integration step that does not do anything.
    /*!
     * Performs a single integration step that does not do anything, except incrementing the
     * step counter.
     * \param stepSize The step size of this step.
     * \return The state history at the end of the interval, which is basically equal to the input.
     */
    virtual Eigen::VectorXd performIntegrationStep( const double stepSize )
    {
        // Increment the step counter.
        numberOfSteps++;

        stepSize_ = stepSize;
        currentIndependentVariable_ += stepSize_;

        return currentState_;
    }

    //! Counter of the number of steps.
    /*!
     * Counter of the number of steps taken.
     */
    int numberOfSteps;

protected:

    //! Last used step size by this integrator.
    /*!
     * Last used step size by this integrator.
     */
    double stepSize_;

    //! Current independent variable value.
    /*!
     * Current value of independent variable stored internally by the integrator.
     */
    double currentIndependentVariable_;

    //! Current state.
    /*!
     * Current state stored internally by the integrator.
     */
    Eigen::VectorXd currentState_;
};

//! Test the amount of steps that NumericalIntegrator::integrateTo takes.
/*!
 * Test the amount of steps that NumericalIntegrator::integrateTo takes.
 * \param intervalStart The start of the integration interval.
 * \param intervalEnd The end of the integration interval.
 * \param initialState The initial state.
 * \param stepSize The step size to take.
 * \return True if actual number of steps is equal to actual number of steps; false otherwise.
 */
bool testIntegrateToFunction( const double intervalStart, const double intervalEnd,
                              const Eigen::VectorXd& initialState,
                              const double stepSize )
{
    DummyNumericalIntegrator integrator( intervalStart, initialState );
    Eigen::VectorXd integratedState = integrator.integrateTo( intervalEnd, stepSize );

    // Calculate expected number of steps.
    int expectedNumberOfSteps = static_cast< int >(
                std::ceil( ( intervalEnd - intervalStart ) / stepSize ) );

    // Check if the integrated state is equal to the initial state.
    // This is an exact comparison, because no arithmetics are performed on the state.
    if ( integratedState != initialState )
    {
        std::cerr << "DummyNumericalIntegrator was not a dummy integrator!" << std::endl;
        return false;
    }

    // Check if the actual number of steps is equal to the expected number of steps.
    if ( integrator.numberOfSteps != expectedNumberOfSteps )
    {
        return false;
    }

    return true;
}

//! Test if the numerical integrator base class is working correctly.
BOOST_AUTO_TEST_CASE( testNumberOfStepsUsingNumericalIntegrator )
{
    // Set random initial state.
    Eigen::VectorXd initialState = ( Eigen::VectorXd( 4 ) << 0.34, 0.24, 0.76, 0.10 ).finished( );

    // Case 1: test number of steps for different start and end times.
    {
        BOOST_CHECK( testIntegrateToFunction( 0.0, 0.0, initialState, 10.0 ) );
        BOOST_CHECK( testIntegrateToFunction( 0.0, 10.0, initialState, 10.0 ) );
        BOOST_CHECK( testIntegrateToFunction( 0.0, 20.0, initialState, 10.0 ) );
        BOOST_CHECK( testIntegrateToFunction( 0.0, 30.0, initialState, 10.0 ) );
    }

    // Case 2: test number of steps for different step sizes.
    {
        BOOST_CHECK( testIntegrateToFunction( 0.0, 10.0, initialState, 2.5 ) );
        BOOST_CHECK( testIntegrateToFunction( 0.0, 10.0, initialState, 3.0 ) );
        BOOST_CHECK( testIntegrateToFunction( 0.0, 10.0, initialState, 3.5 ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
