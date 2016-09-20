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
 *      120207    K. Kumar          Minor comment corrections.
 *      120213    K. Kumar          Updated getCurrentInterval() to getIndependentVariable().
 *      120424    K. Kumar          Added missing this-pointer, to satisfy requirement for
 *                                  accessing base class members.
 *      121128    K. Kumar          Changed base class to ReinitializableNumericalIntegrator, added
 *                                  implementation of modifyCurrentState() function.
 *      121205    D. Dirkx          Migrated namespace to directory-based protocol and added
 *                                  backwards compatibility; added standardized typedefs.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_RUNGE_KUTTA_4_INTEGRATOR_H
#define TUDAT_RUNGE_KUTTA_4_INTEGRATOR_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/NumericalIntegrators/reinitializableNumericalIntegrator.h"

namespace tudat
{
namespace numerical_integrators
{

//! Class that implements the Runge-Kutta 4 integrator.
/*!
 * Class that implements the Runge-Kutta 4, fixed order, fixed step size integrator.
 * \tparam StateType The type of the state. This type should support addition with
 *          StateDerivativeType.
 * \tparam StateDerivativeType The type of the state derivative. This type should support
 *          multiplication with IndependentVariableType and doubles.
 * \tparam IndependentVariableType The type of the independent variable.
 * \sa NumericalIntegrator.
 */
template < typename IndependentVariableType = double, typename StateType = Eigen::VectorXd,
           typename StateDerivativeType = Eigen::VectorXd >
class RungeKutta4Integrator
        : public numerical_integrators::ReinitializableNumericalIntegrator<
        IndependentVariableType, StateType, StateDerivativeType >
{
public:

    //! Typedef for the base class.
    /*!
     * Typedef of the base class with all template parameters filled in.
     */
    typedef numerical_integrators::ReinitializableNumericalIntegrator<
    IndependentVariableType, StateType,
    StateDerivativeType > ReinitializableNumericalIntegratorBase;

    //! Typedef for the state derivative function.
    /*!
     * Typedef to the state derivative function inherited from the base class.
     * \sa NumericalIntegrator::StateDerivativeFunction.
     */
    typedef typename ReinitializableNumericalIntegratorBase::NumericalIntegratorBase::
    StateDerivativeFunction StateDerivativeFunction;

    //! Default constructor.
    /*!
     * Default constructor, taking a state derivative function as argument.
     * \param stateDerivativeFunction State derivative function.
     * \param intervalStart The start of the integration interval.
     * \param initialState The initial state.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    RungeKutta4Integrator( const StateDerivativeFunction& stateDerivativeFunction,
                           const IndependentVariableType intervalStart,
                           const StateType& initialState )
        : ReinitializableNumericalIntegratorBase( stateDerivativeFunction ),
          currentIndependentVariable_( intervalStart ),
          currentState_( initialState ),
          lastIndependentVariable_( intervalStart )
    { }

    //! Get step size of the next step.
    /*!
     * Returns the step size of the next step.
     * \return Step size to be used for the next step.
     */
    virtual IndependentVariableType getNextStepSize( ) const { return stepSize_; }

    //! Get current state.
    /*!
     * Returns the current state of the integrator.
     * \return Current integrated state,
     */
    virtual StateType getCurrentState( ) const { return currentState_; }

    //! Returns the current independent variable.
    /*!
     * Returns the current value of the independent variable of the integrator.
     * \return Current independent variable.
     */
    virtual IndependentVariableType getCurrentIndependentVariable( ) const
    {
        return currentIndependentVariable_;
    }

    //! Perform a single integration step.
    /*!
     * Perform a single integration step.
     * \param stepSize The step size to take.
     * \return The state at the end of the interval,
     */
    virtual StateType performIntegrationStep( const IndependentVariableType stepSize )
    {
        lastIndependentVariable_ = currentIndependentVariable_;
        lastState_ = currentState_;

        // Calculate k1-k4.
        const StateDerivativeType k1 = stepSize * this->stateDerivativeFunction_(
                    currentIndependentVariable_, currentState_ );

        const StateDerivativeType k2 = stepSize * this->stateDerivativeFunction_(
                    currentIndependentVariable_ + stepSize / 2.0,
                    static_cast< StateType >( currentState_ + k1 / 2.0 ) );

        const StateDerivativeType k3 = stepSize * this->stateDerivativeFunction_(
                    currentIndependentVariable_ + stepSize / 2.0,
                    static_cast< StateType >( currentState_ + k2 / 2.0 ) );

        const StateDerivativeType k4 = stepSize * this->stateDerivativeFunction_(
                    currentIndependentVariable_ + stepSize,
                    static_cast< StateType >( currentState_ + k3 ) );

        stepSize_ = stepSize;
        currentIndependentVariable_ += stepSize_;
        currentState_ += ( k1 + 2.0 * k2 + 2.0 * k3 + k4 ) / 6.0;

        // Return the integration result.
        return currentState_;
    }

    //! Rollback internal state to the last state.
    /*!
     * Performs rollback of internal state to the last state. This function can only be called once
     * after calling integrateTo() or performIntegrationStep() unless specified otherwise by
     * implementations, and can not be called before any of these functions have been called. Will
     * return true if the rollback was succesful, and false otherwise.
     * \return True if the rollback was successful.
     */
    virtual bool rollbackToPreviousState( )
    {
        if ( currentIndependentVariable_ == lastIndependentVariable_ )
        {
            return false;
        }

        currentIndependentVariable_ = lastIndependentVariable_;
        currentState_ = lastState_;
        return true;
    }

    //! Modify the state at the current value of the independent variable.
    /*!
     * Modify the state at the current value of the independent variable.
     * \param newState The new state to set the current state to.
     */
    void modifyCurrentState( const StateType& newState )
    {
        this->currentState_ = newState;
        this->lastIndependentVariable_ = currentIndependentVariable_;
    }

protected:

    //! Last used step size.
    /*!
     * Last used step size, passed to either integrateTo() or performIntegrationStep().
     */
    IndependentVariableType stepSize_;

    //! Current independent variable.
    /*!
     * Current independent variable as computed by performIntegrationStep().
     */
    IndependentVariableType currentIndependentVariable_;

    //! Current state.
    /*!
     * Current state as computed by performIntegrationStep().
     */
    StateType currentState_;

    //! Last independent variable.
    /*!
     * Last independent variable value as computed by performIntegrationStep().
     */
    IndependentVariableType lastIndependentVariable_;

    //! Last state.
    /*!
     * Last state as computed by performIntegrationStep().
     */
    StateType lastState_;
};

//! Typedef of RK4 integrator (state/state derivative = VectorXd, independent variable = double).
/*!
 * Typedef of a RK4 integrator with VectorXds as state and state derivative and double as
 * independent variable.
 */
typedef RungeKutta4Integrator< > RungeKutta4IntegratorXd;

//! Typedef of a scalar RK4 integrator.
/*!
 * Typedef of an RK4 integrator with doubles as state and state derivative and independent
 * variable.
 */
typedef RungeKutta4Integrator< double, double, double > RungeKutta4Integratord;

//! Typedef of pointer to default RK4 integrator
/*!
 * Typedef of pointer to a RK4 integrator with VectorXds as state and state derivative and double
 * as independent variable.
 */
typedef boost::shared_ptr< RungeKutta4IntegratorXd > RungeKutta4IntegratorXdPointer;

//! Typedef of pointer to a scalar RK4 integrator.
/*!
 * Typedef of pointer to an RK4 integrator with doubles as state and state derivative and
 * independent variable.
 */
typedef boost::shared_ptr< RungeKutta4Integratord > RungeKutta4IntegratordPointer;

} // namespace integrators
} // namespace tudat

#endif // TUDAT_RUNGE_KUTTA_4_INTEGRATOR_H
