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
           typename StateDerivativeType = Eigen::VectorXd, typename TimeStepType = IndependentVariableType >
class RungeKutta4Integrator
        : public numerical_integrators::ReinitializableNumericalIntegrator<
        IndependentVariableType, StateType, StateDerivativeType, TimeStepType >
{
public:

    //! Typedef for the base class.
    /*!
     * Typedef of the base class with all template parameters filled in.
     */
    typedef numerical_integrators::ReinitializableNumericalIntegrator<
    IndependentVariableType, StateType, StateDerivativeType, TimeStepType > ReinitializableNumericalIntegratorBase;

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
    virtual TimeStepType getNextStepSize( ) const { return stepSize_; }

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
    virtual StateType performIntegrationStep( const TimeStepType stepSize )
    {
        lastIndependentVariable_ = currentIndependentVariable_;
        lastState_ = currentState_;

        // Calculate k1-k4.
        StateDerivativeType k1, k2, k3, k4;
        for ( unsigned int i = 1; i <= 4; i++ )
        {
            IndependentVariableType time;
            StateType state;
            switch ( i ) {
            case 1:
                time = currentIndependentVariable_;
                state = currentState_;
                k1 = stepSize * this->stateDerivativeFunction_( time, state );
                break;
            case 2:
                time = currentIndependentVariable_ + stepSize / 2.0;
                state = static_cast< StateType >( currentState_ + k1 / 2.0 );
                k2 = stepSize * this->stateDerivativeFunction_( time, state );
                break;
            case 3:
                time = currentIndependentVariable_ + stepSize / 2.0;
                state = static_cast< StateType >( currentState_ + k2 / 2.0 );
                k3 = stepSize * this->stateDerivativeFunction_( time, state );
                break;
            case 4:
                time = currentIndependentVariable_ + stepSize;
                state = static_cast< StateType >( currentState_ + k3 );
                k4 = stepSize * this->stateDerivativeFunction_( time, state );
                break;
            }

            // Check if propagation should terminate because the propagation termination condition has been reached
            // while computing k1, k2, k3 or k4. If so, return immediately the current state (not recomputed yet),
            // which will be discarded.
            if ( this->propagationTerminationFunction_( static_cast< double >( time ), TUDAT_NAN ) )
            {
                this->propagationTerminationConditionReachedDuringStep_ = true;
                return currentState_;
            }
        }

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

    //! Get previous independent variable.
    /*!
     * Returns the previoius value of the independent variable of the integrator.
     * \return Previous independent variable.
     */
    IndependentVariableType getPreviousIndependentVariable( )
    {
        return lastIndependentVariable_;
    }

    //! Get previous state value.
    /*!
     * Returns the previous value of the state.
     * \return Previous state
     */
    StateType getPreviousState( )
    {
        return lastState_;
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
    TimeStepType stepSize_;

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

} // namespace numerical_integrators
} // namespace tudat

#endif // TUDAT_RUNGE_KUTTA_4_INTEGRATOR_H
