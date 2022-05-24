/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */
#ifdef NDEBUG
#ifdef TUDAT_BUILD_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif // TUDAT_BUILD_GNU
#endif // NDEBUG

#ifndef TUDAT_RUNGE_KUTTA_FIXED_STEP_INTEGRATOR_H
#define TUDAT_RUNGE_KUTTA_FIXED_STEP_INTEGRATOR_H

#include <memory>

#include <Eigen/Core>

#include "tudat/math/integrators/reinitializableNumericalIntegrator.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"

namespace tudat
{

namespace numerical_integrators
{

//! Class that implements the fixed step Runge Kutta integrator.
/*!
 * Class that implements the fixed order, fixed step size, Runge Kutta integrator.
 * \tparam StateType The type of the state. This type should support addition with
 *          StateDerivativeType.
 * \tparam StateDerivativeType The type of the state derivative. This type should support
 *          multiplication with IndependentVariableType and doubles.
 * \tparam IndependentVariableType The type of the independent variable.
 * \sa NumericalIntegrator.
 */
template< typename IndependentVariableType = double,
          typename StateType = Eigen::VectorXd,
          typename StateDerivativeType = Eigen::VectorXd,
          typename TimeStepType = IndependentVariableType >
class RungeKuttaFixedStepSizeIntegrator
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
     * \param coefficientsSet The fixed step coefficient type to use for integration.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    RungeKuttaFixedStepSizeIntegrator( const StateDerivativeFunction& stateDerivativeFunction,
                                       const IndependentVariableType intervalStart,
                                       const StateType& initialState,
                                       const CoefficientSets& coefficientsSet,
                                       const RungeKuttaCoefficients::OrderEstimateToIntegrate orderToUse = RungeKuttaCoefficients::OrderEstimateToIntegrate::lower ) :
        ReinitializableNumericalIntegratorBase( stateDerivativeFunction ),
        currentIndependentVariable_( intervalStart ),
        currentState_( initialState ),
        lastIndependentVariable_( intervalStart ),
        coefficientsSet_( coefficientsSet ),
        orderToUse_( orderToUse )
    {
        // Load the Butcher tableau coefficients.
        setCoefficients( coefficientsSet );
    }

    // Function to load the Butcher tableau coefficients.
    /*!
     * Function to load the Butcher tableau coefficients.
     * \param coefficientsSet The fixed step coefficient type to use for integration.
     */
    void setCoefficients( const CoefficientSets& coefficientsSet )
    {
        // Load the Butcher tableau coefficients.
        butcherTableau_ = butcherTableau_.get( coefficientsSet );

        // If the coefficient set is not already fixed step size, remove the b column according to the order to use.
        if ( !(butcherTableau_.isFixedStepSize) )
        {
            if (orderToUse_ == RungeKuttaCoefficients::OrderEstimateToIntegrate::lower)
            {
                // Keep the first row of bCoefficients.
                Eigen::MatrixXd oldBCoefficients = butcherTableau_.bCoefficients;
                butcherTableau_.bCoefficients = Eigen::MatrixXd::Zero( 1, oldBCoefficients.cols() );
                butcherTableau_.bCoefficients.row( 0 ) = oldBCoefficients.row( 0 );
            }
            else
            {
                // Keep the second row of bCoefficients.
                Eigen::MatrixXd oldBCoefficients = butcherTableau_.bCoefficients;
                butcherTableau_.bCoefficients = Eigen::MatrixXd::Zero( 1, oldBCoefficients.cols() );
                butcherTableau_.bCoefficients.row( 0 ) = oldBCoefficients.row( 1 );
            }
        }

        currentScaledStateDerivatives_.clear( );
        currentScaledStateDerivatives_.resize( this->butcherTableau_.cCoefficients.rows( ) );

    }

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
        // Save the current state and independent variable.
        lastIndependentVariable_ = currentIndependentVariable_;
        lastState_ = currentState_;

        StateType stateUpdate = StateType::Zero( currentState_.rows( ), currentState_.cols( ) );

        // Compute the k_i state derivatives per stage.
        for ( int stage = 0; stage < this->butcherTableau_.cCoefficients.rows( ); stage++ )
        {
            // Compute the intermediate state.
            stateUpdate.setZero( );
            for ( int column = 0; column < stage; column++ )
            {
                if( this->butcherTableau_.aCoefficients( stage, column ) !=
                        mathematical_constants::getFloatingInteger< typename StateType::Scalar >( 0 ) )
                {
                    stateUpdate += this->butcherTableau_.aCoefficients( stage, column ) *
                            currentScaledStateDerivatives_[ column ];
                }
            }

            // Compute the intermediate state to pass to the state derivative for this stage.
            StateType intermediateState = this->currentState_ + stateUpdate;

            // Compute the state derivative.
            const IndependentVariableType time = this->currentIndependentVariable_ +
                    this->butcherTableau_.cCoefficients( stage ) * stepSize;
            currentScaledStateDerivatives_[ stage ] = stepSize * this->stateDerivativeFunction_( time, intermediateState );

            // Check if propagation should terminate because the propagation termination condition has been reached
            // while computing the intermediate state.
            // If so, return immediately the current state (not recomputed yet), which will be discarded.
            if ( this->propagationTerminationFunction_( static_cast< double >( time ), TUDAT_NAN ) )
            {
                this->propagationTerminationConditionReachedDuringStep_ = true;
                return this->currentState_;
            }
        }

        stateUpdate.setZero( );
        for ( int stage = 0; stage < this->butcherTableau_.cCoefficients.rows( ); stage++ )
        {
            if( this->butcherTableau_.bCoefficients( 0, stage ) !=
                    mathematical_constants::getFloatingInteger< typename StateType::Scalar >( 0 ) )
            {
                // Update the estimate.
                stateUpdate += this->butcherTableau_.bCoefficients( 0, stage ) * currentScaledStateDerivatives_[ stage ];
            }
        }

        // Keep step size unchanged.
        stepSize_ = stepSize;
        // Update the current state and independent variable.
        this->currentIndependentVariable_ += stepSize;
        this->currentState_ += stateUpdate;

        // Return the integration result.
        return currentState_;
    }



    //! Rollback internal state to the last state.
    /*!
     * Performs rollback of internal state to the last state. This function can only be called once
     * after calling integrateTo() or performIntegrationStep() unless specified otherwise by
     * implementations, and can not be called before any of these functions have been called. Will
     * return true if the rollback was successful, and false otherwise.
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
     * Returns the previous value of the independent variable of the integrator.
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

    //! Get the Butcher tableau used.
    /*!
     * Returns the Butcher tableau used by the integrator.
     * \return Butcher tableau used.
     */
    RungeKuttaCoefficients getButcherTableau( )
    {
        return butcherTableau_;
    }

    //! Replace the state with a new value.
    /*!
     * Replace the state with a new value. This allows for discrete jumps in the state, often
     * used in simulations of discrete events. In astro, this relates to simulations of rocket staging,
     * impulsive shots, parachuting, ideal control, etc. The modified state, by default, cannot be rolled back; to do this, either
     * set the flag to true, or store the state before calling this function the first time, and call it again with the initial state
     * as parameter to revert to the state before the discrete change.
     * \param newState The value of the new state.
     * \param allowRollback Boolean denoting whether roll-back should be allowed.
     */
    void modifyCurrentState( const StateType& newState, const bool allowRollback = false )
    {
        currentState_ = newState;
        if ( !allowRollback )
        {
            this->lastIndependentVariable_ = currentIndependentVariable_;
        }
    }

    //! Modify the state and time for the current step.
    /*!
     * Modify the state and time for the current step.
     * \param newState The new state to set the current state to.
     * \param newTime The time to set the current time to.
     * \param allowRollback Boolean denoting whether roll-back should be allowed.
     */
    void modifyCurrentIntegrationVariables( const StateType& newState, const IndependentVariableType newTime,
                                            const bool allowRollback = false )
    {
        currentState_ = newState;
        currentIndependentVariable_ = newTime;
        if ( !allowRollback )
        {
            this->lastIndependentVariable_ = currentIndependentVariable_;
        }
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
     * Current state as computed by performIntegrationStep( ).
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

    //! Type of the coefficient set used for integration.
    CoefficientSets coefficientsSet_;

    //! Butcher tableau used for the integration.
    RungeKuttaCoefficients butcherTableau_;

    //! Vector of state derivatives.
    /*!
     * Vector of state derivatives, i.e. values of k_{i} in Runge-Kutta scheme.
     */
    std::vector< StateDerivativeType > currentScaledStateDerivatives_;

    // Order of Runge-Kutta method to be used.
    RungeKuttaCoefficients::OrderEstimateToIntegrate orderToUse_;
};

extern template class RungeKuttaFixedStepSizeIntegrator < double, Eigen::VectorXd, Eigen::VectorXd >;
extern template class RungeKuttaFixedStepSizeIntegrator < double, Eigen::Vector6d, Eigen::Vector6d >;
extern template class RungeKuttaFixedStepSizeIntegrator < double, Eigen::MatrixXd, Eigen::MatrixXd >;


//! Typedef of RK fixed-step integrator (state/state derivative = VectorXd, independent variable = double).
/*!
 * Typedef of a RK fixed-step integrator with VectorXds as state and state derivative and double as
 * independent variable.
 */
typedef RungeKuttaFixedStepSizeIntegrator< > RungeKuttaFixedStepSizeIntegratorXd;

//! Typedef of a scalar RK fixed-step integrator.
/*!
 * Typedef of a RK fixed-step integrator with doubles as state and state derivative and independent
 * variable.
 */
typedef RungeKuttaFixedStepSizeIntegrator< double, double, double > RungeKuttaFixedStepSizeIntegratord;

//! Typedef of pointer to default RK fixed-step integrator
/*!
 * Typedef of pointer to a RK fixed-step integrator with VectorXds as state and state derivative and double
 * as independent variable.
 */
typedef std::shared_ptr< RungeKuttaFixedStepSizeIntegratorXd > RungeKuttaFixedStepSizeIntegratorXdPointer;

//! Typedef of pointer to a scalar RK fixed-step integrator.
/*!
 * Typedef of pointer to an RK fixed-step integrator with doubles as state and state derivative and
 * independent variable.
 */
typedef std::shared_ptr< RungeKuttaFixedStepSizeIntegratord > RungeKuttaFixedStepSizeIntegratordPointer;

} // namespace numerical_integrators

} // namespace tudat

#endif // TUDAT_RUNGE_KUTTA_FIXED_STEP_INTEGRATOR_H

#ifdef NDEBUG
#ifdef TUDAT_BUILD_GNU
// turn the warnings back on
#pragma GCC diagnostic pop
#endif // TUDAT_BUILD_GNU
#endif // NDEBUG
