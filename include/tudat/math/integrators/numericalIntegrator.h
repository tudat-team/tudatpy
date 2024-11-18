/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_NUMERICAL_INTEGRATOR_H
#define TUDAT_NUMERICAL_INTEGRATOR_H

#include <iostream>
#include <limits>

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/utilityMacros.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace numerical_integrators
{

//! Base class for the numerical integrators.
/*!
 * Base class for numerical integrators.
 * \tparam StateType The type of the state. This type should support addition,
 *          subtraction and assignment operators.
 * \tparam StateDerivativeType The type of the state derivative. This type should support
 *          multiplication with IndependentVariableType and doubles.
 * \tparam IndependentVariableType The type of the independent variable.
 */
template< typename IndependentVariableType = double, typename StateType = Eigen::VectorXd,
          typename StateDerivativeType = StateType, typename TimeStepType = IndependentVariableType >
class NumericalIntegrator
{
public:

    //! Typedef to the state derivative function.
    /*!
     * Typedef to the state derivative function. This should be a pointer to a function or a boost
     * function.
     */
    typedef std::function< StateDerivativeType(
            const IndependentVariableType, const StateType& ) > StateDerivativeFunction;

    //! Default constructor.
    /*!
     * Default constructor, taking a state derivative function as argument.
     * \param stateDerivativeFunction State derivative function.
     */
    NumericalIntegrator( const StateDerivativeFunction& stateDerivativeFunction ) :
        stateDerivativeFunction_( stateDerivativeFunction ),
        propagationTerminationFunction_( [ = ]( const double, const double, const Eigen::MatrixXd& ){ return false; } )
    { }

    //! Default virtual destructor.
    /*!
     * Default virtual destructor.
     */
    virtual ~NumericalIntegrator( ) { }

    //! Get step size of the next step.
    /*!
     * Returns the step size of the next step. Derived classes should override this and provide the
     * last step size that was computed or passed to performIntegrationStep( ).
     * \return Step size to be used for the next step.
     */
    virtual TimeStepType getNextStepSize( ) const = 0;

    //! Get current state.
    /*!
     * Returns the current state of the integrator. Derived classes should override this and
     * provide the computed state by performIntegrationStep( ).
     * \return Current integrated state.
     */
    virtual StateType getCurrentState( ) const = 0;

    //! Get current independent variable.
    /*!
     * Returns the current value of the independent variable of the integrator. Derived classes
     * should override this and provide the computed independent variable by
     * performIntegrationStep().
     * \return Current independent variable.
     */
    virtual IndependentVariableType getCurrentIndependentVariable( ) const = 0;

    //! Rollback internal state to the step performed by performIntegrationStep()
    /*!
     * Performs rollback of internal state to the step performed by performIntegrationStep(). This
     * is not necessarily equal the start of the integration interval after integrateTo() has been.
     * called. This function can only be called once after calling or performIntegrationStep()
     * unless specified otherwise by implementations, and can not be called before integrateTo()
     * has been called.
     * Will return true if the rollback was successful, and false otherwise.
     * \return True if the rollback was successful.
     */
    virtual bool rollbackToPreviousState( ) = 0;

    //! Get previous independent variable.
    /*!
     * Returns the previoius value of the independent variable of the integrator. Derived classes
     * should override this and provide the computed independent variable. If not implemented, throws error.
     * \return Previous independent variable.
     */
    virtual IndependentVariableType getPreviousIndependentVariable( )
    {
        throw std::runtime_error( "Function getPreviousIndependentVariable not implemented in this integrator" );
    }

    //! Get previous state value.
    /*!
     * Returns the previous value of the state. Derived classes
     * should override this and provide the computed state. If not implemented, throws error.
     * \return Previous state
     */
    virtual StateType getPreviousState( )
    {
        throw std::runtime_error( "Function getPreviousState not implemented in this integrator" );
    }

    //! Perform an integration to a specified independent variable value.
    /*!
     * Performs an integration to independentVariableEnd with initial state and initial independent
     * variable value specified by the current state of the integrator and the current independent
     * variable value. This implementation of integrateTo chooses the final step size such that it
     * exactly coincides with the given independentVariableEnd.
     * \param intervalEnd The value of the independent variable at the end of the interval to
     *          integrate over.
     * \param initialStepSize The initial step size to use.
     * \param finalTimeTolerance Tolerance to within which the final time should be reached.
     * \return The state at independentVariableEnd.
     */
    virtual StateType integrateTo(
            const IndependentVariableType intervalEnd,
            const TimeStepType initialStepSize,
            const TimeStepType finalTimeTolerance = std::numeric_limits< TimeStepType >::epsilon( ) );

    //! Perform a single integration step.
    /*!
     * Performs a single integration step from current independent variable and state as specified
     * by initialStateHistory and specified stepSize. This function should determine the next step
     * and make it available to getNextStepSize(), return the current state and store it for
     * getCurrentState(), and store the current independent variable value for
     * getCurrentIndependentVariable().
     * \param stepSize The step size of this step.
     * \return The state at the end of the interval.
     */
    virtual StateType performIntegrationStep( const TimeStepType stepSize ) = 0;

    //! Function to return the function that computes and returns the state derivative
    /*!
     * Function to return the function that computes and returns the state derivative
     * \return Function that returns the state derivative
     */
    StateDerivativeFunction getStateDerivativeFunction( )
    {
        return stateDerivativeFunction_;
    }

    //! Function to return the termination condition was reached during the current step
    /*!
     *  Function to return the termination condition was reached during the current step
     * \return  True if the termination condition was reached during the current step
     */
    bool getPropagationTerminationConditionReached( )
    {
        return propagationTerminationConditionReachedDuringStep_;
    }

    //! Setter for the (optional) propagation termination function
    /*!
     *  Setter for the (optional) propagation termination function, to be evaluated during the intermediate state updates
     *  performed to compute the quantities necessary to integrate the state to a new epoch (e.g. k1-k4 for RK$).
     *  \param terminationFunction Function that returns true if termination condition is reached, false if it has not,
     *  as a function of current time.
     */
    void setPropagationTerminationFunction( std::function< bool( const double, const double, const Eigen::MatrixXd& ) > terminationFunction )
    {
        propagationTerminationFunction_ = terminationFunction;
    }

    //! Function to toggle the use of step-size control
    /*!
     * Function to toggle the use of step-size control. To be implemented in derived classes with variable step sizes
     * \param useStepSizeControl Boolean denoting whether step size control is to be used
     */
    virtual void setStepSizeControl( const bool useStepSizeControl ) { }

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
    virtual void modifyCurrentState( const StateType& newState, const bool allowRollback = false )
    {
        TUDAT_UNUSED_PARAMETER( newState );
        TUDAT_UNUSED_PARAMETER( allowRollback );

        throw std::runtime_error(
                    "Error in numerical integrator. The function to modify the current state has not been implemented" );
    }

    //! Modify the state and time for the current step.
    /*!
     * Modify the state and time for the current step.
     * \param newState The new state to set the current state to.
     * \param newTime The time to set the current time to.
     * \param allowRollback Boolean denoting whether roll-back should be allowed.
     */
    virtual void modifyCurrentIntegrationVariables( const StateType& newState, const IndependentVariableType newTime,
                                                    const bool allowRollback = false )
    {
        TUDAT_UNUSED_PARAMETER( newState );
        TUDAT_UNUSED_PARAMETER( newTime );
        TUDAT_UNUSED_PARAMETER( allowRollback );
        throw std::runtime_error( "Error in numerical integrator. The function to modify the current integration variables has not "
                                  "been implemented in this integrator." );
    }

protected:

    //! Function that returns the state derivative.
    /*!
     * Function that returns the state derivative, as passed to the constructor.
     */
    StateDerivativeFunction stateDerivativeFunction_;

    //! Boolean to denote whether the propagation termination condition was reached during the evaluation of one of the sub-steps
    /*!
     *  Boolean to denote whether the propagation termination condition was reached during the evaluation of one of the sub-steps
     *  necessary to perform the last integration step. Parameter is false by default, and when set to true must be accompanied by
     *  propagationTerminationFunction_ (which is non-active by default)
     */
    bool propagationTerminationConditionReachedDuringStep_ = false;

    //! propagation termination function
    /*!
     *  propagation termination function to be evaluated during the intermediate state updates performed to compute
     *  the quantities necessary to integrate the state to a new epoch.
     *  By default, this function evaluates always to false, so the propagation termination conditions will not be
     *  checked during the integration subteps.
     */
    std::function< bool( const double, const double, const Eigen::MatrixXd& ) > propagationTerminationFunction_;
};


//extern template class NumericalIntegrator < double, Eigen::VectorXd, Eigen::VectorXd >;
//extern template class NumericalIntegrator < double, Eigen::Vector6d, Eigen::Vector6d >;
//extern template class NumericalIntegrator < double, Eigen::MatrixXd, Eigen::MatrixXd >;

//! Perform an integration to a specified independent variable value.
template< typename IndependentVariableType, typename StateType, typename StateDerivativeType, typename TimeStepType >
StateType NumericalIntegrator< IndependentVariableType, StateType, StateDerivativeType, TimeStepType >::integrateTo(
        const IndependentVariableType intervalEnd,
        const TimeStepType initialStepSize,
        const TimeStepType finalTimeTolerance )
{
    TimeStepType stepSize = initialStepSize;

    // Flag to indicate that the integration end value of the independent variable has been
    // reached.
    bool atIntegrationIntervalEnd = static_cast< TimeStepType >( intervalEnd - getCurrentIndependentVariable( ) )
            * stepSize / std::fabs( stepSize )
            <= finalTimeTolerance;

    int loopCounter = 0;
    while ( !atIntegrationIntervalEnd )
    {
//        std::cout<<"Steps "<<stepSize<<" "<<
        // Check if the remaining interval is smaller than the step size.
        if ( std::fabs( static_cast< TimeStepType >( intervalEnd - getCurrentIndependentVariable( ) ) )
             <= std::fabs( stepSize ) *
             ( 1.0 + finalTimeTolerance ) )
        {
            // The next step is beyond the end of the integration interval, so adjust the
            // step size accordingly.
            stepSize = intervalEnd - getCurrentIndependentVariable( );

            // Explicitly flag that the integration interval end is reached. Due to rounding
            // off errors, it may not be possible to use
            // ( currentIndependentVariable >= independentVariableEnd ) // in the while condition.
            atIntegrationIntervalEnd = true;
        }

        // Perform the step.
        performIntegrationStep( stepSize );
        stepSize = getNextStepSize( );

        // Only applicable to adaptive step size methods:
        // Perform additional step(s) to reach intervalEnd exactly, in case the last step was rejected by
        // the variable step size routine and a new step was used that was too small to get to intervalEnd.
        if ( atIntegrationIntervalEnd )
        {
            // As long as intervalEnd is not reached, perform additional steps with the remaining time
            // as suggested step size for the variable step size routine.
            if( std::fabs( static_cast< TimeStepType >( intervalEnd - getCurrentIndependentVariable( ) ) ) >
                    finalTimeTolerance )
            {
                // Ensure that integrateTo function does not get stuck in a loop.
                if( loopCounter < 1000 )
                {
                    atIntegrationIntervalEnd = false;
                    loopCounter++;
                }
                else
                {
                    std::cerr << "Warning, integrateTo function has failed to converge to final time to within tolerances, difference between true and requested final time is " <<
                                 intervalEnd - getCurrentIndependentVariable( ) << ", final time is: " <<
                                 getCurrentIndependentVariable( ) << std::endl;
                }
            }
        }
    }

    return getCurrentState( );
}

//! Typedef for shared-pointer to default numerical integrator.
/*!
 * Typedef for shared-pointer to a default numerical integrator (IndependentVariableType = double,
 * StateType = Eigen::VectorXd, StateDerivativeType = Eigen::VectorXd).
 */
typedef std::shared_ptr< NumericalIntegrator< > > NumericalIntegratorXdPointer;

//! Typedef for a shared-pointer to a scalar numerical integrator.
/*!
 * Typedef for shared-pointer to a scalar numerical integrator (IndependentVariableType = double,
 * StateType = double, StateDerivativeType = double).
 */
typedef std::shared_ptr< NumericalIntegrator< double, double, double > > NumericalIntegratordPointer;

} // namespace numerical_integrators

} // namespace tudat

#endif // TUDAT_NUMERICAL_INTEGRATOR_H
