/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *      Montenbruck, O., Gill, E. Satellite Orbits: Models, Methods, Applications, Springer, 2005.
 *
 */

#ifndef TUDAT_RUNGE_KUTTA_VARIABLE_STEP_SIZE_INTEGRATOR_H
#define TUDAT_RUNGE_KUTTA_VARIABLE_STEP_SIZE_INTEGRATOR_H

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

#include <functional>
#include <memory>

#include <Eigen/Core>

#include <limits>
#include <vector>

#include "tudat/basics/utilityMacros.h"
#include "tudat/math/integrators/reinitializableNumericalIntegrator.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"

namespace tudat
{

namespace numerical_integrators
{

//! Class that implements the Runge-Kutta variable stepsize integrator.
/*!
 * Class that implements the Runge-Kutta variable step size integrator.
 * \tparam StateType The type of the state. This type should be an Eigen::Matrix derived type.
 * \tparam StateDerivativeType The type of the state derivative. This type should be an
 *          Eigen::Matrix derived type.
 * \tparam IndependentVariableType The type of the independent variable. This type should be
 *          either a float or double.
 * \sa NumericalIntegrator.
 */
template< typename IndependentVariableType = double, typename StateType = Eigen::VectorXd,
          typename StateDerivativeType = StateType, typename TimeStepType = IndependentVariableType  >
class RungeKuttaVariableStepSizeIntegrator :
        public ReinitializableNumericalIntegrator<
        IndependentVariableType, StateType, StateDerivativeType, TimeStepType >
{
public:

    //! Typedef to the function used to compute the new step size.
    /*!
     * Typedef to the function used to compute the new step size. This should be a pointer to a
     * function or a boost function.
     */
    typedef std::function< std::pair< TimeStepType, bool >(
            const TimeStepType, const std::pair< TimeStepType, TimeStepType >&, const TimeStepType,
            const std::pair< TimeStepType, TimeStepType >&, const StateType&,
            const StateType&, const StateType&, const StateType& ) > NewStepSizeFunction;

    //! Typedef of the base class.
    /*!
     * Typedef of the base class with all template parameters filled in.
     */
    typedef numerical_integrators::ReinitializableNumericalIntegrator<
    IndependentVariableType, StateType,
    StateDerivativeType, TimeStepType > ReinitializableNumericalIntegratorBase;

    //! Typedef to the state derivative function.
    /*!
     * Typedef to the state derivative function inherited from the base class.
     * \sa NumericalIntegrator::StateDerivativeFunction.
     */
    typedef typename ReinitializableNumericalIntegratorBase::NumericalIntegratorBase::
    StateDerivativeFunction StateDerivativeFunction;

    //! Exception that is thrown if the minimum step size is exceeded.
    /*!
     * Exception thrown by RungeKuttaVariableStepSizeIntegrator< >::
     * computeNextStepSizeAndValidateResult() if the minimum step size is exceeded.
     */
    class MinimumStepSizeExceededError;

    //! Default constructor.
    /*!
     * Default constructor, taking coefficients, a state derivative function, initial conditions,
     * minimum & maximum step size and relative & absolute error tolerance per item in the state
     * vector as argument.
     * \param coefficients Coefficients to use with this integrator.
     * \param stateDerivativeFunction State derivative function.
     * \param intervalStart The start of the integration interval.
     * \param initialState The initial state.
     * \param minimumStepSize The minimum step size to take. If this constraint is violated, an
     *          exception will be thrown and the next state is computed using the minimum value.
     * \param maximumStepSize The maximum step size to take. If this constraint is violated, an
     *          exception will be thrown and the next state is computed using the maximum value.
     * \param relativeErrorTolerance The relative error tolerance, for each individual state
     *          vector element.
     * \param absoluteErrorTolerance The absolute error tolerance, for each individual state
     *          vector element.
     * \param safetyFactorForNextStepSize Safety factor used to scale prediction of next step size.
     * \param maximumFactorIncreaseForNextStepSize Maximum factor increase for next step size.
     * \param minimumFactorDecreaseForNextStepSize Maximum factor decrease for next step size.
     * \param newStepSizeFunction Function that returns the new step size computed, default:
     *          RungeKuttaVariableStepSizeIntegrator::computeNewStepSize.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    RungeKuttaVariableStepSizeIntegrator(
            const RungeKuttaCoefficients& coefficients,
            const StateDerivativeFunction& stateDerivativeFunction,
            const IndependentVariableType intervalStart,
            const StateType& initialState,
            const IndependentVariableType minimumStepSize,
            const IndependentVariableType maximumStepSize,
            const StateType& relativeErrorTolerance,
            const StateType& absoluteErrorTolerance,
            const TimeStepType safetyFactorForNextStepSize = 0.8,
            const TimeStepType maximumFactorIncreaseForNextStepSize = 4.0,
            const TimeStepType minimumFactorDecreaseForNextStepSize = 0.1,
            const NewStepSizeFunction& newStepSizeFunction = 0,
            const bool exceptionIfMinimumStepExceeded = true ) :
        ReinitializableNumericalIntegratorBase( stateDerivativeFunction ),
        currentIndependentVariable_( intervalStart ),
        currentState_( initialState ),
        lastIndependentVariable_( intervalStart ),
        coefficients_( coefficients ),
        minimumStepSize_( std::fabs( static_cast< double >( minimumStepSize ) ) ),
        maximumStepSize_( std::fabs( static_cast< double >( maximumStepSize ) ) ),
        relativeErrorTolerance_( relativeErrorTolerance.array( ).abs( ) ),
        absoluteErrorTolerance_( absoluteErrorTolerance.array( ).abs( ) ),
        safetyFactorForNextStepSize_( std::fabs( static_cast< double >( safetyFactorForNextStepSize ) ) ),
        maximumFactorIncreaseForNextStepSize_( std::fabs( static_cast< double >( maximumFactorIncreaseForNextStepSize ) ) ),
        minimumFactorDecreaseForNextStepSize_( std::fabs( static_cast< double >( minimumFactorDecreaseForNextStepSize ) ) ),
        newStepSizeFunction_( newStepSizeFunction ), exceptionIfMinimumStepExceeded_( exceptionIfMinimumStepExceeded ),
        useStepSizeControl_( true )
    {
        if( !( currentState_.rows( ) == relativeErrorTolerance_.rows( ) ) ||
                !( currentState_.cols( ) == relativeErrorTolerance_.cols( ) ) )
        {
            throw std::runtime_error( "Error when creating variable step-size RK integrator, relative tolerance input size is inconsistent" );
        }

        if( !( currentState_.rows( ) == absoluteErrorTolerance_.rows( ) ) ||
                !( currentState_.cols( ) == absoluteErrorTolerance_.cols( ) ) )
        {
            throw std::runtime_error( "Error when creating variable step-size RK integrator, absolute tolerance input size is inconsistent" );
        }

        // Set default newStepSizeFunction_ to the class method.
        if ( this->newStepSizeFunction_ == 0 )
        {
            this->newStepSizeFunction_ = std::bind(
                        &RungeKuttaVariableStepSizeIntegrator::computeNewStepSize,
                        this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                        std::placeholders::_5, std::placeholders::_6, std::placeholders::_7, std::placeholders::_8 );
        }
    }

    //! Default constructor.
    /*!
     * Default constructor, taking coefficients a state derivative function, initial conditions,
     * minimum & maximum step size and relative & absolute error tolerance for all items in the
     * state vector as argument.
     * \param coefficients Coefficients to use with this integrator.
     * \param stateDerivativeFunction State derivative function.
     * \param intervalStart The start of the integration interval.
     * \param initialState The initial state.
     * \param minimumStepSize The minimum step size to take. If this constraint is violated, an
     *          exception will be thrown and the next state is computed using the minimum value.
     * \param maximumStepSize The maximum step size to take. If this constraint is violated, an
     *          exception will be thrown and the next state is computed using the maximum value.
     * \param relativeErrorTolerance The relative error tolerance, equal for all individual state
     *          vector elements.
     * \param absoluteErrorTolerance The absolute error tolerance, equal for all individual state
     *          vector elements.
     * \param safetyFactorForNextStepSize Safety factor used to scale prediction of next step size.
     * \param maximumFactorIncreaseForNextStepSize Maximum factor increase for next step size.
     * \param minimumFactorDecreaseForNextStepSize Maximum factor decrease for next step size.
     * \param newStepSizeFunction Function that returns the new step size computed, default:
     *          RungeKuttaVariableStepSizeIntegrator::computeNewStepSize.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    RungeKuttaVariableStepSizeIntegrator(
            const RungeKuttaCoefficients& coefficients,
            const StateDerivativeFunction& stateDerivativeFunction,
            const IndependentVariableType intervalStart,
            const StateType& initialState,
            const TimeStepType minimumStepSize,
            const TimeStepType maximumStepSize,
            const typename StateType::Scalar relativeErrorTolerance,
            const typename StateType::Scalar absoluteErrorTolerance,
            const TimeStepType safetyFactorForNextStepSize = 0.8,
            const TimeStepType maximumFactorIncreaseForNextStepSize = 4.0,
            const TimeStepType minimumFactorDecreaseForNextStepSize = 0.1,
            const NewStepSizeFunction& newStepSizeFunction = 0,
            const bool exceptionIfMinimumStepExceeded = true  ) :
        ReinitializableNumericalIntegratorBase( stateDerivativeFunction ),
        currentIndependentVariable_( intervalStart ),
        currentState_( initialState ),
        lastIndependentVariable_( intervalStart ),
        coefficients_( coefficients ),
        minimumStepSize_( std::fabs( static_cast< double >( minimumStepSize ) ) ),
        maximumStepSize_( std::fabs( static_cast< double >( maximumStepSize ) ) ),
        relativeErrorTolerance_( StateType::Constant( initialState.rows( ), initialState.cols( ),
                                                      std::fabs( relativeErrorTolerance ) ) ),
        absoluteErrorTolerance_( StateType::Constant( initialState.rows( ), initialState.cols( ),
                                                      std::fabs( absoluteErrorTolerance ) ) ),
        safetyFactorForNextStepSize_( std::fabs( static_cast< double >( safetyFactorForNextStepSize ) ) ),
        maximumFactorIncreaseForNextStepSize_( std::fabs( static_cast< double >( maximumFactorIncreaseForNextStepSize ) ) ),
        minimumFactorDecreaseForNextStepSize_( std::fabs( static_cast< double >( minimumFactorDecreaseForNextStepSize ) ) ),
        newStepSizeFunction_( newStepSizeFunction ),
        exceptionIfMinimumStepExceeded_( exceptionIfMinimumStepExceeded ),
        useStepSizeControl_( true )
    {
        // Set default newStepSizeFunction_ to the class method.
        if ( newStepSizeFunction_ == 0 )
        {
            this->newStepSizeFunction_ = std::bind(
                        &RungeKuttaVariableStepSizeIntegrator::computeNewStepSize,
                        this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                        std::placeholders::_5, std::placeholders::_6, std::placeholders::_7, std::placeholders::_8 );
        }

        // Raise error if a fixed step coefficient set is used with this variable step integrator.
        if( coefficients_.isFixedStepSize )
        {
            throw std::runtime_error( "Error when creating variable step-size RK integrator, fixed step coefficients are used ("+ coefficients_.name +")." );
        }
    }

    //! Get step size of the next step.
    /*!
     * Returns the step size of the next step.
     * \return Step size to be used for the next step.
     */
    virtual TimeStepType getNextStepSize( ) const { return this->stepSize_; }

    //! Get current state.
    /*!
     * Returns the current state of the integrator.
     * \return Current integrated state.
     */
    virtual StateType getCurrentState( ) const { return this->currentState_; }

    //! Get current independent variable.
    /*!
     * Returns the current value of the independent variable of the integrator.
     * \return Current independent variable.
     */
    virtual IndependentVariableType getCurrentIndependentVariable( ) const
    {
        return this->currentIndependentVariable_;
    }

    //! Get current state derivatives.
    /*!
     * Returns the current state derivatives, i.e., the values of k_{i} (stage evaluations) in
     * Runge-Kutta scheme.
     * \return Current state derivatives evaluated according to stages of Runge-Kutta scheme.
     */
    std::vector< StateDerivativeType > getCurrentStateDerivatives( )
    {
        return currentStateDerivatives_;
    }

    //! Perform a single integration step.
    /*!
     * Perform a single integration step and compute a new step size.
     * \param stepSize The step size to take. If the time step is too large to satisfy the error
     *          constraints, the step is redone until the error constraint is satisfied.
     * \return The state at the end of the interval.
     */
    virtual StateType performIntegrationStep( const TimeStepType stepSize );

    //! Rollback internal state to the last state.
    /*!
     * Performs rollback of the internal state to the last state. This function can only be called
     * once after calling integrateTo( ) or performIntegrationStep( ) unless specified otherwise by
     * implementations, and can not be called before any of these functions have been called. Will
     * return true if the rollback was successful, and false otherwise.
     * \return True if the rollback was successful.
     */
    virtual bool rollbackToPreviousState( )
    {
        if ( this->currentIndependentVariable_ == this->lastIndependentVariable_ )
        {
            return false;
        }

        this->currentIndependentVariable_ = this->lastIndependentVariable_;
        this->currentState_ = this->lastState_;
        return true;
    }

    //! Get previous independent variable.
    /*!
     * Returns the previoius value of the independent variable of the integrator.
     * \return Previous independent variable.
     */
    IndependentVariableType getPreviousIndependentVariable( )
    {
        return this->lastIndependentVariable_;
    }

    //! Get previous state value.
    /*!
     * Returns the previous value of the state.
     * \return Previous state
     */
    StateType getPreviousState( )
    {
        return this->lastState_;
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

    //! Function to toggle the use of step-size control
    /*!
     * Function to toggle the use of step-size control
     * \param useStepSizeControl Boolean denoting whether step size control is to be used
     */
    void setStepSizeControl( const bool useStepSizeControl )
    {
        useStepSizeControl_ = useStepSizeControl;
    }

protected:

    //! Computes the next step size and validates the result.
    /*!
     * Computes the next step size based on a higher and lower order estimate, determines if the
     * error is within bounds and returns a new step size.
     * \param lowerOrderEstimate The integrated result with the lower order coefficients.
     * \param higherOrderEstimate The integrated result with the higher order coefficients.
     * \param stepSize The step size used to obtain these results.
     * \return True if the error was within bounds, false otherwise.
     */
    virtual bool computeNextStepSizeAndValidateResult( const StateType& lowerOrderEstimate,
                                                       const StateType& higherOrderEstimate,
                                                       const TimeStepType stepSize );

    //! Compute new step size.
    /*!
     * Computes the new step size based on a generic definition of the local truncation error.
     * \param stepSize Integration step size of current step.
     * \param orders Pair of lower and higher orders of the two schemes used in variable step size
     *           integration. Note that the order is important (lower first, higher second).
     * \param safetyFactorForNextStepSize Safety factor used to scale prediction of next step size.
     * \param minimumAndMaximumFactorsForNextStepSize Pair of minimum and maximum safety factor
     *          decrease and increase for computing the next step size. Note that the order is
     *          important (minimum first, maximum second).
     * \param relativeErrorTolerance Allowable relative error between integrations using two
     *           schemes.
     * \param absoluteErrorTolerance Allowable relative error between integrations using two
     *           schemes.
     * \param lowerOrderEstimate Numerical integration result using lower order scheme.
     * \param higherOrderEstimate Numerical integration result using higher order scheme.
     * \return Pair with new step size and a boolean denoting whether the step size computation
     *           was succesfull, i.e. whether tolerances etc. are met.
     */
    virtual std::pair< TimeStepType, bool > computeNewStepSize(
            const TimeStepType stepSize,
            const std::pair< TimeStepType, TimeStepType >& orders,
            const TimeStepType safetyFactorForNextStepSize,
            const std::pair< TimeStepType, TimeStepType >& minimumAndMaximumFactorsForNextStepSize,
            const StateType& relativeErrorTolerance, const StateType& absoluteErrorTolerance,
            const StateType& lowerOrderEstimate, const StateType& higherOrderEstimate );

    //! Last used step size.
    /*!
     * Last used step size, passed to either integrateTo( ) or performIntegrationStep( ).
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
     * Last state as computed by performIntegrationStep( ).
     */
    StateType lastState_;

    //! Coefficients for the integrator.
    /*!
     * Coefficients for the integrator, as defined by the Butcher tableau.
     */
    RungeKuttaCoefficients coefficients_;

    //! Minimum step size.
    /*!
     * Minimum step size.
     */
    TimeStepType minimumStepSize_;

    //! Maximum step size.
    /*!
     * Maximum step size.
     */
    TimeStepType maximumStepSize_;

    //! Relative error tolerance.
    /*!
     * Relative error tolerance per element in the state.
     */
    StateType relativeErrorTolerance_;

    //! Absolute error tolerance.
    /*!
     * Absolute error tolerance per element in the state.
     */
    StateType absoluteErrorTolerance_;

    //! Safety factor for next step size.
    /*!
     * Safety factor used to scale prediction of next step size. This is usually picked between
     * 0.8 and 0.9 (Burden and Faires, 2001).
     */
    TimeStepType safetyFactorForNextStepSize_;

    //! Maximum factor increase for next step size.
    /*!
     * The maximum factor by which the next step size can increase compared to the current value.
     * The need for this maximum stems from a need to ensure that the step size changes do not
     * alias with the dynamics of the model being integrated. This is typically set at 4.0, based
     * on numerical experiments (Burden and Faires, 2001).
     */
    TimeStepType maximumFactorIncreaseForNextStepSize_;

    //! Minimum factor decrease for next step size.
    /*!
     * The minimum factor by which the next step size can decrease compared to the current value.
     * The need for this minimum stems from a need to ensure that the step size changes do not
     * alias with the dynamics of the model being integrated. This is typically set at 0.1, based
     * on numerical experiments (Burden and Faires, 2001).
     */
    TimeStepType minimumFactorDecreaseForNextStepSize_;

    //! Function that returns the new step size computed.
    /*!
     * Function that returns the new step size computed, as passed to the constructor.
     */
    NewStepSizeFunction newStepSizeFunction_;

    //! Vector of state derivatives.
    /*!
     * Vector of state derivatives, i.e. values of k_{i} in Runge-Kutta scheme.
     */
    std::vector< StateDerivativeType > currentStateDerivatives_;

    bool exceptionIfMinimumStepExceeded_;

    //! Boolean denoting whether step size control is to be used
    bool useStepSizeControl_;

};

extern template class RungeKuttaVariableStepSizeIntegrator < double, Eigen::VectorXd, Eigen::VectorXd >;
extern template class RungeKuttaVariableStepSizeIntegrator < double, Eigen::Vector6d, Eigen::Vector6d >;
extern template class RungeKuttaVariableStepSizeIntegrator < double, Eigen::MatrixXd, Eigen::MatrixXd >;


//! Perform a single integration step.
template< typename IndependentVariableType, typename StateType, typename StateDerivativeType, typename TimeStepType >
StateType
RungeKuttaVariableStepSizeIntegrator< IndependentVariableType, StateType, StateDerivativeType, TimeStepType >
::performIntegrationStep( const TimeStepType stepSize )
{
    if( !( stepSize == stepSize ) )
    {
        throw std::invalid_argument( "Error in RKF integrator, step size is NaN" );
    }

    // Define and allocated vector for the number of stages.
    currentStateDerivatives_.clear( );
    currentStateDerivatives_.reserve( this->coefficients_.cCoefficients.rows( ) );

    // Define lower and higher order estimates.
    StateType lowerOrderEstimate( this->currentState_ ), higherOrderEstimate( this->currentState_ );

    // Compute the k_i state derivatives per stage.
    for ( int stage = 0; stage < this->coefficients_.cCoefficients.rows( ); stage++ )
    {
        // Compute the intermediate state to pass to the state derivative for this stage.
        StateType intermediateState( this->currentState_ );

        // Compute the intermediate state.
        for ( int column = 0; column < stage; column++ )
        {
            intermediateState += stepSize * this->coefficients_.aCoefficients( stage, column ) *
                    currentStateDerivatives_[ column ];
        }

        // Compute the state derivative.
        const IndependentVariableType time = this->currentIndependentVariable_ +
                this->coefficients_.cCoefficients( stage ) * stepSize;
        currentStateDerivatives_.push_back( this->stateDerivativeFunction_( time, intermediateState ) );

        // Check if propagation should terminate because the propagation termination condition has been reached
        // while computing the intermediate state.
        // If so, return immediately the current state (not recomputed yet), which will be discarded.
        if ( this->propagationTerminationFunction_( static_cast< double >( time ), TUDAT_NAN ) )
        {
            this->propagationTerminationConditionReachedDuringStep_ = true;
            return this->currentState_;
        }

        // Update the estimate.
        lowerOrderEstimate += this->coefficients_.bCoefficients( 0, stage ) * stepSize *
                currentStateDerivatives_[ stage ];
        higherOrderEstimate += this->coefficients_.bCoefficients( 1, stage ) * stepSize *
                currentStateDerivatives_[ stage ];
    }

    // Determine if the error was within bounds and compute a new step size.
    if ( computeNextStepSizeAndValidateResult( lowerOrderEstimate,
                                               higherOrderEstimate, stepSize ) )
    {
        // Accept the current step.
        this->lastIndependentVariable_ = this->currentIndependentVariable_;
        this->lastState_ = this->currentState_;
        this->currentIndependentVariable_ += stepSize;

        switch ( this->coefficients_.orderEstimateToIntegrate )
        {
        case RungeKuttaCoefficients::lower:
            this->currentState_ = lowerOrderEstimate;
            return this->currentState_;

        case RungeKuttaCoefficients::higher:
            this->currentState_ = higherOrderEstimate;
            return this->currentState_;

        default: // The default case will never occur because OrderEstimateToIntegrate is an enum.
            throw std::runtime_error( "Order estimate to integrate is invalid." );
        }
    }
    else
    {
        // Reject current step.
        return performIntegrationStep( this->stepSize_ );
    }
}

//! Compute the next step size and validate the result.
template< typename IndependentVariableType, typename StateType, typename StateDerivativeType, typename TimeStepType >
bool
RungeKuttaVariableStepSizeIntegrator< IndependentVariableType, StateType, StateDerivativeType, TimeStepType >
::computeNextStepSizeAndValidateResult(
        const StateType& lowerOrderEstimate, const StateType& higherOrderEstimate,
        const TimeStepType stepSize )
{
    if( useStepSizeControl_ )
    {
        // Compute new step size using new step size function, which also returns whether the
        // relative error is within bounds or not.
        std::pair< TimeStepType, bool > newStepSizePair = this->newStepSizeFunction_(
                    stepSize, std::make_pair( this->coefficients_.lowerOrder, this->coefficients_.higherOrder ),
                    this->safetyFactorForNextStepSize_, std::make_pair( this->minimumFactorDecreaseForNextStepSize_,
                                                                        this->maximumFactorIncreaseForNextStepSize_ ),
                    this->relativeErrorTolerance_, this->absoluteErrorTolerance_, lowerOrderEstimate, higherOrderEstimate );

        // Check whether change in stepsize does not exceed bounds.
        // If the stepsize is reduced to less than the prescibed minimum factor, set to minimum factor.
        // If the stepsize is increased to more than the prescribed maximum factor, set to maximum
        // factor. These bounds are necessary to prevent the stepsize changes from aliasing
        // with the dynamics of the system of ODEs.
        // Also check if maximum step size is exceeded and step next step size to maximum if necessary.
        // Typically used bounds can be found in (Burden and Faires, 2001).
        if ( newStepSizePair.first / stepSize <= this->minimumFactorDecreaseForNextStepSize_ )
        {
            this->stepSize_ = stepSize * this->minimumFactorDecreaseForNextStepSize_;
        }

        else if ( newStepSizePair.first / stepSize >= maximumFactorIncreaseForNextStepSize_ )
        {
            this->stepSize_ = stepSize * this->maximumFactorIncreaseForNextStepSize_;
        }

        else
        {
            this->stepSize_ = newStepSizePair.first;
        }

        // Check if minimum step size is violated and throw exception if necessary.
        if ( std::fabs( this->stepSize_ ) < std::fabs( this->minimumStepSize_ ) )
        {
            if( exceptionIfMinimumStepExceeded_ )
            {
                throw MinimumStepSizeExceededError( std::fabs( this->minimumStepSize_ ),
                                                    std::fabs( this->stepSize_ ) );
            }
            else
            {
                this->stepSize_ = stepSize / std::fabs( stepSize ) * std::fabs( this->minimumStepSize_ );
                return true;
            }
        }
        else if( std::fabs( this->stepSize_ ) > std::fabs( this->maximumStepSize_ ) )
        {
            this->stepSize_ = stepSize / std::fabs( stepSize ) * std::fabs( this->maximumStepSize_ );
        }

        if( stepSize * this->stepSize_ < 0 )
        {
            throw std::runtime_error( "Error during step size control, step size flipped sign" );
        }

        // Check if computed error in state is too large and reject step if true.
        return newStepSizePair.second;
    }
    else
    {
        this->stepSize_ = stepSize;
        return true;
    }
}

//! Compute new step size.
/*!
 * Computes the new step size based on a generic definition of the local truncation error.
 */
template< typename IndependentVariableType, typename StateType, typename StateDerivativeType, typename TimeStepType >
std::pair< TimeStepType, bool >
RungeKuttaVariableStepSizeIntegrator< IndependentVariableType, StateType, StateDerivativeType, TimeStepType >
::computeNewStepSize(
        const TimeStepType stepSize,
        const std::pair< TimeStepType, TimeStepType >& orders,
        const TimeStepType safetyFactorForNextStepSize,
        const std::pair< TimeStepType, TimeStepType >& minimumAndMaximumFactorsForNextStepSize,
        const StateType& relativeErrorTolerance,
        const StateType& absoluteErrorTolerance,
        const StateType& lowerOrderEstimate,
        const StateType& higherOrderEstimate )
{
    TUDAT_UNUSED_PARAMETER( minimumAndMaximumFactorsForNextStepSize );

    // Compute the truncation error based on the higher and lower order estimates.
    const StateType truncationError_ =
            ( higherOrderEstimate - lowerOrderEstimate ).array( ).abs( );

    // Compute error tolerance based on relative and absolute error tolerances.
    const StateType errorTolerance_ =
            ( higherOrderEstimate.array( ).abs( ) *
              relativeErrorTolerance.array( ) ).matrix( )
            + absoluteErrorTolerance;

    // Compute relative truncation error. This will indicate if the current step satisfies the
    // required tolerances.
    const StateType relativeTruncationError_ = truncationError_.array( ) /
            errorTolerance_.array( );

    // Compute the maximum error based on the largest coefficient in the relative truncation error
    // matrix.
    const typename StateType::Scalar maximumErrorInState_
            = relativeTruncationError_.array( ).abs( ).maxCoeff( );

    // Compute the new step size. This is based off of the equation given in
    // (Montenbruck and Gill, 2005).
    const TimeStepType newStepSize = safetyFactorForNextStepSize * stepSize
            * std::pow( 1.0 / maximumErrorInState_, 1.0 / orders.second );

    // Check if the current state can be accepted.
    const bool isIntegrationStepAccepted = maximumErrorInState_ <= 1.0;

    // Return the computed new step size and whether the current step taken is acceptable. If it
    // isn't, the step will be recomputed with the computed new step size.
    return std::make_pair( newStepSize, isIntegrationStepAccepted );
}

//! Exception that is thrown if the minimum step size is exceeded.
/*!
 * Exception thrown by RungeKuttaVariableStepSizeIntegrator< >::computeNextStepSizeAndValidateResult()
 * if the minimum step size is exceeded.
 */
template< typename IndependentVariableType, typename StateType, typename StateDerivativeType, typename TimeStepType >
class RungeKuttaVariableStepSizeIntegrator< IndependentVariableType, StateType,
        StateDerivativeType, TimeStepType >::MinimumStepSizeExceededError : public std::runtime_error
{
public:

    //! Default constructor.
    /*!
     * Default constructor, initializes the parent runtime_error.
     * \param minimumStepSize_ The minimum step size allowed by the integrator.
     * \param requestedStepSize_ The new calculated step size.
     */
    MinimumStepSizeExceededError( TimeStepType minimumStepSize_,
                                  TimeStepType requestedStepSize_ ) :
        std::runtime_error( "Minimum step size exceeded." ),
        minimumStepSize( minimumStepSize_ ), requestedStepSize( requestedStepSize_ )
    { }

    //! The minimum step size allowed by the integrator.
    /*!
     * The minimum step size allowed by the integrator.
     */
    TimeStepType minimumStepSize;

    //! The new calculated step size.
    /*!
     * The new calculated step size.
     */
    TimeStepType requestedStepSize;

protected:

private:

};

//! Typedef of variable-step size Runge-Kutta integrator (state/state derivative = VectorXd,
//! independent variable = double).
/*!
 * Typedef of a variable-step size Runge-Kutta integrator with VectorXds as state and state
 * derivative and double as independent variable.
 */
typedef RungeKuttaVariableStepSizeIntegrator< > RungeKuttaVariableStepSizeIntegratorXd;

//! Typedef for shared-pointer to RungeKuttaVariableStepSizeIntegratorXd object.
typedef std::shared_ptr< RungeKuttaVariableStepSizeIntegratorXd >
RungeKuttaVariableStepSizeIntegratorXdPointer;

} // namespace numerical_integrators

} // namespace tudat

#endif // TUDAT_RUNGE_KUTTA_VARIABLE_STEP_SIZE_INTEGRATOR_H
