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
 *      Allione, M.S, The N-Body Problem and Special Perturbation Techniques, NASA, 1968.
 *      Bate, R, Fundamentals of astrodynamcs, Dover, 1975.
 *      Shampine, L and Gordon, M., Local Error and Variable Order Adams Codes,
 *          Applied Mathematic and Computation, 1975.
 *
 *    Notes
 *      See Matlab scripts for generating predictor, corrector and error estimation coefficients.
 *
 */

#ifndef TUDAT_ADAMS_BASHFORTH_MOULTON_INTEGRATOR_H
#define TUDAT_ADAMS_BASHFORTH_MOULTON_INTEGRATOR_H

#include <deque>
#include <algorithm>
#include <limits>

#include <memory>

#include <Eigen/Core>

#include "tudat/math/integrators/numericalIntegrator.h"
#include "tudat/math/integrators/reinitializableNumericalIntegrator.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/math/integrators/rungeKuttaVariableStepSizeIntegrator.h"

namespace tudat
{
namespace numerical_integrators
{

//! Adams-Bashforth-Moulton Variable Order and Stepsize integrator.
/*!
 * Class that implements the Adams-Bashforth-Moulton integrator, variable order, variable
 * step size integrator.
 * \tparam StateType The type of the state. This type should support addition with
 *          StateDerivativeType.
 * \tparam StateDerivativeType The type of the state derivative. This type should support
 *          multiplication with IndependentVariableType and doubles.
 * \tparam IndependentVariableType The type of the independent variable.
 * \sa NumericalIntegrator.
 */
template< typename IndependentVariableType = double, typename StateType = Eigen::VectorXd,
          typename StateDerivativeType = Eigen::VectorXd, typename TimeStepType = IndependentVariableType >
class AdamsBashforthMoultonIntegrator
        : public ReinitializableNumericalIntegrator<
        IndependentVariableType, StateType, StateDerivativeType, TimeStepType >
{
    
public:

    //! Typedef for the base class.
    /*!
     * Typedef of the base class with all template parameters filled in.
     */
    typedef ReinitializableNumericalIntegrator< IndependentVariableType, StateType,
    StateDerivativeType, TimeStepType > ReinitializableNumericalIntegratorBase;

    //! Typedef for the state derivative function.
    /*!
     * Typedef to the state derivative function inherited from the base class.
     * \sa NumericalIntegrator::StateDerivativeFunction.
     */
    typedef typename ReinitializableNumericalIntegratorBase::NumericalIntegratorBase::
    StateDerivativeFunction StateDerivativeFunction;

    //! Default constructor.
    /*!
     * Default constructor, taking a state derivative function as
     * argument. Note that with the default settings the integrator
     * has no incentive of changing the stepsize.
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
     * \param bandwidth (optional) bandwidth (between min and max tolerance) (200 default)
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    AdamsBashforthMoultonIntegrator(
            const StateDerivativeFunction& stateDerivativeFunction,
            const IndependentVariableType intervalStart,
            const StateType& initialState,
            const TimeStepType minimumStepSize,
            const TimeStepType maximumStepSize,
            const TimeStepType& initialStepSize,
            const StateType& relativeErrorTolerance,
            const StateType& absoluteErrorTolerance,
            const TimeStepType bandwidth = 200. )
        : ReinitializableNumericalIntegratorBase( stateDerivativeFunction ),
          currentIndependentVariable_( intervalStart ),
          currentState_( initialState ),
          lastIndependentVariable_( intervalStart ),
          minimumStepSize_( std::fabs( static_cast< double >( minimumStepSize ) ) ),
          maximumStepSize_( std::fabs( static_cast< double >( maximumStepSize ) ) ),
          stepSize_( initialStepSize ),
          relativeErrorTolerance_( relativeErrorTolerance.array( ).abs( ) ),
          absoluteErrorTolerance_( absoluteErrorTolerance.array( ).abs( ) ),
          bandwidth_( std::fabs( static_cast< double >( bandwidth ) ) )
    {
        if( !( currentState_.rows( ) == relativeErrorTolerance_.rows( ) ) ||
                !( currentState_.cols( ) == relativeErrorTolerance_.cols( ) ) )
        {
            throw std::runtime_error( "Error when creating ABM integrator, relative tolerance input size is inconsistent" );
        }

        if( !( currentState_.rows( ) == absoluteErrorTolerance_.rows( ) ) ||
                !( currentState_.cols( ) == absoluteErrorTolerance_.cols( ) ) )
        {
            throw std::runtime_error( "Error when creating ABM integrator, absolute tolerance input size is inconsistent" );
        }


        fixedStepSize_ = false;
        strictCompare_ = true;
        fixedOrder_ = false;
        minimumOrder_ = 6;
        maximumOrder_ = 11;
        order_ = minimumOrder_;

        if( ( initialStepSize == minimumStepSize ) && ( initialStepSize == maximumStepSize ) )
        {
            fixedStepSize_ = true;
        }


        
        // Start filling the state and state derivative history deques.
        stateHistory_.push_front( currentState_ );
        derivHistory_.push_front( this->stateDerivativeFunction_(
                                      currentIndependentVariable_, currentState_ ));
    }

    //! Default constructor.
    /*!
     * Default constructor, taking a state derivative function as
     * argument. Note that with the default settings the integrator
     * has no incentive of changing the stepsize.
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
     * \param bandwidth (optional) bandwidth (between min and max tolerance) (200 default)
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    AdamsBashforthMoultonIntegrator(
            const StateDerivativeFunction& stateDerivativeFunction,
            const IndependentVariableType intervalStart,
            const StateType& initialState,
            const TimeStepType minimumStepSize,
            const TimeStepType maximumStepSize,
            const TimeStepType& initialStepSize,
            const typename StateType::Scalar relativeErrorTolerance,
            const typename StateType::Scalar absoluteErrorTolerance,
            const TimeStepType bandwidth = 200. )
        : AdamsBashforthMoultonIntegrator(
              stateDerivativeFunction,
              intervalStart,
              initialState,
              minimumStepSize,
              maximumStepSize,
              initialStepSize,
              StateType::Constant( initialState.rows( ), initialState.cols( ), std::fabs( relativeErrorTolerance ) ),
              StateType::Constant( initialState.rows( ), initialState.cols( ), std::fabs( absoluteErrorTolerance ) ),
              bandwidth ) { }
    
    ~AdamsBashforthMoultonIntegrator( ){ }

    //! Get step size of the next step.
    /*!
     * Returns the step size of the next step. If the
     * performIntegrationStep( stepSize ) function is used during a
     * variable stepsize scheme, this function can be used to
     * find out which stepSize was deemed best for the next
     * integration step by the integrator.
     * \return Step size to be used for the next step.
     */
    virtual TimeStepType getNextStepSize( ) const { return stepSize_; }

    //! Get the order of next step.
    /*!
     * Returns the order of the next step
     * \return Order of the next step
     */
    unsigned int getOrder( ) const { return order_; }

    //! Get step size of the last step.
    /*!
     * Returns the step size of the last step. If the step size is not defined by
     * the user, the step size of last step can be retrieved with this function.
     * \return Step size to be used for the last step.
     */
    virtual TimeStepType getStepSize( ) const { return lastStepSize_; }

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
        if( !( stepSize == stepSize) )
        {
            throw std::runtime_error( "Error in ABM integrator, step size is NaN" );
        }

        // If stepSize is not same as old, clear the step-size dependent histories.
        if ( stepSize != stepSize_ )
        {
            // Pop all values from the history deque (the history is
            // invalid as it is dependent on the stepSize), except for
            // the current state and state derivative.
            while ( stateHistory_.size( ) > 1 )
            {
                stateHistory_.pop_back( );
            }
            while ( derivHistory_.size( ) > 1 )
            {
                derivHistory_.pop_back( );
            }
            stepSize_ = stepSize;
        }
        return performIntegrationStep( );
    }

    //! Perform a single integration step.
    /*!
     * Perform a single integration step. The stepSize as computed by
     * integrator during the last stop.
     * \return The state at the end of the interval.
     */
    StateType performIntegrationStep( )
    {
//        std::cout<<"Orders "<<order_<<" "<<minimumOrder_<<" "<<maximumOrder_<<" "<<fixedStepSize_<<std::endl;
        // Set last* variables for rollback.
        lastStepSize_ = stepSize_;
        lastState_ = stateHistory_.back( );
        lastDerivative_ = derivHistory_.back( );
        lastIndependentVariable_ = currentIndependentVariable_;

        // Remove old elements so enough are left to calculate predicted and corrected.
        // max twice the order, to facilitatie a doubling, halving, and order change.
        while ( stateHistory_.size( ) > order_ * 2 )
        {
            stateHistory_.pop_back( );
        }
        while ( derivHistory_.size( ) > order_ * 2 )
        {
            derivHistory_.pop_back( );
        }
        unsigned int sizeStateHistory = static_cast< unsigned int >( stateHistory_.size( ) );
        unsigned int sizeDerivativeHistory = static_cast< unsigned int >( derivHistory_.size( ) );
        unsigned int possibleOrder = std::min( sizeStateHistory, sizeDerivativeHistory );

        StateType correctedState;
        StateType predictedState;

        // Check if enough history steps are available to perform AM
        // step if not use a single-step method.
        if ( possibleOrder < minimumOrder_ || possibleOrder < order_ )
        {
            correctedState = performSingleStep( );
        }
        else
        {
            StateType predictedState = performPredictorStep( order_, false );
            predictedDerivative_ = this->stateDerivativeFunction_( currentIndependentVariable_ +
                                                                   stepSize_, predictedState );
            correctedState  = performCorrectorStep( predictedState, order_, false );
            absoluteError_ = estimateAbsoluteError( predictedState, correctedState, order_ );
            relativeError_ = estimateRelativeError( predictedState, correctedState, absoluteError_ );
        }

        // Change order to one that gives a higher predicted accuracy
        // Add tolenaces
        StateType predictorAbsoluteError;
        StateType predictorRelativeError;
        
        // If order is not fixed, order is not max yet and enough
        // history is available, then predict the error of an order
        // more.
        if ( !fixedOrder_ && order_ < maximumOrder_ && order_ < possibleOrder )
        {
            predictedState = performPredictorStep( order_ + 1, false );
            correctedState = performCorrectorStep( predictedState, order_ + 1, false );
            predictorAbsoluteError = estimateAbsoluteError( predictedState, correctedState, order_ + 1 );
            predictorRelativeError = estimateRelativeError( predictedState, correctedState, predictorAbsoluteError );

            // If the predicted error is less than the current error,
            // increase the error.
            if ( errorCompare( predictorAbsoluteError, predictorRelativeError, absoluteError_, relativeError_ ) )
            {
                
                order_++;
            }
            // Else if the order is not fixed, the order is not min yet
            // and there is enough history available, then predict the
            // error of an order less.
        }
        else if ( !fixedOrder_ && order_ > minimumOrder_ && order_ - 1 <= possibleOrder )
        {
            predictedState = performPredictorStep( order_ - 1, false );
            correctedState = performCorrectorStep( predictedState, order_ - 1, false );
            predictorAbsoluteError = estimateAbsoluteError( predictedState, correctedState, order_ - 1 );
            predictorRelativeError = estimateRelativeError( predictedState, correctedState, predictorAbsoluteError );
            // If it is less than the current order, lower the order.
            if ( errorCompare( predictorAbsoluteError, predictorRelativeError, absoluteError_, relativeError_ ) )
            {
                order_--;
            } else {
                predictorAbsoluteError = absoluteError_;
                predictorRelativeError = relativeError_;
            }
        }
        else
        {
            predictorAbsoluteError = absoluteError_;
            predictorRelativeError = relativeError_;
        }

        // If the error (after order change) is too big, stepsize
        // isn't fixed and will not become too small, then halve the
        // stepsize.
        if ( errorTooLarge( predictorAbsoluteError, predictorRelativeError )
             && std::fabs( stepSize_ / 2.0 )> minimumStepSize_ && !fixedStepSize_ )
        {
            // Set up new data for halving
            std::deque< StateType > tempStateHistory;
            std::deque< StateType > tempDerivativeHistory;
            StateType midState = currentState_;
            StateDerivativeType midDerivative = midState;
            unsigned int interpolationStateIndex;
            unsigned int interpolationDerivativeIndex;

            // Only if there is enough historical data, is it possible to create the intermediate
            // points. Redefine order based on available data, it could drop, but can't be lower
            // that minimum order
            unsigned int possibleHalvingOrder = std::min( 2 * possibleOrder - 1, maximumOrder_ );
            order_ = std::min( possibleHalvingOrder, order_ );
            
            // FIXME: make this a setting?
            // Reduce order. Too much backlog aversly affects the accuracy when halving.
            if( order_ > minimumOrder_ )
            {
                order_ = minimumOrder_;
            }
            
            for( unsigned int i = 0; i < possibleHalvingOrder ; i++ )
            {
                // If states are even, they already exist, no need to interpolate
                if ( i % 2 == 0 )
                {
                    tempStateHistory.push_back( stateHistory_.at( i / 2 ));
                    tempDerivativeHistory.push_back( derivHistory_.at( i / 2 ));
                } else {
                    // Reset midpoint state and deriv to zero
                    midState = 0 * midState;
                    midDerivative = midState;
                    interpolationDerivativeIndex = ( order_ - 1 ) * ( order_ - 1 ) + ( i - 1 ) / 2;
                    interpolationStateIndex = interpolationDerivativeIndex - order_ + 1;
                    for( unsigned int j = 0; j < order_; j++ )
                    {
                        midState += interpolationCoefficients[ interpolationStateIndex ][ j ] *
                                stateHistory_.at( j ) + interpolationCoefficients[ interpolationStateIndex ][ order_ + j ] *
                                derivHistory_.at( j ) * stepSize_;
                        midDerivative += interpolationCoefficients[ interpolationDerivativeIndex ][ j ] *
                                stateHistory_.at( j ) / stepSize_
                                + interpolationCoefficients[ interpolationDerivativeIndex ][ order_ + j ] *
                                derivHistory_.at( j );
                    }
                    tempStateHistory.push_back( midState );
                    tempDerivativeHistory.push_back( midDerivative );
                }
            }
            
            // Set the new history and stepsize
            stateHistory_ = tempStateHistory;
            derivHistory_ = tempDerivativeHistory;
            stepSize_ = stepSize_ / 2.0;
            
            // Temporarily turn halving off.
            fixedStepSize_ = true;
            StateType halvedState = performIntegrationStep( );
            fixedStepSize_ = false;
            return halvedState;
        } // end if ( errorTooLarge( ...
        
        // If the error (after order change ) is too small, the
        // stepsize isn't fixed and the and will not become too big,
        // then double the stepsize.
        if ( errorTooSmall( predictorAbsoluteError, predictorRelativeError )
             && sizeDerivativeHistory >= 2 * order_
             && std::fabs( stepSize_ * 2.0 ) <= maximumStepSize_ && !fixedStepSize_ )
        {
            
            // Predict error after doubling, to prevent error from becoming too big
            // This prevents fluttering and throwing away states from the history that will be
            // needed again later.
            // Note that two important assumptions are made here:
            // 1. The error made at the current step is representative of the error that
            //    would be made if we blindly double and only compute the error next step.
            // 2. The difference in the derivative of the predicted state (predictedDerivative_)
            //    at the normal stepsize (already computed) with the doubled stepsize (not computed)
            //    is neglibile. This assumption saves one function evaluation.
            // It's possible to reuse previously defined variables here except for correctedState
            // which is still used below.
            predictedState = performPredictorStep( order_ , true );
            StateType doubleStepCorrectedState = performCorrectorStep( predictedState, order_, true );
            predictorAbsoluteError = estimateAbsoluteError( predictedState, doubleStepCorrectedState, order_ );
            predictorRelativeError = estimateRelativeError( predictedState, doubleStepCorrectedState,
                                                            predictorAbsoluteError );

            // Only update the history if the error will not be too large
            if ( !errorTooLarge( predictorAbsoluteError, predictorRelativeError ) )
            {
                // Note that the history should be at least 7 to allow successful
                // continuation of the AM scheme.
                std::deque< StateType > tempStateHistory;
                std::deque< StateType > tempDerivativeHistory;
                
                // Use old history to fill new history, skipping every other entry starting at 1
                for( unsigned int i = 1; i < sizeStateHistory; i += 2 )
                {
                    tempStateHistory.push_back( stateHistory_.at( i ));
                }
                for( unsigned int i = 1; i < sizeDerivativeHistory; i += 2 )
                {
                    tempDerivativeHistory.push_back( derivHistory_.at( i ));
                }
                
                // Set the new history and stepsize
                stateHistory_ = tempStateHistory;
                derivHistory_ = tempDerivativeHistory;
                stepSize_ = stepSize_ * 2.0;
            }
        } // end if ( errorTooSmall( ...

        // Move computed state to history
        currentIndependentVariable_ += lastStepSize_;
        currentState_ = correctedState;
        stateHistory_.push_front( currentState_ );
        derivHistory_.push_front( this->stateDerivativeFunction_(
                                      currentIndependentVariable_, currentState_ ) );
        return currentState_;
    }

    //! Rollback internal state to the last state.
    /*!
     * Performs rollback of internal state to the last state. This function can only be called once
     * after calling integrateTo( ) or performIntegrationStep( ) unless specified otherwise by
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
        stepSize_ = lastStepSize_;
        stateHistory_.push_back( lastState_ );
        derivHistory_.push_back( lastDerivative_ );
        stateHistory_.pop_front( );
        derivHistory_.pop_front( );
        derivHistory_.pop_front( );
        currentState_ = stateHistory_.front( );
        // Recalculate the derivative in order to make sure that all
        // update functions inside state derivative model get reactivated
        derivHistory_.push_front( this->stateDerivativeFunction_(
                                      currentIndependentVariable_, currentState_ ) );
        return true;
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

        // Clear the history and initiate with new state and derivative.
        stateHistory_[ 0 ] = currentState_;
        derivHistory_[ 0 ] = this->stateDerivativeFunction_(
                    currentIndependentVariable_, currentState_ );
        if ( !allowRollback )
        {
            lastIndependentVariable_ = currentIndependentVariable_;
        }

    }

    //! Return maximum truncation error.
    /*!
     * Return the truncation error to be estimated by computError( ).
     * \return Double value of the truncation error.
     */
    double getMaximumError( ){ return absoluteError_.cwiseMin( relativeError_ ).array( ).maxCoeff( ); }

    //! Return absolute truncation error.
    /*!
     *
     * \return Vector of the absolute truncation error.
     */
    StateType getAbsoluteError( ){ return absoluteError_; }

    //! Return relative truncation error.
    /*!
     *
     * \return Vector of the relative truncation error.
     */
    StateType getRelativeError( ){ return relativeError_; }

    //! Return list derivative .
    /*!
     * Return the last value of the derivative
     * \return last derivative value
     */
    StateDerivativeType getLastDerivative( ){ return derivHistory_.front( ); }

    //! (Un )set fixed step size.
    /*!
     * Fixes the step if set to true.
     * \param fixedStepSize false or true to ( un)set fixed step size (default false).
     */
    void setFixedStepSize( bool fixedStepSize ) { fixedStepSize_ = fixedStepSize; }

    //! (Un)set fixed order.
    /*!
     * Fixes the order if set to true.
     * \param fixedOrder false or true to ( un)set fixed order (default false).
     */
    void setFixedOrder( bool fixedOrder ) { fixedOrder_ = fixedOrder; }

    //! Set order.
    /*!
     * Sets the order to specified.
     * \param order desired.
     */
    void setOrder( unsigned int order ) { order_ = order; }

    //! (Un)set fixed order.
    /*!
     * For a error vector to be better than another all its elements
     * must be equal or small. Not strict compares the norm of the
     * vectors instead.
     * \param strictCompare false or true to ( un)set strict compare (default true).
     */
    void setStrictCompare( bool strictCompare ) { strictCompare_ = strictCompare; }

    //! Set the minimum and maximum step size.
    /*!
     * Change the bounds of the stepsize.
     * \param minimumStepSize minimum stepsize
     * \param maximumStepSize maximum stepsize
     */
    void setStepSizeBounds( IndependentVariableType minimumStepSize, IndependentVariableType maximumStepSize )
    {
        minimumStepSize_ = minimumStepSize;
        maximumStepSize_ = maximumStepSize;
    }

    //! Set the minimum and maximum truncation error.
    /*!
     * Change the bounds of the error.
     * \param absoluteErrorTolerance absolute error tolerance.
     * \param relativeErrorTolerance relative error tolerance.
     */
    void setErrorTolerances( StateType& absoluteErrorTolerance, StateType& relativeErrorTolerance )
    {
        absoluteErrorTolerance_ = absoluteErrorTolerance;
        relativeErrorTolerance_ = relativeErrorTolerance;
    }
    
    //! Set step size.
    /*!
     * Set step size, to be used with empty performIntegrationStep( ).
     * \param stepSize of the new integration step.
     */
    void setStepSize( IndependentVariableType stepSize ) { stepSize_ = stepSize; }

    //! Function to reset the minimum order of the integrator
    /*!
     * Function to reset the minimum order of the integrator
     * \param minimumOrder New minimum order of the integrator
     */
    void setMinimumOrder( unsigned int minimumOrder )
    {
        minimumOrder_ = minimumOrder;
        if( order_ < minimumOrder )
        {
            order_ = minimumOrder_;
        }
    }

    //! Function to reset the maximum order of the integrator
    /*!
     * Function to reset the maximum order of the integrator
     * \param maximumOrder New maximum order of the integrator
     */
    void setMaximumOrder( unsigned int maximumOrder )
    {
        maximumOrder_ = maximumOrder;
        if( order_ > maximumOrder_ )
        {
            order_ = maximumOrder_;
        }
    }

    IndependentVariableType getPreviousIndependentVariable( )
    {
        return lastIndependentVariable_;
    }

    StateType getPreviousState( )
    {
        return lastState_;
    }

protected:


    //! Current independent variable.
    /*!
     * Current independent variable as computed by performIntegrationStep( ).
     */
    IndependentVariableType currentIndependentVariable_;

    //! Current state.
    /*!
     * Current state as computed by performIntegrationStep( ).
     */
    StateType currentState_;

    //! Last independent variable.
    /*!
     * Last independent variable value as computed by performIntegrationStep( ).
     */
    IndependentVariableType lastIndependentVariable_;

    //! Minimum step size.
    /*!
     * Parameter that limits the value of the step size.
     * Important to have to prevent infinite halving.
     * Not initialised, constructor does so.
     */
    TimeStepType minimumStepSize_;

    //! Maximum step size.
    /*!
     * Parameter that limits the value of the step size. Not really important, can
     * be set pretty high.
     * Not initialised, constructor does so.
     */
    TimeStepType maximumStepSize_;


    //! Last used step size.
    /*!
     * Last used step size, passed to either integrateTo( ) or performIntegrationStep( ).
     */
    TimeStepType stepSize_;

    //! Maximum relative error.
    /*!
     * Parameter that limits the truncation error. Which drives stepsize control.
     * Not initialised, constructor does so.
     */
    StateType relativeErrorTolerance_;
    
    //! Maximum absolute error.
    /*!
     * Parameter that limits the truncation error. Which drives stepsize control.
     * Not initialised, constructor does so.
     */
    StateType absoluteErrorTolerance_;

    //! Bandwidth.
    /*!
     * Parameter that increases performance. Which drives stepsize control.
     * Advisable to set at least 200 to prevent flutter.
     * Not initialised, constructor does so.
     */
    TimeStepType bandwidth_;

    //! Fixed step size.
    /*!
     * Fixes to step if set to true, can be manipulated through setFixedStepSize( ).
     */
    bool fixedStepSize_;

    //! Strict compare
    /*!
     * For a error vector to be better than another all its elements
     * must be equal or small. Not strict compares the norm of the
     * vectors instead.
     */
    bool strictCompare_;

    //! Fixed order.
    /*!
     * Fixes to order if set to true, can be manipulated through setFixedOrder( ).
     */
    bool fixedOrder_;

    //! Absolute truncation error.
    /*!
     * Truncation error to be estimated by estimateAbsoluteError( ).
     */
    StateType absoluteError_;

    //! Relative truncation error.
    /*!
     * Truncation error to be estimated by estimateRelativeError( ).
     */
    StateType relativeError_;


    //! Minimum order.
    /*!
     * Minimum order allowed. True mininum is 1, but very imprecise.
     * Advised to be at least three.
     */
    unsigned int minimumOrder_;

    //! Maximum order.
    /*!
     * Maximum order allowed. True maximum is 12, based on predetermined coefficients.
     * However, double precision has shown to be incapable of estimation twelfth-order errors.
     */
    unsigned int maximumOrder_;

    //! Order.
    /*!
     * Dynamic order. Starts at minimumOrder_ and is changed using order-control.
     * Advisable to change for fixth-order integration runs.
     */
    unsigned int order_;

    //! Numerical Single Step Integrator.
    /*!
     * Single step integrator to perform intitialization of multi-step history.
     */
    //RungeKuttaVariableStepSizeIntegrator< IndependentVariableType, StateType, StateType, TimeStepType > singleStepIntegrator_;

    //! Truncation error coefficients.
    /*!
     * Coefficients for estimating the truncation error.
     */
    const static double truncationErrorCoefficients[ 12 ];

    //! Extrapolation coefficients.
    /*!
     * Coefficients of predictor and corrector extrapolation polynomials.
     */
    const static double extrapolationCoefficients[ 24 ][ 12 ];

    //! Interpolation coefficients.
    /*!
     * Coefficients of state and derivative mid-points used during halving.
     */
    const static double interpolationCoefficients[ 132 ][ 24 ];

    //! Perform integration step.
    /*!
     * Perform integration step using built-in Runge-Kutta fourth
     * order method. Ideally this could be replaced by a pointer to
     * single-step method of choice by the user.
     * \return new integrated state
     */
    StateType performSingleStep( )
    {
        typename StateType::Scalar maximumTolerance = 1.0;
        
        // FIXME: define these integrators in the private scope of class and initialize in constructor.
        // Perhaps with an user-selectable RK variant. The modifyCurrentState statements below should
        // are already implemented to ensure that this will work.
        
        // Single step integrator with fixed stepsize and no error control
        RungeKuttaVariableStepSizeIntegrator< IndependentVariableType, StateType, StateType, TimeStepType >
                singleFixedStepIntegrator_( RungeKuttaCoefficients::get(
                                                rungeKutta87DormandPrince ),
                                            this->stateDerivativeFunction_, currentIndependentVariable_, currentState_,
                                            stepSize_, stepSize_, stepSize_, maximumTolerance, maximumTolerance );
        
        // Single step integrator with variable stepsize and error control
        RungeKuttaVariableStepSizeIntegrator< IndependentVariableType, StateType, StateType, TimeStepType >
                singleVariableStepIntegrator_(
                    RungeKuttaCoefficients::get( CoefficientSets::rungeKutta87DormandPrince ),
                    this->stateDerivativeFunction_, currentIndependentVariable_, currentState_,
                    minimumStepSize_, maximumStepSize_, stepSize_, relativeErrorTolerance_, absoluteErrorTolerance_ );
        
        StateType currentState;
        if( fixedStepSize_ )
        {
            singleFixedStepIntegrator_.modifyCurrentState( currentState_ );
            currentState = singleFixedStepIntegrator_.performIntegrationStep( stepSize_ );
            
            // Find out the step size that was used during integration
            lastStepSize_ = singleFixedStepIntegrator_.getCurrentIndependentVariable( ) - currentIndependentVariable_;
        }
        else
        {
            singleVariableStepIntegrator_.modifyCurrentState( currentState_ );
            currentState = singleVariableStepIntegrator_.performIntegrationStep( stepSize_ );
            
            // Find out the step size that was used during integration
            lastStepSize_ = singleVariableStepIntegrator_.getCurrentIndependentVariable( ) - currentIndependentVariable_;
        }
        
        // Even if a different step size is suggested, let's stick with the old one, since the goal is to start
        // filling up the deque at a constant stepsize interval
        stepSize_ = lastStepSize_; // singleStepIntegrator_.getNextStepSize( );

        // Disregard the ABAM error control in the performIntegrationStep function when using single steps.
        // The lines below will make sure to trigger neither errorTooLarge or errorTooSmall.
        absoluteError_ = absoluteErrorTolerance_ / bandwidth_;
        relativeError_ = relativeErrorTolerance_ / bandwidth_;
        return currentState;
    }
    
    //! Perform predictor step.
    /*!
     * Using the order find predicted estimate using the Adams-Bashforth predictor
     * \param order Order of the integration.
     * \param doubleStep Boolean if stepsize should be considered double, true for estimating doubling error.
     * \return state after predictor step
     */
    StateType performPredictorStep( unsigned int order, bool doubleStep )
    {
        // Calculate predicted state
        unsigned int stepsToSkip = static_cast< unsigned int>( doubleStep );
        TimeStepType stepSize = stepSize_ * static_cast< double >( stepsToSkip + 1 );
        StateType predictedState = stateHistory_.at( stepsToSkip );
        for ( unsigned int i = 0; i < order; i++ )
        {
            predictedState += extrapolationCoefficients[ order * 2 - 2 ][ i ] * stepSize *
                    derivHistory_.at( i * ( stepsToSkip + 1 ) + stepsToSkip );
        }
        return predictedState;
    }

    //! Perform correcter step.
    /*!
     * Using the order and predicted state, find corrected estimate using the Adams-Moulton corrector
     * \param predictedState by the predictor.
     * \param order of the integration.
     * \param doubleStep boolean if stepsize should be considered double, true for estimating doubling error.
     * \return  after corrector step
     */
    StateType performCorrectorStep( StateType predictedState, unsigned int order, bool doubleStep )
    {
        unsigned int stepsToSkip = static_cast< unsigned int>( doubleStep );
        TimeStepType stepSize = stepSize_ * static_cast< double >( stepsToSkip + 1 );
        StateType correctedState = stateHistory_.at( stepsToSkip ) + extrapolationCoefficients[ order * 2 - 1 ][ 0 ] *
                stepSize * predictedDerivative_;
        for ( unsigned int i = 1; i < order; i++ )
        {
            correctedState += stepSize * extrapolationCoefficients[ order * 2 - 1 ][ i ] *
                    derivHistory_.at( ( i - 1 ) * ( stepsToSkip + 1 ) + stepsToSkip );
        }
        return correctedState;
    }

    //! Estimate the absolute error
    /*!
     * Based on order, predicted and corrected step find the error.
     * \param predictedState by the predictor.
     * \param correctedState by the corrector.
     * \param order of the integration.
     * \return absolute error vector.
     */
    StateType estimateAbsoluteError( StateType predictedState, StateType correctedState, unsigned int order )
    {
        // Estimate the maximum truncation error
        return truncationErrorCoefficients[ order ] * ( predictedState - correctedState ).cwiseAbs( ).array( );
    }

    //! Estimate the relative error
    /*!
     * Based on order, predicted and corrected step find the error.
     * \param predictedState by the predictor.
     * \param correctedState by the corrector.
     * \param absoluteError
     * \return relative error vector.
     */
    StateType estimateRelativeError( StateType predictedState, StateType correctedState, StateType absoluteError )
    {
        // Estimate the maximum truncation error
        return absoluteError.cwiseQuotient( ( correctedState.cwiseAbs( ) ).cwiseMax( predictedState.cwiseAbs( ) ) );
    }

    //! Compare two errors
    /*!
     * Compares two errors (absolute and relative form) and return
     * true if first one is better than second (smaller is better).
     * \param absoluteError1 absolute error one.
     * \param relativeError1 relative error one.
     * \param absoluteError2 absolute error two.
     * \param relativeError2 relative error two.
     * \return true if one is better than two, false otherwise.
     */
    bool errorCompare( StateType absoluteError1, StateType relativeError1,
                       StateType absoluteError2, StateType relativeError2 )
    {
        // Find compound error
        StateType c1 = absoluteError1.cwiseMin( relativeError1 );
        StateType c2 = absoluteError2.cwiseMin( relativeError2 );
        bool oneBetter = true;
        if( strictCompare_ )
        {
            // Needs to be better or equal for each component
            for( int i = 0; i < c1.size( ); ++i )
            {
                oneBetter = oneBetter && ( c1( i ) <= c2( i ) );
            }
        } else {
            // Needs to be overal better
            oneBetter = ( c1.norm( ) <= c2.norm( ) );
        }
        return oneBetter;
    }

    //! Checks if error is too large
    /*!
     * Checks error (absolute and relative form) and return true if
     * the exceed tolerance limits.
     * \param absoluteError absolute error.
     * \param relativeError relative error.
     * \return true if one error is too big, false if within limits
     */
    bool errorTooLarge( StateType absoluteError, StateType relativeError )
    {
        bool belowLimit = true;
        // All components needs to be below the upper limit (tol)
        for( int i = 0; i < absoluteError.size( ); ++i )
        {
            belowLimit = belowLimit &&
                    ( absoluteError( i ) < absoluteErrorTolerance_( i )
                      || relativeError( i ) < relativeErrorTolerance_( i ) );
        }
        return !belowLimit;
    }

    //! Checks if error is too small
    /*!
     * Checks error (absolute and relative form) and return true if
     * the exceed tolerance limits.
     * \param absoluteError absolute error.
     * \param relativeError relative error.
     * \return true if one error is too small, false if within limits
     */
    bool errorTooSmall( StateType absoluteError, StateType relativeError )
    {
        bool belowLimit = true;
        // All components need to be above lower limit ( tol / bw )
        for( int i = 0; i < absoluteError.size( ); ++i )
        {
            belowLimit = belowLimit &&
                    ( absoluteError( i ) <= absoluteErrorTolerance_( i ) / bandwidth_
                      || relativeError( i ) <= relativeErrorTolerance_( i ) / bandwidth_ );
        }
        return belowLimit;
    }



    //! Predicted derivative state.
    /*!
     * Current state as computed by performIntegrationStep( ), after execution of
     * performPredictorStep( ).
     */
    StateDerivativeType predictedDerivative_;

    //! State history.
    /*!
     * History of states, dynamic deque size depends on order.
     */
    std::deque<StateType> stateHistory_;

    //! Derivative history.
    /*!
     * History of derivatives, dynamic deque size depends on order.
     */
    std::deque<StateType> derivHistory_;

    //! Last state.
    /*!
     * Last state as computed by performIntegrationStep( ).
     */
    StateType lastState_;

    //! Last step size.
    /*!
     * Last step size as computed by performIntegrationStep( ).
     */
    IndependentVariableType lastStepSize_;

    //! Last state derivative.
    /*!
     * Last state derivative as computed by performIntegrationStep( ).
     */
    StateDerivativeType lastDerivative_;
};

//extern template class AdamsBashforthMoultonIntegrator < double, Eigen::VectorXd, Eigen::VectorXd >;
//extern template class AdamsBashforthMoultonIntegrator < double, Eigen::Vector6d, Eigen::Vector6d >;
//extern template class AdamsBashforthMoultonIntegrator < double, Eigen::MatrixXd, Eigen::MatrixXd >;


//! Typedef of Adam-Bashforh-Moulton integrator (state/state derivative = VectorXd, independent variable = double).
/*!
 * Typedef of a Adams-Bashforth-Moulton integrator with VectorXds as state and state derivative and double as
 * independent variable.
 */
typedef AdamsBashforthMoultonIntegrator< > AdamsBashforthMoultonIntegratorXd;

//! Typedef of a scalar Adams-Bashforth-Moulton integrator.
/*!
 * Typedef of an Adams-Bashforth-Moulton integrator with doubles as state and state derivative and independent
 * variable.
 */
typedef AdamsBashforthMoultonIntegrator< double, double, double, double > AdamsBashforthMoultonIntegratord;

//! Typedef of pointer to default Adams-Bashforth-Moulton integrator
/*!
 * Typedef of pointer to a Adams-Bashforth-Moulton integrator with VectorXds as state and state derivative and double
 * as independent variable.
 */
typedef std::shared_ptr< AdamsBashforthMoultonIntegratorXd > AdamsBashforthMoultonIntegratorXdPointer;

//! Typedef of pointer to a scalar Adams-Bashforth-Moulton integrator.
/*!
 * Typedef of pointer to an Adams-Bashforth-Moulton integrator with doubles as state and state derivative and
 * independent variable.
 */
typedef std::shared_ptr< AdamsBashforthMoultonIntegratord > AdamsBashforthMoultonIntegratordPointer;

//! Truncation error coefficients ( size: 1 x o )
/*!
 * truncationErrorCoefficients( o - 1 )
 */
template< typename IndependentVariableType, typename StateType, typename StateDerivativeType, typename TimeStepType>
const double AdamsBashforthMoultonIntegrator< IndependentVariableType, StateType, StateDerivativeType, TimeStepType >::truncationErrorCoefficients[ 12 ] = {
    +2.00000000000000000e+00, +6.00000000000000000e+00, +1.00000000000000000e+01, +1.42105263157894743e+01,
    +1.85925925925925917e+01, +2.31170336037079949e+01, +2.77629090909090905e+01, +3.25146526080169664e+01,
    +3.73602765314851339e+01, +4.22902727113587602e+01, +4.72969156506348156e+01, +5.23738028029830431e+01
};

//! Extrapolation coefficients    ( size: o * 2 X o )
/*! 
 * extrapolationCoefficients( o * 2 - 2, i )                 <--- predictor coefficients
 * extrapolationCoefficients( o * 2 - 1, i )                 <--- corrector coefficients
 */
template< typename IndependentVariableType, typename StateType, typename StateDerivativeType, typename TimeStepType >
const double AdamsBashforthMoultonIntegrator< IndependentVariableType, StateType, StateDerivativeType, TimeStepType >::extrapolationCoefficients[ 24 ][ 12 ] = {
    {+1.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.50000000000000000e+00, -5.00000000000000000e-01, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+5.00000000000000000e-01, +5.00000000000000000e-01, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.91666666666666674e+00, -1.33333333333333326e+00, +4.16666666666666685e-01, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+4.16666666666666685e-01, +6.66666666666666630e-01, -8.33333333333333287e-02, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.29166666666666652e+00, -2.45833333333333348e+00, +1.54166666666666674e+00, -3.75000000000000000e-01,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.75000000000000000e-01, +7.91666666666666630e-01, -2.08333333333333343e-01, +4.16666666666666644e-02,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.64027777777777795e+00, -3.85277777777777786e+00, +3.63333333333333330e+00, -1.76944444444444438e+00,
     +3.48611111111111094e-01, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.48611111111111094e-01, +8.97222222222222254e-01, -3.66666666666666641e-01, +1.47222222222222227e-01,
     -2.63888888888888888e-02, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.97013888888888911e+00, -5.50208333333333321e+00, +6.93194444444444446e+00, -5.06805555555555554e+00,
     +1.99791666666666656e+00, -3.29861111111111105e-01, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.29861111111111105e-01, +9.90972222222222254e-01, -5.54166666666666696e-01, +3.34722222222222199e-01,
     -1.20138888888888892e-01, +1.87499999999999993e-02, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.28573082010582018e+00, -7.39563492063492056e+00, +1.16658234126984119e+01, -1.13798941798941797e+01,
     +6.73179563492063515e+00, -2.22341269841269851e+00, +3.15591931216931243e-01, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.15591931216931243e-01, +1.07658730158730154e+00, -7.68204365079365070e-01, +6.20105820105820160e-01,
     -3.34176587301587280e-01, +1.04365079365079369e-01, -1.42691798941798949e-02, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.58995535714285730e+00, -9.52520667989417902e+00, +1.80545386904761891e+01, -2.20277529761904773e+01,
     +1.73796544312169310e+01, -8.61212797619047699e+00, +2.44516369047619042e+00, -3.04224537037037057e-01,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.04224537037037057e-01, +1.15615906084656084e+00, -1.00691964285714275e+00, +1.01796461640211633e+00,
     -7.32035383597883560e-01, +3.43080357142857117e-01, -9.38409391534391485e-02, +1.13673941798941806e-02,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.88482335758377406e+00, -1.18841506834215167e+01, +2.63108427028218692e+01, -3.85403610008818376e+01,
     +3.80204144620811277e+01, -2.51247360008818355e+01, +1.07014677028218692e+01, -2.66316854056437391e+00,
     +2.94868000440917100e-01, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.94868000440917100e-01, +1.23101135361552028e+00, -1.26890266754850090e+00, +1.54193066578483240e+00,
     -1.38699294532627859e+00, +8.67046406525573188e-01, -3.55823963844797198e-01, +8.62196869488536105e-02,
     -9.35653659611992983e-03, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+4.17179880401234549e+00, -1.44669297012786604e+01, +3.66419587742504405e+01, -6.26462985008818336e+01,
     +7.41793207120811218e+01, -6.12836422508818330e+01, +3.48074052028218688e+01, -1.29942846119929456e+01,
     +2.87764701829806002e+00, -2.86975446428571423e-01, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.86975446428571423e-01, +1.30204433972663147e+00, -1.55303461199294524e+00, +2.20490520282186964e+00,
     -2.38145475088183423e+00, +1.86150821208112882e+00, -1.01879850088183432e+00, +3.70351631393298075e-01,
     -8.03895227072310425e-02, +7.89255401234567958e-03, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+4.45198840045628241e+00, -1.72688256657180261e+01, +4.92504906142275942e+01, -9.62690500741542365e+01,
     +1.33019135965307839e+02, -1.31891420554753893e+02, +9.36472204560485864e+01, -4.66170361852653485e+01,
     +1.54861788582752116e+01, -3.08887141086793848e+00, +2.80189596443936706e-01, +0.00000000000000000e+00},
    {+2.80189596443936706e-01, +1.36990283957297843e+00, -1.85839786130150708e+00, +3.01920720097803441e+00,
     -3.80648324765512269e+00, +3.57154240820907498e+00, -2.44382699765512257e+00, +1.18465362954946296e+00,
     -3.85752772015792833e-01, +7.57510538586927407e-02, -6.78584998463470715e-03, +0.00000000000000000e+00},
    {+4.72625394048788117e+00, -2.02857466060656151e+01, +6.43350953159655461e+01, -1.41522864179368099e+02,
     +2.23526764175735536e+02, -2.58602100049352657e+02, +2.20357899950647351e+02, -1.37124664395693031e+02,
     +6.07399929634890583e+01, -1.81734761126058864e+01, +3.29711053679152633e+00, -2.74265540031599087e-01},
    {+2.74265540031599087e-01, +1.43506746010869279e+00, -2.18422096398007870e+00, +3.99667650901374838e+00,
     -5.76142186372655107e+00, +6.30845647070907489e+00, -5.18074106015512292e+00, +3.13959224562089156e+00,
     -1.36322208005150713e+00, +4.01574156537264193e-01, -7.19504705203489886e-02, +5.92405641233766274e-03}
};


//! Interpolation coefficients    ( size: o^2-o X 2 * o )
/*! 
 * interpolationCoefficients( ( o - 1 )^2 - o + 1 + j, i )     <--- state coefficients for state history
 * interpolationCoefficients( ( o - 1 )^2 - o + 1 + j, o + i ) <--- state coefficients for deriv history
 * interpolationCoefficients( ( o - 1 )^2 + j, i )             <--- deriv coefficients for state history
 * interpolationCoefficients( ( o - 1 )^2 + j, o + i )         <--- deriv coefficients for deriv history
 *
 * Order k occupies : registers ( k - 1 )^2 - k + 1 to k^2 - k - 1
 */
template< typename IndependentVariableType, typename StateType, typename StateDerivativeType, typename TimeStepType >
const double AdamsBashforthMoultonIntegrator< IndependentVariableType, StateType, StateDerivativeType, TimeStepType >::interpolationCoefficients[ 132 ][ 24 ] = {
    {+5.00000000000000000e-01, +5.00000000000000000e-01, -1.25000000000000000e-01, +1.25000000000000000e-01,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.50000000000000000e+00, -1.50000000000000000e+00, -2.50000000000000000e-01, -2.50000000000000000e-01,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.51562500000000000e-01, +5.62500000000000000e-01, +8.59375000000000000e-02, -7.03125000000000000e-02,
     +2.81250000000000000e-01, +2.34375000000000000e-02, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+8.59375000000000000e-02, +5.62500000000000000e-01, +3.51562500000000000e-01, -2.34375000000000000e-02,
     -2.81250000000000000e-01, +7.03125000000000000e-02, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.45312500000000000e+00, -1.50000000000000000e+00, +4.68750000000000000e-02, -2.34375000000000000e-01,
     -1.87500000000000000e-01, +1.56250000000000000e-02, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-4.68750000000000000e-02, +1.50000000000000000e+00, -1.45312500000000000e+00, +1.56250000000000000e-02,
     -1.87500000000000000e-01, -2.34375000000000000e-01, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.76692708333333315e-01, +4.39453125000000000e-01, +2.44140625000000000e-01, +3.97135416666666644e-02,
     -4.88281250000000000e-02, +4.39453125000000000e-01, +1.46484375000000000e-01, +9.76562500000000000e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.53906250000000000e-02, +4.74609375000000000e-01, +4.74609375000000000e-01, +2.53906250000000000e-02,
     -5.85937500000000000e-03, -1.58203125000000000e-01, +1.58203125000000000e-01, +5.85937500000000000e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.97135416666666644e-02, +2.44140625000000000e-01, +4.39453125000000000e-01, +2.76692708333333315e-01,
     -9.76562500000000000e-03, -1.46484375000000000e-01, -4.39453125000000000e-01, +4.88281250000000000e-02,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.33897569444444442e+00, -1.69921875000000000e+00, +2.92968750000000000e-01, +6.72743055555555525e-02,
     -2.01822916666666657e-01, +5.85937500000000000e-02, +2.14843750000000000e-01, +1.69270833333333322e-02,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.95312500000000000e-02, +1.58203125000000000e+00, -1.58203125000000000e+00, -1.95312500000000000e-02,
     -3.90625000000000000e-03, -3.16406250000000000e-01, -3.16406250000000000e-01, -3.90625000000000000e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-6.72743055555555525e-02, -2.92968750000000000e-01, +1.69921875000000000e+00, -1.33897569444444442e+00,
     +1.69270833333333322e-02, +2.14843750000000000e-01, +5.85937500000000000e-02, -2.01822916666666657e-01,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.30534871419270843e-01, +1.99381510416666657e-01, +2.99072265625000000e-01, +2.47233072916666657e-01,
     +2.37782796223958322e-02, -3.73840332031250000e-02, +5.98144531250000000e-01, +4.48608398437500000e-01,
     +1.19628906250000000e-01, +5.34057617187500000e-03, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.10626220703125000e-02, +4.02832031250000000e-01, +4.94384765625000000e-01, +8.54492187500000000e-02,
     +6.27136230468750000e-03, -2.28881835937500000e-03, -1.09863281250000000e-01, +2.47192382812500000e-01,
     +3.66210937500000000e-02, +1.37329101562500000e-03, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+6.27136230468750000e-03, +8.54492187500000000e-02, +4.94384765625000000e-01, +4.02832031250000000e-01,
     +1.10626220703125000e-02, -1.37329101562500000e-03, -3.66210937500000000e-02, -2.47192382812500000e-01,
     +1.09863281250000000e-01, +2.28881835937500000e-03, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.37782796223958322e-02, +2.47233072916666657e-01, +2.99072265625000000e-01, +1.99381510416666657e-01,
     +2.30534871419270843e-01, -5.34057617187500000e-03, -1.19628906250000000e-01, -4.48608398437500000e-01,
     -5.98144531250000000e-01, +3.73840332031250000e-02, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.23414781358506942e+00, -2.25206163194444464e+00, +4.10156250000000000e-01, +5.50672743055555580e-01,
     +5.70848253038194475e-02, -1.75882975260416657e-01, +4.21549479166666685e-01, +9.14306640625000000e-01,
     +2.75716145833333315e-01, +1.29191080729166661e-02, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.72424316406250000e-02, +1.56738281250000000e+00, -1.58203125000000000e+00, -4.88281250000000000e-03,
     +2.28881835937500000e-03, -3.35693359375000000e-03, -3.07617187500000000e-01, -2.96630859375000000e-01,
     +4.88281250000000000e-03, +5.49316406250000000e-04, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-2.28881835937500000e-03, +4.88281250000000000e-03, +1.58203125000000000e+00, -1.56738281250000000e+00,
     -1.72424316406250000e-02, +5.49316406250000000e-04, +4.88281250000000000e-03, -2.96630859375000000e-01,
     -3.07617187500000000e-01, -3.35693359375000000e-03, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-5.70848253038194475e-02, -5.50672743055555580e-01, -4.10156250000000000e-01, +2.25206163194444464e+00,
     -1.23414781358506942e+00, +1.29191080729166661e-02, +2.75716145833333315e-01, +9.14306640625000000e-01,
     +4.21549479166666685e-01, -1.75882975260416657e-01, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.98845672607421864e-01, -1.26171112060546875e-01, +0.00000000000000000e+00, +6.45996093750000000e-01,
     +2.65216827392578125e-01, +1.61125183105468757e-02, -3.02810668945312500e-02, +7.57026672363281250e-01,
     +1.00936889648437500e+00, +6.05621337890625000e-01, +1.08146667480468750e-01, +3.36456298828125000e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+5.86929321289062465e-03, +3.50475311279296875e-01, +4.48608398437500000e-01, +1.49536132812500000e-01,
     +4.31785583496093750e-02, +2.33230590820312491e-03, -1.12152099609375000e-03, -8.41140747070312500e-02,
     +3.36456298828125000e-01, +1.12152099609375000e-01, +1.68228149414062500e-02, +4.80651855468750000e-04,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.70516967773437500e-03, +4.05311584472656250e-02, +4.57763671875000000e-01, +4.57763671875000000e-01,
     +4.05311584472656250e-02, +1.70516967773437500e-03, -3.43322753906250000e-04, -1.43051147460937500e-02,
     -1.71661376953125000e-01, +1.71661376953125000e-01, +1.43051147460937500e-02, +3.43322753906250000e-04,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.33230590820312491e-03, +4.31785583496093750e-02, +1.49536132812500000e-01, +4.48608398437500000e-01,
     +3.50475311279296875e-01, +5.86929321289062465e-03, -4.80651855468750000e-04, -1.68228149414062500e-02,
     -1.12152099609375000e-01, -3.36456298828125000e-01, +8.41140747070312500e-02, +1.12152099609375000e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.61125183105468757e-02, +2.65216827392578125e-01, +6.45996093750000000e-01, +0.00000000000000000e+00,
     -1.26171112060546875e-01, +1.98845672607421864e-01, -3.36456298828125000e-03, -1.08146667480468750e-01,
     -6.05621337890625000e-01, -1.00936889648437500e+00, -7.57026672363281250e-01, +3.02810668945312500e-02,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.14502166748046874e+00, -3.17310333251953125e+00, -4.48608398437500000e-01, +1.67907714843750000e+00,
     +7.50617980957031250e-01, +4.69949340820312519e-02, -1.55923461914062506e-01, +8.69979858398437500e-01,
     +2.50579833984375000e+00, +1.66497802734375000e+00, +3.09677124023437500e-01, +9.84802246093749965e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.24606323242187501e-02, +1.51805877685546875e+00, -1.62780761718750000e+00, +5.55419921875000000e-02,
     +3.92532348632812500e-02, +2.49298095703124993e-03, -2.28576660156250017e-03, -2.83584594726562500e-01,
     -2.11486816406250000e-01, +7.90405273437500000e-02, +1.63421630859375000e-02, +5.21850585937500043e-04,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+7.36999511718749957e-04, +3.33786010742187500e-02, +1.60217285156250000e+00, -1.60217285156250000e+00,
     -3.33786010742187500e-02, -7.36999511718749957e-04, -1.37329101562500000e-04, -9.53674316406250000e-03,
     -3.43322753906250000e-01, -3.43322753906250000e-01, -9.53674316406250000e-03, -1.37329101562500000e-04,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-2.49298095703124993e-03, -3.92532348632812500e-02, -5.55419921875000000e-02, +1.62780761718750000e+00,
     -1.51805877685546875e+00, -1.24606323242187501e-02, +5.21850585937500043e-04, +1.63421630859375000e-02,
     +7.90405273437500000e-02, -2.11486816406250000e-01, -2.83584594726562500e-01, -2.28576660156250017e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-4.69949340820312519e-02, -7.50617980957031250e-01, -1.67907714843750000e+00, +4.48608398437500000e-01,
     +3.17310333251953125e+00, -1.14502166748046874e+00, +9.84802246093749965e-03, +3.09677124023437500e-01,
     +1.66497802734375000e+00, +2.50579833984375000e+00, +8.69979858398437500e-01, -1.55923461914062506e-01,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.75567102432250988e-01, -5.19067955017089799e-01, -9.54169034957885742e-01, +8.14224243164062500e-01,
     +1.18784308433532715e+00, +2.83847618103027333e-01, +1.17549419403076179e-02, -2.54445075988769531e-02,
     +9.16002273559570312e-01, +1.90833806991577148e+00, +2.03556060791015625e+00, +8.17859172821044922e-01,
     +1.01778030395507812e-01, +2.31313705444335938e-03, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.51176261901855460e-03, +3.11137962341308583e-01, +3.54856252670288086e-01, +1.68228149414062500e-01,
     +1.33425951004028320e-01, +2.77627944946289076e-02, +1.07712745666503911e-03, -6.30855560302734375e-04,
     -6.81324005126953125e-02, +4.25827503204345703e-01, +2.52342224121093750e-01, +8.51655006408691406e-02,
     +9.73320007324218750e-03, +2.10285186767578125e-04, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+6.19173049926757812e-04, +2.26640701293945312e-02, +4.16189432144165039e-01, +4.67300415039062500e-01,
     +8.03172588348388672e-02, +1.24769210815429688e-02, +4.32729721069335938e-04, -1.16825103759765625e-04,
     -7.00950622558593750e-03, -1.31428241729736328e-01, +2.33650207519531250e-01, +4.38094139099121094e-02,
     +4.20570373535156250e-03, +8.34465026855468750e-05, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+4.32729721069335938e-04, +1.24769210815429688e-02, +8.03172588348388672e-02, +4.67300415039062500e-01,
     +4.16189432144165039e-01, +2.26640701293945312e-02, +6.19173049926757812e-04, -8.34465026855468750e-05,
     -4.20570373535156250e-03, -4.38094139099121094e-02, -2.33650207519531250e-01, +1.31428241729736328e-01,
     +7.00950622558593750e-03, +1.16825103759765625e-04, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.07712745666503911e-03, +2.77627944946289076e-02, +1.33425951004028320e-01, +1.68228149414062500e-01,
     +3.54856252670288086e-01, +3.11137962341308583e-01, +3.51176261901855460e-03, -2.10285186767578125e-04,
     -9.73320007324218750e-03, -8.51655006408691406e-02, -2.52342224121093750e-01, -4.25827503204345703e-01,
     +6.81324005126953125e-02, +6.30855560302734375e-04, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.17549419403076179e-02, +2.83847618103027333e-01, +1.18784308433532715e+00, +8.14224243164062500e-01,
     -9.54169034957885742e-01, -5.19067955017089799e-01, +1.75567102432250988e-01, -2.31313705444335938e-03,
     -1.01778030395507812e-01, -8.17859172821044922e-01, -2.03556060791015625e+00, -1.90833806991577148e+00,
     -9.16002273559570312e-01, +2.54445075988769531e-02, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.06965185165405274e+00, -4.44927726745605501e+00, -3.56388330459594727e+00, +2.20886230468750000e+00,
     +3.76655817031860352e+00, +9.29008712768554679e-01, +3.90795326232910162e-02, -1.40271568298339838e-01,
     +1.38576736450195304e+00, +5.43146610260009766e+00, +6.33638000488281250e+00, +2.63933658599853516e+00,
     +3.34912872314453103e-01, +7.70511627197265590e-03, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+8.99847030639648431e-03, +1.45979255676269526e+00, -1.76864862442016602e+00, +8.11767578125000000e-02,
     +1.75287723541259766e-01, +4.16869354248046858e-02, +1.70618057250976557e-03, -1.56612396240234379e-03,
     -2.59984588623046853e-01, -7.84063339233398438e-02, +2.89993286132812500e-01, +1.20583534240722656e-01,
     +1.48933410644531243e-02, +3.35121154785156228e-04, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+6.20174407958984332e-04, +3.11756134033203125e-02, +1.59591436386108398e+00, -1.60217285156250000e+00,
     -2.71201133728027344e-02, +1.46598815917968754e-03, +1.16825103759765625e-04, -1.13487243652343750e-04,
     -8.67843627929687500e-03, -3.37958335876464844e-01, -3.33786010742187500e-01, -4.17232513427734375e-03,
     +7.20977783203125000e-04, +2.38418579101562500e-05, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-1.16825103759765625e-04, -1.46598815917968754e-03, +2.71201133728027344e-02, +1.60217285156250000e+00,
     -1.59591436386108398e+00, -3.11756134033203125e-02, -6.20174407958984332e-04, +2.38418579101562500e-05,
     +7.20977783203125000e-04, -4.17232513427734375e-03, -3.33786010742187500e-01, -3.37958335876464844e-01,
     -8.67843627929687500e-03, -1.13487243652343750e-04, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-1.70618057250976557e-03, -4.16869354248046858e-02, -1.75287723541259766e-01, -8.11767578125000000e-02,
     +1.76864862442016602e+00, -1.45979255676269526e+00, -8.99847030639648431e-03, +3.35121154785156228e-04,
     +1.48933410644531243e-02, +1.20583534240722656e-01, +2.89993286132812500e-01, -7.84063339233398438e-02,
     -2.59984588623046853e-01, -1.56612396240234379e-03, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-3.90795326232910162e-02, -9.29008712768554679e-01, -3.76655817031860352e+00, -2.20886230468750000e+00,
     +3.56388330459594727e+00, +4.44927726745605501e+00, -1.06965185165405274e+00, +7.70511627197265590e-03,
     +3.34912872314453103e-01, +2.63933658599853516e+00, +6.33638000488281250e+00, +5.43146610260009766e+00,
     +1.38576736450195304e+00, -1.40271568298339838e-01, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.57650237424033030e-01, -9.67527401447296098e-01, -2.90258220434188852e+00, -5.37515223026275635e-01,
     +3.01666706800460815e+00, +1.92311002016067500e+00, +3.01186215877532970e-01, +9.01128734861101482e-03,
     -2.19393968582153320e-02, +1.07503044605255127e+00, +3.22509133815765381e+00, +5.37515223026275635e+00,
     +3.83939445018768311e+00, +1.07503044605255127e+00, +9.77300405502319336e-02, +1.68764591217041016e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.27924925940377375e-03, +2.80525696277618419e-01, +2.23275554180145275e-01, +7.95140862464904785e-02,
     +2.57625639438629150e-01, +1.36348807811737055e-01, +1.98608517646789544e-02, +5.70115021296909927e-04,
     -3.89456748962402344e-04, -5.72501420974731445e-02, +5.15251278877258301e-01, +4.77084517478942871e-01,
     +2.86250710487365723e-01, +7.36073255538940430e-02, +6.36112689971923828e-03, +1.06215476989746094e-04,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.69676957811628077e-04, +1.40628218650817871e-02, +3.79696190357208252e-01, +4.43570315837860107e-01,
     +1.14999711513519287e-01, +4.18730378150939941e-02, +5.38319349288940430e-03, +1.45052160535539892e-04,
     -4.82797622680664062e-05, -3.94284725189208984e-03, -1.06456875801086426e-01, +2.95713543891906738e-01,
     +9.85711812973022461e-02, +2.12913751602172852e-02, +1.68979167938232422e-03, +2.68220901489257812e-05,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.14142894744873047e-04, +4.72265481948852539e-03, +4.79421019554138184e-02, +4.47221100330352783e-01,
     +4.47221100330352783e-01, +4.79421019554138184e-02, +4.72265481948852539e-03, +1.14142894744873047e-04,
     -2.08616256713867188e-05, -1.43110752105712891e-03, -2.14666128158569336e-02, -1.78888440132141113e-01,
     +1.78888440132141113e-01, +2.14666128158569336e-02, +1.43110752105712891e-03, +2.08616256713867188e-05,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.45052160535539892e-04, +5.38319349288940430e-03, +4.18730378150939941e-02, +1.14999711513519287e-01,
     +4.43570315837860107e-01, +3.79696190357208252e-01, +1.40628218650817871e-02, +2.69676957811628077e-04,
     -2.68220901489257812e-05, -1.68979167938232422e-03, -2.12913751602172852e-02, -9.85711812973022461e-02,
     -2.95713543891906738e-01, +1.06456875801086426e-01, +3.94284725189208984e-03, +4.82797622680664062e-05,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+5.70115021296909927e-04, +1.98608517646789544e-02, +1.36348807811737055e-01, +2.57625639438629150e-01,
     +7.95140862464904785e-02, +2.23275554180145275e-01, +2.80525696277618419e-01, +2.27924925940377375e-03,
     -1.06215476989746094e-04, -6.36112689971923828e-03, -7.36073255538940430e-02, -2.86250710487365723e-01,
     -4.77084517478942871e-01, -5.15251278877258301e-01, +5.72501420974731445e-02, +3.89456748962402344e-04,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+9.01128731746033036e-03, +3.01186215023935555e-01, +1.92311001601042575e+00, +3.01666706432530907e+00,
     -5.37515219346976547e-01, -2.90258220019163904e+00, -9.67527400593698794e-01, +1.57650237455183723e-01,
     -1.68764590616339127e-03, -9.77300402558880094e-02, -1.07503044340345588e+00, -3.83939444282908493e+00,
     -5.37515222290415817e+00, -3.22509133550855864e+00, -1.07503044575820739e+00, +2.19393968642223498e-02,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.00536631535510623e+00, -6.06153930187225320e+00, -1.05877360868453980e+01, -2.69861400127410889e+00,
     +1.03499573469161987e+01, +6.86686347484588655e+00, +1.09270061016082765e+00, +3.30016427137413831e-02,
     -1.27699027742658344e-01, +1.95713057518005362e+00, +1.01715135097503655e+01, +1.83858964443206787e+01,
     +1.35715711116790771e+01, +3.86829581260681143e+00, +3.55611944198608421e-01, +6.18807247706822019e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+6.66023653380724834e-03, +1.40146034717559820e+00, -2.02092505931854260e+00, -9.17452573776245117e-02,
     +4.12647128105163574e-01, +2.52207942008972175e-01, +3.85592603683471702e-02, +1.13540250427868910e-03,
     -1.10846246991838735e-03, -2.39277505874633784e-01, +9.24924373626709040e-02, +7.21753835678100586e-01,
     +5.09385824203491211e-01, +1.39397192001342762e-01, +1.24505519866943363e-02, +2.12185723440987714e-04,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.89553089531100526e-04, +2.54136323928833008e-02, +1.57094299793243408e+00, -1.61938369274139404e+00,
     -3.65078449249267578e-03, +2.23818540573120124e-02, +3.79264354705810547e-03, +1.13796214668118221e-04,
     -6.83580126081194216e-05, -6.63399696350097656e-03, -3.21060419082641602e-01, -2.91019678115844727e-01,
     +3.44216823577880859e-02, +1.31127834320068359e-02, +1.23381614685058594e-03, +2.12873731340680818e-05,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.43152454921177426e-05, +2.11803913116455095e-03, +4.15021181106567383e-02, +1.60999596118927002e+00,
     -1.60999596118927002e+00, -4.15021181106567383e-02, -2.11803913116455095e-03, -3.43152454921177426e-05,
     -5.96046447753906250e-06, -5.72443008422851562e-04, -1.43110752105712891e-02, -3.57776880264282227e-01,
     -3.57776880264282227e-01, -1.43110752105712891e-02, -5.72443008422851562e-04, -5.96046447753906250e-06,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-1.13796214668118221e-04, -3.79264354705810547e-03, -2.23818540573120124e-02, +3.65078449249267578e-03,
     +1.61938369274139404e+00, -1.57094299793243408e+00, -2.54136323928833008e-02, -3.89553089531100526e-04,
     +2.12873731340680818e-05, +1.23381614685058594e-03, +1.31127834320068359e-02, +3.44216823577880859e-02,
     -2.91019678115844727e-01, -3.21060419082641602e-01, -6.63399696350097656e-03, -6.83580126081194216e-05,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-1.13540250427868910e-03, -3.85592603683471702e-02, -2.52207942008972175e-01, -4.12647128105163574e-01,
     +9.17452573776245117e-02, +2.02092505931854260e+00, -1.40146034717559820e+00, -6.66023653380724834e-03,
     +2.12185723440987714e-04, +1.24505519866943363e-02, +1.39397192001342762e-01, +5.09385824203491211e-01,
     +7.21753835678100586e-01, +9.24924373626709040e-02, -2.39277505874633784e-01, -1.10846246991838735e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-3.30016427137413831e-02, -1.09270061016082765e+00, -6.86686347484588655e+00, -1.03499573469161987e+01,
     +2.69861400127410889e+00, +1.05877360868453980e+01, +6.06153930187225320e+00, -1.00536631535510623e+00,
     +6.18807247706822019e-03, +3.55611944198608421e-01, +3.86829581260681143e+00, +1.35715711116790771e+01,
     +1.83858964443206787e+01, +1.01715135097503655e+01, +1.95713057518005362e+00, -1.27699027742658344e-01,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.43380447240945480e-01, -1.46327941545418327e+00, -6.21501976624131203e+00, -6.04704625904560089e+00,
     +3.85653460398316383e+00, +7.54014410078525543e+00, +2.86110246554017067e+00, +3.17024748240198384e-01,
     +7.15907495136239736e-03, -1.92826730199158192e-02, +1.23409107327461243e+00, +5.03920521587133408e+00,
     +1.20940925180912018e+01, +1.34978711139410734e+01, +6.71894028782844543e+00, +1.37432869523763657e+00,
     +9.49300825595855713e-02, +1.28551153466105461e-03, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.56893622063632512e-03, +2.55985748342105301e-01, +6.04704625904560061e-02, -1.88130328059196467e-01,
     +3.02352312952280045e-01, +4.09718236327171303e-01, +1.42590843886137020e-01, +1.51121518441608982e-02,
     +3.31635896249541208e-04, -2.57102306932210922e-04, -4.93636429309844971e-02, +6.04704625904560089e-01,
     +8.06272834539413452e-01, +7.55880782380700111e-01, +3.45545500516891479e-01, +6.71894028782844543e-02,
     +4.48760390281677246e-03, +5.93313015997409821e-05, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.33169432436781274e-04, +9.37709850924355569e-03, +3.48868053406476974e-01, +3.93594726920127869e-01,
     +1.24240759760141373e-01, +9.30314809083938599e-02, +2.79313512146472931e-02, +2.76509140219007238e-03,
     +5.82684463422213263e-05, -2.28197313845157623e-05, -2.43410468101501465e-03, -8.94533470273017883e-02,
     +3.57813388109207153e-01, +1.86361139640212059e-01, +7.15626776218414307e-02, +1.27790495753288269e-02,
     +8.11368227005004883e-04, +1.03726051747798920e-05, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.77657124772667885e-05, +2.12068855762481689e-03, +3.09924222528934479e-02, +4.20208945870399475e-01,
     +4.52811364084482193e-01, +7.56698101758956909e-02, +1.66634581983089447e-02, +1.46649777889251709e-03,
     +2.90473690256476402e-05, -6.60074874758720398e-06, -5.91427087783813477e-04, -1.20749697089195251e-02,
     -1.44899636507034302e-01, +2.26405682042241096e-01, +4.82998788356781006e-02, +7.24498182535171509e-03,
     +4.22447919845581055e-04, +5.13391569256782532e-06, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.90473690001379217e-05, +1.46649777793568482e-03, +1.66634581913182762e-02, +7.56698101626502168e-02,
     +4.52811364084482193e-01, +4.20208945883644935e-01, +3.09924222598841163e-02, +2.12068855858164916e-03,
     +3.77657125027765071e-05, -5.13391568787484181e-06, -4.22447919545230110e-04, -7.24498182167241603e-03,
     -4.82998788209609009e-02, -2.26405682019245491e-01, +1.44899636521751501e-01, +1.20749697125988251e-02,
     +5.91427088084164367e-04, +6.60074875228018748e-06, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+5.82684497256277760e-05, +2.76509152693440053e-03, +2.79313521090848935e-02, +9.30314825611349983e-02,
     +1.24240759714150134e-01, +3.93594725261499856e-01, +3.48868052553247543e-01, +9.37709839477981165e-03,
     +1.33169429442758357e-04, -1.03726057989466970e-05, -8.11368266350978519e-04, -1.27790500499584056e-02,
     -7.15626794909253600e-02, -1.86361142514664457e-01, -3.57813389919422298e-01, +8.94533465821066021e-02,
     +2.43410464527325223e-03, +2.28197308354366933e-05, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.31634927390439809e-04, +1.51121168888332189e-02, +1.42590598977642286e-01, +4.09717797255810423e-01,
     +3.02352339397242043e-01, -1.88129885602880409e-01, +6.04706838042647771e-02, +2.55985777386097346e-01,
     +1.56893696559988840e-03, -5.93311222292188100e-05, -4.48759281055607136e-03, -6.71892715972146892e-02,
     -3.45544993200419870e-01, -7.55880016603587612e-01, -8.06272361072493227e-01, -6.04704511548265988e-01,
     +4.93636519508236851e-02, +2.57102443127284913e-04, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+7.15908796620545740e-03, +3.17025211999806211e-01, +2.86110567416673822e+00, +7.54014976423758032e+00,
     +3.85653417345918248e+00, -6.04705197809949357e+00, -6.21502258893001258e+00, -1.46327978269364678e+00,
     +1.43380437893639523e-01, -1.28551394926207157e-03, -9.49302302298302009e-02, -1.37433042467376332e+00,
     -6.71894690501837122e+00, -1.34978810099907864e+01, -1.20940985839749562e+01, -5.03920666930116212e+00,
     -1.23409118707608334e+00, +1.92826713130167145e-02, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+9.49916392436479073e-01, -7.99052355788191981e+00, -2.34983395263552666e+01, -2.42317339718341813e+01,
     +1.35587006807327271e+01, +2.88107025921344757e+01, +1.11282700225710869e+01, +1.24472418184183087e+00,
     +2.82831863547694309e-02, -1.17377519740590033e-01, +2.57579697029931198e+00, +1.72367779165506363e+01,
     +4.45933583378791809e+01, +5.13119869865477085e+01, +2.59685662388801575e+01, +5.36728061735630035e+00,
     +3.73393765517643528e-01, +5.08274337542908520e-03, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+5.06248598747259526e-03, +1.34611425358421943e+00, -2.38935719221830389e+00, -7.00841332674026485e-01,
     +5.10129332542419434e-01, +8.73781739473342922e-01, +3.19647338092327093e-01, +3.46909550501375771e-02,
     +7.72420162411064531e-04, -8.10866829540048345e-04, -2.21504621846335265e-01, +2.94613113999366738e-01,
     +1.46784793138504033e+00, +1.57767564430832863e+00, +7.60714066028594926e-01, +1.52182617783546442e-01,
     +1.03456480162484304e-02, +1.38441020888941616e-04, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.40627411106715387e-04, +2.02690186549206175e-02, +1.53681624680757523e+00, -1.67547538876533508e+00,
     +5.73694705963134766e-03, +7.96751797199249240e-02, +2.95079872012138367e-02, +3.15993720171402898e-03,
     +6.94447092483846536e-05, -4.06079260366303590e-05, -4.98060669217790901e-03, -3.02308425307273865e-01,
     -2.22019851207733154e-01, +1.32846180349588394e-01, +7.00963139533996582e-02, +1.39776617288589478e-02,
     +9.38986028943743071e-04, +1.24231779149600433e-05, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.81138012983969285e-05, +1.88543200492858882e-03, +3.98026779294013977e-02, +1.60677596926689148e+00,
     -1.60999596118927002e+00, -3.82821261882781982e-02, -4.18598949909210216e-04, +1.98291880743844165e-04,
     +6.20144419372081757e-06, -4.81959432363510132e-06, -4.99427318572998047e-04, -1.34166330099105835e-02,
     -3.54199111461639404e-01, -3.52186616510152817e-01, -1.07333064079284668e-02, +3.21999192237854004e-04,
     +6.70552253723144531e-05, +1.14087015390396118e-06, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-6.20144419372081757e-06, -1.98291880743844165e-04, +4.18598949909210216e-04, +3.82821261882781982e-02,
     +1.60999596118927002e+00, -1.60677596926689148e+00, -3.98026779294013977e-02, -1.88543200492858882e-03,
     -2.81138012983969285e-05, +1.14087015390396118e-06, +6.70552253723144531e-05, +3.21999192237854004e-04,
     -1.07333064079284668e-02, -3.52186616510152817e-01, -3.54199111461639404e-01, -1.34166330099105835e-02,
     -4.99427318572998047e-04, -4.81959432363510132e-06, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-6.94447083810542626e-05, -3.15993716918173137e-03, -2.95079869635311159e-02, -7.96751792695787192e-02,
     -5.73694705963134766e-03, +1.67547538831498888e+00, -1.53681624704525799e+00, -2.02690186874529138e-02,
     -2.40627411974045791e-04, +1.24231777553986042e-05, +9.38986018731810968e-04, +1.39776616037627800e-02,
     +7.00963134530149801e-02, +1.32846179567737338e-01, -2.22019851708117832e-01, -3.02308425432370020e-01,
     -4.98060670238984100e-03, -4.06079261961917964e-05, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-7.72420692276723801e-04, -3.46909745786555881e-02, -3.19647478056542644e-01, -8.73781997951460498e-01,
     -5.10129325183821258e-01, +7.00841592094044619e-01, +2.38935732558921510e+00, -1.34611423570059485e+00,
     -5.06248551990830329e-03, +1.38441118643787826e-04, +1.03456541764462857e-02, +1.52182692068594472e-01,
     +7.60714358459284257e-01, +1.57767609387268259e+00, +1.46784821439672397e+00, +2.94613183574912008e-01,
     -2.21504616262811233e-01, -8.10866743799239903e-04, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-2.82831773131981754e-02, -1.24472385303626520e+00, -1.11282676989943887e+01, -2.88106983785747026e+01,
     -1.35587008823583162e+01, +2.42317297324663343e+01, +2.34983373834350964e+01, +7.99052327414643315e+00,
     -9.49916399770992959e-01, +5.08274170366634654e-03, +3.73393661410600097e-01, +5.36727937651269205e+00,
     +2.59685614098442912e+01, +5.13119796455642359e+01, +4.45933537669240678e+01, +1.72367768047474037e+01,
     +2.57579688199313095e+00, -1.17377521083168138e-01, +0.00000000000000000e+00, +0.00000000000000000e+00,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.31713967173940516e-01, -2.00019855902480348e+00, -1.12869137039940277e+01, -2.02268098248168826e+01,
     -4.45815400220453739e+00, +1.88783558364957571e+01, +1.56130912201479077e+01, +4.01158489885606961e+00,
     +3.31486748477410775e-01, +5.84341888916530862e-03, -1.71996682183817029e-02, +1.39317312568891793e+00,
     +7.43025667034089565e+00, +2.42721717897802591e+01, +3.90088475192897022e+01, +3.03402147372253239e+01,
     +1.10328053589910269e+01, +1.71467461623251438e+00, +9.28782083792611957e-02, +1.01174518931657076e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.12921532404254240e-03, +2.35833977702506165e-01, -1.28918442476008616e-01, -7.13887405581772327e-01,
     +0.00000000000000000e+00, +8.33011474460363388e-01, +6.11236928962171078e-01, +1.49407705425151749e-01,
     +1.19797477340658331e-02, +2.06798449480196560e-04, -1.78543268702924252e-04, -4.33860142948105931e-02,
     +6.94176228716969490e-01, +1.25980130396783352e+00, +1.70073176035657525e+00, +1.21480840025469661e+00,
     +4.19933767989277840e-01, +6.31069298833608627e-02, +3.33738571498543024e-03, +3.57086537405848503e-05,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+7.21069190481175067e-05, +6.59207254232439636e-03, +3.22847039704876271e-01, +3.21949222125113010e-01,
     +8.39867535978555679e-02, +1.51176156476140022e-01, +9.11284843459725380e-02, +2.06362010378922725e-02,
     +1.58535672068995021e-03, +2.66065300878373866e-05, -1.19028845801949501e-05, -1.60688941832631826e-03,
     -7.71306920796632767e-02, +4.19933767989277840e-01, +3.14950325991958380e-01, +1.88970195595175028e-01,
     +5.99905382841825485e-02, +8.57007689774036407e-03, +4.38242568634450436e-04, +4.57803253084421158e-06,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.46516605229419104e-05, +1.07224212843396523e-03, +2.12629185989909389e-02, +3.93677554093427362e-01,
     +4.38321400433870745e-01, +9.74047556518588542e-02, +3.97736085577538728e-02, +7.89562705900142409e-03,
     +5.68081013621456730e-04, +9.16080251845303859e-06, -2.46509443969982938e-06, -2.79541709461491358e-04,
     -7.45444558562725467e-03, -1.21755944565040741e-01, +2.73950875270881744e-01, +9.13169584234739540e-02,
     +2.43511889128855061e-02, +3.19476239381343485e-03, +1.55300949699003474e-04, +1.56869646160346215e-06,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+7.54705833776709430e-06, +4.97425283443287194e-04, +7.74188549709179351e-03, +5.16204934361104151e-02,
     +4.40132644468251566e-01, +4.40132647442229019e-01, +5.16204975136459049e-02, +7.74188682571564681e-03,
     +4.97425414173448916e-04, +7.54706100116315168e-06, -1.28347866460781303e-06, -1.33665141650248187e-04,
     -2.99409931118706360e-03, -2.71686801367499384e-02, -1.83388598693741689e-01, +1.83388606133284393e-01,
     +2.71686834448486049e-02, +2.99409991939772469e-03, +1.33665179719730297e-04, +1.28347913552842859e-06,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+9.16085703518856217e-06, +5.68083579031833591e-04, +7.89565200382255118e-03, +3.97736812713868973e-02,
     +9.74048033643651889e-02, +4.38321344923090050e-01, +3.93677485056865339e-01, +2.12628972100434056e-02,
     +1.07224011323401843e-03, +1.46516211255000619e-05, -1.56870616562147584e-06, -1.55301705081626842e-04,
     -3.19477401412863928e-03, -2.43512497713536391e-02, -9.13170902155063802e-02, -2.73951002190142512e-01,
     +1.21755890208507273e-01, +7.45443595870371942e-03, +2.79541128894849613e-04, +2.46508751888174960e-06,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.65995379717739779e-05, +1.58503712011004482e-03, +2.06331825569399142e-02, +9.11199656328655627e-02,
     +1.51170926673912714e-01, +8.39934038262308419e-02, +3.21957073649479231e-01, +3.22849408864297027e-01,
     +6.59229102844777117e-03, +7.21111097450849963e-05, -4.57678197187676538e-06, -4.38147744308884008e-04,
     -8.56865445081999612e-03, -5.99832655994326294e-02, -1.88954802442480030e-01, -3.14935819852499299e-01,
     -4.19927681202097058e-01, +7.71317495294168709e-02, +1.60695204719183988e-03, +1.19036186353955091e-05,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.06882836166528652e-04, +1.19835852943434556e-02, +1.49443798569260239e-01, +6.11338426155645043e-01,
     +8.33073447839713377e-01, -7.93216140181930018e-05, -7.13980773757024068e-01, -1.28946592804079746e-01,
     +2.35831381960086239e-01, +1.12916551990715000e-03, -3.57237607264315583e-05, -3.33852575839872334e-03,
     -6.31239659301457978e-02, -4.20020610432565944e-01, -1.21499180396397910e+00, -1.70090433473940239e+00,
     -1.25987365043702937e+00, -6.94188792328895854e-01, +4.33852701847605507e-02, +1.78534544096153842e-04,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+5.84232327853112781e-03, +3.31438515087110308e-01, +4.01114489050686895e+00, +1.56118927500698188e+01,
     +1.88776699814148898e+01, -4.45720040693400854e+00, -2.02257358301221153e+01, -1.12865968965815924e+01,
     -2.00016983703712814e+00, +1.31714510317624184e-01, -1.01154795709039610e-03, -9.28637583993429799e-02,
     -1.71446440219759011e+00, -1.10317593994716230e+01, -3.03380532882243017e+01, -3.90068529972618023e+01,
     -2.42713501565512004e+01, -7.43011621233977859e+00, -1.39316492370256828e+00, +1.71997631628774857e-02,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+9.01559677770751189e-01, -1.02181452649552842e+01, -4.45653642971372719e+01, -8.32231241696824640e+01,
     -2.11809894014149904e+01, +7.59081415105611086e+01, +6.42840566032876524e+01, +1.66823079739982383e+01,
     +1.38700067837727348e+00, +2.45566891949830536e-02, -1.08744830155508621e-01, +3.23563873984052686e+00,
     +2.71637488396040041e+01, +9.52074920199811459e+01, +1.57470194748602808e+02, +1.24403180978260934e+02,
     +4.56832902692258358e+01, +7.14787774054067437e+00, +3.89081904963989367e-01, +4.25423549059482809e-03,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.94187312548028353e-03, +1.29466774054608136e+00, -2.87389244888053863e+00, -2.04993868805468082e+00,
     -2.72117081657052062e-01, +1.95683262683451176e+00, +1.52669813670217991e+00, +3.82197934944106588e-01,
     +3.10686186813673333e-02, +5.41287758544253556e-04, -6.10713503279146689e-04, -2.06251400356580100e-01,
     +5.23317490837403754e-01, +2.62945940718054771e+00, +4.00329866912215948e+00, +2.99833429511636496e+00,
     +1.06312369927763939e+00, +1.62314189864056463e-01, +8.67727694899908183e-03, +9.35757776633614521e-05,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.52406144628977527e-04, +1.62310440185935989e-02, +1.49892491782654314e+00, -1.78044024544457602e+00,
     -5.42664993554353714e-02, +1.64358972385525715e-01, +1.22818235928813621e-01, +2.98083970817376170e-02,
     +2.37208249983472804e-03, +4.06889143340293343e-05, -2.48437207783498462e-05, -3.78240614996424749e-03,
     -2.84396417971168269e-01, -1.31354574114084244e-01, +3.21417837403714657e-01, +2.43242754600942135e-01,
     +8.40759836137294769e-02, +1.25549866684845513e-02, +6.59723134179200491e-04, +7.01975151305160790e-06,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.62271465950134676e-05, +1.34008181048557158e-03, +3.46705224364995956e-02, +1.59250216248134779e+00,
     -1.61824719049036503e+00, -2.68109049648046494e-02, +1.23346560945113495e-02, +3.86361431862626776e-03,
     +3.25145956594496965e-04, +5.68521050966264331e-06, -2.69630820386939560e-06, -3.37708974257111549e-04,
     -1.09934248030185699e-02, -3.41900531202554703e-01, -3.26527305878698826e-01, +1.29135092720389366e-02,
     +9.93725284934043884e-03, +1.66883692145347595e-03, +9.09843947738409042e-05, +9.82415965861744327e-07,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.74050549735878477e-06, +1.53032943303075836e-04, +3.57581013553073812e-03, +4.64886325252726607e-02,
     +1.61381970099166594e+00, -1.61381970096517491e+00, -4.64886325377822721e-02, -3.57581014759841007e-03,
     -1.53032945165720966e-04, -1.74050554837821960e-06, -2.85217436360145729e-07, -3.81900381384162586e-05,
     -1.19763972312324630e-03, -1.81124538649256837e-02, -3.66777203353466807e-01, -3.66777203370023674e-01,
     -1.81124538870014763e-02, -1.19763972988114241e-03, -3.81900387297321817e-05, -2.85217445746112735e-07,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-5.68517357936472920e-06, -3.25144190684626539e-04, -3.86359686798811320e-03, -1.23346042988449669e-02,
     +2.68109402067867320e-02, +1.61824715148666765e+00, -1.59250221250064850e+00, -3.46705381735283452e-02,
     -1.34008331189983099e-03, -1.62271762805516308e-05, +9.82409409387204218e-07, +9.09838769711604782e-05,
     +1.66882884385274031e-03, +9.93720997181121897e-03, +1.29134152112183754e-02, -3.26527397588447288e-01,
     -3.41900570945321114e-01, -1.09934319209980414e-02, -3.37709408092931303e-04, -2.69631342753662335e-06,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-4.06962434149808983e-05, -2.37242564361984801e-03, -2.98117158476178673e-02, -1.22827849483402779e-01,
     -1.64365196245687989e-01, +5.42738757015746309e-02, +1.78044931564126108e+00, -1.49892212461581664e+00,
     -1.62307822120609929e-02, -1.52401051214476650e-04, +7.02105714349959747e-06, +6.59824307978859039e-04,
     +1.25565359820103764e-02, +8.40840607251977096e-02, +2.43260166129948191e-01, +3.21434528812950038e-01,
     -1.31347457947227353e-01, -2.84395163304921705e-01, -3.78233082087522324e-03, -2.48428267162970709e-05,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-5.41195184368196774e-04, -3.10643215750732310e-02, -3.82156714368028727e-01, -1.52657975656450962e+00,
     -1.95675724713596288e+00, +2.72025763935512976e-01, +2.04982780770945894e+00, +2.87385850549513977e+00,
     -1.29467090779643335e+00, -3.94193451573538953e-03, +9.35592621117027070e-05, +8.67600715649941880e-03,
     +1.62294884636288606e-01, +1.06302371332511325e+00, +2.99812003671387117e+00, +4.00309436639150018e+00,
     +2.62937271931861494e+00, +5.23302271062764945e-01, -2.06252310754490548e-01, -6.10724273574036352e-04,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-2.45582736633880375e-02, -1.38707234934187862e+00, -1.66829785881952120e+01, -6.42859320283444049e+01,
     -7.59092737422990211e+01, +2.11824600327574188e+01, +8.32248412039895840e+01, +4.45658798652320485e+01,
     +1.02181926510425480e+01, -9.01558771177692009e-01, +4.25451939943239949e-03, +3.89103226258010471e-01,
     +7.14819492006097601e+00, +4.56849004053709749e+01, +1.24406568459581450e+02, +1.57473371027016356e+02,
     +9.52088192846091488e+01, +2.71639786476775349e+01, +3.23565231356250749e+00, -1.08744671407890112e-01,
     +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.21976395487899208e-01, -2.57356519719443932e+00, -1.85357676235941788e+01, -5.00274463990769718e+01,
     -4.37740155991923530e+01, +2.43395944891963154e+01, +5.69516337996174116e+01, +2.87643645605671097e+01,
     +5.38360518596450000e+00, +3.44748548075528127e-01, +4.87184014918244732e-03, -1.55227005670894869e-02,
     +1.55227005670894869e+00, +1.04778228827854036e+01, +4.47053776332177222e+01, +9.77930135726637673e+01,
     +1.09528175201383419e+02, +6.22319177280587610e+01, +1.71943760127760470e+01, +2.09556457655708073e+00,
     +9.13100033358205110e-02, +8.16984240373130888e-04, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+8.41657573948309447e-04, +2.18957919158648512e-01, -3.41451806251435674e-01, -1.58335231427502410e+00,
     -1.13775805264594965e+00, +1.00308465049602091e+00, +1.81198504680651240e+00, +8.61672917879851785e-01,
     +1.56102663969639760e-01, +9.78128750942968903e-03, +1.36029778358144019e-04, -1.28997511637862772e-04,
     -3.86992534913588315e-02, +7.83659883200016338e-01, +1.85756416758522391e+00, +3.41327415793784894e+00,
     +3.51079627673607320e+00, +1.89626342107658274e+00, +5.06608409341424704e-01, +6.02815294769243337e-02,
     +2.57995023275725543e-03, +2.27642667596228421e-05, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+4.18992164798380294e-05, +4.82580491647933638e-03, +3.00699942824524336e-01, +2.31851187501368766e-01,
     -3.28073256241623312e-02, +1.70073176035657525e-01, +2.14921459701145068e-01, +9.32676027462418650e-02,
     +1.61312927632804651e-02, +9.81588657720879314e-04, +1.33712612642695956e-05, -6.69537257635965943e-06,
     -1.11589542939327657e-03, -6.77906473356415518e-02, +4.82066825497895479e-01, +4.92109884362434968e-01,
     +4.25182940089143813e-01, +2.10904236155329272e-01, +5.35629806108772755e-02, +6.16278612142195925e-03,
     +2.57514329859986901e-04, +2.23179085878655314e-06, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+6.39862705590442580e-06, +5.91689086420686314e-04, +1.52671233211488094e-02, +3.69441731616994062e-01,
     +4.07248268772108235e-01, +1.02883773156448644e-01, +7.28760059619570216e-02, +2.70671595193145709e-02,
     +4.35998942167824168e-03, +2.54489055145051010e-04, +3.37146172874670349e-06, -1.04150240159647788e-06,
     -1.45810336188744041e-04, -4.92109884524027583e-03, -1.04983442008590444e-01, +3.21511791082934595e-01,
     +1.54325659688240219e-01, +6.43023581907858355e-02, +1.49976345606131350e-02, +1.64036627976905719e-03,
     +6.62774254336478806e-05, +5.60808984344565668e-07, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+2.36061738476230663e-06, +1.96858685272322670e-04, +4.01303149291108024e-03, +3.66598377272726178e-02,
     +4.21199386635278283e-01, +4.43800424180060271e-01, +7.19121911754456028e-02, +1.93123936439546368e-02,
     +2.75094830524146190e-03, +1.50648929418847085e-04, +1.91860776014516654e-06, -3.88249163353392979e-07,
     -4.99178540623333505e-05, -1.41517390391671344e-03, -1.67724615348980276e-02, -1.54097246404619692e-01,
     +2.21900376691788803e-01, +5.13659013758716240e-02, +1.00635368996339816e-02, +1.01084758868052885e-03,
     +3.88254664232258700e-05, +3.17663247501702670e-07, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.91900226354104274e-06, +1.50671984842940578e-04, +2.75123921494359932e-03, +1.93135955139930905e-02,
     +7.19138061808269752e-02, +4.43800250776314464e-01, +4.21197706630486668e-01, +3.66587961513697505e-02,
     +4.01281094311281648e-03, +1.96843219336181717e-04, +2.36038250997812851e-06, -3.17731302900352492e-07,
     -3.88318735292556176e-05, -1.01096993326787901e-03, -1.00643585870976195e-02, -5.13682816951281443e-02,
     -2.21903623774022229e-01, +1.54095107266073866e-01, +1.67717980243888634e-02, +1.41508516130706597e-03,
     +4.99136813133241661e-05, +3.88209390747013144e-07, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.34439005721081967e-06, +2.52967194483324077e-04, +4.34147720176341657e-03, +2.69934187938876656e-02,
     +7.27814027289576809e-02, +1.02898545790307314e-01, +4.07348502025573711e-01, +3.69501815198740369e-01,
     +1.52795683796704031e-02, +5.92546826665598498e-04, +6.41146989329350400e-06, -5.56108838771992010e-07,
     -6.58504328314824894e-05, -1.63247643052240148e-03, -1.49462132746126554e-02, -6.41574005061052860e-02,
     -1.54132703028196177e-01, -3.21387420882495167e-01, +1.05021280798385172e-01, +4.92607484477072116e-03,
     +1.46040915448190348e-04, +1.04367294979719001e-06, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.37273487260711734e-05, +1.00173896086833743e-03, +1.63782594373876531e-02, +9.42600881433504900e-02,
     +2.16211377032417995e-01, +1.69890986620349654e-01, -3.41667369369998281e-02, +2.31027130319305496e-01,
     +3.00527851225456111e-01, +4.81385868259462044e-03, +4.17191665434350604e-05, -2.29355602169530475e-06,
     -2.63158204645883725e-04, -6.26774453024670202e-03, -5.42517227077433384e-02, -2.12859521580576416e-01,
     -4.27804156341670461e-01, -4.93811240601359724e-01, -4.82587984635245826e-01, +6.77216622030841625e-02,
     +1.11267876408270680e-03, +6.66491303265767214e-06, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.33036469182410457e-04, +9.62009696620871359e-03, +1.54216633394945152e-01, +8.54436580711991689e-01,
     +1.80313041003844909e+00, +1.00491093669114284e+00, -1.12819205740302686e+00, -1.57780427343165686e+00,
     -3.40328094329817865e-01, +2.19033953150615940e-01, +8.42777741965422544e-04, -2.22407351893920721e-05,
     -2.53424454828443581e-03, -5.94662155792526348e-02, -5.01457406563253105e-01, -1.88213631245988555e+00,
     -3.49244389578619119e+00, -3.40169893582171357e+00, -1.85411022271952741e+00, -7.83213550277818360e-01,
     +3.87196103016010373e-02, +1.29186387996816797e-04, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+4.87296979405118982e-03, +3.44840572547135904e-01, +5.38503900208132880e+00, +2.87712901470889335e+01,
     +5.69624808452245475e+01, +2.43399748353381220e+01, -4.37846648425588754e+01, -5.00347235511693853e+01,
     -1.85373997612231634e+01, -2.57368473097423944e+00, +1.21974513851543412e-01, -8.17165333535844140e-04,
     -9.13338662109639959e-02, -2.09612694181650161e+00, -1.71987437737919500e+01, -6.22459921490213972e+01,
     -1.09549000605801012e+02, -9.78076412648248805e+01, -4.47101566546576237e+01, -1.04784902216484692e+01,
     -1.55230260523224861e+00, +1.55223803269937989e-02, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+8.58965292186519402e-01, -1.27279778002339725e+01, -7.63226764151595347e+01, -2.13917403027742836e+02,
     -1.93904822569321055e+02, +9.95143226761370840e+01, +2.45750919832518179e+02, +1.25557246504088567e+02,
     +2.36489328906034331e+01, +1.52093033549496548e+00, +2.15622814286390851e-02, -1.01410146176839050e-01,
     +3.93193439084811081e+00, +4.05109876486052869e+01, +1.84768314669573954e+02, +4.15357032747997437e+02,
     +4.72154046531813208e+02, +2.70783765538653824e+02, +7.52971981991348542e+01, +9.21983197055149972e+00,
     +4.03167911623198039e-01, +3.61740929273230609e-03, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+3.13428851060647317e-03, +1.24721364707287008e+00, -3.47239084537513731e+00, -4.50291318198361878e+00,
     -3.49196431684540576e+00, +2.42519176416099080e+00, +4.92001171904848889e+00, +2.40226516016397840e+00,
     +4.41173648621993264e-01, +2.78878200964500431e-02, +3.90296528783936768e-04, -4.71592861798680602e-04,
     -1.93076863194749287e-01, +7.75166946893607767e-01, +4.31418498645403581e+00, +8.83752135472605005e+00,
     +9.49125582505948806e+00, +5.24684758303919807e+00, +1.42222422507724588e+00, +1.70917204179983989e-01,
     +7.36789704976780792e-03, +6.53679428588865055e-05, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+9.96964993359876944e-05, +1.31421443100789709e-02, +1.46009057728647629e+00, -1.93896297892859693e+00,
     -2.60945149057079107e-01, +1.96242055669426929e-01, +3.41305434994865209e-01, +1.58861398041172303e-01,
     +2.83817316497389782e-02, +1.76078636146307851e-03, +2.43031731183921879e-05, -1.57600256248510314e-05,
     -2.92424305198004568e-03, -2.68035295188643163e-01, -2.22385362056749213e-02, +6.33444673439953476e-01,
     +6.60678315209224820e-01, +3.51820759533438832e-01, +9.27521281742623854e-02, +1.09207662067741951e-02,
     +4.63531510571816153e-04, +4.06305341693084824e-06, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+9.05261580489646631e-06, +9.20196063537080623e-04, +2.93998830072255378e-02, +1.57103038316028010e+00,
     -1.64614703223909076e+00, -2.23834160717181208e-02, +4.18743085443801510e-02, +2.12430380706603306e-02,
     +3.81486670301714979e-03, +2.35490319548944959e-04, +3.22982635469912836e-06, -1.45965249361409477e-06,
     -2.21015387505149304e-04, -8.77156235283108390e-03, -3.27104586181059698e-01, -2.84289369315013751e-01,
     +6.93086490306589337e-02, +4.60258992837386804e-02, +1.24489190013734180e-02, +1.46575075572283374e-03,
     +6.19001285292811990e-05, +5.39457282898515956e-07, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.40215793686603336e-06, +1.31905148851615907e-04, +3.29092553802271285e-03, +4.52252029514151471e-02,
     +1.61195178830599772e+00, -1.61381969798347091e+00, -4.46207187477793907e-02, -2.31238332974844632e-03,
     +1.31850452977965907e-04, +1.93871608927796554e-05, +3.38344903921246952e-07, -2.27458598075135604e-07,
     -3.24141590052162660e-05, -1.08067826575579200e-03, -1.72807286219915784e-02, -3.64230046867608848e-01,
     -3.63109301012092156e-01, -1.55653015423990251e-02, -3.65917191440422813e-04, +7.87708481586484098e-05,
     +5.49062412504573782e-06, +5.77583687549083253e-08, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-3.37930497541136063e-07, -1.93624669560720610e-05, -1.31532774886517758e-04, +2.31372319920852137e-03,
     +4.46225707747645167e-02, +1.61381955827312695e+00, -1.61195369240697461e+00, -4.52264136785025805e-02,
     -3.29118676015982027e-03, -1.31923783421593918e-04, -1.40244570131481778e-06, +5.76871014621317277e-08,
     +5.48379472605266305e-06, +7.86381864741455882e-05, -3.66823089833269087e-04, -1.55679685306647185e-02,
     -3.63112996963218415e-01, -3.64232519668397925e-01, -1.72815074193404562e-02, -1.08078401696249478e-03,
     -3.24192072934489627e-05, -2.27507453185129059e-07, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-3.27071501428337479e-06, -2.37851468050804464e-04, -3.84431701785634822e-03, -2.13632915294433498e-02,
     -4.20334983631035081e-02, +2.24030088495325158e-02, +1.64631359299866098e+00, -1.57092824176927071e+00,
     -2.93784036045260551e-02, -9.18697440380675614e-04, -9.02994054768786228e-06, +5.46524660886749544e-07,
     +6.25582398196813591e-05, +1.47818849885469721e-03, +1.25316547796829741e-02, +4.62634693468385594e-02,
     +6.96301522973142989e-02, -2.84079061880278483e-01, -3.27039753511622733e-01, -8.76293592162488526e-03,
     -2.20611454058617669e-04, -1.45581447240213470e-06, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-2.37981646191077112e-05, -1.73166217927467922e-03, -2.80186365375615261e-02, -1.57378269941534310e-01,
     -3.39339098261707606e-01, -1.96479324863575111e-01, +2.58889255576282007e-01, +1.93769919143670299e+00,
     -1.46035691136961421e+00, -1.31607667979984425e-02, -9.99788971000597490e-05, +3.97573944909288326e-06,
     +4.55411706815484526e-04, +1.07674023158741278e-02, +9.17319317448162946e-02, +3.48889711734689301e-01,
     +6.56707787220566330e-01, +6.30843906560388912e-01, -2.30415841691528241e-02, -2.68142339860793777e-01,
     -2.92926519064159313e-03, -1.58078410049050934e-05, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-3.96137489453936818e-04, -2.82142412895799929e-02, -4.45125729016738692e-01, -2.41794628678634593e+00,
     -4.94004509410320747e+00, -2.42198255920410155e+00, +3.51322641062284102e+00, +4.51562747722010283e+00,
     +3.47502037162993682e+00, -1.24703263056991132e+00, -3.13158101354254429e-03, +6.63831449849208564e-05,
     +7.45960993388735533e-03, +1.72604374190134024e-01, +1.43318228143716730e+00, +5.27765572722965448e+00,
     +9.53217928439128492e+00, +8.86385459664836084e+00, +4.32218577722356301e+00, +7.76217875847568073e-01,
     -1.93028214767015482e-01, -4.71135342385006653e-04, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {-2.15661071403811348e-02, -1.52111427474358751e+00, -2.36508091510973664e+01, -1.25563257395346270e+02,
     -2.45756161868294612e+02, -9.95108267856362403e+01, +1.93911442644113095e+02, +2.13920181710441256e+02,
     +7.63230806470995304e+01, +1.27279957367076175e+01, -8.58965156103019711e-01, +3.61808935432900419e-03,
     +4.03221686658249756e-01, +9.22068807264593282e+00, +7.53019469190844575e+01, +2.70794986410164711e+02,
     +4.72166321027351898e+02, +4.15363361557544692e+02, +1.84769796077367261e+02, +4.05111285884645298e+01,
     +3.93193861736329398e+00, -1.01410126434251610e-01, +0.00000000000000000e+00, +0.00000000000000000e+00},
    {+1.13711216456042832e-01, -3.17963117741119294e+00, -2.83984572978950673e+01, -1.05433092758028124e+02,
     -1.61358472394895216e+02, -3.72700040615818580e+01, +1.41379602184017131e+02, +1.38964450943640259e+02,
     +4.88353861528489688e+01, +6.98539988938819079e+00, +3.56975351445089462e-01, +4.13195201577819336e-03,
     -1.41436176654679002e-02, +1.71137773752161593e+00, +1.42614811460134661e+01, +7.70119981884727167e+01,
     +2.20034280538493476e+02, +3.35430036554236722e+02, +2.74442757180739136e+02, +1.18479997213034949e+02,
     +2.56706660628242389e+01, +2.51673196694355283e+00, +9.00725125011376804e-02, +6.73505603117519058e-04},
    {+6.45258908630968000e-04, +2.04594737295870549e-01, -5.74478516523088567e-01, -2.88763809210097744e+00,
     -3.89176278503457862e+00, -4.19112915311416145e-01, +3.80306163893692428e+00, +3.43157977294530214e+00,
     +1.16236963700722695e+00, +1.62480172304921788e-01, +8.16773734671057475e-03, +9.33542244732867477e-05,
     -9.62150861596455798e-05, -3.49260762759513455e-02, +8.73151906898783636e-01, +2.61945572069635091e+00,
     +6.28669372967124218e+00, +8.80137122153973905e+00, +6.84551095008646371e+00, +2.85758805894147372e+00,
     +6.04489781699157902e-01, +5.82101271265855758e-02, +2.05447507505596150e-03, +1.51918557094177231e-05},
    {+2.57452168554735078e-05, +3.64789094170364907e-03, +2.81654165105140131e-01, +1.25686191584658685e-01,
     -2.55691496282192288e-01, +5.68879026322974754e-02, +3.77317721540748607e-01, +3.00226351490502086e-01,
     +9.65161152553914020e-02, +1.30785925469873686e-02, +6.43579847248821257e-04, +7.24012065861584432e-06,
     -3.99785676563624293e-06, -8.06234447736642323e-04, -6.04675835802481743e-02, +5.44208252222233568e-01,
     +7.25611002962978091e-01, +8.53318539484462235e-01, +6.09513242488901597e-01, +2.41870334320992697e-01,
     +4.94734774747485062e-02, +4.65135258309601427e-03, +1.61246889547328465e-04, +1.17584022518713027e-06},
    {+3.06259924382005836e-06, +3.49228778536175670e-04, +1.13615978833297299e-02, +3.47772298198451346e-01,
     +3.62052273065046448e-01, +8.03779483552925084e-02, +1.06098890388201861e-01, +6.87949522573243072e-02,
     +2.04150740023692594e-02, +2.64672598108818156e-03, +1.26553708042481585e-04, +1.39478307381784322e-06,
     -4.84169515200541841e-07, -8.20183153107163115e-05, -3.41742978268659526e-03, -9.22706035801804242e-02,
     +3.69082412258282200e-01, +2.41133841416868278e-01, +1.44680304144360700e-01, +5.27260581210833579e-02,
     +1.02522890353134999e-02, +9.32026272208239860e-04, +3.15455044794477353e-05, +2.25945761315889778e-07},
    {+8.48795617232961390e-07, +8.73822893102994158e-05, +2.25711430587143339e-03, +2.69685526856586685e-02,
     +4.01138085031401470e-01, +4.34039946462657533e-01, +8.68094208164374853e-02, +3.77113644504017201e-02,
     +9.75124176880897212e-03, +1.18123771984057855e-03, +5.42243026484026918e-05, +5.81371346211903380e-07,
     -1.35546616617780277e-07, -2.10880404391572603e-05, -7.38108426633137579e-04, -1.10719910090082800e-02,
     -1.32867839379710445e-01, +2.60427963306457400e-01, +8.68114368363867434e-02, +2.65755186885939902e-02,
     +4.74572426625809434e-03, +4.10131900636017961e-04, +1.34227260841803371e-05, +9.38666802053354043e-08},
    {+4.96198919738902950e-07, +4.80713262163315637e-05, +1.11367000223126962e-03, +1.03307315888892821e-02,
     +5.34834553014765235e-02, +4.35083533312191573e-01, +4.35048066535202338e-01, +5.34278225091289471e-02,
     +1.03063122558433931e-02, +1.10952310594628696e-03, +4.78249477926716453e-05, +4.92916161630190124e-07,
     -7.97689516175147179e-08, -1.17845659085876208e-05, -3.78440664865063356e-04, -4.76447597049723792e-03,
     -3.17402744248264507e-02, -1.86512914988586603e-01, +1.86405630587153814e-01, +3.16852357063910420e-02,
     +4.75056378278265051e-03, +3.76868994918949239e-04, +1.17202850504166403e-05, +7.92226345050764720e-08},
    {+5.06207266281322842e-07, +4.91325797200047718e-05, +1.10382767947292426e-03, +9.34119907913816792e-03,
     +3.68871437983726153e-02, +8.64182929658595389e-02, +4.34648555320953756e-01, +4.01894622664971279e-01,
     +2.72636163826979229e-02, +2.30240044912537927e-03, +8.98245444298377321e-05, +8.78327992312841269e-07,
     -8.11561545371664361e-08, -1.20613622808822854e-05, -3.79754772802738466e-04, -4.49976430954376980e-03,
     -2.56837884634620560e-02, -8.52163037071747981e-02, -2.58962901844033144e-01, +1.33558341388131413e-01,
     +1.12322940377992050e-02, +7.54726830850457496e-04, +2.17109053240741850e-05, +1.40388662251762854e-07},
    {+2.59017745901164656e-06, +2.09211011827320752e-04, +3.92873691146099387e-03, +2.73425915414492147e-02,
     +8.30329415733500609e-02, +1.13224006832725224e-01, +7.00648533028884746e-02, +3.48812936739885671e-01,
     +3.42527662869219873e-01, +1.05470820474377249e-02, +3.04865481760785970e-04, +2.52151053564038583e-06,
     -4.27440930459477635e-07, -5.35367175356901091e-05, -1.43166832572539748e-03, -1.43671856261534288e-02,
     -6.78819682906240673e-02, -1.72186918164722230e-01, -2.66732768695769840e-01, -3.81292216454515442e-01,
     +8.94054554891649395e-02, +3.11752815596925666e-03, +7.06804497787781946e-05, +3.95347915841149361e-07},
    {+2.19711522579101417e-06, +3.33089867636882813e-04, +8.78150981671081514e-03, +7.58603390667829253e-02,
     +2.63258485325438552e-01, +3.64957093133586263e-01, +8.69963492296687152e-02, -2.24086974086985496e-01,
     +1.36899865821035138e-01, +2.83243825924258896e-01, +3.72757508969368857e-03, +2.66436969478141871e-05,
     -3.09983172781922205e-07, -7.62815689260157922e-05, -2.90859146157732527e-03, -3.64618428007420226e-02,
     -1.98246270596108359e-01, -5.37160151578271350e-01, -7.91555929458877205e-01, -6.98496820414551189e-01,
     -5.38333252998141054e-01, +6.10370818552799857e-02, +8.26227235492566978e-04, +4.14365977557307720e-06},
    {+1.10556682886181037e-04, +9.46345129378198778e-03, +1.83967554373186293e-01, +1.28509223214697443e+00,
     +3.69737021006557809e+00, +3.95066455531207517e+00, -6.03862463320473553e-01, -4.14528528630945647e+00,
     -2.99098679257047495e+00, -5.90856380608582921e-01, +2.03688308315558675e-01, +6.34054618947242632e-04,
     -1.80466283190310879e-05, -2.39271107674061435e-03, -6.64032683861923156e-02, -6.75477557246641380e-01,
     -3.12991466015007491e+00, -7.35638398714062269e+00, -9.28991990569190484e+00, -6.52504129771601171e+00,
     -2.67645822674054390e+00, -8.79215366278865318e-01, +3.46936616881181256e-02, +9.43725054885398304e-05},
    {+4.00709773728189639e-03, +3.48766295218042865e-01, +6.86396848134057880e+00, +4.82088447516690906e+01,
     +1.37740829201086427e+02, +1.40839296870673309e+02, -3.63447384718574042e+01, -1.60253334836371067e+02,
     -1.05010562919787958e+02, -2.83345963331938258e+01, -3.17623198217635672e+00, +1.13751845661890838e-01,
     -6.52287765769117372e-04, -8.78621001243751337e-02, -2.46863035015036303e+00, -2.52899208570967566e+01,
     -1.17127618399058079e+02, -2.72068118691447921e+02, -3.33285436426285628e+02, -2.19038822886903063e+02,
     -7.67840815492022557e+01, -1.42381493793694851e+01, -1.71051323246679643e+00, +1.41502682163508232e-02},
    {+8.21111642458840718e-01, -1.55052310917325702e+01, -1.21546196284493547e+02, -4.68161718185059158e+02,
     -7.34036033241807900e+02, -1.84327024364905071e+02, +6.33028442086270616e+02, +6.32214942749138572e+02,
     +2.23706111798549870e+02, +3.21388791678490549e+01, +1.64759808349714554e+00, +1.91176402341768116e-02,
     -9.50945894868213681e-02, +4.66093437781892117e+00, +5.78564280098422969e+01, +3.32961244103407807e+02,
     +9.76464615214135847e+02, +1.50986320303781827e+03, +1.24643121691703277e+03, +5.41412374857682153e+02,
     +1.17832592317889066e+02, +1.15916930816274544e+01, +4.15976043571572529e-01, +3.11715918807859652e-03},
    {+2.53756256112447299e-03, +1.20353326656273230e+00, -4.18180306241542876e+00, -8.47865370622103143e+00,
     -1.19015956288591447e+01, -1.93976363189576650e+00, +1.09854389514221111e+01, +1.02615980061976249e+01,
     +3.52608754578014505e+00, +4.97187115804961999e-01, +2.51448716886603467e-02, +2.88709374012662045e-04,
     -3.72001723127688593e-04, -1.81604727196619409e-01, +1.04751055232035073e+00, +6.63513928455618718e+00,
     +1.76007859441805117e+01, +2.56469713186001158e+01, +2.03822799747897072e+01, +8.62383100566928285e+00,
     +1.84118074858945735e+00, +1.78492940790983129e-01, +6.33197791144734113e-03, +4.70100676655231130e-05},
    {+6.72982063857198529e-05, +1.07760796654440324e-02, +1.42176539272545144e+00, -2.15305789825342053e+00,
     -7.11783009519793453e-01, -3.47492658445844455e-02, +6.68448549477034248e-01, +5.78583080604570310e-01,
     +1.92105585963000342e-01, +2.65089816376322222e-02, +1.32022359135116826e-03, +1.49817469288787256e-05,
     -1.03511007040793491e-05, -2.30246782805243987e-03, -2.53308531877597221e-01, +1.02943778009440781e-01,
     +1.10473970796322507e+00, +1.52672550709394272e+00, +1.16017687563726213e+00, +4.75744495685960722e-01,
     +9.93103025004009876e-02, +9.46697381160837258e-03, +3.31496053972625223e-04, +2.43577114933823496e-06},
    {+5.14841136836391782e-06, +6.35602193281654872e-04, +2.47999810593790351e-02, +1.54540126873257133e+00,
     -1.69991873223179168e+00, -4.96382224655285559e-02, +8.10869219871312025e-02, +7.11450104177856757e-02,
     +2.31796416143802766e-02, +3.14690562391077696e-03, +1.54736441772425950e-04, +1.73821573935104832e-06,
     -8.07669336362795916e-07, -1.46192706786905599e-04, -7.00267736760407794e-03, -3.12099759608578387e-01,
     -2.27930616203335795e-01, +1.72597120900183337e-01, +1.42139686923898723e-01, +5.78260125783099863e-02,
     +1.18948858553951948e-02, +1.11901088465267895e-03, +3.87566062186709513e-05, +2.82230310963710485e-07},
    {+7.65116265329614628e-07, +8.53492285729881859e-05, +2.53621552159675936e-03, +4.10050938606389831e-02,
     +1.60305310158863357e+00, -1.61839704233881521e+00, -3.81753548446124580e-02, +5.98178394211032298e-03,
     +3.37335965308532045e-03, +5.10325898575241183e-04, +2.61036977661833955e-05, +2.98676182769257051e-07,
     -1.21114419294838318e-07, -2.01818941713794623e-05, -7.90756164874351140e-04, -1.48143903937389236e-02,
     -3.54936536673319480e-01, -3.46014832759978685e-01, +4.08354981020808358e-04, +7.21206484668517663e-03,
     +1.83030567800973490e-03, +1.84221085255131632e-04, +6.57159601537917382e-06, +4.85822880003726783e-08},
    {+9.50545469760443839e-08, +1.14171462390892460e-05, +3.51363630585618597e-04, +4.90508507614678940e-03,
     +4.97668507973926533e-02, +1.61598783190450179e+00, -1.61599336435096852e+00, -4.97650808361139937e-02,
     -4.90217871202878411e-03, -3.50571816047030881e-04, -1.13538867144820165e-05, -9.40075399908022440e-08,
     -1.48730679251375957e-08, -2.65945573150105783e-06, -1.09029824007050787e-04, -1.91261592554521631e-03,
     -2.11782628390372284e-02, -3.72986334820056997e-01, -3.72983018537664912e-01, -2.11731706242871236e-02,
     -1.91048055958034191e-03, -1.08694461071118663e-04, -2.64198756598605008e-06, -1.46938022685541688e-08},
    {-4.41373783980161750e-07, -3.60368019865108164e-05, -6.64972236264215985e-04, -4.21025609743261482e-03,
     -7.70086065565173550e-03, +3.73197759540019949e-02, +1.61964566696922319e+00, -1.60145802629692224e+00,
     -4.03753563271191315e-02, -2.43873505895644046e-03, -8.00572911715536473e-05, -7.00783936637902197e-07,
     +7.26023948438688329e-08, +9.21063846011469051e-06, +2.44423285656077799e-04, +2.32709124801032657e-03,
     +9.04250876860767348e-03, +3.72776491896278754e-03, -3.42930753903648966e-01, -3.53468833444474495e-01,
     -1.44708825177229002e-02, -7.54902921723405395e-04, -1.88304803473717496e-05, -1.10559281349981141e-07},
    {+2.41415721956815448e-07, -1.56469532542195267e-05, -9.60292449100800562e-04, -1.12229370629323474e-02,
     -4.62793688051977847e-02, -6.83289916714086981e-02, +3.17944674436393057e-02, +1.67665938763148459e+00,
     -1.55467818667644697e+00, -2.62477873879196980e-02, -7.14768301734897984e-04, -6.11718285022672012e-06,
     -5.05362170808449568e-08, +1.88692771882941662e-06, +2.70647251216903290e-04, +4.83341587773362773e-03,
     +3.15847671256881402e-02, +9.41567687512098184e-02, +1.27659561055318144e-01, -2.49479172200165605e-01,
     -3.17179457153714595e-01, -7.53646936598475516e-03, -1.66441431508265508e-04, -9.66768710685903563e-07},
    {-2.94806602448759522e-05, -2.28657191248979406e-03, -4.09686100618275248e-02, -2.67448192901593162e-01,
     -7.27085598557017931e-01, -7.35301393947734305e-01, +1.46450249824100931e-01, +8.46372193059346301e-01,
     +2.20467408390121955e+00, -1.41395411358018164e+00, -1.03602339314986091e-02, -6.23312320789096754e-05,
     +4.89395121806636411e-06, +5.90878370935503101e-04, +1.51726071872762547e-02, +1.44875077049894441e-01,
     +6.38741704234350505e-01, +1.44796479664591282e+00, +1.78769710622841105e+00, +1.22621613788962147e+00,
     +1.30804164652529586e-01, -2.50454327273303468e-01, -2.19671895979898319e-03, -9.53820441798706342e-06},
    {-2.76642537617777728e-04, -2.41282162030401100e-02, -4.78872451212265349e-01, -3.41426289662884708e+00,
     -1.00038635110294809e+01, -1.08250108131177250e+01, +1.77012576990871051e+00, +1.16476518665496673e+01,
     +8.37134655150653018e+00, +4.16435893266242285e+00, -1.20451863536605730e+00, -2.54995453229817527e-03,
     +4.50498111905334346e-05, +6.07284782342196183e-03, +1.71690748898417400e-01, +1.77848382946543171e+00,
     +8.37118534812304205e+00, +1.98889714867617684e+01, +2.51593125579660644e+01, +1.73561489046289594e+01,
     +6.57522961938856909e+00, +1.04100682347803275e+00, -1.81858452304734880e-01, -3.74044462690546384e-04},
    {-1.93937539473814750e-02, -1.66646788951276847e+00, -3.24278501580332374e+01, -2.25245878087926116e+02,
     -6.35325689051146469e+02, -6.34517587253520901e+02, +1.86618670801984678e+02, +7.36895475196541952e+02,
     +4.69277811202846920e+02, +1.21717450729002238e+02, +1.55144584799105925e+01, -8.21000216199603949e-01,
     +3.16378010236818483e-03, +4.21011006290568202e-01, +1.17048170798283220e+01, +1.18753563626098710e+02,
     +5.44765196408638531e+02, +1.25244632865145763e+03, +1.51539811336286630e+03, +9.79075771894341528e+02,
     +3.33567552160577122e+02, +5.79192561925530072e+01, +4.66328701295883352e+00, -9.50763248527119881e-02}

};

} // namespace integrators
} // namespace tudat

#endif // TUDAT_ADAMS_BASHFORTH_MOULTON_INTEGRATOR_H
