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

#ifndef TUDAT_BULIRSCH_STOER_VARIABLE_STEP_SIZE_INTEGRATOR_H
#define TUDAT_BULIRSCH_STOER_VARIABLE_STEP_SIZE_INTEGRATOR_H

#include <cmath>

#include <boost/assign/std/vector.hpp>

#include <Eigen/Core>

#include <tudat/math/integrators/stepSizeController.h>
#include <tudat/math/integrators/numericalIntegrator.h>
#include <tudat/math/basic/mathematicalConstants.h>

namespace tudat
{

namespace numerical_integrators
{

// Types of sequences available for extrapolation integration methods
//! @get_docstring(ExtrapolationMethodStepSequences.__docstring__)
enum ExtrapolationMethodStepSequences
{
    bulirsch_stoer_sequence = 0,
    deufelhard_sequence = 1
};

// Function to retrieve the sequence of number of steps to used for Bulirsch-Stoer integration
/*
 * Function to retrieve the sequence of number of steps to used for Bulirsch-Stoer integration
 * \param extrapolationMethodStepSequenceType Type of sequence that is to be retrieved
 * \param lengthOfSequence Length of the sequence that is to be retrieved (default 12)
 * \return Function to retrieve the sequence of number of steps to used for Bulirsch-Stoer integration
 */
std::vector< unsigned int > getBulirschStoerStepSequence(
        const ExtrapolationMethodStepSequences& extrapolationMethodStepSequenceType = bulirsch_stoer_sequence,
        const unsigned int lengthOfSequence = 12 );

// Class that implements the Bulirsch-Stoer variable stepsize integrator.
/*
 * Class that implements the Bulirsch-Stoer variable step size integrator.
 * \tparam StateType The type of the state. This type should be an Eigen::Matrix derived type.
 * \tparam StateDerivativeType The type of the state derivative. This type should be an
 *          Eigen::Matrix derived type.
 * \tparam IndependentVariableType The type of the independent variable. This type should be
 *          either a float or double.
 * \sa NumericalIntegrator.
 */
template< typename IndependentVariableType = double, typename StateType = Eigen::VectorXd,
           typename StateDerivativeType = Eigen::VectorXd, typename TimeStepType = double >
class BulirschStoerVariableStepSizeIntegrator :
        public NumericalIntegrator< IndependentVariableType, StateType, StateDerivativeType, TimeStepType >
{
public:

    // Typedef of the base class.
    /*
     * Typedef of the base class with all template parameters filled in.
     */
    typedef NumericalIntegrator< IndependentVariableType, StateType, StateDerivativeType, TimeStepType > Base;

    // Typedef to the state derivative function.
    /*
     * Typedef to the state derivative function inherited from the base class.
     * \sa NumericalIntegrator::StateDerivativeFunction.
     */
    typedef typename Base::StateDerivativeFunction StateDerivativeFunction;



    BulirschStoerVariableStepSizeIntegrator(
        const std::vector< unsigned int >& sequence,
        const StateDerivativeFunction& stateDerivativeFunction,
        const IndependentVariableType intervalStart,
        const StateType& initialState,
        const TimeStepType initialStepSize,
        const std::shared_ptr< IntegratorStepSizeController< TimeStepType, StateType > > stepSizeController,
        const std::shared_ptr< IntegratorStepSizeValidator< TimeStepType > > stepSizeValidator ):
            Base( stateDerivativeFunction ), currentIndependentVariable_( intervalStart ),
            currentState_( initialState ), lastIndependentVariable_( intervalStart ),
            sequence_( sequence ),
            stepSize_( initialStepSize ),
            stepSizeController_( stepSizeController ),
            stepSizeValidator_( stepSizeValidator )
        {
            maximumStepIndex_ = sequence_.size( ) - 1;
            subSteps_.resize( maximumStepIndex_ + 1 );

            integratedStates_.resize( maximumStepIndex_ + 1  );
            for( unsigned int i = 0; i < maximumStepIndex_ + 1 ; i++ )
            {
                integratedStates_[ i ].resize( maximumStepIndex_ + 1  );
            }

            useFixedStep_ = false;
            stepSizeController_->initialize( initialState );
        }

    BulirschStoerVariableStepSizeIntegrator(
        const std::vector< unsigned int >& sequence,
        const StateDerivativeFunction& stateDerivativeFunction,
        const IndependentVariableType intervalStart,
        const StateType& initialState,
        const TimeStepType initialStepSize ):
        Base( stateDerivativeFunction ), currentIndependentVariable_( intervalStart ),
        currentState_( initialState ), lastIndependentVariable_( intervalStart ),
        sequence_( sequence ),
        stepSize_( initialStepSize ),
        stepSizeController_( nullptr ),
        stepSizeValidator_( nullptr )
    {
        if( sequence_.size( ) <= 0 )
        {
            throw std::runtime_error( "Error when creating BS integrator, sequence is empty." );
        }
        maximumStepIndex_ = static_cast< unsigned int >( sequence_.size( ) ) - 1;
        subSteps_.resize( maximumStepIndex_ + 1 );

        integratedStates_.resize( maximumStepIndex_ + 1  );
        for( unsigned int i = 0; i < maximumStepIndex_ + 1 ; i++ )
        {
            integratedStates_[ i ].resize( maximumStepIndex_ + 1  );
        }

        useFixedStep_ = true;
    }

// Default constructor.
    /*
     * Default constructor, taking sequence, a state derivative function, initial conditions,
     * minimum step size and relative error tolerance per item in the state vector as argument.
     * \param sequence Rational function sequence used by algorithm.
     * \param stateDerivativeFunction State derivative function.
     * \param intervalStart The start of the integration interval.
     * \param initialState The initial state.
     * \param minimumStepSize The minimum step size to take. If this constraint is violated, a
     *          flag will be set that can be retrieved with isMinimumStepSizeViolated( ).
     * \param maximumStepSize The maximum step size to take.
     * \param relativeErrorTolerance The relative error tolerance, for each individual state
     *          vector element.
     * \param absoluteErrorTolerance The absolute error tolerance, for each individual state
     *          vector element.
     * \param safetyFactorForNextStepSize Safety factor used to scale prediction of next step size.
     * \param maximumFactorIncreaseForNextStepSize Maximum factor increase for next step size.
     * \param minimumFactorDecreaseForNextStepSize Maximum factor decrease for next step size.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    BulirschStoerVariableStepSizeIntegrator(
        const std::vector< unsigned int >& sequence,
        const StateDerivativeFunction& stateDerivativeFunction,
        const IndependentVariableType intervalStart,  const StateType& initialState,
        const TimeStepType minimumStepSize,
        const TimeStepType maximumStepSize,
        const TimeStepType initialStepSize,
        const StateType& relativeErrorTolerance,
        const StateType& absoluteErrorTolerance,
        const TimeStepType safetyFactorForNextStepSize = 0.6,
        const TimeStepType maximumFactorIncreaseForNextStepSize = 4.0,
        const TimeStepType minimumFactorDecreaseForNextStepSize = 0.1 ):
        Base( stateDerivativeFunction ), currentIndependentVariable_( intervalStart ),
        currentState_( initialState ), lastIndependentVariable_( intervalStart ),
        sequence_( sequence ),
        stepSize_( initialStepSize )
    {
        maximumStepIndex_ = sequence_.size( ) - 1;
        subSteps_.resize( maximumStepIndex_ + 1 );

        integratedStates_.resize( maximumStepIndex_ + 1  );
        for( unsigned int i = 0; i < maximumStepIndex_ + 1 ; i++ )
        {
            integratedStates_[ i ].resize( maximumStepIndex_ + 1  );
        }

        useFixedStep_ = false;
        if( ( initialStepSize == minimumStepSize ) && ( initialStepSize == maximumStepSize ) &&
            !relativeErrorTolerance.allFinite( ) && !absoluteErrorTolerance.allFinite( ) )
        {
            useFixedStep_ = true;
        }

        stepSizeController_ = std::make_shared< PerElementIntegratorStepSizeController< TimeStepType, StateType > >(
            relativeErrorTolerance, absoluteErrorTolerance,
            static_cast< double >( safetyFactorForNextStepSize ),
            static_cast< double >( 2 * maximumStepIndex_ - 1 ),
            static_cast< double >( minimumFactorDecreaseForNextStepSize ),
            static_cast< double >( maximumFactorIncreaseForNextStepSize ) );
        stepSizeController_->initialize( initialState );

        stepSizeValidator_=
            std::make_shared< BasicIntegratorStepSizeValidator< TimeStepType > >( minimumStepSize, maximumStepSize);

    }

    // Default constructor.
    /*
     * Default constructor, taking coefficients, a state derivative function, initial conditions,
     * minimum step size and relative error tolerance for all items in the state vector as argument.
     * \param sequence Rational function sequence used by algorithm.
     * \param stateDerivativeFunction State derivative function.
     * \param intervalStart The start of the integration interval.
     * \param initialState The initial state.
     * \param minimumStepSize The minimum step size to take. If this constraint is violated, a
     *          flag will be set that can be retrieved with isMinimumStepSizeViolated( ).
     * \param maximumStepSize The maximum step size to take.
     * \param relativeErrorTolerance The relative error tolerance, equal for all individual state
     *          vector elements.
     * \param absoluteErrorTolerance The absolute error tolerance, for each individual state
     *          vector element.
     * \param safetyFactorForNextStepSize Safety factor used to scale prediction of next step size.
     * \param maximumFactorIncreaseForNextStepSize Maximum factor increase for next step size.
     * \param minimumFactorDecreaseForNextStepSize Maximum factor decrease for next step size.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    BulirschStoerVariableStepSizeIntegrator(
            const std::vector< unsigned int >& sequence,
            const StateDerivativeFunction& stateDerivativeFunction,
            const IndependentVariableType intervalStart, const StateType& initialState,
            const TimeStepType minimumStepSize,
            const TimeStepType maximumStepSize,
            const TimeStepType stepSize,
            const typename StateType::Scalar relativeErrorTolerance = 1.0e-12,
            const typename StateType::Scalar absoluteErrorTolerance = 1.0e-12,
            const TimeStepType safetyFactorForNextStepSize = 0.75,
            const TimeStepType maximumFactorIncreaseForNextStepSize = 4.0,
            const TimeStepType minimumFactorDecreaseForNextStepSize = 0.1 ):
        Base( stateDerivativeFunction ), currentIndependentVariable_( intervalStart ),
        currentState_( initialState ), lastIndependentVariable_( intervalStart ),
        sequence_( sequence ), stepSize_( stepSize )
    {
        maximumStepIndex_ = sequence_.size( ) - 1;
        subSteps_.resize( maximumStepIndex_ + 1 );

        integratedStates_.resize( maximumStepIndex_ + 1  );
        for( unsigned int i = 0; i < maximumStepIndex_ + 1 ; i++ )
        {
            integratedStates_[ i ].resize( maximumStepIndex_ + 1  );
        }

        useFixedStep_ = false;
        if( ( stepSize == minimumStepSize ) && ( stepSize == maximumStepSize ) &&
            std::isinf( relativeErrorTolerance ) && std::isinf( absoluteErrorTolerance ) )
        {
            useFixedStep_ = true;
        }

        stepSizeController_ = std::make_shared< PerElementIntegratorStepSizeController< TimeStepType, StateType > >(
            StateType::Constant( initialState.rows( ), initialState.cols( ),
                                 relativeErrorTolerance ),
            StateType::Constant( initialState.rows( ), initialState.cols( ),
                                 absoluteErrorTolerance ),
            static_cast< double >( safetyFactorForNextStepSize ), static_cast< double >( 2 * maximumStepIndex_ - 1 ),
            static_cast< double >( minimumFactorDecreaseForNextStepSize ),
            static_cast< double >( maximumFactorIncreaseForNextStepSize ) );
        stepSizeController_->initialize( initialState );

        stepSizeValidator_=
            std::make_shared< BasicIntegratorStepSizeValidator< TimeStepType > >( minimumStepSize, maximumStepSize);
    }

    ~BulirschStoerVariableStepSizeIntegrator( ){ }

    // Get step size of the next step.
    /*
     * Returns the step size of the next step.
     * \return Step size to be used for the next step.
     */
    virtual TimeStepType getNextStepSize( ) const { return stepSize_; }

    // Get current state.
    /*
     * Returns the current state of the integrator.
     * \return Current integrated state.
     */
    virtual StateType getCurrentState( ) const { return currentState_; }

    // Returns the current independent variable.
    /*
     * Returns the current value of the independent variable of the integrator.
     * \return Current independent variable.
     */
    virtual IndependentVariableType getCurrentIndependentVariable( ) const
    {
        return currentIndependentVariable_;
    }

    // Perform a single integration step.
    /*
     * Perform a single integration step and compute a new step size.
     * \param stepSize The step size to take. If the time step is too large to satisfy the error
     *          constraints, the step is redone until the error constraint is satisfied.
     * \return The state at the end of the interval.
     */
    virtual StateType performIntegrationStep( const TimeStepType stepSize )
    {
        if( !( stepSize == stepSize) )
        {
            throw std::runtime_error( "Error in BS integrator, step size is NaN" );
        }

        StateType stateAtFirstPoint_;
        StateType stateAtCenterPoint_;
        StateType stateAtLastPoint_;

        bool stepSuccessful = 0;

        // Compute sub steps to take.
        for ( unsigned int p = 0; p <  subSteps_.size( ); p++ )
        {
            subSteps_.at( p ) = stepSize / static_cast< double >(
                        sequence_.at( p ) );
        }

        double errorScaleTerm = TUDAT_NAN;
        for( unsigned int i = 0; i <= maximumStepIndex_; i++ )
        {
            // Compute Euler step and set as state at center point for use with mid-point method.
            stateAtCenterPoint_ = currentState_ + subSteps_.at( i )
                    * this->stateDerivativeFunction_( currentIndependentVariable_, currentState_ );

            // Apply modified mid-point rule.
            stateAtFirstPoint_ = currentState_;
            IndependentVariableType independentVariableAtFirstPoint_ = currentIndependentVariable_;
            for ( unsigned int j = 0; j < sequence_.at( i ) - 1; j++ )
            {
                stateAtLastPoint_ = executeMidPointMethod(
                            stateAtFirstPoint_, stateAtCenterPoint_,
                            independentVariableAtFirstPoint_, subSteps_.at( i ) );

                if ( j < sequence_.at( i ) - 2 )
                {
                    stateAtFirstPoint_ = stateAtCenterPoint_;
                    stateAtCenterPoint_ = stateAtLastPoint_;
                    independentVariableAtFirstPoint_ += subSteps_.at( i );
                }
            }

            // Apply end-point correction.
            integratedStates_[ i ][ 0 ]
                    = 0.5 * ( stateAtLastPoint_ + stateAtCenterPoint_+ subSteps_.at( i ) * this->stateDerivativeFunction_(
                                  currentIndependentVariable_ + stepSize, stateAtLastPoint_ ) );

            for ( unsigned int k = 1; k < i + 1; k++ )
            {
                integratedStates_[ i ][ k ] =
                        integratedStates_[ i ][ k - 1 ] + 1.0 /
                        ( pow( subSteps_.at( i - k ), 2.0 ) / std::pow( subSteps_.at( i ), 2.0 ) - 1.0 )
                        * ( integratedStates_[ i ][ k - 1 ] - integratedStates_[ i - 1 ][ k - 1 ] );
            }
        }

        if( computeNextStepSizeAndValidateResult(
            integratedStates_.at( maximumStepIndex_ ).at( maximumStepIndex_ - 1 ),
            integratedStates_.at( maximumStepIndex_ ).at( maximumStepIndex_ ), stepSize ) )
        {
            this->lastIndependentVariable_ = this->currentIndependentVariable_;
            this->lastState_ = this->currentState_;
            this->currentIndependentVariable_ += stepSize;
            currentState_ = integratedStates_[ maximumStepIndex_ ][ maximumStepIndex_ ];
        }
        else
        {
            performIntegrationStep( this->stepSize_ );
        }

        return currentState_;
    }

    bool computeNextStepSizeAndValidateResult(
        const StateType& lowerOrderEstimate,
        const StateType& higherOrderEstimate,
        const TimeStepType stepSize )
    {
        if( !useFixedStep_ )
        {
            // Compute new step size using new step size function, which also returns whether the
            // relative error is within bounds or not.
            std::pair< TimeStepType, bool > recommendedNewStepSizePair = stepSizeController_->computeNewStepSize(
                lowerOrderEstimate, higherOrderEstimate, stepSize );
            std::pair< TimeStepType, bool > validatedNewStepSizePair = stepSizeValidator_->validateStep(
                recommendedNewStepSizePair, stepSize );

            this->stepSize_ = validatedNewStepSizePair.first;
            return validatedNewStepSizePair.second;
        }
        else
        {
            this->stepSize_ = stepSize;
            return true;
        }
    }

    // Rollback internal state to the last state.
    /*
     * Performs rollback of the internal state to the last state. This function can only be called
     * once after calling integrateTo( ) or performIntegrationStep( ) unless specified otherwise by
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

    IndependentVariableType getPreviousIndependentVariable( )
    {
        return lastIndependentVariable_;
    }

    // Get previous state value.
    /*
     * Returns the previous value of the state.
     * \return Previous state
     */
    StateType getPreviousState( )
    {
        return lastState_;
    }

    void modifyCurrentState( const StateType& newState, const bool allowRollback = false )
    {
        currentState_ = newState;
        if ( !allowRollback )
        {
            this->lastIndependentVariable_ = currentIndependentVariable_;
        }
    }

private:

    // Current independent variable.
    /*
     * Current independent variable as computed by performIntegrationStep().
     */
    IndependentVariableType currentIndependentVariable_;

    // Current state.
    /*
     * Current state as computed by performIntegrationStep( ).
     */
    StateType currentState_;

    // Last independent variable.
    /*
     * Last independent variable value as computed by performIntegrationStep().
     */
    IndependentVariableType lastIndependentVariable_;

    // Last state.
    /*
     * Last state as computed by performIntegrationStep( ).
     */
    StateType lastState_;

    // Sequence for the integrator.
    /*
     * Rational function sequence for the integrator.
     */
    std::vector< unsigned int > sequence_;


    // Last used step size.
    /*
     * Last used step size, passed to either integrateTo( ) or performIntegrationStep( ).
     */
    TimeStepType stepSize_;


    // Execute mid-point method.
    /*
     * Executes mid-point method, given a known state and state derivative.
     * \param stateAtFirstPoint State at first point.
     * \param stateAtCenterPoint State at center point.
     * \param independentVariableAtFirstPoint Independent variable at first point.
     * \param subStepSize Sub step size between successive states used by mid-point method.
     * \return Result of midpoint method
     */
    StateType executeMidPointMethod( StateType stateAtFirstPoint, StateType stateAtCenterPoint,
                                     const IndependentVariableType independentVariableAtFirstPoint,
                                     const IndependentVariableType subStepSize )
    {
        return stateAtFirstPoint + 2.0 * subStepSize
                * this->stateDerivativeFunction_( independentVariableAtFirstPoint + subStepSize,
                                                  stateAtCenterPoint );
    }

    std::vector< std::vector< StateType > > integratedStates_;

    unsigned int maximumStepIndex_;

    std::vector< double > subSteps_;


    std::shared_ptr< IntegratorStepSizeController< TimeStepType, StateType > > stepSizeController_;

    std::shared_ptr< IntegratorStepSizeValidator< TimeStepType > > stepSizeValidator_;

    bool useFixedStep_;

};

extern template class BulirschStoerVariableStepSizeIntegrator < double, Eigen::VectorXd, Eigen::VectorXd >;
extern template class BulirschStoerVariableStepSizeIntegrator < double, Eigen::Vector6d, Eigen::Vector6d >;
extern template class BulirschStoerVariableStepSizeIntegrator < double, Eigen::MatrixXd, Eigen::MatrixXd >;


// Typedef of variable-step size Bulirsch-Stoer integrator (state/state derivative = VectorXd,
// independent variable = double).
/*
 * Typedef of a variable-step size Bulirsch-Stoer integrator with VectorXds as state and state
 * derivative and double as independent variable.
 */
typedef BulirschStoerVariableStepSizeIntegrator< > BulirschStoerVariableStepSizeIntegratorXd;

} // namespace numerical_integrators

} // namespace tudat

#endif // TUDAT_BULIRSCH_STOER_VARIABLE_STEP_SIZE_INTEGRATOR_H
