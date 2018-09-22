/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INTEGRATEEQUATIONS_H
#define TUDAT_INTEGRATEEQUATIONS_H

#include <Eigen/Core>
#include <boost/lambda/lambda.hpp>
#include <chrono>
#include <limits>

#include <map>

#include "Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Basics/timeType.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/Mathematics/RootFinders/createRootFinder.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTermination.h"

namespace tudat
{

namespace propagators
{

//! Function to determine, for a given time step of the numerical integrator, the error in termination dependent variable
/*!
 *  Function to determine, for a given time step of the numerical integrator, the error in termination dependent variable. This
 *  function is used as input for the root finder when the propagation must terminate exactly on a dependent variable value
 *  \param timeStep Time step to take with the numerical integrator
 *  \param integrator Numerical integrator used for propagation
 *  \param dependentVariableTerminationCondition Settings used to determine value/type of dependent variable at which propagation
 *  is to terminate
 *  \return The difference between the reached and required value of the termination dependent variable
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double, typename TimeStepType = TimeType  >
TimeStepType getTerminationDependentVariableErrorForGivenTimeStep(
        TimeStepType timeStep,
        const std::shared_ptr< numerical_integrators::NumericalIntegrator< TimeType, StateType, StateType, TimeStepType > >
        integrator,
        const std::shared_ptr< SingleVariableLimitPropagationTerminationCondition > dependentVariableTerminationCondition )
{
    // Perform integration step
    integrator->performIntegrationStep( timeStep );

    // Retrieve value of dependent variable after time step
    integrator->getStateDerivativeFunction( )(
                integrator->getCurrentIndependentVariable( ), integrator->getCurrentState( ) );
    TimeStepType dependentVariableError =
            static_cast< TimeStepType >( dependentVariableTerminationCondition->getStopConditionError( ) );

    // Undo step
    integrator->rollbackToPreviousState( );

    return dependentVariableError;
}

//! Function that propagates to an exact final condition (within tolerance) for dependent variable termination condition
/*!
 * Function that propagates to an exact final condition (within tolerance) for dependent variable termination condition.
 * Determines the time step that is to be taken by using a root finder, and returns (by reference) the converged final time
 * and state.
 * \param integrator Numerical integrator that is used for propagation. Upon input to this function, the integrator is rolled
 * back to the secondToLastTime/secondToLastState
 * \param dependentVariableTerminationCondition Termination condition that is to be used
 * \param secondToLastTime Second to last time (e.g. last time at which integration did not exceed termination condition)
 * \param lastTime Time at which integration first exceeded termination condition
 * \param secondToLastState Second to last state (e.g. state at last time where integration did not exceed termination condition)
 * \param lastState State at time where integration first exceeded termination condition
 * \param endTime Time at which exact termination condition is met (returned by reference).
 * \param endState State at time where exact termination condition is met (returned by reference).
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double, typename TimeStepType = TimeType  >
void getFinalStateForExactDependentVariableTerminationCondition(
        const std::shared_ptr< numerical_integrators::NumericalIntegrator< TimeType, StateType, StateType, TimeStepType > >
        integrator,
        const std::shared_ptr< SingleVariableLimitPropagationTerminationCondition > dependentVariableTerminationCondition,
        const TimeType secondToLastTime,
        const TimeType lastTime,
        const StateType& secondToLastState,
        const StateType& lastState,
        TimeType& endTime,
        StateType& endState )
{
    TUDAT_UNUSED_PARAMETER( secondToLastState );

    // Function for which the root (zero value) occurs at the required end time/state
    std::function< TimeStepType( TimeStepType ) > dependentVariableErrorFunction =
            std::bind( &getTerminationDependentVariableErrorForGivenTimeStep< StateType, TimeType, TimeStepType >, std::placeholders::_1,
                       integrator, dependentVariableTerminationCondition );

    // Create root finder.
    bool increasingTime = static_cast< double >( lastTime - secondToLastTime ) > 0.0;
    std::shared_ptr< root_finders::RootFinderCore< TimeStepType > > finalConditionRootFinder;
    if( increasingTime )
    {
        finalConditionRootFinder = root_finders::createRootFinder< TimeStepType >(
                    dependentVariableTerminationCondition->getTerminationRootFinderSettings( ),
                    static_cast< TimeStepType >( std::numeric_limits< double >::min( ) ),
                    static_cast< TimeStepType >( lastTime - secondToLastTime ),
                    static_cast< TimeStepType >( std::numeric_limits< double >::min( ) ) );
    }
    else
    {
        finalConditionRootFinder = root_finders::createRootFinder< TimeStepType >(
                    dependentVariableTerminationCondition->getTerminationRootFinderSettings( ),
                    static_cast< TimeStepType >( lastTime - secondToLastTime ),
                    static_cast< TimeStepType >( -std::numeric_limits< double >::min( ) ),
                    static_cast< TimeStepType >( lastTime - secondToLastTime ) );

    }

    // Solve root-finding problem.
    TimeStepType finalTimeStep;
    try
    {
        finalTimeStep = finalConditionRootFinder->execute(
                    std::make_shared< basic_mathematics::FunctionProxy< TimeStepType, TimeStepType > >(
                        dependentVariableErrorFunction ), ( lastTime - secondToLastTime ) / 2.0 );

        endState = integrator->performIntegrationStep( finalTimeStep );
        endTime = integrator->getCurrentIndependentVariable( );
    }
    // If dependent variable has no root in given interval, set end time and state at NaN
    catch( std::runtime_error& caughtException )
    {
        std::cerr << "Warning in propagation to exact dependent variable value. Root finder could not find a "
                     "root to the function. Returning time and state as NaNs. Caught exception: "
                  << caughtException.what( ) << std::endl;
        endTime = TUDAT_NAN;
        endState = StateType::Constant( lastState.rows( ), lastState.cols( ), TUDAT_NAN );
    }
}

//! Function that propagates to an exact final condition (within tolerance) for hybrid termination condition
/*!
 * Function that propagates to an exact final condition (within tolerance) for hybrid termination condition. Determines
 * the termination time/state for each of the constituent termination condition, and chooses the hybrid termiantion time/state
 * accordingly
 * \param integrator Numerical integrator that is used for propagation. Upon input to this function, the integrator is rolled
 * back to the secondToLastTime/secondToLastState
 * \param hyrbidTerminationCondition Termination condition that is to be used
 * \param secondToLastTime Second to last time (e.g. last time at which integration did not exceed termination condition)
 * \param lastTime Time at which integration first exceeded termination condition
 * \param secondToLastState Second to last state (e.g. state at last time where integration did not exceed termination condition)
 * \param lastState State at time where integration first exceeded termination condition
 * \param endTime Time at which exact termination condition is met (returned by reference).
 * \param endState State at time where exact termination condition is met (returned by reference).
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double, typename TimeStepType = TimeType  >
void getFinalStateForExactHybridVariableTerminationCondition(
        const std::shared_ptr< numerical_integrators::NumericalIntegrator< TimeType, StateType, StateType, TimeStepType > >
        integrator,
        const std::shared_ptr< HybridPropagationTerminationCondition > hyrbidTerminationCondition,
        const TimeType secondToLastTime,
        const TimeType lastTime,
        const StateType& secondToLastState,
        const StateType& lastState,
        TimeType& endTime,
        StateType& endState )
{

    std::vector< std::shared_ptr< PropagationTerminationCondition > > terminationConditionList =
            hyrbidTerminationCondition->getPropagationTerminationConditions( );

    // Create list of converged times/states for each constituent condition
    std::vector< TimeType > endTimes;
    endTimes.resize( terminationConditionList.size( ) );
    std::vector< StateType > endStates;
    endStates.resize( terminationConditionList.size( ) );

    // Iterate over all constituent termination conditions, and determine separate end times/states
    unsigned int minimumTimeIndex = 0, maximumTimeIndex = 0;
    TimeStepType minimumTimeStep, maximumTimeStep;
    bool timesAreSet = false;
    for( unsigned int i = 0; i < terminationConditionList.size( ); i++ )
    {
        // Determine single termination condition
        getFinalStateForExactTerminationCondition(
                    integrator, terminationConditionList.at( i ),secondToLastTime, lastTime, secondToLastState, lastState,
                    endTimes[ i ], endStates[ i ] );

        // If converged time is found, check if it is smallest/highest converged time
        if( endTimes[ i ] == endTimes[ i ] )
        {
            TimeStepType currentFinalTimeStep = endTimes[ i ] - secondToLastTime;

            if( !timesAreSet )
            {
                maximumTimeStep = currentFinalTimeStep;
                maximumTimeIndex = i;
                minimumTimeStep = currentFinalTimeStep;
                minimumTimeIndex = i;
                timesAreSet = true;
            }
            else
            {
                if( currentFinalTimeStep < minimumTimeStep )
                {
                    minimumTimeStep = currentFinalTimeStep;
                    minimumTimeIndex = i;
                }

                if( currentFinalTimeStep > maximumTimeStep )
                {
                    maximumTimeStep = currentFinalTimeStep;
                    maximumTimeIndex = i;
                }
            }
        }
    }

    // Check if any of the conditions were valid
    if( !timesAreSet )
    {
        throw std::runtime_error( "Error when propagating exactly to hybrid condition, no conditions met" );
    }

    // Set converged end time/state
    bool propagationIsForwards = ( ( lastTime - secondToLastTime ) > 0.0 ) ? true : false;
    if( ( propagationIsForwards && !hyrbidTerminationCondition->getFulfillSingleCondition( ) ) ||
            ( !propagationIsForwards && hyrbidTerminationCondition->getFulfillSingleCondition( ) ) )
    {
        endState = endStates[ maximumTimeIndex ];
        endTime = endTimes[ maximumTimeIndex ];
    }
    else if( ( propagationIsForwards && hyrbidTerminationCondition->getFulfillSingleCondition( ) ) ||
             ( !propagationIsForwards && !hyrbidTerminationCondition->getFulfillSingleCondition( ) ) )
    {
        endState = endStates[ minimumTimeIndex ];
        endTime = endTimes[ minimumTimeIndex ];
    }
    else
    {
        throw std::runtime_error( "Error when propagating to exact final hybrid condition, case not recognized" );
    }

}

//! Function that propagates to an exact final condition (within tolerance) for arbitrary termination condition
/*!
 * Function that propagates to an exact final condition (within tolerance) for arbitrary termination condition.
 * \param integrator Numerical integrator that is used for propagation. Upon input to this function, the integrator is at
 * the lastTime/lastState
 * \param terminationCondition Termination condition that is to be used
 * \param secondToLastTime Second to last time (e.g. last time at which integration did not exceed termination condition)
 * \param lastTime Time at which integration first exceeded termination condition
 * \param secondToLastState Second to last state (e.g. state at last time where integration did not exceed termination condition)
 * \param lastState State at time where integration first exceeded termination condition
 * \param endTime Time at which exact termination condition is met (returned by reference).
 * \param endState State at time where exact termination condition is met (returned by reference).
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double, typename TimeStepType = TimeType  >
void getFinalStateForExactTerminationCondition(
        const std::shared_ptr< numerical_integrators::NumericalIntegrator< TimeType, StateType, StateType, TimeStepType > > integrator,
        const std::shared_ptr< PropagationTerminationCondition > terminationCondition,
        const TimeType secondToLastTime,
        const TimeType lastTime,
        const StateType& secondToLastState,
        const StateType& lastState,
        TimeType& endTime,
        StateType& endState )
{
    // Check type of termination condition
    switch( terminationCondition->getTerminationType( ) )
    {
    case time_stopping_condition:
    {
        // Check input consistency
        std::shared_ptr< FixedTimePropagationTerminationCondition > timeTerminationCondition =
                std::dynamic_pointer_cast< FixedTimePropagationTerminationCondition >( terminationCondition );

        // Determine final time step and propagate
        TimeStepType finalTimeStep = timeTerminationCondition->getStopTime( ) - secondToLastTime;

        integrator->rollbackToPreviousState( );
        endState = integrator->performIntegrationStep( finalTimeStep );
        endTime = integrator->getCurrentIndependentVariable( );

        break;
    }
    case cpu_time_stopping_condition:
    {
        // No exact final condition on CPU time is possible
        std::cerr << "Error, cannot propagate to exact CPU time, returning state after condition violation." << std::endl;

        endTime = lastTime;
        endState = lastState;
        break;
    }
    case dependent_variable_stopping_condition:
    {
        integrator->rollbackToPreviousState( );

        std::shared_ptr< SingleVariableLimitPropagationTerminationCondition > dependentVariableTerminationCondition =
                std::dynamic_pointer_cast< SingleVariableLimitPropagationTerminationCondition >( terminationCondition );
        getFinalStateForExactDependentVariableTerminationCondition(
                    integrator, dependentVariableTerminationCondition, secondToLastTime, lastTime,
                    secondToLastState, lastState, endTime, endState );

        break;
    }
    case custom_stopping_condition:
    {
        // No exact final condition on custom condition is possible
        std::cerr << "Error, cannot propagate to exact custom condition, returning state after condition violation." << std::endl;

        endTime = lastTime;
        endState = lastState;
        break;
    }
    case hybrid_stopping_condition:
    {
        std::shared_ptr< HybridPropagationTerminationCondition > hyrbidTerminationCondition =
                std::dynamic_pointer_cast< HybridPropagationTerminationCondition >( terminationCondition );

        getFinalStateForExactHybridVariableTerminationCondition(
                    integrator, hyrbidTerminationCondition, secondToLastTime, lastTime,
                    secondToLastState, lastState, endTime, endState );
        break;
    }
    default:
        throw std::runtime_error( "Error when propagating to exact final condition, did not recognize termination time" );
    }
}

//! Function that propagates to an exact final condition (within tolerance) for arbitrary termination condition
/*!
 * Function that propagates to an exact final condition (within tolerance) for arbitrary termination condition.
 * \param integrator Numerical integrator that is used for propagation. Upon input to this function, the integrator is at
 * the final time/state encountered by the propagation
 * \param propagationTerminationCondition Termination condition that is to be used
 * \param timeStep Last time step taken by integrator.
 * \param dependentVariableFunction Function returning dependent variables (obtained from environment and state
 * derivative model).
 * \param solutionHistory History of state variables that are to be saved given as map
 * (time as key; returned by reference)
 * \param dependentVariableHistory History of dependent variables that are to be saved given as map
 * (time as key; returned by reference)
 * \param currentCpuTime Current run time of propagation.
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double, typename TimeStepType = TimeType  >
void propagateToExactTerminationCondition(
        const std::shared_ptr< numerical_integrators::NumericalIntegrator< TimeType, StateType, StateType, TimeStepType > > integrator,
        const std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition,
        const TimeStepType timeStep,
        const std::function< Eigen::VectorXd( ) > dependentVariableFunction,
        std::map< TimeType, StateType >& solutionHistory,
        std::map< TimeType, Eigen::VectorXd >& dependentVariableHistory,
        const double currentCpuTime )
{
    // Turn off step size control
    integrator->setStepSizeControl( false );

    // Determine exact final time/state
    TimeType endTime;
    StateType endState;
    getFinalStateForExactTerminationCondition(
                integrator, propagationTerminationCondition,
                integrator->getPreviousIndependentVariable( ),
                integrator->getCurrentIndependentVariable( ),
                integrator->getPreviousState( ),
                integrator->getCurrentState( ),
                endTime, endState );

    // Check if any dependent variables are saved. If so, remove last entry
    bool recomputeDependentVariables = false;
    if( dependentVariableHistory.size( ) > 0 )
    {
        if( dependentVariableHistory.rbegin( )->first == solutionHistory.rbegin( )->first )
        {
            if( timeStep > 0 )
            {
                dependentVariableHistory.erase( std::prev( dependentVariableHistory.end() ) );
            }
            else
            {
                dependentVariableHistory.erase(  dependentVariableHistory.begin( ) );
            }
            recomputeDependentVariables = true;
        }
    }

    // Remove state entry last added, and enter converged final state
    if( timeStep > 0 )
    {
        solutionHistory.erase( std::prev( solutionHistory.end() ) );
        solutionHistory[ endTime ] = endState;
    }
    else
    {
        solutionHistory.erase( solutionHistory.begin( ) );
        solutionHistory[ endTime ] = endState;
    }

    // Recompute final dependent variables, if required
    if( recomputeDependentVariables )
    {
        integrator->getStateDerivativeFunction( )( endTime, endState );
        dependentVariableHistory[ endTime ] = dependentVariableFunction( );

        // Check stopping conditions to be able to save details
        propagationTerminationCondition->checkStopCondition( endTime, currentCpuTime );
    }

    // Turn step size control back on
    integrator->setStepSizeControl( true );
}

//! Function to numerically integrate a given first order differential equation
/*!
 *  Function to numerically integrate a given first order differential equation, with the state derivative a function of
 *  a single independent variable and the current state
 *  \param integrator Numerical integrator used for propagation
 *  \param initialTimeStep Time step to use for first step of numerical integration
 *  \param propagationTerminationCondition Object to determine when/how the propagation is to be stopped at the current time.
 *  \param solutionHistory History of state variables that are to be saved given as map
 *  (time as key; returned by reference)
 *  \param dependentVariableHistory History of dependent variables that are to be saved given as map
 *  (time as key; returned by reference)
 *  \param cumulativeComputationTimeHistory History of cumulative computation times that are to be saved given
 *  as map (time as key; returned by reference)
 *  \param dependentVariableFunction Function returning dependent variables (obtained from environment and state
 *  derivative model).
 *  \param statePostProcessingFunction Function to post-process state after numerical integration (obtained from state derivative model).
 *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration time
 *  steps, with n = saveFrequency).
 *  \param printInterval Frequency with which to print progress to console (nan = never).
 *  \param initialClockTime Initial clock time from which to determine cumulative computation time.
 *  By default now(), i.e. the moment at which this function is called.
 *  \return Event that triggered the termination of the propagation
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double, typename TimeStepType = TimeType  >
std::shared_ptr< PropagationTerminationDetails > integrateEquationsFromIntegrator(
        const std::shared_ptr< numerical_integrators::NumericalIntegrator< TimeType, StateType, StateType, TimeStepType > > integrator,
        const TimeStepType initialTimeStep,
        const std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition,
        std::map< TimeType, StateType >& solutionHistory,
        std::map< TimeType, Eigen::VectorXd >& dependentVariableHistory,
        std::map< TimeType, double >& cumulativeComputationTimeHistory,
        const std::function< Eigen::VectorXd( ) > dependentVariableFunction = std::function< Eigen::VectorXd( ) >( ),
        const std::function< void( StateType& ) > statePostProcessingFunction = std::function< void( StateType& ) >( ),
        const int saveFrequency = TUDAT_NAN,
        const TimeType printInterval = TUDAT_NAN,
        const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) )
{
    std::shared_ptr< PropagationTerminationDetails > propagationTerminationReason;

    // Get Initial state and time.
    TimeType currentTime = integrator->getCurrentIndependentVariable( );
    TimeType initialTime = currentTime;
    StateType newState = integrator->getCurrentState( );

    // Initialization of numerical solutions for variational equations
    solutionHistory.clear( );
    solutionHistory[ currentTime ] = newState;

    dependentVariableHistory.clear( );
    if( !( dependentVariableFunction == nullptr ) )
    {
        integrator->getStateDerivativeFunction( )( currentTime, newState );
        dependentVariableHistory[ currentTime ] = dependentVariableFunction( );
    }

    // CPU time
    cumulativeComputationTimeHistory.clear( );
    double currentCPUTime = std::chrono::duration_cast< std::chrono::nanoseconds >(
                std::chrono::steady_clock::now( ) - initialClockTime ).count( ) * 1.0e-9;
    cumulativeComputationTimeHistory[ currentTime ] = currentCPUTime;

    // Set initial time step and total integration time.
    TimeStepType timeStep = initialTimeStep;
    TimeType previousTime = currentTime;
    TimeType previousPrintTime = TUDAT_NAN;

    int saveIndex = 0;

    propagationTerminationReason = std::make_shared< PropagationTerminationDetails >(
                unknown_propagation_termination_reason );
    bool breakPropagation = 0;

    // Perform numerical integration steps until end time reached.
    do
    {
        try
        {
            if( ( newState.allFinite( ) == true ) && ( !newState.hasNaN( ) ) )
            {
                previousTime = currentTime;

                // Perform integration step.
                newState = integrator->performIntegrationStep( timeStep );
                if( statePostProcessingFunction != nullptr )
                {
                    statePostProcessingFunction( newState );
                    integrator->modifyCurrentState( newState, true );
                }

                // Check if the termination condition was reached during evaluation of integration sub-steps.
                // If evaluation of the termination condition during integration sub-steps is disabled,
                // this function returns always `false`.
                // If the termination condition was reached, the last step could not be computed correctly because some
                // of the integrator sub-steps were not computed. Thus, return immediately without saving the `newState`.
                if( integrator->getPropagationTerminationConditionReached( ) )
                {
                    propagationTerminationReason = std::make_shared< PropagationTerminationDetails >(
                                termination_condition_reached );
                    break;
                }

                // Update epoch and step-size
                currentTime = integrator->getCurrentIndependentVariable( );
                timeStep = integrator->getNextStepSize( );

                // Save integration result in map
                saveIndex++;
                saveIndex = saveIndex % saveFrequency;
                if( saveIndex == 0 )
                {
                    solutionHistory[ currentTime ] = newState;

                    if( !( dependentVariableFunction == nullptr ) )
                    {
                        integrator->getStateDerivativeFunction( )( currentTime, newState );
                        dependentVariableHistory[ currentTime ] = dependentVariableFunction( );
                    }
                }
            }
            else
            {
                std::cerr << "Error, propagation terminated at t=" + std::to_string( static_cast< double >( currentTime ) ) +
                             ", found NaN/Inf entry, returning propagation data up to current time" << std::endl;
                breakPropagation = true;
                propagationTerminationReason = std::make_shared< PropagationTerminationDetails >(
                            nan_or_inf_detected_in_state );
            }

            currentCPUTime = std::chrono::duration_cast< std::chrono::nanoseconds >(
                        std::chrono::steady_clock::now( ) - initialClockTime ).count( ) * 1.0e-9;
            cumulativeComputationTimeHistory[ currentTime ] = currentCPUTime;

            // Print solutions
            if( printInterval == printInterval )
            {
                if( !( previousPrintTime == previousPrintTime ) ||
                        std::fabs( static_cast< double >( currentTime - previousPrintTime ) ) > printInterval )
                {
                    previousPrintTime = currentTime;
                    std::cout << "Current time: "<<currentTime<<"; time since initial time: "<<currentTime - initialTime<<std::endl;
                    std::cout << "Current time step: "<<timeStep<<std::endl;
                    std::cout << "Current state (transposed): "<<newState.transpose( ) <<std::endl<<std::endl;
                }
            }

            if( propagationTerminationCondition->checkStopCondition( static_cast< double >( currentTime ), currentCPUTime ) )
            {
                if( propagationTerminationCondition->getTerminateExactlyOnFinalCondition( ) )
                {
                    propagateToExactTerminationCondition(
                                integrator, propagationTerminationCondition,
                                timeStep, dependentVariableFunction,
                                solutionHistory, dependentVariableHistory, currentCPUTime );
                }

                // Set termination details
                if( propagationTerminationCondition->getTerminationType( ) != hybrid_stopping_condition )
                {
                    propagationTerminationReason = std::make_shared< PropagationTerminationDetails >(
                                termination_condition_reached,
                                propagationTerminationCondition->getTerminateExactlyOnFinalCondition( ) );
                }
                else
                {
                    if( std::dynamic_pointer_cast< HybridPropagationTerminationCondition >( propagationTerminationCondition )
                            == nullptr )
                    {
                        throw std::runtime_error( "Error when saving termination reason, type is hybrid, but class is not." );
                    }
                    propagationTerminationReason = std::make_shared< PropagationTerminationDetailsFromHybridCondition >(
                                propagationTerminationCondition->getTerminateExactlyOnFinalCondition( ),
                                std::dynamic_pointer_cast< HybridPropagationTerminationCondition >(
                                    propagationTerminationCondition ) );
                }
                breakPropagation = true;
            }
        }
        catch( const std::exception& caughtException )
        {
            std::cerr << caughtException.what( ) << std::endl;
            std::cerr << "Error, propagation terminated at t=" + std::to_string( static_cast< double >( currentTime ) ) +
                         ", returning propagation data up to current time." << std::endl;
            breakPropagation = true;
            propagationTerminationReason = std::make_shared< PropagationTerminationDetails >(
                        runtime_error_caught_in_propagation );
        }
    }
    while( !breakPropagation );

    return propagationTerminationReason;
}

extern template std::shared_ptr< PropagationTerminationDetails > integrateEquationsFromIntegrator<
Eigen::MatrixXd, double, double >(
        const std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::MatrixXd, Eigen::MatrixXd, double > > integrator,
        const double initialTimeStep,
        const std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition,
        std::map< double, Eigen::MatrixXd >& solutionHistory,
        std::map< double, Eigen::VectorXd >& dependentVariableHistory,
        std::map< double, double >& cumulativeComputationTimeHistory,
        const std::function< Eigen::VectorXd( ) > dependentVariableFunction,
        const std::function< void( Eigen::MatrixXd& ) > statePostProcessingFunction,
        const int saveFrequency,
        const double printInterval,
        const std::chrono::steady_clock::time_point initialClockTime );


extern template std::shared_ptr< PropagationTerminationDetails > integrateEquationsFromIntegrator<
Eigen::VectorXd, double, double >(
        const std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::VectorXd, Eigen::VectorXd, double > > integrator,
        const double initialTimeStep,
        const std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition,
        std::map< double, Eigen::VectorXd >& solutionHistory,
        std::map< double, Eigen::VectorXd >& dependentVariableHistory,
        std::map< double, double >& cumulativeComputationTimeHistory,
        const std::function< Eigen::VectorXd( ) > dependentVariableFunction,
        const std::function< void( Eigen::VectorXd& ) > statePostProcessingFunction,
        const int saveFrequency,
        const double printInterval,
        const std::chrono::steady_clock::time_point initialClockTime );


//! Interface class for integrating some state derivative function.
/*!
 *  Interface class for integrating some state derivative function.. This class is used instead of a single templated free
 *  function to allow ObservationModel the integrator etc. to adapt its time step variable to long double if the Time
 *  object is used as TimeType. This class has template specializations for double/Time TimeType, and contains a single
 *  integrateEquations function that performs the required operation.
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double >
class EquationIntegrationInterface
{
public:

    //! Function to numerically integrate a given first order differential equation
    /*!
     *  Function to numerically integrate a given first order differential equation, with the state derivative a function of
     *  a single independent variable and the current state
     *  \param stateDerivativeFunction Function returning the state derivative from current time and state.
     *  \param solutionHistory History of numerical states given as map (time as key; returned by reference)
     *  \param initialState Initial state
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagationTerminationCondition Object to determine when/how the propagation is to be stopped at the current time.
     *  \param dependentVariableHistory History of dependent variables that are to be saved given as map
     *  (time as key; returned by reference)
     *  \param cumulativeComputationTimeHistory History of cumulative computation times that are to be saved given
     *  as map (time as key; returned by reference)
     *  \param dependentVariableFunction Function returning dependent variables (obtained from environment and state
     *  derivative model).
     *  \param statePostProcessingFunction Function to post-process state after numerical integration (obtained from state derivative model).
     *  \param printInterval Frequency with which to print progress to console (nan = never).
     *  \param initialClockTime Initial clock time from which to determine cumulative computation time.
     *  By default now(), i.e. the moment at which this function is called.
     *  \return Event that triggered the termination of the propagation
     */
    static std::shared_ptr< PropagationTerminationDetails > integrateEquations(
            std::function< StateType( const TimeType, const StateType& ) > stateDerivativeFunction,
            std::map< TimeType, StateType >& solutionHistory,
            const StateType initialState,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition,
            std::map< TimeType, Eigen::VectorXd >& dependentVariableHistory,
            std::map< TimeType, double >& cumulativeComputationTimeHistory,
            const std::function< Eigen::VectorXd( ) > dependentVariableFunction = std::function< Eigen::VectorXd( ) >( ),
            const std::function< void( StateType& ) > statePostProcessingFunction = std::function< void( StateType& ) >( ),
            const TimeType printInterval = TUDAT_NAN,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) );

};

//! Interface class for integrating some state derivative function.
template< typename StateType >
class EquationIntegrationInterface< StateType, double >
{
public:

    //! Function to numerically integrate a given first order differential equation
    /*!
     *  Function to numerically integrate a given first order differential equation, with the state derivative a function of
     *  a single independent variable and the current state
     *  \param stateDerivativeFunction Function returning the state derivative from current time and state.
     *  \param solutionHistory History of numerical states given as map (time as key; returned by reference)
     *  \param initialState Initial state
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagationTerminationCondition Object to determine when/how the propagation is to be stopped at the current time.
     *  \param dependentVariableHistory History of dependent variables that are to be saved given as map
     *  (time as key; returned by reference)
     *  \param cumulativeComputationTimeHistory History of cumulative computation times that are to be saved given
     *  as map (time as key; returned by reference)
     *  \param dependentVariableFunction Function returning dependent variables (obtained from environment and state
     *  derivative model).
     *  \param statePostProcessingFunction Function to post-process state after numerical integration (obtained from state derivative model).
     *  \param printInterval Frequency with which to print progress to console (nan = never).
     *  \param initialClockTime Initial clock time from which to determine cumulative computation time.
     *  By default now(), i.e. the moment at which this function is called.
     *  \return Event that triggered the termination of the propagation
     */
    static std::shared_ptr< PropagationTerminationDetails > integrateEquations(
            std::function< StateType( const double, const StateType& ) > stateDerivativeFunction,
            std::map< double, StateType >& solutionHistory,
            const StateType initialState,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
            const std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition,
            std::map< double, Eigen::VectorXd >& dependentVariableHistory,
            std::map< double, double >& cumulativeComputationTimeHistory,
            const std::function< Eigen::VectorXd( ) > dependentVariableFunction = std::function< Eigen::VectorXd( ) >( ),
            const std::function< void( StateType& ) > statePostProcessingFunction = std::function< void( StateType& ) >( ),
            const double printInterval = TUDAT_NAN,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) )
    {
        std::function< bool( const double, const double ) > stopPropagationFunction =
                std::bind( &PropagationTerminationCondition::checkStopCondition, propagationTerminationCondition, std::placeholders::_1, std::placeholders::_2 );

        // Create numerical integrator.
        std::shared_ptr< numerical_integrators::NumericalIntegrator< double, StateType, StateType > > integrator =
                numerical_integrators::createIntegrator< double, StateType >(
                    stateDerivativeFunction, initialState, integratorSettings );

        if( integratorSettings->assessPropagationTerminationConditionDuringIntegrationSubsteps_ )
        {
            integrator->setPropagationTerminationFunction( stopPropagationFunction );
        }

        return integrateEquationsFromIntegrator< StateType, double >(
                    integrator, integratorSettings->initialTimeStep_, propagationTerminationCondition, solutionHistory,
                    dependentVariableHistory,
                    cumulativeComputationTimeHistory,
                    dependentVariableFunction,
                    statePostProcessingFunction,
                    integratorSettings->saveFrequency_,
                    printInterval,
                    initialClockTime );
    }

};

//! Interface class for integrating some state derivative function.
template< typename StateType >
class EquationIntegrationInterface< StateType, Time >
{
public:

    //! Function to numerically integrate a given first order differential equation
    /*!
     *  Function to numerically integrate a given first order differential equation, with the state derivative a function of
     *  a single independent variable and the current state
     *  \param stateDerivativeFunction Function returning the state derivative from current time and state.
     *  \param solutionHistory History of numerical states given as map (time as key; returned by reference)
     *  \param initialState Initial state
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagationTerminationCondition Object to determine when/how the propagation is to be stopped at the current time.
     *  \param dependentVariableHistory History of dependent variables that are to be saved given as map
     *  (time as key; returned by reference)
     *  \param cumulativeComputationTimeHistory History of cumulative computation times that are to be saved given
     *  as map (time as key; returned by reference)
     *  \param dependentVariableFunction Function returning dependent variables (obtained from environment and state
     *  derivative model).
     *  \param statePostProcessingFunction Function to post-process state after numerical integration (obtained from state derivative model).
     *  \param printInterval Frequency with which to print progress to console (nan = never).
     *  \param initialClockTime Initial clock time from which to determine cumulative computation time.
     *  By default now(), i.e. the moment at which this function is called.
     *  \return Event that triggered the termination of the propagation
     */
    static std::shared_ptr< PropagationTerminationDetails > integrateEquations(
            std::function< StateType( const Time, const StateType& ) > stateDerivativeFunction,
            std::map< Time, StateType >& solutionHistory,
            const StateType initialState,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
            const std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition,
            std::map< Time, Eigen::VectorXd >& dependentVariableHistory,
            std::map< Time, double >& cumulativeComputationTimeHistory,
            const std::function< Eigen::VectorXd( ) > dependentVariableFunction = std::function< Eigen::VectorXd( ) >( ),
            const std::function< void( StateType& ) > statePostProcessingFunction = std::function< void( StateType& ) >( ),
            const Time printInterval = TUDAT_NAN,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) )
    {
        std::function< bool( const double, const double ) > stopPropagationFunction =
                std::bind( &PropagationTerminationCondition::checkStopCondition, propagationTerminationCondition, std::placeholders::_1, std::placeholders::_2 );

        // Create numerical integrator.
        std::shared_ptr< numerical_integrators::NumericalIntegrator< Time, StateType, StateType, long double > > integrator =
                numerical_integrators::createIntegrator< Time, StateType, long double  >(
                    stateDerivativeFunction, initialState, integratorSettings );

        if( integratorSettings->assessPropagationTerminationConditionDuringIntegrationSubsteps_ )
        {
            integrator->setPropagationTerminationFunction( stopPropagationFunction );
        }

        return integrateEquationsFromIntegrator< StateType, Time, long double >(
                    integrator, integratorSettings->initialTimeStep_, propagationTerminationCondition, solutionHistory,
                    dependentVariableHistory,
                    cumulativeComputationTimeHistory,
                    dependentVariableFunction,
                    statePostProcessingFunction,
                    integratorSettings->saveFrequency_,
                    printInterval,
                    initialClockTime );
    }

};

} // namespace propagators

} // namespace tudat

#endif // TUDAT_INTEGRATEEQUATIONS_H
