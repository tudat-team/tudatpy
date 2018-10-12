/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/propagationTermination.h"

namespace tudat
{

namespace propagators
{

//! Function to check whether the propagation is to be be stopped
bool FixedTimePropagationTerminationCondition::checkStopCondition( const double time, const double cpuTime )
{
    bool stopPropagation = false;

    // Check whether stop time has been reached
    if( propagationDirectionIsPositive_ && ( time >= stopTime_ ) )
    {
        stopPropagation = true;
    }
    else if( !propagationDirectionIsPositive_ && ( time <= stopTime_ ) )
    {
        stopPropagation = true;
    }
    return stopPropagation;
}

//! Function to check whether the propagation is to be be stopped
bool FixedCPUTimePropagationTerminationCondition::checkStopCondition( const double time, const double cpuTime )
{
    return cpuTime >= cpuStopTime_;
}

//! Function to check whether the propagation is to be be stopped
bool SingleVariableLimitPropagationTerminationCondition::checkStopCondition( const double time, const double cpuTime  )
{
    bool stopPropagation = false;
    double currentVariable = variableRetrievalFunction_( );

    if( useAsLowerBound_ && ( currentVariable < limitingValue_ ) )
    {
        stopPropagation = true;
    }
    else if( !useAsLowerBound_ && ( currentVariable > limitingValue_ ) )
    {
        stopPropagation = true;
    }
    return stopPropagation;
}

//! Function to check whether the propagation is to be be stopped
bool HybridPropagationTerminationCondition::checkStopCondition( const double time, const double cpuTime )
{
    // Check if single condition is fulfilled.
    bool stopPropagation = -1;
    unsigned int stopIndex = 0;
    if( fulfillSingleCondition_ )
    {
        stopPropagation = false;
        for( unsigned int i = 0; i < propagationTerminationCondition_.size( ); i++ )
        {
            if( propagationTerminationCondition_.at( i )->checkStopCondition( time, cpuTime ) )
            {
                stopIndex = i;
                stopPropagation = true;
                isConditionMetWhenStopping_[ i ] = true;
                break;
            }
            else
            {
                isConditionMetWhenStopping_[ i ] = false;
            }
        }
    }
    // Check all conditions are fulfilled.
    else
    {
        stopPropagation = true;
        for( unsigned int i = 0; i < propagationTerminationCondition_.size( ); i++ )
        {
            if( !propagationTerminationCondition_.at( i )->checkStopCondition( time, cpuTime ) )
            {
                stopIndex = i;
                stopPropagation = false;
                isConditionMetWhenStopping_[ i ] = false;
                break;
            }
            else
            {
                isConditionMetWhenStopping_[ i ] = true;
            }
        }
    }

    // Save if conditions were met
    if( stopPropagation )
    {
        for( unsigned int i = ( stopIndex + 1 ); i < propagationTerminationCondition_.size( ); i++ )
        {
            isConditionMetWhenStopping_[ i ] = propagationTerminationCondition_.at( i )->checkStopCondition( time, cpuTime );
        }
    }

    return stopPropagation;
}


//! Function to create propagation termination conditions from associated settings
std::shared_ptr< PropagationTerminationCondition > createPropagationTerminationConditions(
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double initialTimeStep )
{
    std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition;

    // Check termination type.
    switch( terminationSettings->terminationType_ )
    {
    case time_stopping_condition:
    {
        // Create stopping time termination condition.
        std::shared_ptr< PropagationTimeTerminationSettings > timeTerminationSettings =
                std::dynamic_pointer_cast< PropagationTimeTerminationSettings >( terminationSettings );
        propagationTerminationCondition = std::make_shared< FixedTimePropagationTerminationCondition >(
                    timeTerminationSettings->terminationTime_, ( initialTimeStep > 0 ),
                    timeTerminationSettings->terminateExactlyOnFinalCondition_ );
        break;
    }
    case cpu_time_stopping_condition:
    {
        // Create stopping time termination condition.
        std::shared_ptr< PropagationCPUTimeTerminationSettings > cpuTimeTerminationSettings =
                std::dynamic_pointer_cast< PropagationCPUTimeTerminationSettings >( terminationSettings );
        propagationTerminationCondition = std::make_shared< FixedCPUTimePropagationTerminationCondition >(
                    cpuTimeTerminationSettings->cpuTerminationTime_ );
        break;
    }
    case dependent_variable_stopping_condition:
    {
        std::shared_ptr< PropagationDependentVariableTerminationSettings > dependentVariableTerminationSettings =
                std::dynamic_pointer_cast< PropagationDependentVariableTerminationSettings >( terminationSettings );

        // Get dependent variable function
        std::function< double( ) > dependentVariableFunction;
        if( getDependentVariableSaveSize( dependentVariableTerminationSettings->dependentVariableSettings_ ) == 1 )
        {
            dependentVariableFunction =
                    getDoubleDependentVariableFunction(
                        dependentVariableTerminationSettings->dependentVariableSettings_, bodyMap );
        }
        else
        {
            throw std::runtime_error( "Error, cannot make stopping condition from vector dependent variable" );
        }

        // Create dependent variable termination condition.
        propagationTerminationCondition = std::make_shared< SingleVariableLimitPropagationTerminationCondition >(
                    dependentVariableTerminationSettings->dependentVariableSettings_,
                    dependentVariableFunction, dependentVariableTerminationSettings->limitValue_,
                    dependentVariableTerminationSettings->useAsLowerLimit_,
                    dependentVariableTerminationSettings->terminateExactlyOnFinalCondition_,
                    dependentVariableTerminationSettings->terminationRootFinderSettings_ );
        break;
    }
    case custom_stopping_condition:
    {
        std::shared_ptr< PropagationCustomTerminationSettings > customTerminationSettings =
                std::dynamic_pointer_cast< PropagationCustomTerminationSettings >( terminationSettings );

        // Create dependent variable termination condition.
        propagationTerminationCondition = std::make_shared< CustomTerminationCondition >(
                    customTerminationSettings->checkStopCondition_,
                    customTerminationSettings->terminateExactlyOnFinalCondition_ );
        break;
    }
    case hybrid_stopping_condition:
    {
        std::shared_ptr< PropagationHybridTerminationSettings > hybridTerminationSettings =
                std::dynamic_pointer_cast< PropagationHybridTerminationSettings >( terminationSettings );

        // Recursively create termination condition list.
        std::vector< std::shared_ptr< PropagationTerminationCondition > > propagationTerminationConditionList;
        for( unsigned int i = 0; i < hybridTerminationSettings->terminationSettings_.size( ); i++ )
        {
            propagationTerminationConditionList.push_back(
                        createPropagationTerminationConditions(
                            hybridTerminationSettings->terminationSettings_.at( i ),
                            bodyMap, initialTimeStep ) );
        }
        propagationTerminationCondition = std::make_shared< HybridPropagationTerminationCondition >(
                    propagationTerminationConditionList, hybridTerminationSettings->fulfillSingleCondition_,
                    hybridTerminationSettings->terminateExactlyOnFinalCondition_ );
        break;
    }
    default:
        std::string errorMessage = "Error, stopping condition type " + std::to_string(
                    terminationSettings->terminationType_ ) + " not recognized when making stopping conditions object";
        throw std::runtime_error( errorMessage );
    }
    return propagationTerminationCondition;

} // namespace propagators

} // namespace tudat

}
