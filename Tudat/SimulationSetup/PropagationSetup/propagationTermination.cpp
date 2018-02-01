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
    bool stopPropagation = 0;

    // Check whether stop time has been reached
    if( propagationDirectionIsPositive_ && ( time >= stopTime_ ) )
    {
        stopPropagation = 1;
    }
    else if( !propagationDirectionIsPositive_ && ( time <= stopTime_ ) )
    {
        stopPropagation = 1;
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
    bool stopPropagation = 0;
    double currentVariable = variableRetrievalFuntion_( );
    if( useAsLowerBound_ && ( currentVariable < limitingValue_ ) )
    {
        stopPropagation = 1;
    }
    else if( !useAsLowerBound_ && ( currentVariable > limitingValue_ ) )
    {
        stopPropagation = 1;
    }
    return stopPropagation;
}

//! Function to check whether the propagation is to be be stopped
bool HybridPropagationTerminationCondition::checkStopCondition( const double time, const double cpuTime )
{
    // Check if single condition is fulfilled.
    bool stopPropagation = -1;
    if( fulFillSingleCondition_ )
    {
        stopPropagation = 0;
        for( unsigned int i = 0; i < propagationTerminationCondition_.size( ); i++ )
        {
            if( propagationTerminationCondition_.at( i )->checkStopCondition( time, cpuTime ) )
            {
                stopPropagation = 1;
                break;
            }
        }
    }
    // Check all conditions are fulfilled.
    else
    {
        stopPropagation = 1;
        for( unsigned int i = 0; i < propagationTerminationCondition_.size( ); i++ )
        {
            if( !propagationTerminationCondition_.at( i )->checkStopCondition( time, cpuTime ) )
            {
                stopPropagation = 0;
                break;
            }
        }
    }

    // Save if conditions were met
    if( stopPropagation )
    {
        for( unsigned int i = 0; i < propagationTerminationCondition_.size( ); i++ )
        {
            if( propagationTerminationCondition_.at( i )->checkStopCondition( time, cpuTime ) )
            {
                isConditionMetWhenStopping_[ i ] = false;
            }
            else
            {
                isConditionMetWhenStopping_[ i ] = true;
            }
        }
    }

    return stopPropagation;

}


//! Function to create propagation termination conditions from associated settings
boost::shared_ptr< PropagationTerminationCondition > createPropagationTerminationConditions(
        const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double initialTimeStep )
{
    boost::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition;

    // Check termination type.
    switch( terminationSettings->terminationType_ )
    {
    case time_stopping_condition:
    {
        // Create stopping time termination condition.
        boost::shared_ptr< PropagationTimeTerminationSettings > timeTerminationSettings =
                boost::dynamic_pointer_cast< PropagationTimeTerminationSettings >( terminationSettings );
        propagationTerminationCondition = boost::make_shared< FixedTimePropagationTerminationCondition >(
                    timeTerminationSettings->terminationTime_, ( initialTimeStep > 0 ),
                    timeTerminationSettings->terminateExactlyOnFinalCondition_ );
        break;
    }
    case cpu_time_stopping_condition:
    {
        // Create stopping time termination condition.
        boost::shared_ptr< PropagationCPUTimeTerminationSettings > cpuTimeTerminationSettings =
                boost::dynamic_pointer_cast< PropagationCPUTimeTerminationSettings >( terminationSettings );
        propagationTerminationCondition = boost::make_shared< FixedCPUTimePropagationTerminationCondition >(
                    cpuTimeTerminationSettings->cpuTerminationTime_ );
        break;
    }
    case dependent_variable_stopping_condition:
    {
        boost::shared_ptr< PropagationDependentVariableTerminationSettings > dependentVariableTerminationSettings =
                boost::dynamic_pointer_cast< PropagationDependentVariableTerminationSettings >( terminationSettings );

        // Get dependent variable function
        boost::function< double( ) > dependentVariableFunction;
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
        propagationTerminationCondition = boost::make_shared< SingleVariableLimitPropagationTerminationCondition >(
                    dependentVariableTerminationSettings->dependentVariableSettings_,
                    dependentVariableFunction, dependentVariableTerminationSettings->limitValue_,
                    dependentVariableTerminationSettings->useAsLowerLimit_,
                    dependentVariableTerminationSettings->terminateExactlyOnFinalCondition_,
                    dependentVariableTerminationSettings->terminationRootFinderSettings_ );
        break;
    }
    case hybrid_stopping_condition:
    {
        boost::shared_ptr< PropagationHybridTerminationSettings > hybridTerminationSettings =
                boost::dynamic_pointer_cast< PropagationHybridTerminationSettings >( terminationSettings );

        // Recursively create termination condition list.
        std::vector< boost::shared_ptr< PropagationTerminationCondition > > propagationTerminationConditionList;
        for( unsigned int i = 0; i < hybridTerminationSettings->terminationSettings_.size( ); i++ )
        {
            propagationTerminationConditionList.push_back(
                        createPropagationTerminationConditions(
                            hybridTerminationSettings->terminationSettings_.at( i ),
                            bodyMap, initialTimeStep ) );
        }
        propagationTerminationCondition = boost::make_shared< HybridPropagationTerminationCondition >(
                    propagationTerminationConditionList, hybridTerminationSettings->fulFillSingleCondition_,
                    hybridTerminationSettings->terminateExactlyOnFinalCondition_ );
        break;
    }
    default:
        std::string errorMessage = "Error, stopping condition type " + std::to_string(
                    terminationSettings->terminationType_ ) + "not recognized when making stopping conditions object";
        throw std::runtime_error( errorMessage );
        break;
    }
    return propagationTerminationCondition;

} // namespace propagators

} // namespace tudat

}
