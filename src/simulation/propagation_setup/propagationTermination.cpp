/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/propagation_setup/propagationTermination.h"

namespace tudat
{

namespace propagators
{

//! Function to check whether the propagation is to be be stopped
bool FixedTimePropagationTerminationCondition::checkStopCondition( const double time, const double cpuTime, const Eigen::MatrixXd& state )
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
bool FixedCPUTimePropagationTerminationCondition::checkStopCondition( const double time, const double cpuTime, const Eigen::MatrixXd& state )
{
    return cpuTime >= cpuStopTime_;
}

//! Function to check whether the propagation is to be be stopped
bool SingleVariableLimitPropagationTerminationCondition::checkStopCondition( const double time, const double cpuTime, const Eigen::MatrixXd& state )
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
bool HybridPropagationTerminationCondition::checkStopCondition( const double time, const double cpuTime, const Eigen::MatrixXd& state )
{
    // Check if single condition is fulfilled.
    bool stopPropagation = -1;
    unsigned int stopIndex = 0;
    if( fulfillSingleCondition_ )
    {
        stopPropagation = false;
        for( unsigned int i = 0; i < propagationTerminationCondition_.size( ); i++ )
        {
            if( propagationTerminationCondition_.at( i )->checkStopCondition( time, cpuTime, state ) )
            {
                stopIndex = i;
                stopPropagation = true;
                isConditionMetWhenStopping_[ i ] = true;
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
            if( !propagationTerminationCondition_.at( i )->checkStopCondition( time, cpuTime, state ) )
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
            isConditionMetWhenStopping_[ i ] = propagationTerminationCondition_.at( i )->checkStopCondition( time, cpuTime, state );
        }
    }

    return stopPropagation;
}

bool HybridPropagationTerminationCondition::iterateToExactTermination( )
{
    bool iterateToExactCondition = 0;
    if( checkTerminationToExactCondition_ )
    {
        for( unsigned int i = 0; i < propagationTerminationCondition_.size( ); i++ )
        {
            if( isConditionMetWhenStopping_.at( i ) && propagationTerminationCondition_.at( i )->getcheckTerminationToExactCondition( ) )
            {
                iterateToExactCondition = true;
            }
        }
    }
    return iterateToExactCondition;
}

//! Function to check whether the propagation is to be be stopped
bool NonSequentialPropagationTerminationCondition::checkStopCondition( const double time, const double cpuTime, const Eigen::MatrixXd& state )
{
    // Check if single condition is fulfilled.
    bool stopPropagation = false;
    if ( forwardPropagationTerminationCondition_->checkStopCondition( time, cpuTime, state ) ||
    backwardPropagationTerminationCondition_->checkStopCondition( time, cpuTime, state ) )
    {
        stopPropagation = true;
    }

    return stopPropagation;
}

} // namespace propagators

} // namespace tudat

