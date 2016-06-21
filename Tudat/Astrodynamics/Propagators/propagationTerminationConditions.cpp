
#include "Tudat/Astrodynamics/Propagators/propagationTerminationConditions.h"

namespace tudat
{

namespace propagators
{


bool FixedTimePropagationTerminationCondition::checkStopCondition( const double time  )
{
    bool stopPropagation = 0;

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


bool SingleVariableLimitPropagationTerminationCondition::checkStopCondition( const double time  )
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


bool HybridPropagationTerminationCondition::checkStopCondition( const double time )
{
    if( fulFillSingleCondition_ )
    {
        bool stopPropagation = 0;
        for( unsigned int i = 0; i < propagationTerminationCondition_.size( ); i++ )
        {
            if( propagationTerminationCondition_.at( i )->checkStopCondition( time ) )
            {
                stopPropagation = 1;
                break;
            }
        }
        return stopPropagation;
    }
    else
    {
        bool stopPropagation = 1;
        for( unsigned int i = 0; i < propagationTerminationCondition_.size( ); i++ )
        {
            if( !propagationTerminationCondition_.at( i )->checkStopCondition( time ) )
            {
                stopPropagation = 0;
                break;
            }
        }
        return stopPropagation;
    }
}



boost::shared_ptr< PropagationTerminationCondition > createPropagationTerminationConditions(
        const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double initialTimeStep )
{
    boost::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition;
    switch( terminationSettings->terminationType_ )
    {
    case time_stopping_condition:
    {
        boost::shared_ptr< PropagationTimeTerminationSettings > timeTerminationSettings =
                boost::dynamic_pointer_cast< PropagationTimeTerminationSettings >( terminationSettings );
        propagationTerminationCondition = boost::make_shared< FixedTimePropagationTerminationCondition >(
                    timeTerminationSettings->terminationTime_, ( initialTimeStep > 0 ) ? true : false );
        break;
    }
    case dependent_variable_stopping_condition:
    {
        boost::shared_ptr< PropagationDependentVariableTerminationSettings > dependentVariableTerminationSettings =
                boost::dynamic_pointer_cast< PropagationDependentVariableTerminationSettings >( terminationSettings );

        boost::function< double( ) > dependentVariableFunction;

        if( getDependentVariableSize(
                    dependentVariableTerminationSettings->dependentVariableSettings_->variableType_ ) == 1 )
        {
            dependentVariableFunction =
                    getDoubleDependentVariableFunction(
                        dependentVariableTerminationSettings->dependentVariableSettings_, bodyMap );
        }
        else
        {
            throw std::runtime_error( "Error, cannot make stopping condition from vector dependent variable" );
        }

        propagationTerminationCondition = boost::make_shared< SingleVariableLimitPropagationTerminationCondition >(
                    dependentVariableTerminationSettings->dependentVariableSettings_,
                    dependentVariableFunction, dependentVariableTerminationSettings->limitValue_,
                    dependentVariableTerminationSettings->useAsLowerLimit_ );
        break;
    }
    case hybrid_stopping_condition:
    {
        boost::shared_ptr< PropagationHybridTerminationSettings > hybridTerminationSettings =
                boost::dynamic_pointer_cast< PropagationHybridTerminationSettings >( terminationSettings );

        std::vector< boost::shared_ptr< PropagationTerminationCondition > > propagationTerminationConditionList;
        for( unsigned int i = 0; i < hybridTerminationSettings->terminationSettings_.size( ); i++ )
        {
            propagationTerminationConditionList.push_back(
                        createPropagationTerminationConditions(
                            hybridTerminationSettings->terminationSettings_.at( i ),
                            bodyMap, initialTimeStep ) );
        }
        propagationTerminationCondition = boost::make_shared< HybridPropagationTerminationCondition >(
                    propagationTerminationConditionList, hybridTerminationSettings->fulFillSingleCondition_ );
        break;
    }
    default:

        break;
    }
    return propagationTerminationCondition;
}

}

}
