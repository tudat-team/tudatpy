/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_CONDITION_H
#define TUDAT_OBSERVATION_CONDITION_H

#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/simulation/estimation_setup/observationDatasetRows.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"

namespace tudat
{

namespace observation_models
{

//! Inspectable condition kinds used by ObservationCondition.
enum class ObservationConditionType {
    all,
    observable_type,
    link_definition,
    link_end_type,
    link_end,
    time_bounds,
    active,
    residual_absolute_value_greater_than,
    observation_absolute_value_greater_than,
    dependent_variable_greater_than,
    and_condition,
    or_condition,
    not_condition,
    custom
};

inline std::string getObservationConditionTypeString( const ObservationConditionType conditionType )
{
    switch( conditionType )
    {
        case ObservationConditionType::all:
            return "all";
        case ObservationConditionType::observable_type:
            return "observable_type";
        case ObservationConditionType::link_definition:
            return "link_definition";
        case ObservationConditionType::link_end_type:
            return "link_end_type";
        case ObservationConditionType::link_end:
            return "link_end";
        case ObservationConditionType::time_bounds:
            return "time_bounds";
        case ObservationConditionType::active:
            return "active";
        case ObservationConditionType::residual_absolute_value_greater_than:
            return "residual_absolute_value_greater_than";
        case ObservationConditionType::observation_absolute_value_greater_than:
            return "observation_absolute_value_greater_than";
        case ObservationConditionType::dependent_variable_greater_than:
            return "dependent_variable_greater_than";
        case ObservationConditionType::and_condition:
            return "and";
        case ObservationConditionType::or_condition:
            return "or";
        case ObservationConditionType::not_condition:
            return "not";
        case ObservationConditionType::custom:
            return "custom";
        default:
            throw std::runtime_error( "Error when getting observation condition type string, type not recognized." );
    }
}

//! Row-level condition used by the new observation-selection API.
/*!
 * Conditions are evaluated for individual observation events, not complete
 * observation sets. Public factory functions build an inspectable query tree,
 * so a selection can be logged, exposed to Python and combined predictably.
 * Legacy/custom callables remain supported as opaque custom leaves. Legacy
 * parsers and filters can be adapted to conditions, but are not part of this
 * type.
 */
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
class ObservationCondition
{
public:
    using Evaluator = std::function< bool( const ObservationDataset< ObservationScalarType, TimeType >&, const ObservationId ) >;

    ObservationCondition( ): type_( ObservationConditionType::all ) {}

    explicit ObservationCondition( const Evaluator& evaluator ): type_( ObservationConditionType::custom ), customEvaluator_( evaluator ) {}

    bool operator( )( const ObservationDataset< ObservationScalarType, TimeType >& dataset, const ObservationId observationId ) const
    {
        switch( type_ )
        {
            case ObservationConditionType::all:
                return true;
            case ObservationConditionType::observable_type:
                return evaluateObservableType( dataset, observationId );
            case ObservationConditionType::link_definition:
                return evaluateLinkDefinition( dataset, observationId );
            case ObservationConditionType::link_end_type:
                return evaluateLinkEndType( dataset, observationId );
            case ObservationConditionType::link_end:
                return evaluateLinkEnd( dataset, observationId );
            case ObservationConditionType::time_bounds:
                return evaluateTimeBounds( dataset, observationId );
            case ObservationConditionType::active:
                return dataset.getObservationRow( observationId ).isActive_;
            case ObservationConditionType::residual_absolute_value_greater_than:
                return evaluateVectorLimit( dataset.getResidualValue( observationId ), vectorLimit_, "residual" );
            case ObservationConditionType::observation_absolute_value_greater_than:
                return evaluateVectorLimit( dataset.getObservationValue( observationId ), vectorLimit_, "observation" );
            case ObservationConditionType::dependent_variable_greater_than:
                return evaluateDependentVariable( dataset, observationId );
            case ObservationConditionType::and_condition:
                return children_.at( 0 )( dataset, observationId ) && children_.at( 1 )( dataset, observationId );
            case ObservationConditionType::or_condition:
                return children_.at( 0 )( dataset, observationId ) || children_.at( 1 )( dataset, observationId );
            case ObservationConditionType::not_condition:
                return !children_.at( 0 )( dataset, observationId );
            case ObservationConditionType::custom:
                return customEvaluator_( dataset, observationId );
            default:
                throw std::runtime_error( "Error when evaluating observation condition, condition type not recognized." );
        }
    }

    ObservationCondition operator&&( const ObservationCondition& other ) const
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::and_condition;
        condition.children_ = { *this, other };
        return condition;
    }

    ObservationCondition operator||( const ObservationCondition& other ) const
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::or_condition;
        condition.children_ = { *this, other };
        return condition;
    }

    ObservationCondition operator!( ) const
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::not_condition;
        condition.children_ = { *this };
        return condition;
    }

    ObservationConditionType getConditionType( ) const
    {
        return type_;
    }

    std::string getConditionTypeString( ) const
    {
        return getObservationConditionTypeString( type_ );
    }

    const std::vector< ObservationCondition >& getChildConditions( ) const
    {
        return children_;
    }

    ObservableType getObservableType( ) const
    {
        return observableType_;
    }

    LinkDefinition getLinkDefinition( ) const
    {
        return linkDefinition_;
    }

    LinkEndType getLinkEndType( ) const
    {
        return linkEndType_;
    }

    LinkEndId getLinkEndId( ) const
    {
        return linkEndId_;
    }

    std::pair< TimeType, TimeType > getTimeBounds( ) const
    {
        return std::make_pair( startTime_, endTime_ );
    }

    const Eigen::VectorXd& getVectorLimit( ) const
    {
        return vectorLimit_;
    }

    std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > getDependentVariableSettings( ) const
    {
        return dependentVariableSettings_;
    }

    bool getReturnFirstCompatibleDependentVariableSettings( ) const
    {
        return returnFirstCompatibleSettings_;
    }

    static ObservationCondition all( )
    {
        return ObservationCondition( );
    }

    static ObservationCondition observableType( const ObservableType observableType )
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::observable_type;
        condition.observableType_ = observableType;
        return condition;
    }

    static ObservationCondition linkDefinition( const LinkDefinition& linkDefinition )
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::link_definition;
        condition.linkDefinition_ = linkDefinition;
        return condition;
    }

    static ObservationCondition linkEndType( const LinkEndType linkEndType )
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::link_end_type;
        condition.linkEndType_ = linkEndType;
        return condition;
    }

    static ObservationCondition linkEnd( const LinkEndType linkEndType, const LinkEndId& linkEndId )
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::link_end;
        condition.linkEndType_ = linkEndType;
        condition.linkEndId_ = linkEndId;
        return condition;
    }

    static ObservationCondition timeBounds( const TimeType startTime, const TimeType endTime )
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::time_bounds;
        condition.startTime_ = startTime;
        condition.endTime_ = endTime;
        return condition;
    }

    static ObservationCondition active( )
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::active;
        return condition;
    }

    static ObservationCondition rejected( )
    {
        return !active( );
    }

    static ObservationCondition residualAbsoluteValueGreaterThan( const Eigen::VectorXd& residualLimit )
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::residual_absolute_value_greater_than;
        condition.vectorLimit_ = residualLimit;
        return condition;
    }

    static ObservationCondition observationAbsoluteValueGreaterThan( const Eigen::VectorXd& observationLimit )
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::observation_absolute_value_greater_than;
        condition.vectorLimit_ = observationLimit;
        return condition;
    }

    static ObservationCondition dependentVariableGreaterThan(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
            const Eigen::VectorXd& dependentVariableLimit,
            const bool returnFirstCompatibleSettings = false )
    {
        ObservationCondition condition;
        condition.type_ = ObservationConditionType::dependent_variable_greater_than;
        condition.dependentVariableSettings_ = dependentVariableSettings;
        condition.vectorLimit_ = dependentVariableLimit;
        condition.returnFirstCompatibleSettings_ = returnFirstCompatibleSettings;
        return condition;
    }

private:
    bool evaluateObservableType( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                 const ObservationId observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
        return dataset.getObservationSetMetadata( row.setId_ ).observableType_ == observableType_;
    }

    bool evaluateLinkDefinition( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                 const ObservationId observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = dataset.getObservationSetMetadata( row.setId_ );
        return dataset.getLinkDefinition( metadata.linkDefinitionId_ ) == linkDefinition_;
    }

    bool evaluateLinkEndType( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                              const ObservationId observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = dataset.getObservationSetMetadata( row.setId_ );
        return dataset.getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_.count( linkEndType_ ) > 0;
    }

    bool evaluateLinkEnd( const ObservationDataset< ObservationScalarType, TimeType >& dataset, const ObservationId observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = dataset.getObservationSetMetadata( row.setId_ );
        const LinkEnds& linkEnds = dataset.getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        return linkEnds.count( linkEndType_ ) > 0 && linkEnds.at( linkEndType_ ) == linkEndId_;
    }

    bool evaluateTimeBounds( const ObservationDataset< ObservationScalarType, TimeType >& dataset, const ObservationId observationId ) const
    {
        const TimeType observationTime = dataset.getObservationTime( observationId );
        return observationTime >= startTime_ && observationTime <= endTime_;
    }

    bool evaluateDependentVariable( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                    const ObservationId observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = dataset.getObservationSetMetadata( row.setId_ );
        const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping =
                dataset.getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
        if( dependentVariableBookkeeping == nullptr )
        {
            return false;
        }

        std::vector< std::pair< int, int > > compatibleIndicesAndSizes;
        for( const auto& settingsIterator : dependentVariableBookkeeping->getSettingsIndicesAndSizes( ) )
        {
            if( dependentVariableSettings_->areSettingsCompatible( settingsIterator.second ) )
            {
                compatibleIndicesAndSizes.push_back( settingsIterator.first );
            }
        }

        if( compatibleIndicesAndSizes.empty( ) )
        {
            return false;
        }
        if( compatibleIndicesAndSizes.size( ) > 1 && !returnFirstCompatibleSettings_ )
        {
            throw std::runtime_error(
                    "Error when evaluating dependent-variable condition, multiple compatible dependent variables found." );
        }

        const std::pair< int, int >& indexAndSize = compatibleIndicesAndSizes.at( 0 );
        const Eigen::VectorXd dependentVariables = dataset.getDependentVariables( observationId );
        if( indexAndSize.first + indexAndSize.second > dependentVariables.size( ) || vectorLimit_.size( ) != indexAndSize.second )
        {
            throw std::runtime_error( "Error when evaluating dependent-variable condition, dependent-variable size is inconsistent." );
        }

        return evaluateGreaterThanLimit( dependentVariables.segment( indexAndSize.first, indexAndSize.second ), vectorLimit_ );
    }

    static bool evaluateGreaterThanLimit( const Eigen::VectorXd& values, const Eigen::VectorXd& limit )
    {
        for( int i = 0; i < values.size( ); ++i )
        {
            if( values( i ) > limit( i ) )
            {
                return true;
            }
        }
        return false;
    }

    template< typename VectorType >
    static bool evaluateVectorLimit( const VectorType& values, const Eigen::VectorXd& limit, const std::string& valueName )
    {
        if( limit.size( ) != values.size( ) )
        {
            throw std::runtime_error( "Error when evaluating " + valueName +
                                      " condition, limit size is incompatible with observation size." );
        }
        for( int i = 0; i < values.size( ); ++i )
        {
            if( std::fabs( values( i ) ) > limit( i ) )
            {
                return true;
            }
        }
        return false;
    }

    //! Node type used to evaluate and inspect this condition.
    ObservationConditionType type_ = ObservationConditionType::all;

    //! Child nodes used by logical AND/OR/NOT conditions.
    std::vector< ObservationCondition > children_;

    //! Opaque evaluator used only for custom callable leaves.
    Evaluator customEvaluator_;

    //! Observable type parameter for observable-type leaves.
    ObservableType observableType_ = undefined_observation_model;

    //! Full link definition parameter for link-definition leaves.
    LinkDefinition linkDefinition_;

    //! Link-end type parameter for link-end-type and link-end leaves.
    LinkEndType linkEndType_ = unidentified_link_end;

    //! Link-end id parameter for link-end leaves.
    LinkEndId linkEndId_;

    //! Inclusive start time for time-bound leaves.
    TimeType startTime_ = TimeType( );

    //! Inclusive end time for time-bound leaves.
    TimeType endTime_ = TimeType( );

    //! Component limit used by residual, observation and dependent-variable leaves.
    Eigen::VectorXd vectorLimit_;

    //! Dependent-variable setting matched by dependent-variable leaves.
    std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > dependentVariableSettings_;

    //! Whether dependent-variable leaves accept the first match when multiple settings are compatible.
    bool returnFirstCompatibleSettings_ = false;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_CONDITION_H
