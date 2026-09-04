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

//! Inspectable condition kinds used by ObservationSelectionCondition.
enum class ObservationSelectionConditionType {
    all,
    observable_type,
    link_definition,
    link_end_type,
    link_end,
    set_id,
    time_bounds,
    time_greater_equal,
    time_greater_than,
    time_less_equal,
    time_less_than,
    active,
    residual_absolute_value_greater_than,
    observation_absolute_value_greater_than,
    dependent_variable_greater_than,
    and_condition,
    or_condition,
    not_condition,
    custom
};

inline std::string getObservationSelectionConditionTypeString( const ObservationSelectionConditionType conditionType )
{
    switch( conditionType )
    {
        case ObservationSelectionConditionType::all:
            return "all";
        case ObservationSelectionConditionType::observable_type:
            return "observable_type";
        case ObservationSelectionConditionType::link_definition:
            return "link_definition";
        case ObservationSelectionConditionType::link_end_type:
            return "link_end_type";
        case ObservationSelectionConditionType::link_end:
            return "link_end";
        case ObservationSelectionConditionType::set_id:
            return "set_id";
        case ObservationSelectionConditionType::time_bounds:
            return "time_bounds";
        case ObservationSelectionConditionType::time_greater_equal:
            return "time_greater_equal";
        case ObservationSelectionConditionType::time_greater_than:
            return "time_greater_than";
        case ObservationSelectionConditionType::time_less_equal:
            return "time_less_equal";
        case ObservationSelectionConditionType::time_less_than:
            return "time_less_than";
        case ObservationSelectionConditionType::active:
            return "active";
        case ObservationSelectionConditionType::residual_absolute_value_greater_than:
            return "residual_absolute_value_greater_than";
        case ObservationSelectionConditionType::observation_absolute_value_greater_than:
            return "observation_absolute_value_greater_than";
        case ObservationSelectionConditionType::dependent_variable_greater_than:
            return "dependent_variable_greater_than";
        case ObservationSelectionConditionType::and_condition:
            return "and";
        case ObservationSelectionConditionType::or_condition:
            return "or";
        case ObservationSelectionConditionType::not_condition:
            return "not";
        case ObservationSelectionConditionType::custom:
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
class ObservationSelectionCondition
{
public:
    using Evaluator = std::function< bool( const ObservationDataset< ObservationScalarType, TimeType >&, const int ) >;

    ObservationSelectionCondition( ): type_( ObservationSelectionConditionType::all ) {}

    explicit ObservationSelectionCondition( const Evaluator& evaluator ):
        type_( ObservationSelectionConditionType::custom ), customEvaluator_( evaluator )
    {}

    bool operator( )( const ObservationDataset< ObservationScalarType, TimeType >& dataset, const unsigned int observationId ) const
    {
        switch( type_ )
        {
            case ObservationSelectionConditionType::all:
                return true;
            case ObservationSelectionConditionType::observable_type:
                return evaluateObservableType( dataset, observationId );
            case ObservationSelectionConditionType::link_definition:
                return evaluateLinkDefinition( dataset, observationId );
            case ObservationSelectionConditionType::link_end_type:
                return evaluateLinkEndType( dataset, observationId );
            case ObservationSelectionConditionType::link_end:
                return evaluateLinkEnd( dataset, observationId );
            case ObservationSelectionConditionType::set_id:
                return dataset.getObservationRow( observationId ).setId_ == setId_;
            case ObservationSelectionConditionType::time_bounds:
                return evaluateTimeBounds( dataset, observationId );
            case ObservationSelectionConditionType::time_greater_equal:
                return dataset.getObservationTime( observationId ) >= timeValue_;
            case ObservationSelectionConditionType::time_greater_than:
                return dataset.getObservationTime( observationId ) > timeValue_;
            case ObservationSelectionConditionType::time_less_equal:
                return dataset.getObservationTime( observationId ) <= timeValue_;
            case ObservationSelectionConditionType::time_less_than:
                return dataset.getObservationTime( observationId ) < timeValue_;
            case ObservationSelectionConditionType::active:
                return dataset.getObservationRow( observationId ).isActive_;
            case ObservationSelectionConditionType::residual_absolute_value_greater_than:
                return evaluateVectorLimit( dataset.getResidualValue( observationId ), vectorLimit_, "residual" );
            case ObservationSelectionConditionType::observation_absolute_value_greater_than:
                return evaluateVectorLimit( dataset.getObservationValue( observationId ), vectorLimit_, "observation" );
            case ObservationSelectionConditionType::dependent_variable_greater_than:
                return evaluateDependentVariable( dataset, observationId );
            case ObservationSelectionConditionType::and_condition:
                return children_.at( 0 )( dataset, observationId ) && children_.at( 1 )( dataset, observationId );
            case ObservationSelectionConditionType::or_condition:
                return children_.at( 0 )( dataset, observationId ) || children_.at( 1 )( dataset, observationId );
            case ObservationSelectionConditionType::not_condition:
                return !children_.at( 0 )( dataset, observationId );
            case ObservationSelectionConditionType::custom:
                return customEvaluator_( dataset, observationId );
            default:
                throw std::runtime_error( "Error when evaluating observation condition, condition type not recognized." );
        }
    }

    ObservationSelectionCondition operator&&( const ObservationSelectionCondition& other ) const
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::and_condition;
        condition.children_ = { *this, other };
        return condition;
    }

    ObservationSelectionCondition operator||( const ObservationSelectionCondition& other ) const
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::or_condition;
        condition.children_ = { *this, other };
        return condition;
    }

    ObservationSelectionCondition operator!( ) const
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::not_condition;
        condition.children_ = { *this };
        return condition;
    }

    ObservationSelectionConditionType getConditionType( ) const
    {
        return type_;
    }

    std::string getConditionTypeString( ) const
    {
        return getObservationSelectionConditionTypeString( type_ );
    }

    const std::vector< ObservationSelectionCondition >& getChildConditions( ) const
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

    unsigned int getSetId( ) const
    {
        return setId_;
    }

    std::pair< TimeType, TimeType > getTimeBounds( ) const
    {
        return std::make_pair( startTime_, endTime_ );
    }

    TimeType getTimeValue( ) const
    {
        return timeValue_;
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

    static ObservationSelectionCondition all( )
    {
        return ObservationSelectionCondition( );
    }

    static ObservationSelectionCondition observableType( const ObservableType observableType )
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::observable_type;
        condition.observableType_ = observableType;
        return condition;
    }

    static ObservationSelectionCondition linkDefinition( const LinkDefinition& linkDefinition )
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::link_definition;
        condition.linkDefinition_ = linkDefinition;
        return condition;
    }

    static ObservationSelectionCondition linkEndType( const LinkEndType linkEndType )
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::link_end_type;
        condition.linkEndType_ = linkEndType;
        return condition;
    }

    static ObservationSelectionCondition linkEnd( const LinkEndType linkEndType, const LinkEndId& linkEndId )
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::link_end;
        condition.linkEndType_ = linkEndType;
        condition.linkEndId_ = linkEndId;
        return condition;
    }

    static ObservationSelectionCondition setId( const unsigned int setId )
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::set_id;
        condition.setId_ = setId;
        return condition;
    }

    static ObservationSelectionCondition timeBounds( const TimeType startTime, const TimeType endTime )
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::time_bounds;
        condition.startTime_ = startTime;
        condition.endTime_ = endTime;
        return condition;
    }

    static ObservationSelectionCondition timeGreaterEqual( const TimeType time )
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::time_greater_equal;
        condition.timeValue_ = time;
        return condition;
    }

    static ObservationSelectionCondition timeGreaterThan( const TimeType time )
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::time_greater_than;
        condition.timeValue_ = time;
        return condition;
    }

    static ObservationSelectionCondition timeLessEqual( const TimeType time )
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::time_less_equal;
        condition.timeValue_ = time;
        return condition;
    }

    static ObservationSelectionCondition timeLessThan( const TimeType time )
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::time_less_than;
        condition.timeValue_ = time;
        return condition;
    }

    static ObservationSelectionCondition active( )
    {
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::active;
        return condition;
    }

    static ObservationSelectionCondition rejected( )
    {
        return !active( );
    }

    static ObservationSelectionCondition residualAbsoluteValueGreaterThan( const Eigen::VectorXd& residualLimit )
    {
        // Selects a row when any scalar residual component exceeds the absolute limit.
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::residual_absolute_value_greater_than;
        condition.vectorLimit_ = residualLimit;
        return condition;
    }

    static ObservationSelectionCondition residualAbsoluteValueGreaterThan( const double residualLimit )
    {
        return residualAbsoluteValueGreaterThan( ( Eigen::VectorXd( 1 ) << residualLimit ).finished( ) );
    }

    static ObservationSelectionCondition observationAbsoluteValueGreaterThan( const Eigen::VectorXd& observationLimit )
    {
        // Selects a row when any scalar observation component exceeds the absolute limit.
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::observation_absolute_value_greater_than;
        condition.vectorLimit_ = observationLimit;
        return condition;
    }

    static ObservationSelectionCondition observationAbsoluteValueGreaterThan( const double observationLimit )
    {
        return observationAbsoluteValueGreaterThan( ( Eigen::VectorXd( 1 ) << observationLimit ).finished( ) );
    }

    static ObservationSelectionCondition dependentVariableGreaterThan(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
            const Eigen::VectorXd& dependentVariableLimit,
            const bool returnFirstCompatibleSettings = false )
    {
        if( dependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating dependent-variable observation condition, settings are null." );
        }
        // Selects a row when any compatible dependent-variable component is greater than the signed limit.
        ObservationSelectionCondition condition;
        condition.type_ = ObservationSelectionConditionType::dependent_variable_greater_than;
        condition.dependentVariableSettings_ = dependentVariableSettings;
        condition.vectorLimit_ = dependentVariableLimit;
        condition.returnFirstCompatibleSettings_ = returnFirstCompatibleSettings;
        return condition;
    }

    static ObservationSelectionCondition dependentVariableGreaterThan(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
            const double dependentVariableLimit,
            const bool returnFirstCompatibleSettings = false )
    {
        return dependentVariableGreaterThan(
                dependentVariableSettings, ( Eigen::VectorXd( 1 ) << dependentVariableLimit ).finished( ), returnFirstCompatibleSettings );
    }

private:
    bool evaluateObservableType( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                 const unsigned int observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
        return dataset.getObservationSetMetadata( row.setId_ ).observableType_ == observableType_;
    }

    bool evaluateLinkDefinition( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                 const unsigned int observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = dataset.getObservationSetMetadata( row.setId_ );
        return dataset.getLinkDefinition( metadata.linkDefinitionId_ ) == linkDefinition_;
    }

    bool evaluateLinkEndType( const ObservationDataset< ObservationScalarType, TimeType >& dataset, const unsigned int observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = dataset.getObservationSetMetadata( row.setId_ );
        return dataset.getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_.count( linkEndType_ ) > 0;
    }

    bool evaluateLinkEnd( const ObservationDataset< ObservationScalarType, TimeType >& dataset, const unsigned int observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = dataset.getObservationSetMetadata( row.setId_ );
        const LinkEnds& linkEnds = dataset.getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        return linkEnds.count( linkEndType_ ) > 0 && linkEnds.at( linkEndType_ ) == linkEndId_;
    }

    bool evaluateTimeBounds( const ObservationDataset< ObservationScalarType, TimeType >& dataset, const unsigned int observationId ) const
    {
        const TimeType observationTime = dataset.getObservationTime( observationId );
        return observationTime >= startTime_ && observationTime <= endTime_;
    }

    bool evaluateDependentVariable( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                    const unsigned int observationId ) const
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
        if( indexAndSize.first + indexAndSize.second > dependentVariables.size( ) )
        {
            throw std::runtime_error( "Error when evaluating dependent-variable condition, dependent-variable size is inconsistent." );
        }

        return evaluateGreaterThanLimit(
                dependentVariables.segment( indexAndSize.first, indexAndSize.second ), vectorLimit_, "dependent-variable" );
    }

    static double limitForComponent( const Eigen::VectorXd& limit,
                                     const int componentIndex,
                                     const int expectedSize,
                                     const std::string& valueName )
    {
        if( limit.size( ) == 1 )
        {
            return limit( 0 );
        }
        if( limit.size( ) == expectedSize )
        {
            return limit( componentIndex );
        }
        throw std::runtime_error( "Error when evaluating " + valueName + " condition, limit size is " + std::to_string( limit.size( ) ) +
                                  " but expected either 1 or " + std::to_string( expectedSize ) + "." );
    }

    static bool evaluateGreaterThanLimit( const Eigen::VectorXd& values, const Eigen::VectorXd& limit, const std::string& valueName )
    {
        for( int i = 0; i < values.size( ); ++i )
        {
            if( values( i ) > limitForComponent( limit, i, values.size( ), valueName ) )
            {
                return true;
            }
        }
        return false;
    }

    template< typename VectorType >
    static bool evaluateVectorLimit( const VectorType& values, const Eigen::VectorXd& limit, const std::string& valueName )
    {
        for( int i = 0; i < values.size( ); ++i )
        {
            if( std::fabs( values( i ) ) > limitForComponent( limit, i, values.size( ), valueName ) )
            {
                return true;
            }
        }
        return false;
    }

    //! Node type used to evaluate and inspect this condition.
    ObservationSelectionConditionType type_ = ObservationSelectionConditionType::all;

    //! Child nodes used by logical AND/OR/NOT conditions.
    std::vector< ObservationSelectionCondition > children_;

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

    //! Observation set id parameter for set-id leaves.
    unsigned int setId_ = 0;

    //! Inclusive start time for time-bound leaves.
    TimeType startTime_ = TimeType( );

    //! Inclusive end time for time-bound leaves.
    TimeType endTime_ = TimeType( );

    //! Single time parameter for one-sided time comparison leaves.
    TimeType timeValue_ = TimeType( );

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
