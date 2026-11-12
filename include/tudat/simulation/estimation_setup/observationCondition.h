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

//! Row-level condition used to select observations in an ObservationDataset.
/*!
 * Conditions are evaluated for individual observation events, not complete observation sets. They are created with the
 * static factory functions or from a custom row predicate, and combined with &&, || and !.
 */
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
class ObservationSelectionCondition
{
    //! Condition kinds, used to dispatch the evaluation.
    enum class Type {
        all,
        observable_type,
        link_definition,
        link_end_type,
        link_end,
        set_id,
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

public:
    using Evaluator = std::function< bool( const ObservationDataset< ObservationScalarType, TimeType >&, const unsigned int ) >;

    //! Create a condition that selects every observation row.
    ObservationSelectionCondition( ) = default;

    //! Create a condition evaluated by a custom function of the dataset and observation id.
    explicit ObservationSelectionCondition( const Evaluator& evaluator ): type_( Type::custom ), customEvaluator_( evaluator ) {}

    //! Evaluate this condition for one observation identity in a dataset.
    bool operator( )( const ObservationDataset< ObservationScalarType, TimeType >& dataset, const unsigned int observationId ) const
    {
        switch( type_ )
        {
            case Type::all:
                return true;
            case Type::observable_type:
                return getSetMetadata( dataset, observationId ).observableType_ == observableType_;
            case Type::link_definition:
                return dataset.getLinkDefinition( getSetMetadata( dataset, observationId ).linkDefinitionId_ ) == linkDefinition_;
            case Type::link_end_type:
                return dataset.getLinkDefinition( getSetMetadata( dataset, observationId ).linkDefinitionId_ )
                               .linkEnds_.count( linkEndType_ ) > 0;
            case Type::link_end: {
                const LinkEnds& linkEnds =
                        dataset.getLinkDefinition( getSetMetadata( dataset, observationId ).linkDefinitionId_ ).linkEnds_;
                return linkEnds.count( linkEndType_ ) > 0 && linkEnds.at( linkEndType_ ) == linkEndId_;
            }
            case Type::set_id:
                return dataset.getObservationRow( observationId ).setId_ == setId_;
            case Type::time_greater_equal:
                return dataset.getObservationTime( observationId ) >= timeValue_;
            case Type::time_greater_than:
                return dataset.getObservationTime( observationId ) > timeValue_;
            case Type::time_less_equal:
                return dataset.getObservationTime( observationId ) <= timeValue_;
            case Type::time_less_than:
                return dataset.getObservationTime( observationId ) < timeValue_;
            case Type::active:
                return dataset.getObservationRow( observationId ).isActive_;
            case Type::residual_absolute_value_greater_than:
                return isAnyComponentGreaterThanLimit( dataset.getResidualValue( observationId ), true, "residual" );
            case Type::observation_absolute_value_greater_than:
                return isAnyComponentGreaterThanLimit( dataset.getObservationValue( observationId ), true, "observation" );
            case Type::dependent_variable_greater_than:
                return evaluateDependentVariable( dataset, observationId );
            case Type::and_condition:
                return children_.at( 0 )( dataset, observationId ) && children_.at( 1 )( dataset, observationId );
            case Type::or_condition:
                return children_.at( 0 )( dataset, observationId ) || children_.at( 1 )( dataset, observationId );
            case Type::not_condition:
                return !children_.at( 0 )( dataset, observationId );
            case Type::custom:
                return customEvaluator_( dataset, observationId );
            default:
                throw std::runtime_error( "Error when evaluating observation condition, condition type not recognized." );
        }
    }

    //! Combine two conditions with logical AND.
    ObservationSelectionCondition operator&&( const ObservationSelectionCondition& other ) const
    {
        return ObservationSelectionCondition( Type::and_condition, { *this, other } );
    }

    //! Combine two conditions with logical OR.
    ObservationSelectionCondition operator||( const ObservationSelectionCondition& other ) const
    {
        return ObservationSelectionCondition( Type::or_condition, { *this, other } );
    }

    //! Negate this condition.
    ObservationSelectionCondition operator!( ) const
    {
        return ObservationSelectionCondition( Type::not_condition, { *this } );
    }

    //! Create a condition that selects every observation row.
    static ObservationSelectionCondition all( )
    {
        return ObservationSelectionCondition( );
    }

    //! Create a condition matching one observable type.
    static ObservationSelectionCondition observableType( const ObservableType observableType )
    {
        ObservationSelectionCondition condition( Type::observable_type );
        condition.observableType_ = observableType;
        return condition;
    }

    //! Create a condition matching one complete link definition.
    static ObservationSelectionCondition linkDefinition( const LinkDefinition& linkDefinition )
    {
        ObservationSelectionCondition condition( Type::link_definition );
        condition.linkDefinition_ = linkDefinition;
        return condition;
    }

    //! Create a condition matching observations that have a link end with the given role (for any link-end identity).
    static ObservationSelectionCondition linkEndType( const LinkEndType linkEndType )
    {
        ObservationSelectionCondition condition( Type::link_end_type );
        condition.linkEndType_ = linkEndType;
        return condition;
    }

    //! Create a condition matching one link-end role and identity.
    static ObservationSelectionCondition linkEnd( const LinkEndType linkEndType, const LinkEndId& linkEndId )
    {
        ObservationSelectionCondition condition( Type::link_end );
        condition.linkEndType_ = linkEndType;
        condition.linkEndId_ = linkEndId;
        return condition;
    }

    //! Create a condition matching one logical observation set.
    static ObservationSelectionCondition setId( const unsigned int setId )
    {
        ObservationSelectionCondition condition( Type::set_id );
        condition.setId_ = setId;
        return condition;
    }

    //! Create a condition matching epochs within inclusive bounds.
    static ObservationSelectionCondition timeBounds( const TimeType startTime, const TimeType endTime )
    {
        return timeGreaterEqual( startTime ) && timeLessEqual( endTime );
    }

    //! Create a condition matching epochs greater than or equal to a limit.
    static ObservationSelectionCondition timeGreaterEqual( const TimeType time )
    {
        ObservationSelectionCondition condition( Type::time_greater_equal );
        condition.timeValue_ = time;
        return condition;
    }

    //! Create a condition matching epochs strictly greater than a limit.
    static ObservationSelectionCondition timeGreaterThan( const TimeType time )
    {
        ObservationSelectionCondition condition( Type::time_greater_than );
        condition.timeValue_ = time;
        return condition;
    }

    //! Create a condition matching epochs less than or equal to a limit.
    static ObservationSelectionCondition timeLessEqual( const TimeType time )
    {
        ObservationSelectionCondition condition( Type::time_less_equal );
        condition.timeValue_ = time;
        return condition;
    }

    //! Create a condition matching epochs strictly less than a limit.
    static ObservationSelectionCondition timeLessThan( const TimeType time )
    {
        ObservationSelectionCondition condition( Type::time_less_than );
        condition.timeValue_ = time;
        return condition;
    }

    //! Create a condition matching active observation rows.
    static ObservationSelectionCondition active( )
    {
        return ObservationSelectionCondition( Type::active );
    }

    //! Create a condition matching rejected observation rows.
    static ObservationSelectionCondition rejected( )
    {
        return !active( );
    }

    //! Create a condition matching rows where any absolute residual component exceeds its limit (a length-one limit applies to
    //! all components).
    static ObservationSelectionCondition residualAbsoluteValueGreaterThan( const Eigen::VectorXd& residualLimit )
    {
        ObservationSelectionCondition condition( Type::residual_absolute_value_greater_than );
        condition.vectorLimit_ = residualLimit;
        return condition;
    }

    //! Create a condition matching rows where any absolute observation component exceeds its limit (a length-one limit applies
    //! to all components).
    static ObservationSelectionCondition observationAbsoluteValueGreaterThan( const Eigen::VectorXd& observationLimit )
    {
        ObservationSelectionCondition condition( Type::observation_absolute_value_greater_than );
        condition.vectorLimit_ = observationLimit;
        return condition;
    }

    //! Create a condition matching rows where any compatible dependent-variable component is greater than its signed limit.
    static ObservationSelectionCondition dependentVariableGreaterThan(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
            const Eigen::VectorXd& dependentVariableLimit,
            const bool returnFirstCompatibleSettings = false )
    {
        if( dependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating dependent-variable observation condition, settings are null." );
        }
        ObservationSelectionCondition condition( Type::dependent_variable_greater_than );
        condition.dependentVariableSettings_ = dependentVariableSettings;
        condition.vectorLimit_ = dependentVariableLimit;
        condition.returnFirstCompatibleSettings_ = returnFirstCompatibleSettings;
        return condition;
    }

private:
    //! Create a condition of the given kind, with child conditions for logical operators.
    explicit ObservationSelectionCondition( const Type type, const std::vector< ObservationSelectionCondition >& children = {} ):
        type_( type ), children_( children )
    {}

    //! Return the metadata of the observation set containing one observation.
    static const ObservationSetMetadata< ObservationScalarType, TimeType >& getSetMetadata(
            const ObservationDataset< ObservationScalarType, TimeType >& dataset,
            const unsigned int observationId )
    {
        return dataset.getObservationSetMetadata( dataset.getObservationRow( observationId ).setId_ );
    }

    //! Evaluate a dependent-variable threshold condition for one row.
    bool evaluateDependentVariable( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                    const unsigned int observationId ) const
    {
        const Eigen::VectorXd dependentVariable =
                dataset.getSingleDependentVariable( observationId, dependentVariableSettings_, returnFirstCompatibleSettings_ );
        return dependentVariable.size( ) > 0 && isAnyComponentGreaterThanLimit( dependentVariable, false, "dependent-variable" );
    }

    //! Return whether any (absolute) value is greater than vectorLimit_, which has length one or the size of values.
    template< typename VectorType >
    bool isAnyComponentGreaterThanLimit( const VectorType& values, const bool useAbsoluteValues, const std::string& valueName ) const
    {
        for( int i = 0; i < values.size( ); ++i )
        {
            if( vectorLimit_.size( ) != 1 && vectorLimit_.size( ) != values.size( ) )
            {
                throw std::runtime_error( "Error when evaluating " + valueName + " condition, limit size is " +
                                          std::to_string( vectorLimit_.size( ) ) + " but expected either 1 or " +
                                          std::to_string( values.size( ) ) + "." );
            }
            if( ( useAbsoluteValues ? std::fabs( values( i ) ) : values( i ) ) > vectorLimit_( vectorLimit_.size( ) == 1 ? 0 : i ) )
            {
                return true;
            }
        }
        return false;
    }

    //! Condition kind used to dispatch the evaluation.
    Type type_ = Type::all;

    //! Child conditions of logical AND/OR/NOT conditions.
    std::vector< ObservationSelectionCondition > children_;

    //! Row predicate of custom conditions.
    Evaluator customEvaluator_;

    //! Parameters of the condition kinds above; only those belonging to type_ are used.
    ObservableType observableType_ = undefined_observation_model;
    LinkDefinition linkDefinition_;
    LinkEndType linkEndType_ = unidentified_link_end;
    LinkEndId linkEndId_;
    unsigned int setId_ = 0;
    TimeType timeValue_ = TimeType( );
    Eigen::VectorXd vectorLimit_;
    std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > dependentVariableSettings_;
    bool returnFirstCompatibleSettings_ = false;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_CONDITION_H
