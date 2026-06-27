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
#include <vector>

#include <Eigen/Core>

#include "tudat/simulation/estimation_setup/observationDatasetRows.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"

namespace tudat
{

namespace observation_models
{

//! Row-level condition used by the new observation-selection API.
/*!
 * Conditions are evaluated for individual observation events, not complete
 * observation sets. They are the core selection mechanism for viewers,
 * active/rejected status changes and reduced datasets. Legacy parsers and
 * filters can be adapted to conditions, but are not part of this type.
 */
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
class ObservationCondition
{
public:
    using Evaluator = std::function< bool( const ObservationDataset< ObservationScalarType, TimeType >&, const ObservationId ) >;

    ObservationCondition( ):
        evaluator_( []( const ObservationDataset< ObservationScalarType, TimeType >&, const ObservationId ) { return true; } )
    {}

    explicit ObservationCondition( const Evaluator& evaluator ): evaluator_( evaluator ) {}

    bool operator( )( const ObservationDataset< ObservationScalarType, TimeType >& dataset, const ObservationId observationId ) const
    {
        return evaluator_( dataset, observationId );
    }

    ObservationCondition operator&&( const ObservationCondition& other ) const
    {
        const Evaluator leftEvaluator = evaluator_;
        const Evaluator rightEvaluator = other.evaluator_;
        return ObservationCondition(
                [ leftEvaluator, rightEvaluator ]( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                                   const ObservationId observationId ) {
                    return leftEvaluator( dataset, observationId ) && rightEvaluator( dataset, observationId );
                } );
    }

    ObservationCondition operator||( const ObservationCondition& other ) const
    {
        const Evaluator leftEvaluator = evaluator_;
        const Evaluator rightEvaluator = other.evaluator_;
        return ObservationCondition(
                [ leftEvaluator, rightEvaluator ]( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                                   const ObservationId observationId ) {
                    return leftEvaluator( dataset, observationId ) || rightEvaluator( dataset, observationId );
                } );
    }

    ObservationCondition operator!( ) const
    {
        const Evaluator evaluator = evaluator_;
        return ObservationCondition( [ evaluator ]( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                                    const ObservationId observationId ) { return !evaluator( dataset, observationId ); } );
    }

    static ObservationCondition all( )
    {
        return ObservationCondition( );
    }

    static ObservationCondition observableType( const ObservableType observableType )
    {
        return ObservationCondition( [ observableType ]( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                                         const ObservationId observationId ) {
            const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
            return dataset.getObservationSetMetadata( row.setId_ ).observableType_ == observableType;
        } );
    }

    static ObservationCondition linkDefinition( const LinkDefinition& linkDefinition )
    {
        return ObservationCondition( [ linkDefinition ]( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                                         const ObservationId observationId ) {
            const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
            const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = dataset.getObservationSetMetadata( row.setId_ );
            return dataset.getLinkDefinition( metadata.linkDefinitionId_ ) == linkDefinition;
        } );
    }

    static ObservationCondition linkEndType( const LinkEndType linkEndType )
    {
        return ObservationCondition( [ linkEndType ]( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                                      const ObservationId observationId ) {
            const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
            const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = dataset.getObservationSetMetadata( row.setId_ );
            return dataset.getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_.count( linkEndType ) > 0;
        } );
    }

    static ObservationCondition linkEnd( const LinkEndType linkEndType, const LinkEndId& linkEndId )
    {
        return ObservationCondition( [ linkEndType, linkEndId ]( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                                                 const ObservationId observationId ) {
            const ObservationDatasetRow< TimeType >& row = dataset.getObservationRow( observationId );
            const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = dataset.getObservationSetMetadata( row.setId_ );
            const LinkEnds& linkEnds = dataset.getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
            return linkEnds.count( linkEndType ) > 0 && linkEnds.at( linkEndType ) == linkEndId;
        } );
    }

    static ObservationCondition timeBounds( const TimeType startTime, const TimeType endTime )
    {
        return ObservationCondition( [ startTime, endTime ]( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                                             const ObservationId observationId ) {
            const TimeType observationTime = dataset.getObservationTime( observationId );
            return observationTime >= startTime && observationTime <= endTime;
        } );
    }

    static ObservationCondition active( )
    {
        return ObservationCondition(
                []( const ObservationDataset< ObservationScalarType, TimeType >& dataset, const ObservationId observationId ) {
                    return dataset.getObservationRow( observationId ).isActive_;
                } );
    }

    static ObservationCondition rejected( )
    {
        return !active( );
    }

    static ObservationCondition residualAbsoluteValueGreaterThan( const Eigen::VectorXd& residualLimit )
    {
        return ObservationCondition( [ residualLimit ]( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                                        const ObservationId observationId ) {
            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residual = dataset.getResidualValue( observationId );
            if( residualLimit.size( ) != residual.size( ) )
            {
                throw std::runtime_error( "Error when evaluating residual condition, limit size is incompatible with observation size." );
            }
            for( int i = 0; i < residual.size( ); ++i )
            {
                if( std::fabs( residual( i ) ) > residualLimit( i ) )
                {
                    return true;
                }
            }
            return false;
        } );
    }

    static ObservationCondition observationAbsoluteValueGreaterThan( const Eigen::VectorXd& observationLimit )
    {
        return ObservationCondition( [ observationLimit ]( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                                           const ObservationId observationId ) {
            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observation = dataset.getObservationValue( observationId );
            if( observationLimit.size( ) != observation.size( ) )
            {
                throw std::runtime_error(
                        "Error when evaluating observation condition, limit size is incompatible with observation size." );
            }
            for( int i = 0; i < observation.size( ); ++i )
            {
                if( std::fabs( observation( i ) ) > observationLimit( i ) )
                {
                    return true;
                }
            }
            return false;
        } );
    }

    static ObservationCondition dependentVariableGreaterThan(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
            const Eigen::VectorXd& dependentVariableLimit,
            const bool returnFirstCompatibleSettings = false )
    {
        return ObservationCondition( [ dependentVariableSettings, dependentVariableLimit, returnFirstCompatibleSettings ](
                                             const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                             const ObservationId observationId ) {
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
                if( dependentVariableSettings->areSettingsCompatible( settingsIterator.second ) )
                {
                    compatibleIndicesAndSizes.push_back( settingsIterator.first );
                }
            }

            if( compatibleIndicesAndSizes.empty( ) )
            {
                return false;
            }
            if( compatibleIndicesAndSizes.size( ) > 1 && !returnFirstCompatibleSettings )
            {
                throw std::runtime_error(
                        "Error when evaluating dependent-variable condition, multiple compatible dependent variables found." );
            }

            const std::pair< int, int >& indexAndSize = compatibleIndicesAndSizes.at( 0 );
            const Eigen::VectorXd dependentVariables = dataset.getDependentVariables( observationId );
            if( indexAndSize.first + indexAndSize.second > dependentVariables.size( ) ||
                dependentVariableLimit.size( ) != indexAndSize.second )
            {
                throw std::runtime_error( "Error when evaluating dependent-variable condition, dependent-variable size is inconsistent." );
            }

            const Eigen::VectorXd selectedDependentVariables = dependentVariables.segment( indexAndSize.first, indexAndSize.second );
            for( int i = 0; i < selectedDependentVariables.size( ); ++i )
            {
                if( selectedDependentVariables( i ) > dependentVariableLimit( i ) )
                {
                    return true;
                }
            }
            return false;
        } );
    }

private:
    //! Callable that evaluates one observation row in the context of a dataset.
    Evaluator evaluator_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_CONDITION_H
