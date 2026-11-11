/* Implementation of the ObservationDataset functions that retrieve the data of selected observations. */
#ifndef TUDAT_OBSERVATION_DATASET_INSPECTION_IMPLEMENTATION_H
#define TUDAT_OBSERVATION_DATASET_INSPECTION_IMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{
namespace observation_models
{
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< double > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getNumericalObservationMetadata(
        const std::string& key,
        const std::vector< unsigned int >& observationIds ) const
{
    const auto& metadata = numericalObservationMetadata_.at( key );
    std::vector< double > values;
    values.reserve( observationIds.size( ) );
    for( const auto id : observationIds )
    {
        values.push_back( metadata.at( id ) );
    }
    return values;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::resolveObservationIds(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering,
        const bool includeRejected ) const
{
    if( ordering != ObservationOrdering::internal && ordering != ObservationOrdering::estimation )
    {
        throw std::invalid_argument( "Observation ordering must be internal or estimation." );
    }
    // Select the observations in internal order; the requested order does not affect the selection.
    std::vector< unsigned int > selected;
    selected.reserve( observationRows_.size( ) );
    for( const auto& row : observationRows_ )
    {
        if( ( includeRejected || row.isActive_ ) && condition( *this, row.observationId_ ) )
        {
            selected.push_back( row.observationId_ );
        }
    }
    if( ordering == ObservationOrdering::internal )
    {
        return selected;
    }
    const std::unordered_set< unsigned int > membership( selected.begin( ), selected.end( ) );
    selected.clear( );
    // Estimation order: by observable type, link ends, set, and order within the set.
    for( const auto id : getObservationIdsInObservationVectorOrder( ) )
    {
        if( membership.count( id ) != 0 )
        {
            selected.push_back( id );
        }
    }
    return selected;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< TimeType > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getTimes(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    std::vector< TimeType > result;
    for( const auto id : resolveObservationIds( condition, ordering ) )
    {
        result.push_back( getObservationTime( id ) );
    }
    return result;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservations(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > result;
    for( const auto id : resolveObservationIds( condition, ordering ) )
    {
        result.push_back( getObservationValue( id ) );
    }
    return result;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getResiduals(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > result;
    for( const auto id : resolveObservationIds( condition, ordering ) )
    {
        result.push_back( getResidualValue( id ) );
    }
    return result;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getSetIds(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    std::vector< unsigned int > result;
    for( const auto id : resolveObservationIds( condition, ordering ) )
    {
        result.push_back( getObservationRow( id ).setId_ );
    }
    return result;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< ObservationDatasetRow< TimeType > > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getRows(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    std::vector< ObservationDatasetRow< TimeType > > result;
    for( const auto id : resolveObservationIds( condition, ordering ) )
    {
        result.push_back( getObservationRow( id ) );
    }
    return result;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::VectorXd > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getDependentVariableValues(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    std::vector< Eigen::VectorXd > result;
    for( const auto id : resolveObservationIds( condition, ordering ) )
    {
        result.push_back( getDependentVariables( id ) );
    }
    return result;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< std::pair< unsigned int, unsigned int > > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getScalarComponents(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    std::vector< std::pair< unsigned int, unsigned int > > result;
    for( const auto id : resolveObservationIds( condition, ordering ) )
    {
        for( unsigned int component = 0; component < getObservationRow( id ).scalarSize_; ++component )
        {
            result.emplace_back( id, component );
        }
    }
    return result;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
typename ObservationDataset< ObservationScalarType, TimeType, Dummy >::InspectionMetadata
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getMetadata(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    InspectionMetadata result;
    for( const auto id : resolveObservationIds( condition, ordering ) )
    {
        const auto setId = getObservationRow( id ).setId_;
        if( result.count( setId ) != 0 )
        {
            continue;
        }
        const auto& metadata = getObservationSetMetadata( setId );
        const auto& ancillary = getAncillarySettings( metadata.ancillarySettingsId_ );
        const auto& bookkeeping = getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
        std::map< std::pair< int, int >, std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > layout;
        if( bookkeeping )
        {
            for( const auto& entry : bookkeeping->getSettingsIndicesAndSizes( ) )
            {
                layout.emplace( entry.first, entry.second ? entry.second->clone( ) : nullptr );
            }
        }
        result.emplace( setId,
                        std::make_tuple( metadata,
                                         getLinkDefinition( metadata.linkDefinitionId_ ),
                                         ancillary ? std::make_shared< ObservationAncillarySimulationSettings >( *ancillary ) : nullptr,
                                         std::move( layout ) ) );
    }
    return result;
}

}  // namespace observation_models
}  // namespace tudat

#endif
