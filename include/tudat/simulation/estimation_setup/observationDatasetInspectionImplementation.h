/* Detached observation inspection: implementation details, not a public view API. */
#ifndef TUDAT_OBSERVATION_DATASET_INSPECTION_IMPLEMENTATION_H
#define TUDAT_OBSERVATION_DATASET_INSPECTION_IMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{
namespace observation_models
{
namespace detail
{
//! Transient event indices for a single extraction. No source ownership, payload or revision state.
//! Used only while its source dataset is unchanged, within the existing mutation model.
template< typename ObservationScalarType, typename TimeType >
struct ObservationSelectionIndices {
    using Dataset = ObservationDataset< ObservationScalarType, TimeType >;
    using Vector = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >;

    //! Resolve selected observation identities once for the requested ordering.
    ObservationSelectionIndices( const Dataset& dataset,
                                 const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
                                 const ObservationOrdering ordering ): ids( dataset.resolveObservationIds( condition, ordering ) )
    {}

    //! Gather one detached value per selected observation identity.
    template< typename Value, typename Getter >
    std::vector< Value > gather( Getter getter ) const
    {
        std::vector< Value > result;
        result.reserve( ids.size( ) );
        for( const auto id : ids )
        {
            result.push_back( getter( id ) );
        }
        return result;
    }

    //! Gather event times for the current selection.
    std::vector< TimeType > getTimes( const Dataset& dataset ) const
    {
        return gather< TimeType >( [ & ]( const auto id ) { return dataset.getObservationTime( id ); } );
    }

    //! Gather observed vectors for the current selection.
    std::vector< Vector > getObservations( const Dataset& dataset ) const
    {
        return gather< Vector >( [ & ]( const auto id ) { return dataset.getObservationValue( id ); } );
    }

    //! Gather residual vectors for the current selection.
    std::vector< Vector > getResiduals( const Dataset& dataset ) const
    {
        return gather< Vector >( [ & ]( const auto id ) { return dataset.getResidualValue( id ); } );
    }

    //! Return the selected stable observation identities.
    std::vector< unsigned int > getObservationIds( const Dataset& ) const
    {
        return ids;
    }

    //! Gather owning set identities for the current selection.
    std::vector< unsigned int > getSetIds( const Dataset& dataset ) const
    {
        return gather< unsigned int >( [ & ]( const auto id ) { return dataset.getObservationRow( id ).setId_; } );
    }

    //! Gather detached event rows for the current selection.
    std::vector< ObservationDatasetRow< TimeType > > getRows( const Dataset& dataset ) const
    {
        return gather< ObservationDatasetRow< TimeType > >( [ & ]( const auto id ) { return dataset.getObservationRow( id ); } );
    }

    //! Gather dependent-variable vectors for the current selection.
    std::vector< Eigen::VectorXd > getDependentVariableValues( const Dataset& dataset ) const
    {
        return gather< Eigen::VectorXd >( [ & ]( const auto id ) { return dataset.getDependentVariables( id ); } );
    }

    //! Expand selected events to stable observation/component pairs.
    std::vector< std::pair< unsigned int, unsigned int > > getScalarComponents( const Dataset& dataset ) const
    {
        std::vector< std::pair< unsigned int, unsigned int > > result;
        for( const auto id : ids )
        {
            for( unsigned int component = 0; component < dataset.getObservationRow( id ).scalarSize_; ++component )
            {
                result.emplace_back( id, component );
            }
        }
        return result;
    }

    //! Gather the effective diagonal weights for the current selection.
    Eigen::VectorXd getWeightDiagonal( const Dataset& dataset ) const
    {
        return dataset.observationWeights_.getDiagonal( dataset.getScalarComponentIdsForObservationSelection( ids, {} ) );
    }

    //! Gather the effective principal weight matrix for the current selection.
    Eigen::SparseMatrix< double > getWeightMatrix( const Dataset& dataset ) const
    {
        return dataset.observationWeights_.restricted( dataset.getScalarComponentIdsForObservationSelection( ids, {} ) ).sparseMatrix( );
    }

    //! Clone metadata for sets represented in the current selection.
    typename Dataset::InspectionMetadata getMetadata( const Dataset& dataset ) const
    {
        typename Dataset::InspectionMetadata result;
        for( const auto id : ids )
        {
            const auto setId = dataset.getObservationRow( id ).setId_;
            if( result.count( setId ) != 0 )
            {
                continue;
            }
            const auto& metadata = dataset.getObservationSetMetadata( setId );
            const auto& ancillary = dataset.getAncillarySettings( metadata.ancillarySettingsId_ );
            const auto& bookkeeping = dataset.getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
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
                                             dataset.getLinkDefinition( metadata.linkDefinitionId_ ),
                                             ancillary ? std::make_shared< ObservationAncillarySimulationSettings >( *ancillary ) : nullptr,
                                             std::move( layout ) ) );
        }
        return result;
    }

    //! Stable observation identities in the requested extraction order.
    std::vector< unsigned int > ids;
};
}  // namespace detail

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
    // Resolve membership in storage order once, independently of the output sequence.
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
    // This is the existing authoritative observable/link/set/event ordering, including equal-time ties.
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
    return detail::ObservationSelectionIndices< ObservationScalarType, TimeType >( *this, condition, ordering ).getTimes( *this );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservations(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    return detail::ObservationSelectionIndices< ObservationScalarType, TimeType >( *this, condition, ordering ).getObservations( *this );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getResiduals(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    return detail::ObservationSelectionIndices< ObservationScalarType, TimeType >( *this, condition, ordering ).getResiduals( *this );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationIds(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    return detail::ObservationSelectionIndices< ObservationScalarType, TimeType >( *this, condition, ordering ).getObservationIds( *this );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getSetIds(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    return detail::ObservationSelectionIndices< ObservationScalarType, TimeType >( *this, condition, ordering ).getSetIds( *this );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< ObservationDatasetRow< TimeType > > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getRows(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    return detail::ObservationSelectionIndices< ObservationScalarType, TimeType >( *this, condition, ordering ).getRows( *this );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::VectorXd > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getDependentVariableValues(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    return detail::ObservationSelectionIndices< ObservationScalarType, TimeType >( *this, condition, ordering )
            .getDependentVariableValues( *this );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< std::pair< unsigned int, unsigned int > > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getScalarComponents(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    return detail::ObservationSelectionIndices< ObservationScalarType, TimeType >( *this, condition, ordering )
            .getScalarComponents( *this );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::VectorXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getWeightDiagonal(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    return detail::ObservationSelectionIndices< ObservationScalarType, TimeType >( *this, condition, ordering ).getWeightDiagonal( *this );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::SparseMatrix< double > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getWeightMatrix(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    return detail::ObservationSelectionIndices< ObservationScalarType, TimeType >( *this, condition, ordering ).getWeightMatrix( *this );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
typename ObservationDataset< ObservationScalarType, TimeType, Dummy >::InspectionMetadata
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getMetadata(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const ObservationOrdering ordering ) const
{
    return detail::ObservationSelectionIndices< ObservationScalarType, TimeType >( *this, condition, ordering ).getMetadata( *this );
}

}  // namespace observation_models
}  // namespace tudat

#endif
