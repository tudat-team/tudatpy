#ifndef TUDAT_SINGLE_OBSERVATION_SET_UTILITIES_H
#define TUDAT_SINGLE_OBSERVATION_SET_UTILITIES_H

#include "tudat/simulation/estimation_setup/singleObservationSet.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
//! Creates a filtered copy of a single observation set, preserving weights and residuals.
std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > filterObservations(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > singleObservationSet,
        const std::shared_ptr< ObservationFilterBase > observationFilter,
        const bool saveFilteredObservations )
{
    if( !observationFilter->filterOut( ) )
    {
        throw std::runtime_error(
                "Error when creating new single observation set post-filtering, the filterOut "
                "option should be set to true" );
    }

    // Clone metadata, values and weights before applying the filter in-place on the clone.
    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > newObservationSet =
            std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                    singleObservationSet->getObservableType( ),
                    singleObservationSet->getLinkEnds( ),
                    singleObservationSet->getObservationsReference( ),
                    singleObservationSet->getObservationTimesReference( ),
                    singleObservationSet->getReferenceLinkEnd( ),
                    singleObservationSet->getObservationsDependentVariablesReference( ),
                    singleObservationSet->getDependentVariableBookkeeping( ),
                    singleObservationSet->getAncillarySettings( ) );
    if( singleObservationSet->getWeightsMatrixType( ) == diagonal_weights_matrix )
    {
        newObservationSet->setTabulatedWeights( singleObservationSet->getBaseWeightsDiagonalVector( ) );
    }
    else if( singleObservationSet->getWeightsMatrixType( ) == block_diagonal_weights_matrix )
    {
        newObservationSet->setBlockDiagonalWeights( singleObservationSet->getBlockDiagonalWeightMatrices( ) );
    }
    if( singleObservationSet->hasFullWeightMatrixContribution( ) )
    {
        newObservationSet->setFullWeightMatrix( singleObservationSet->getFullWeightMatrix( ) );
    }
    newObservationSet->setResiduals( singleObservationSet->getResidualsReference( ) );

    newObservationSet->filterObservations( observationFilter, saveFilteredObservations );
    return newObservationSet;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
//! Splits one observation set into multiple sets according to a splitter rule.
std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > splitObservationSet(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > observationSet,
        const std::shared_ptr< ObservationSetSplitterBase > observationSetSplitter,
        const bool printWarning )
{
    if( observationSet->hasFullWeightMatrixContribution( ) )
    {
        throw std::runtime_error(
                "Error when splitting single observation set: this operation is not supported when full weights are defined." );
    }

    if( printWarning && observationSet->getFilteredObservationSet( ) != nullptr )
    {
        std::cerr << "Warning when splitting single observation set, the filtered observation set "
                     "pointer is not empty and "
                     " any filtered observation will be lost after splitting."
                  << std::endl;
    }

    // Gather source data that will be partitioned into new observation sets.
    std::vector< int > rawStartIndicesNewSets = { 0 };
    std::vector< TimeType > observationTimes = observationSet->getObservationTimes( );
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations = observationSet->getObservations( );
    std::vector< Eigen::VectorXd > dependentVariables = observationSet->getObservationsDependentVariables( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > weightsVector = observationSet->getBaseWeightsDiagonalVector( );
    ObservationWeightsMatrixType weightsMatrixType = observationSet->getWeightsMatrixType( );
    std::vector< Eigen::MatrixXd > blockWeights;
    if( weightsMatrixType == block_diagonal_weights_matrix )
    {
        blockWeights = observationSet->getBlockDiagonalWeightMatrices( );
    }

    // Compute raw split boundaries according to the selected splitter mode.
    switch( observationSetSplitter->getSplitterType( ) )
    {
        case time_tags_splitter: {
            std::vector< double > timeTags =
                    std::dynamic_pointer_cast< ObservationSetSplitter< std::vector< double > > >( observationSetSplitter )
                            ->getSplitterValue( );
            for( auto currentTimeTag : timeTags )
            {
                if( currentTimeTag > observationSet->getTimeBounds( ).first )
                {
                    bool detectedStartSet = false;
                    int indexObs = rawStartIndicesNewSets.at( rawStartIndicesNewSets.size( ) - 1 );
                    while( !detectedStartSet && indexObs < static_cast< int >( observationTimes.size( ) ) )
                    {
                        if( observationTimes.at( indexObs ) > currentTimeTag )
                        {
                            rawStartIndicesNewSets.push_back( indexObs );
                            detectedStartSet = true;
                        }
                        indexObs++;
                    }
                }
            }
            rawStartIndicesNewSets.push_back( static_cast< int >( observationTimes.size( ) ) );
            break;
        }
        case time_interval_splitter: {
            double maxTimeInterval =
                    std::dynamic_pointer_cast< ObservationSetSplitter< double > >( observationSetSplitter )->getSplitterValue( );
            for( unsigned int i = 1; i < observationTimes.size( ); i++ )
            {
                if( ( observationTimes.at( i ) - observationTimes.at( i - 1 ) ) > maxTimeInterval )
                {
                    rawStartIndicesNewSets.push_back( i );
                }
            }
            rawStartIndicesNewSets.push_back( static_cast< int >( observationTimes.size( ) ) );
            break;
        }
        case time_span_splitter: {
            double maxTimeSpan =
                    std::dynamic_pointer_cast< ObservationSetSplitter< double > >( observationSetSplitter )->getSplitterValue( );
            if( observationSet->getNumberOfObservables( ) > 0 )
            {
                double referenceEpoch = observationTimes.at( 0 );
                for( unsigned int i = 1; i < observationTimes.size( ); i++ )
                {
                    if( ( observationTimes.at( i ) - referenceEpoch ) > maxTimeSpan )
                    {
                        rawStartIndicesNewSets.push_back( i );
                        referenceEpoch = observationTimes.at( i );
                    }
                }
                rawStartIndicesNewSets.push_back( static_cast< int >( observationTimes.size( ) ) );
            }
            break;
        }
        case nb_observations_splitter: {
            int maxNbObs = std::dynamic_pointer_cast< ObservationSetSplitter< int > >( observationSetSplitter )->getSplitterValue( );
            if( maxNbObs < observationSetSplitter->getMinNumberObservations( ) )
            {
                throw std::runtime_error(
                        "Error when splitting observation sets, the maximum number of observations "
                        "cannot be smaller than the minimum number of observations." );
            }
            for( int ind = maxNbObs; ind < static_cast< int >( observationSet->getNumberOfObservables( ) ); ind += maxNbObs )
            {
                rawStartIndicesNewSets.push_back( ind );
            }
            rawStartIndicesNewSets.push_back( static_cast< int >( observationTimes.size( ) ) );
            break;
        }
        default:
            throw std::runtime_error( "Observation set splitter type not recognised." );
    }

    // Convert raw boundaries to [start, size] pairs while enforcing minimum set size.
    std::vector< std::pair< int, int > > indicesNewSets;
    for( unsigned int j = 1; j < rawStartIndicesNewSets.size( ); j++ )
    {
        if( ( rawStartIndicesNewSets.at( j ) - rawStartIndicesNewSets.at( j - 1 ) ) >= observationSetSplitter->getMinNumberObservations( ) )
        {
            indicesNewSets.push_back( std::make_pair( rawStartIndicesNewSets.at( j - 1 ),
                                                      rawStartIndicesNewSets.at( j ) - rawStartIndicesNewSets.at( j - 1 ) ) );
        }
    }

    // Materialize each split as a new single observation set with consistent weights/residuals.
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > newObsSets;
    for( unsigned int k = 0; k < indicesNewSets.size( ); k++ )
    {
        int startIndex = indicesNewSets.at( k ).first;
        int sizeCurrentSet = indicesNewSets.at( k ).second;

        std::vector< Eigen::VectorXd > newDependentVariables;
        if( !dependentVariables.empty( ) )
        {
            newDependentVariables = utilities::getStlVectorSegment(
                    observationSet->getObservationsDependentVariablesReference( ), startIndex, sizeCurrentSet );
        }

        std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > newSet =
                std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                        observationSet->getObservableType( ),
                        observationSet->getLinkEnds( ),
                        utilities::getStlVectorSegment( observationSet->getObservationsReference( ), startIndex, sizeCurrentSet ),
                        utilities::getStlVectorSegment( observationSet->getObservationTimesReference( ), startIndex, sizeCurrentSet ),
                        observationSet->getReferenceLinkEnd( ),
                        newDependentVariables,
                        observationSet->getDependentVariableBookkeeping( ),
                        observationSet->getAncillarySettings( ) );

        if( weightsMatrixType == diagonal_weights_matrix )
        {
            Eigen::Matrix< double, Eigen::Dynamic, 1 > newWeightsVector =
                    weightsVector.segment( startIndex * observationSet->getSingleObservableSize( ),
                                           sizeCurrentSet * observationSet->getSingleObservableSize( ) );
            newSet->setTabulatedWeights( newWeightsVector );
        }
        else if( weightsMatrixType == block_diagonal_weights_matrix )
        {
            std::vector< Eigen::MatrixXd > newBlockWeights;
            newBlockWeights.reserve( sizeCurrentSet );
            for( int i = 0; i < sizeCurrentSet; i++ )
            {
                newBlockWeights.push_back( blockWeights.at( startIndex + i ) );
            }
            newSet->setBlockDiagonalWeights( newBlockWeights );
        }

        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > newResiduals =
                utilities::getStlVectorSegment( observationSet->getResidualsReference( ), startIndex, sizeCurrentSet );
        newSet->setResiduals( newResiduals );

        newObsSets.push_back( newSet );
    }

    return newObsSets;
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_SINGLE_OBSERVATION_SET_UTILITIES_H
