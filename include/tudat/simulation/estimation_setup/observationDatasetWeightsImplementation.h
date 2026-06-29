/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONDATASETWEIGHTSIMPLEMENTATION_H
#define TUDAT_OBSERVATIONDATASETWEIGHTSIMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantWeightForSet( const unsigned int setId,
                                                                                            const Eigen::VectorXd& weight )
{
    if( weight.size( ) != static_cast< int >( getObservationSetMetadata( setId ).observableSize_ ) )
    {
        throw std::runtime_error( "Error when setting dataset weights, weight size is inconsistent." );
    }
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        setWeightValue( observationId, weight );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightVectorForSet( const unsigned int setId,
                                                                                          const Eigen::VectorXd& weightVector )
{
    const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
    const std::vector< unsigned int >& observationIds = observationIdsBySet_.at( setId );
    if( weightVector.size( ) != static_cast< int >( observationIds.size( ) * observableSize ) )
    {
        throw std::runtime_error( "Error when setting dataset weights, vector size is inconsistent." );
    }
    for( std::size_t i = 0; i < observationIds.size( ); ++i )
    {
        setWeightValue( observationIds.at( i ), weightVector.segment( i * observableSize, observableSize ) );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightMatrixForSet( const unsigned int setId,
                                                                                          const Eigen::MatrixXd& weightMatrix )
{
    if( weightMatrix.rows( ) != static_cast< int >( getTotalScalarSizeForSet( setId ) ) ||
        weightMatrix.cols( ) != static_cast< int >( getTotalScalarSizeForSet( setId ) ) )
    {
        throw std::runtime_error( "Error when setting dataset set weight matrix, matrix size is inconsistent." );
    }
    observationWeights_.setSetWeightBlock( setId, weightMatrix );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
bool ObservationDataset< ObservationScalarType, TimeType, Dummy >::hasWeightMatrixForSet( const unsigned int setId ) const
{
    return observationWeights_.hasSetWeightBlock( setId );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightMatrixForObservation( const unsigned int observationId,
                                                                                                  const Eigen::MatrixXd& weightMatrix )
{
    const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
    if( weightMatrix.rows( ) != static_cast< int >( row.scalarSize_ ) || weightMatrix.cols( ) != static_cast< int >( row.scalarSize_ ) )
    {
        throw std::runtime_error( "Error when setting dataset observation weight matrix, matrix size is inconsistent." );
    }
    observationWeights_.setWeightBlock( observationId, weightMatrix );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
bool ObservationDataset< ObservationScalarType, TimeType, Dummy >::hasWeightMatrixForObservation( const unsigned int observationId ) const
{
    return observationWeights_.hasObservationWeightBlock( observationId );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::addExtraWeightBlock( const ObservationWeightBlock& weightBlock )
{
    observationWeights_.addExtraWeightBlock( weightBlock );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightBlock( const std::vector< unsigned int >& rowObservationIds,
                                                                                   const std::vector< unsigned int >& columnObservationIds,
                                                                                   const Eigen::MatrixXd& weightBlock,
                                                                                   const std::vector< unsigned int >& rowComponents,
                                                                                   const std::vector< unsigned int >& columnComponents,
                                                                                   const bool makeSymmetric )
{
    const std::vector< unsigned int > rowScalarComponentIds =
            getScalarComponentIdsForObservationSelection( rowObservationIds, rowComponents );
    const std::vector< unsigned int > columnScalarComponentIds =
            getScalarComponentIdsForObservationSelection( columnObservationIds, columnComponents );

    if( weightBlock.rows( ) != static_cast< int >( rowScalarComponentIds.size( ) ) ||
        weightBlock.cols( ) != static_cast< int >( columnScalarComponentIds.size( ) ) )
    {
        throw std::runtime_error( "Error when setting dataset weight block, matrix size is inconsistent with selected observations." );
    }

    ObservationWeightBlock datasetWeightBlock;
    datasetWeightBlock.rowScalarComponentIds_ = rowScalarComponentIds;
    datasetWeightBlock.columnScalarComponentIds_ = columnScalarComponentIds;
    datasetWeightBlock.weightBlock_ = weightBlock;
    addExtraWeightBlock( datasetWeightBlock );

    if( makeSymmetric )
    {
        if( rowScalarComponentIds == columnScalarComponentIds )
        {
            if( weightBlock.rows( ) != weightBlock.cols( ) || !weightBlock.isApprox( weightBlock.transpose( ) ) )
            {
                throw std::runtime_error(
                        "Error when setting symmetric dataset weight block, block with identical row and column selection is not "
                        "symmetric." );
            }
        }
        else
        {
            ObservationWeightBlock transposedWeightBlock;
            transposedWeightBlock.rowScalarComponentIds_ = columnScalarComponentIds;
            transposedWeightBlock.columnScalarComponentIds_ = rowScalarComponentIds;
            transposedWeightBlock.weightBlock_ = weightBlock.transpose( );
            addExtraWeightBlock( transposedWeightBlock );
        }
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const std::vector< ObservationWeightBlock >& ObservationDataset< ObservationScalarType, TimeType, Dummy >::getExtraWeightBlocks( ) const
{
    return observationWeights_.getExtraWeightBlocks( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
bool ObservationDataset< ObservationScalarType, TimeType, Dummy >::hasExtraWeightBlocks( ) const
{
    return observationWeights_.hasExtraWeightBlocks( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantWeight(
        const double weight,
        const std::shared_ptr< ObservationCollectionParser > observationParser )
{
    std::vector< unsigned int > setIds = getObservationSetIds( observationParser );
    if( setIds.empty( ) )
    {
        std::cerr << "Warning when setting constant weights, no observation dataset set found for specified observation parser. "
                     "Weights not set";
    }
    for( const unsigned int setId : setIds )
    {
        setConstantWeightForSet( setId, weight );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantWeight(
        const Eigen::VectorXd weight,
        const std::shared_ptr< ObservationCollectionParser > observationParser )
{
    std::vector< unsigned int > setIds = getObservationSetIds( observationParser );
    if( setIds.empty( ) )
    {
        std::cerr << "Warning when setting constant weights, no observation dataset set found for specified observation parser. "
                     "Weights not set";
    }
    for( const unsigned int setId : setIds )
    {
        setConstantWeightForSet( setId, weight );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantWeightForObservableType( const ObservableType observableType,
                                                                                                       const double weight )
{
    setConstantWeight( weight, observationParser( observableType ) );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantWeightForObservableType( const ObservableType observableType,
                                                                                                       const Eigen::VectorXd& weight )
{
    setConstantWeight( weight, observationParser( observableType ) );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantWeightPerObservableType(
        const std::map< ObservableType, double >& weightsPerObservableType )
{
    for( const auto& weightIterator : weightsPerObservableType )
    {
        setConstantWeightForObservableType( weightIterator.first, weightIterator.second );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantWeightPerObservableType(
        const std::initializer_list< std::pair< ObservableType, double > > weightsPerObservableType )
{
    for( const auto& weightIterator : weightsPerObservableType )
    {
        setConstantWeightForObservableType( weightIterator.first, weightIterator.second );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantWeightPerObservableType(
        const std::map< ObservableType, Eigen::VectorXd >& weightsPerObservableType )
{
    for( const auto& weightIterator : weightsPerObservableType )
    {
        setConstantWeightForObservableType( weightIterator.first, weightIterator.second );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantWeightPerObservable(
        const std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser )
{
    for( auto parserIt : weightsPerObservationParser )
    {
        setConstantWeight( parserIt.second, parserIt.first );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantWeightPerObservable(
        const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd > weightsPerObservationParser )
{
    for( auto parserIt : weightsPerObservationParser )
    {
        setConstantWeight( parserIt.second, parserIt.first );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setTabulatedWeights(
        const Eigen::VectorXd tabulatedWeights,
        const std::shared_ptr< ObservationCollectionParser > observationParser )
{
    std::vector< unsigned int > setIds = getObservationSetIds( observationParser );
    if( setIds.empty( ) )
    {
        std::cerr << "Warning when setting tabulated weights, no observation dataset set found for specified observation parser. "
                     "Weights not set";
        return;
    }

    bool areSetsSameSize = true;
    int totalSizeAllSets = static_cast< int >( getTotalScalarSizeForSet( setIds.at( 0 ) ) );
    for( unsigned int i = 1; i < setIds.size( ); ++i )
    {
        const int currentSetSize = static_cast< int >( getTotalScalarSizeForSet( setIds.at( i ) ) );
        totalSizeAllSets += currentSetSize;
        if( currentSetSize != static_cast< int >( getTotalScalarSizeForSet( setIds.at( 0 ) ) ) )
        {
            areSetsSameSize = false;
        }
    }

    int startSet = 0;
    for( const unsigned int setId : setIds )
    {
        const int currentSetSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
        if( tabulatedWeights.size( ) == totalSizeAllSets )
        {
            setWeightVectorForSet( setId, tabulatedWeights.segment( startSet, currentSetSize ) );
            startSet += currentSetSize;
        }
        else if( areSetsSameSize && tabulatedWeights.size( ) == static_cast< int >( getTotalScalarSizeForSet( setIds.at( 0 ) ) ) )
        {
            setWeightVectorForSet( setId, tabulatedWeights );
        }
        else
        {
            throw std::runtime_error(
                    "Error when setting tabulated weights, the size of the input weight vector should be consistent with either the "
                    "size of each individual observation set, or the combined size of all selected observation sets." );
        }
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setTabulatedWeights(
        const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd > weightsPerObservationParser )
{
    for( auto parserIt : weightsPerObservationParser )
    {
        setTabulatedWeights( parserIt.second, parserIt.first );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::replaceObservationSetData(
        const unsigned int setId,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType >& times,
        const std::vector< Eigen::VectorXd >& dependentVariables,
        const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals )
{
    replaceObservationSetDataWithSourceRows(
            setId, observations, times, dependentVariables, weights, residuals, std::vector< unsigned int >( ) );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETWEIGHTSIMPLEMENTATION_H
