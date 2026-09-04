/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONDATASETLEGACYIMPLEMENTATION_H
#define TUDAT_OBSERVATIONDATASETLEGACYIMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setObservationVectorForSet(
        const unsigned int setId,
        const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationVector )
{
    const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
    const std::vector< unsigned int >& observationIds = observationIdsBySet_.at( setId );
    if( observationVector.size( ) != static_cast< int >( observationIds.size( ) * observableSize ) )
    {
        throw std::runtime_error( "Error when setting dataset observations, vector size is inconsistent." );
    }
    for( std::size_t i = 0; i < observationIds.size( ); ++i )
    {
        setObservationValue( observationIds.at( i ), observationVector.segment( i * observableSize, observableSize ) );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setResidualVectorForSet(
        const unsigned int setId,
        const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualVector )
{
    const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
    const std::vector< unsigned int >& observationIds = observationIdsBySet_.at( setId );
    if( residualVector.size( ) != static_cast< int >( observationIds.size( ) * observableSize ) )
    {
        throw std::runtime_error( "Error when setting dataset residuals, vector size is inconsistent." );
    }
    for( std::size_t i = 0; i < observationIds.size( ); ++i )
    {
        setResidualValue( observationIds.at( i ), residualVector.segment( i * observableSize, observableSize ) );
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
        setConstantSingleObservationScalarWeightForSet( setId, weight );
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
        setConstantSingleObservationDiagonalWeightForSet( setId, weight );
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
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getFilteredObservationIndices(
        const unsigned int setId,
        const std::shared_ptr< ObservationFilterBase >& observationFilter ) const
{
    const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
    const unsigned int numberOfObservations = static_cast< unsigned int >( getNumberOfObservationsForSet( setId ) );
    const bool useOppositeCondition = observationFilter->useOppositeCondition( );

    std::vector< unsigned int > indicesToRemove;
    switch( observationFilter->getFilterType( ) )
    {
        case residual_filtering: {
            Eigen::VectorXd residualCutOff = Eigen::VectorXd::Zero( observableSize );
            if( std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter ) != nullptr )
            {
                residualCutOff = std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter )->getFilterValue( ) *
                        Eigen::VectorXd::Ones( observableSize );
            }
            else if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter ) != nullptr )
            {
                if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( ).size( ) !=
                    static_cast< int >( observableSize ) )
                {
                    throw std::runtime_error(
                            "Error when performing residual filtering, size of the residual cut off vector inconsistent with "
                            "observable size." );
                }
                residualCutOff = std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( );
            }

            for( unsigned int j = 0; j < numberOfObservations; ++j )
            {
                const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservationResidual =
                        getResidualValue( observationIdsBySet_.at( setId ).at( j ) );
                bool removeObservation = false;
                for( unsigned int k = 0; k < observableSize; ++k )
                {
                    if( ( !useOppositeCondition && ( std::fabs( singleObservationResidual[ k ] ) > residualCutOff[ k ] ) ) ||
                        ( useOppositeCondition && ( std::fabs( singleObservationResidual[ k ] ) <= residualCutOff[ k ] ) ) )
                    {
                        removeObservation = true;
                    }
                }
                if( removeObservation )
                {
                    indicesToRemove.push_back( j );
                }
            }
            break;
        }
        case absolute_value_filtering: {
            Eigen::VectorXd absoluteValueCutOff = Eigen::VectorXd::Zero( observableSize );
            if( std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter ) != nullptr )
            {
                absoluteValueCutOff = std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter )->getFilterValue( ) *
                        Eigen::VectorXd::Ones( observableSize );
            }
            else if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter ) != nullptr )
            {
                if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( ).size( ) !=
                    static_cast< int >( observableSize ) )
                {
                    throw std::runtime_error(
                            "Error when performing observation value filtering, size of the filter value inconsistent with observable "
                            "size." );
                }
                absoluteValueCutOff =
                        std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( );
            }

            for( unsigned int j = 0; j < numberOfObservations; ++j )
            {
                const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservation =
                        getObservationValue( observationIdsBySet_.at( setId ).at( j ) );
                bool removeObservation = false;
                for( unsigned int k = 0; k < observableSize; ++k )
                {
                    if( ( !useOppositeCondition && ( singleObservation[ k ] > absoluteValueCutOff[ k ] ) ) ||
                        ( useOppositeCondition && ( singleObservation[ k ] <= absoluteValueCutOff[ k ] ) ) )
                    {
                        removeObservation = true;
                    }
                }
                if( removeObservation )
                {
                    indicesToRemove.push_back( j );
                }
            }
            break;
        }
        case epochs_filtering: {
            const std::vector< double > filterEpochs =
                    std::dynamic_pointer_cast< ObservationFilter< std::vector< double > > >( observationFilter )->getFilterValue( );
            for( unsigned int j = 0; j < numberOfObservations; ++j )
            {
                const TimeType singleObservationTime = getObservationTime( observationIdsBySet_.at( setId ).at( j ) );
                if( ( !useOppositeCondition && ( std::count( filterEpochs.begin( ), filterEpochs.end( ), singleObservationTime ) > 0 ) ) ||
                    ( useOppositeCondition && ( std::count( filterEpochs.begin( ), filterEpochs.end( ), singleObservationTime ) == 0 ) ) )
                {
                    indicesToRemove.push_back( j );
                }
            }
            break;
        }
        case time_bounds_filtering: {
            const std::pair< double, double > timeBounds =
                    std::dynamic_pointer_cast< ObservationFilter< std::pair< double, double > > >( observationFilter )->getFilterValue( );
            for( unsigned int j = 0; j < numberOfObservations; ++j )
            {
                const TimeType singleObservationTime = getObservationTime( observationIdsBySet_.at( setId ).at( j ) );
                if( ( !useOppositeCondition &&
                      ( ( singleObservationTime >= timeBounds.first ) && ( singleObservationTime <= timeBounds.second ) ) ) ||
                    ( useOppositeCondition &&
                      ( ( singleObservationTime < timeBounds.first ) || ( singleObservationTime > timeBounds.second ) ) ) )
                {
                    indicesToRemove.push_back( j );
                }
            }
            break;
        }
        case dependent_variable_filtering: {
            if( std::dynamic_pointer_cast< ObservationDependentVariableFilter >( observationFilter ) == nullptr )
            {
                throw std::runtime_error(
                        "Error when performing dependent variable filtering, inconsistent filter input (should be "
                        "ObservationDependentVariableFilter object)." );
            }

            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > settings =
                    std::dynamic_pointer_cast< ObservationDependentVariableFilter >( observationFilter )->getDependentVariableSettings( );
            const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
            const unsigned int dependentVariableSize =
                    getObservationDependentVariableSize( settings, getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_ );

            const Eigen::VectorXd dependentVariableCutOff =
                    std::dynamic_pointer_cast< ObservationDependentVariableFilter >( observationFilter )->getFilterValue( );
            if( dependentVariableCutOff.size( ) != static_cast< int >( dependentVariableSize ) )
            {
                throw std::runtime_error(
                        "Error when performing dependent variable filtering, size of the dependent variable cut off vector "
                        "inconsistent with dependent variable size." );
            }

            const Eigen::MatrixXd singleDependentVariableValues = getSingleDependentVariableForSet( setId, settings );
            if( ( singleDependentVariableValues.rows( ) != numberOfObservations ) ||
                ( singleDependentVariableValues.cols( ) != dependentVariableSize ) )
            {
                throw std::runtime_error(
                        "Error when performing dependent variable filtering, size of observation dependent variables is inconsistent "
                        "with the number of observations and presupposed dependent variable size." );
            }

            for( unsigned int j = 0; j < numberOfObservations; ++j )
            {
                bool removeObservation = false;
                for( unsigned int k = 0; k < dependentVariableSize; ++k )
                {
                    if( ( !useOppositeCondition && ( singleDependentVariableValues( j, k ) ) > dependentVariableCutOff[ k ] ) ||
                        ( useOppositeCondition && ( singleDependentVariableValues( j, k ) ) <= dependentVariableCutOff[ k ] ) )
                    {
                        removeObservation = true;
                    }
                }
                if( removeObservation )
                {
                    indicesToRemove.push_back( j );
                }
            }
            break;
        }
        default:
            throw std::runtime_error( "Observation filter type not recognised." );
    }

    return indicesToRemove;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::moveObservationsToSet(
        const unsigned int sourceSetId,
        ObservationDataset< ObservationScalarType, TimeType >& targetDataset,
        const unsigned int targetSetId,
        const std::vector< unsigned int >& indices,
        const bool removeFromSource )
{
    if( this == &targetDataset && sourceSetId == targetSetId )
    {
        throw std::runtime_error( "Error when moving observations in dataset, source and target sets are identical." );
    }
    const ObservationSetMetadata< ObservationScalarType, TimeType >& sourceMetadata = getObservationSetMetadata( sourceSetId );
    const ObservationSetMetadata< ObservationScalarType, TimeType >& targetMetadata =
            targetDataset.getObservationSetMetadata( targetSetId );
    if( sourceMetadata.observableType_ != targetMetadata.observableType_ ||
        sourceMetadata.referenceLinkEnd_ != targetMetadata.referenceLinkEnd_ ||
        !( getLinkDefinition( sourceMetadata.linkDefinitionId_ ) == targetDataset.getLinkDefinition( targetMetadata.linkDefinitionId_ ) ) )
    {
        throw std::runtime_error( "Error when moving observations in dataset, source and target set metadata are incompatible." );
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
    std::vector< TimeType > times;
    std::vector< Eigen::VectorXd > dependentVariables;
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;
    std::vector< unsigned int > sourceObservationIds;
    std::vector< std::size_t > selectedSourceScalarRows;
    const FlattenedObservationData< ObservationScalarType, TimeType > sourceFlattenedData =
            createComputationFlattenedObservationData( true );

    const bool hasDependentVariables = !getDependentVariablesForSet( sourceSetId ).empty( );
    for( const unsigned int index : indices )
    {
        if( index >= getNumberOfObservationsForSet( sourceSetId ) )
        {
            throw std::runtime_error( "Error when moving observations in dataset, index is out of bounds." );
        }
        const unsigned int observationId = observationIdsBySet_.at( sourceSetId ).at( index );
        sourceObservationIds.push_back( observationId );
        observations.push_back( getObservationValue( observationId ) );
        times.push_back( getObservationTime( observationId ) );
        residuals.push_back( getResidualValue( observationId ) );
        for( unsigned int componentIndex = 0; componentIndex < sourceMetadata.observableSize_; ++componentIndex )
        {
            selectedSourceScalarRows.push_back( sourceFlattenedData.getFlattenedRow( observationId, componentIndex ) );
        }
        if( hasDependentVariables )
        {
            dependentVariables.push_back( getDependentVariables( observationId ) );
        }
    }

    const Eigen::MatrixXd selectedWeightMatrix =
            selectSubmatrix( sourceFlattenedData.getSparseWeightMatrix( ).toDense( ), selectedSourceScalarRows, selectedSourceScalarRows );
    std::vector< TimeType > combinedTimes = targetDataset.getObservationTimesForSet( targetSetId );
    std::vector< int > sourceIndicesAfterSorting( combinedTimes.size( ), -1 );
    combinedTimes.insert( combinedTimes.end( ), times.begin( ), times.end( ) );
    for( std::size_t i = 0; i < times.size( ); ++i )
    {
        sourceIndicesAfterSorting.push_back( static_cast< int >( i ) );
    }
    const std::vector< std::size_t > sortingPermutation = getTimeSortingPermutation( combinedTimes );
    reorderVector( sourceIndicesAfterSorting, sortingPermutation );

    targetDataset.addObservationsToSet( targetSetId, observations, times, dependentVariables, {}, residuals, true );
    const std::vector< unsigned int >& targetObservationIds = targetDataset.getObservationIdsForSet( targetSetId );
    std::vector< unsigned int > addedTargetObservationIds( observations.size( ) );
    for( std::size_t i = 0; i < sourceIndicesAfterSorting.size( ); ++i )
    {
        if( sourceIndicesAfterSorting.at( i ) >= 0 )
        {
            addedTargetObservationIds.at( sourceIndicesAfterSorting.at( i ) ) = targetObservationIds.at( i );
        }
    }

    bool hasCrossObservationWeight = false;
    for( std::size_t i = 0; i < sourceObservationIds.size( ); ++i )
    {
        const unsigned int targetObservationId = addedTargetObservationIds.at( i );
        targetDataset.copyObservationStateAndWeightFrom( *this, sourceObservationIds.at( i ), targetObservationId );
        targetDataset.setWeightMatrixForObservation( targetObservationId,
                                                     selectedWeightMatrix.block( i * sourceMetadata.observableSize_,
                                                                                 i * sourceMetadata.observableSize_,
                                                                                 sourceMetadata.observableSize_,
                                                                                 sourceMetadata.observableSize_ ) );
        for( std::size_t j = 0; j < sourceObservationIds.size( ); ++j )
        {
            if( i != j &&
                !selectedWeightMatrix
                         .block( i * sourceMetadata.observableSize_,
                                 j * sourceMetadata.observableSize_,
                                 sourceMetadata.observableSize_,
                                 sourceMetadata.observableSize_ )
                         .isZero( 0.0 ) )
            {
                hasCrossObservationWeight = true;
            }
        }
    }
    if( hasCrossObservationWeight )
    {
        targetDataset.setWeightBlock( addedTargetObservationIds, addedTargetObservationIds, selectedWeightMatrix );
    }
    if( removeFromSource )
    {
        removeObservationsFromSet( sourceSetId, indices );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::eraseDuplicateObservationsFromSet( const unsigned int setId,
                                                                                                      const bool printWarning )
{
    const std::vector< TimeType > observationTimes = getObservationTimesForSet( setId );
    std::set< TimeType > retainedTimes;
    std::vector< unsigned int > indicesToRemove;
    for( unsigned int i = 0; i < observationTimes.size( ); ++i )
    {
        if( !retainedTimes.insert( observationTimes.at( i ) ).second )
        {
            indicesToRemove.push_back( i );
        }
    }

    if( !indicesToRemove.empty( ) )
    {
        const std::size_t beforeCount = getNumberOfObservationsForSet( setId );
        removeObservationsFromSet( setId, indicesToRemove );
        if( printWarning )
        {
            std::cerr << "[WARNING] Detected and removed " << beforeCount - getNumberOfObservationsForSet( setId )
                      << "duplicate observations when creating observation dataset" << std::endl;
        }
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getComputedObservationVectorForSet( const unsigned int setId ) const
{
    return getObservationVectorForSet( setId ) - getResidualVectorForSet( setId );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationVectorForSet( const unsigned int setId ) const
{
    return createSetVector( setId, observedValues_ );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::VectorXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getWeightVectorForSet( const unsigned int setId ) const
{
    Eigen::VectorXd weights = Eigen::VectorXd::Zero( getTotalScalarSizeForSet( setId ) );
    std::size_t currentIndex = 0;
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
        weights.segment( currentIndex, row.scalarSize_ ) = observationWeights_.getObservationWeightVector( observationId, row.scalarSize_ );
        currentIndex += row.scalarSize_;
    }
    if( observationWeights_.hasSetWeightBlock( setId ) )
    {
        weights = observationWeights_.getSetWeightBlock( setId ).diagonal( );
    }
    return weights;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getResidualVectorForSet( const unsigned int setId ) const
{
    return createSetVector( setId, residualValues_ );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationSetIds(
        const std::shared_ptr< ObservationCollectionParser >& observationParser ) const
{
    std::vector< unsigned int > setIds;
    for( const unsigned int setId : getSetIdsInOrderedFlattenedDataOrder( ) )
    {
        if( isObservationSetSelectedByLegacyParser( *this, setId, observationParser ) )
        {
            setIds.push_back( setId );
        }
    }
    return setIds;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< std::pair< int, int > > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationSetStartAndSize( ) const
{
    std::vector< std::pair< int, int > > startAndSize;
    startAndSize.reserve( getNumberOfObservationSets( ) );

    int currentIndex = 0;
    for( const unsigned int setId : getSetIdsInOrderedFlattenedDataOrder( ) )
    {
        const int currentSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
        startAndSize.push_back( std::make_pair( currentIndex, currentSize ) );
        currentIndex += currentSize;
    }
    return startAndSize;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationSetStartAndSizeByLink( ) const
{
    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > startAndSizeByLink;

    int currentIndex = 0;
    for( const unsigned int setId : getSetIdsInOrderedFlattenedDataOrder( ) )
    {
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
        const LinkEnds linkEnds = getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        const int currentSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
        startAndSizeByLink[ metadata.observableType_ ][ linkEnds ].push_back( std::make_pair( currentIndex, currentSize ) );
        currentIndex += currentSize;
    }
    return startAndSizeByLink;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::map< ObservableType, std::map< LinkEnds, std::pair< int, int > > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationTypeAndLinkEndStartAndSize( ) const
{
    std::map< ObservableType, std::map< LinkEnds, std::pair< int, int > > > startAndSize;

    int currentIndex = 0;
    for( const unsigned int setId : getSetIdsInOrderedFlattenedDataOrder( ) )
    {
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
        const LinkEnds linkEnds = getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        const int currentSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
        if( startAndSize[ metadata.observableType_ ].count( linkEnds ) == 0 )
        {
            startAndSize[ metadata.observableType_ ][ linkEnds ] = std::make_pair( currentIndex, currentSize );
        }
        else
        {
            startAndSize[ metadata.observableType_ ][ linkEnds ].second += currentSize;
        }
        currentIndex += currentSize;
    }
    return startAndSize;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::map< ObservableType, std::pair< int, int > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservableTypeStartAndSize( ) const
{
    std::map< ObservableType, std::pair< int, int > > startAndSize;
    int currentIndex = 0;

    for( const unsigned int setId : getSetIdsInOrderedFlattenedDataOrder( ) )
    {
        const ObservableType observableType = getObservationSetMetadata( setId ).observableType_;
        const int currentSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
        if( startAndSize.count( observableType ) == 0 )
        {
            startAndSize[ observableType ] = std::make_pair( currentIndex, currentSize );
        }
        else
        {
            startAndSize[ observableType ].second += currentSize;
        }
        currentIndex += currentSize;
    }
    return startAndSize;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::map< ObservableType, std::vector< LinkEnds > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getLinkEndsPerObservableType( ) const
{
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservableType;
    for( const unsigned int setId : getSetIdsInOrderedFlattenedDataOrder( ) )
    {
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
        const LinkEnds linkEnds = getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        std::vector< LinkEnds >& linkEndsForType = linkEndsPerObservableType[ metadata.observableType_ ];
        if( std::count( linkEndsForType.begin( ), linkEndsForType.end( ), linkEnds ) == 0 )
        {
            linkEndsForType.push_back( linkEnds );
        }
    }
    return linkEndsPerObservableType;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationVectorForObservableType(
        const ObservableType observableType ) const
{
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observations =
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( getTotalScalarSizeForObservableType( observableType ) );

    int currentIndex = 0;
    for( const unsigned int setId : getObservationSetIdsForObservableType( observableType ) )
    {
        const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > setObservations = getObservationVectorForSet( setId );
        observations.segment( currentIndex, setObservations.size( ) ) = setObservations;
        currentIndex += setObservations.size( );
    }
    return observations;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setObservationVectorForObservableType(
        const ObservableType observableType,
        const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observations )
{
    if( observations.size( ) != static_cast< int >( getTotalScalarSizeForObservableType( observableType ) ) )
    {
        throw std::runtime_error( "Error when setting observable-type observation vector, input size is inconsistent." );
    }

    int currentIndex = 0;
    for( const unsigned int setId : getObservationSetIdsForObservableType( observableType ) )
    {
        const int setSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
        setObservationVectorForSet( setId, observations.segment( currentIndex, setSize ) );
        currentIndex += setSize;
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightVectorForObservableType( const ObservableType observableType,
                                                                                                     const Eigen::VectorXd& weights )
{
    setTabulatedWeights( weights, observationParser( observableType ) );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setResidualVector(
        const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualVector )
{
    if( residualVector.size( ) != static_cast< int >( getTotalScalarSize( ) ) )
    {
        throw std::runtime_error(
                "Error when setting dataset residual vector, input size is inconsistent with total scalar observation size." );
    }

    int currentIndex = 0;
    for( unsigned int setId = 0; setId < getNumberOfObservationSets( ); ++setId )
    {
        const int currentSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > setResiduals = residualVector.segment( currentIndex, currentSize );
        setResidualVectorForSet( setId, setResiduals );
        currentIndex += currentSize;
    }
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETLEGACYIMPLEMENTATION_H
