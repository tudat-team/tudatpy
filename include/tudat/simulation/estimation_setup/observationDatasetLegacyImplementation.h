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

namespace observation_legacy_detail
{

//! Preserve legacy broadcasting of either one vector per equal-sized set or one concatenated vector.
template< typename Sets, typename SizeFunction, typename WeightSetter >
void setTabulatedWeights( const Sets& sets, const Eigen::VectorXd& weights, SizeFunction size, WeightSetter assign )
{
    if( sets.empty( ) )
    {
        return;
    }
    const int firstSize = size( sets.front( ) );
    int totalSize = 0;
    bool equalSizes = true;
    for( const auto& set : sets )
    {
        const int currentSize = size( set );
        totalSize += currentSize;
        equalSizes = equalSizes && currentSize == firstSize;
    }
    const bool concatenated = weights.size( ) == totalSize;
    if( !concatenated && !( equalSizes && weights.size( ) == firstSize ) )
    {
        throw std::runtime_error( "Tabulated weights must match each equal-sized observation set or all selected sets combined." );
    }
    int offset = 0;
    for( const auto& set : sets )
    {
        const int currentSize = size( set );
        assign( set, weights.segment( concatenated ? offset : 0, currentSize ) );
        offset += currentSize;
    }
}

}  // namespace observation_legacy_detail

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
    completeWeightData_.reset( );
    observationWeights_.setDiagonal( getScalarComponentIdsForObservationSelection( observationIdsBySet_.at( setId ), {} ), weightVector );
    auto& structure = setMetadata_.at( setId ).weightStructure_;
    if( structure != ObservationWeightStructure::inter_set_weights )
    {
        structure = ObservationWeightStructure::diagonal;
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

    observation_legacy_detail::setTabulatedWeights(
            setIds,
            tabulatedWeights,
            [ this ]( unsigned int setId ) { return getTotalScalarSizeForSet( setId ); },
            [ this ]( unsigned int setId, const Eigen::VectorXd& weights ) { setWeightVectorForSet( setId, weights ); } );
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
    return observationWeights_.getDiagonal( getScalarComponentIdsForObservationSelection( observationIdsBySet_.at( setId ), {} ) );
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
    for( const unsigned int setId : getSetIdsInObservationVectorOrder( ) )
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
    for( const unsigned int setId : getSetIdsInObservationVectorOrder( ) )
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
    for( const unsigned int setId : getSetIdsInObservationVectorOrder( ) )
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
std::map< ObservableType, std::pair< int, int > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservableTypeStartAndSize( ) const
{
    std::map< ObservableType, std::pair< int, int > > startAndSize;
    int currentIndex = 0;

    for( const unsigned int setId : getSetIdsInObservationVectorOrder( ) )
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
