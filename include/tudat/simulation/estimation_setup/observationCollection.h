/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_COLLECTION_H
#define TUDAT_OBSERVATION_COLLECTION_H

#include <Eigen/Core>
#include <functional>
#include <memory>
#include <vector>

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/simulation/estimation_setup/observationsProcessing.h"
#include "tudat/simulation/estimation_setup/singleObservationSet.h"

namespace tudat
{

namespace observation_models
{

using namespace simulation_setup;


template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if<
                  is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value,
                  int >::type = 0 >
std::map< ObservableType,
          std::map< LinkEnds,
                    std::vector< std::shared_ptr<
                            SingleObservationSet< ObservationScalarType, TimeType > > > > >
createSortedObservationSetList(
        const std::vector<
                std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                observationSetList )
{
    std::map< ObservableType,
              std::map< LinkEnds,
                        std::vector< std::shared_ptr<
                                SingleObservationSet< ObservationScalarType, TimeType > > > > >
            sortedObservations;
    for( unsigned int i = 0; i < observationSetList.size( ); i++ )
    {
        sortedObservations[ observationSetList.at( i )->getObservableType( ) ]
                          [ observationSetList.at( i )->getLinkEnds( ).linkEnds_ ]
                                  .push_back( observationSetList.at( i ) );
    }
    return sortedObservations;
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if<
                  is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value,
                  int >::type = 0 >
class ObservationCollection
{
   public:
    typedef std::map<
            ObservableType,
            std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > >
            SortedObservationSets;

    ObservationCollection(
            const SortedObservationSets& observationSetList = SortedObservationSets( ) ) :
        observationSetList_( observationSetList )
    {
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    ObservationCollection(
            const std::vector<
                    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                    observationSetList ) :
        observationSetList_( createSortedObservationSetList< ObservationScalarType, TimeType >(
                observationSetList ) )
    {
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationVector( )
    {
        return concatenatedObservations_;
    }

    const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >&
    getObservationVectorReference( )
    {
        return concatenatedObservations_;
    }

    void setObservations(
            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& newObservations )
    {
        if( newObservations.size( ) != totalObservableSize_ )
        {
            throw std::runtime_error(
                    "Error when resetting observations in ObservationCollection, size of input "
                    "observation vector is inconsistent." );
        }

        unsigned int startIndexObsSet = 0;
        for( auto observableIt: observationSetList_ )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto set: linkEndsIt.second )
                {
                    unsigned int sizeCurrentSet = set->getTotalObservationSetSize( );
                    set->setObservations(
                            newObservations.segment( startIndexObsSet, sizeCurrentSet ) );
                    startIndexObsSet += sizeCurrentSet;
                }
            }
        }
    }

    void setResiduals(
            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& newResiduals )
    {
        if( newResiduals.size( ) != totalObservableSize_ )
        {
            throw std::runtime_error(
                    "Error when resetting observations in ObservationCollection, size of input "
                    "observation vector is inconsistent." );
        }

        unsigned int startIndexObsSet = 0;
        for( auto observableIt: observationSetList_ )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto set: linkEndsIt.second )
                {
                    unsigned int sizeCurrentSet = set->getTotalObservationSetSize( );
                    set->setResiduals( newResiduals.segment( startIndexObsSet, sizeCurrentSet ) );
                    startIndexObsSet += sizeCurrentSet;
                }
            }
        }
    }

    void setObservations(
            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observations,
            const std::shared_ptr< ObservationCollectionParser > observationParser )
    {
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObsSets = getSingleObservationSets( observationParser );
        if( singleObsSets.empty( ) )
        {
            std::cerr << "Warning when setting observations, no single observation set found for "
                         "specified observation parser. Weights not set";
        }
        int totalSizeAllObsSets = 0;
        for( unsigned int k = 0; k < singleObsSets.size( ); k++ )
        {
            totalSizeAllObsSets += singleObsSets.at( k )->getTotalObservationSetSize( );
        }

        if( observations.size( ) == totalSizeAllObsSets )
        {
            unsigned int startObsSet = 0;
            for( auto obsSet: singleObsSets )
            {
                obsSet->setObservations( observations.segment(
                        startObsSet, obsSet->getTotalObservationSetSize( ) ) );
                startObsSet += obsSet->getTotalObservationSetSize( );
            }
        }
        else
        {
            throw std::runtime_error(
                    "Error when setting observations, the size of the input observation vector "
                    "should be consistent with "
                    "the combined size of all required observation sets." );
        }
    }

    void setObservations(
            const std::map< std::shared_ptr< ObservationCollectionParser >,
                            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
                    observationsPerParser )
    {
        for( auto parserIt: observationsPerParser )
        {
            setObservations( parserIt.second, parserIt.first );
        }
    }

    void setResiduals( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals,
                       const std::shared_ptr< ObservationCollectionParser > observationParser )
    {
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObsSets = getSingleObservationSets( observationParser );
        if( singleObsSets.empty( ) )
        {
            std::cerr << "Warning when setting residuals, no single observation set found for "
                         "specified observation parser. Weights not set";
        }
        int totalSizeAllObsSets = 0;
        for( unsigned int k = 0; k < singleObsSets.size( ); k++ )
        {
            totalSizeAllObsSets += singleObsSets.at( k )->getTotalObservationSetSize( );
        }

        if( residuals.size( ) == totalSizeAllObsSets )
        {
            unsigned int startObsSet = 0;
            for( auto obsSet: singleObsSets )
            {
                obsSet->setResiduals(
                        residuals.segment( startObsSet, obsSet->getTotalObservationSetSize( ) ) );
                startObsSet += obsSet->getTotalObservationSetSize( );
            }
        }
        else
        {
            throw std::runtime_error(
                    "Error when setting residuals, the size of the input residual vector should be "
                    "consistent with "
                    "the combined size of all required observation sets." );
        }
    }

    void setResiduals( const std::map< std::shared_ptr< ObservationCollectionParser >,
                                       Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
                               residualsPerParser )
    {
        for( auto parserIt: residualsPerParser )
        {
            setResiduals( parserIt.second, parserIt.first );
        }
    }

    std::vector< TimeType > getConcatenatedTimeVector( )
    {
        return concatenatedTimes_;
    }

    std::vector< double > getConcatenatedDoubleTimeVector( )
    {
        return utilities::staticCastVector< double, TimeType >( concatenatedTimes_ );
    }

    std::pair< TimeType, TimeType > getTimeBounds( )
    {
        return std::make_pair(
                *std::min_element( concatenatedTimes_.begin( ), concatenatedTimes_.end( ) ),
                *std::max_element( concatenatedTimes_.begin( ), concatenatedTimes_.end( ) ) );
    }

    std::vector< std::pair< TimeType, TimeType > > getTimeBoundsPerSet(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< std::pair< TimeType, TimeType > > timeBounds;
        for( auto observableIt: getSingleObservationSetsIndices( observationParser ) )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto index: linkEndsIt.second )
                {
                    timeBounds.push_back( observationSetList_.at( observableIt.first )
                                                  .at( linkEndsIt.first )
                                                  .at( index )
                                                  ->getTimeBounds( ) );
                }
            }
        }
        return timeBounds;
    }

    std::vector< int > getConcatenatedLinkEndIds( )
    {
        return concatenatedLinkEndIds_;
    }

    std::vector< int > getConcatenatedLinkEndIds(
            std::shared_ptr< ObservationCollectionParser > observationParser )
    {
        std::vector< int > subsetConcatenatedLinkEndIds;
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );

        for( auto obsIt: observationSetsIndices )
        {
            for( auto linkEndsIt: obsIt.second )
            {
                std::vector< unsigned int > listRelevantSetIndices = linkEndsIt.second;
                for( auto index: listRelevantSetIndices )
                {
                    std::pair< int, int > singleSetStartAndSize =
                            observationSetStartAndSize_.at( obsIt.first )
                                    .at( linkEndsIt.first )
                                    .at( index );
                    for( unsigned int k = singleSetStartAndSize.first;
                         k < ( singleSetStartAndSize.first + singleSetStartAndSize.second );
                         k++ )
                    {
                        subsetConcatenatedLinkEndIds.push_back( concatenatedLinkEndIds_[ k ] );
                    }
                }
            }
        }
        return subsetConcatenatedLinkEndIds;
    }

    std::map< observation_models::LinkEnds, int > getLinkEndIdentifierMap( )
    {
        return linkEndIds_;
    }

    std::map< int, observation_models::LinkEnds > getInverseLinkEndIdentifierMap( )
    {
        return inverseLinkEndIds_;
    }

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > >
    getObservationSetStartAndSize( )
    {
        return observationSetStartAndSize_;
    }

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > >&
    getObservationSetStartAndSizeReference( )
    {
        return observationSetStartAndSize_;
    }

    std::vector< std::pair< int, int > > getConcatenatedObservationSetStartAndSize( )
    {
        return concatenatedObservationSetStartAndSize_;
    }

    std::map< ObservableType, std::map< LinkEnds, std::pair< int, int > > >
    getObservationTypeAndLinkEndStartAndSize( )
    {
        return observationTypeAndLinkEndStartAndSize_;
    }

    std::map< ObservableType, std::map< int, std::vector< std::pair< int, int > > > >
    getObservationSetStartAndSizePerLinkEndIndex( )
    {
        return observationSetStartAndSizePerLinkEndIndex_;
    }

    std::map< ObservableType, std::pair< int, int > > getObservationTypeStartAndSize( )
    {
        return observationTypeStartAndSize_;
    }

    int getTotalObservableSize( )
    {
        return totalObservableSize_;
    }

    SortedObservationSets getObservationsSets( )
    {
        return observationSetList_;
    }

    const SortedObservationSets& getObservationsReference( ) const
    {
        return observationSetList_;
    }

    std::vector< LinkEnds > getConcatenatedLinkEndIdNames( )
    {
        return concatenatedLinkEndIdNames_;
    }

    std::map< ObservableType, std::vector< LinkDefinition > > getLinkDefinitionsPerObservable( )
    {
        return linkDefinitionsPerObservable_;
    }

    std::vector< LinkDefinition > getLinkDefinitionsForSingleObservable(
            const ObservableType observableType )
    {
        if( linkDefinitionsPerObservable_.count( observableType ) > 0 )
        {
            return linkDefinitionsPerObservable_.at( observableType );
        }
        else
        {
            return std::vector< LinkDefinition >( );
        }
    }

    std::map< ObservableType,
              std::map< int,
                        std::vector< std::shared_ptr<
                                SingleObservationSet< ObservationScalarType, TimeType > > > > >
    getSortedObservationSets( )
    {
        std::map< ObservableType,
                  std::map< int,
                            std::vector< std::shared_ptr<
                                    SingleObservationSet< ObservationScalarType, TimeType > > > > >
                observationSetListIndexSorted;
        for( auto it1: observationSetList_ )
        {
            for( auto it2: it1.second )
            {
                observationSetListIndexSorted[ it1.first ][ linkEndIds_[ it2.first ] ] = it2.second;
            }
        }
        return observationSetListIndexSorted;
    }

    std::map< ObservableType, std::map< int, std::vector< std::pair< double, double > > > >
    getSortedObservationSetsTimeBounds( )
    {
        std::map< ObservableType, std::map< int, std::vector< std::pair< double, double > > > >
                observationSetTimeBounds;
        for( auto it1: observationSetList_ )
        {
            for( auto it2: it1.second )
            {
                int currentLinkEndId = linkEndIds_[ it2.first ];
                for( unsigned int i = 0; i < it2.second.size( ); i++ )
                {
                    std::pair< TimeType, TimeType > timeBounds =
                            it2.second.at( i )->getTimeBounds( );
                    observationSetTimeBounds[ it1.first ][ currentLinkEndId ].push_back(
                            { static_cast< double >( timeBounds.first ),
                              static_cast< double >( timeBounds.second ) } );
                }
            }
        }
        return observationSetTimeBounds;
    }

    std::map< ObservableType, std::vector< LinkEnds > > getLinkEndsPerObservableType( )
    {
        std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservableType;

        for( auto observableTypeIt = observationSetList_.begin( );
             observableTypeIt != observationSetList_.end( );
             ++observableTypeIt )
        {
            std::vector< LinkEnds > linkEndsVector;

            for( auto linkEndsIt = observableTypeIt->second.begin( );
                 linkEndsIt != observableTypeIt->second.end( );
                 ++linkEndsIt )
            {
                linkEndsVector.push_back( linkEndsIt->first );
            }

            linkEndsPerObservableType[ observableTypeIt->first ] = linkEndsVector;
        }

        return linkEndsPerObservableType;
    }

    std::vector< ObservableType > getObservableTypes(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        std::vector< ObservableType > observableTypes;
        for( auto observableIt: observationSetsIndices )
        {
            if( std::count( observableTypes.begin( ),
                            observableTypes.end( ),
                            observableIt.first ) == 0 )
            {
                observableTypes.push_back( observableIt.first );
            }
        }
        return observableTypes;
    }

    std::vector< LinkEnds > getLinkEnds(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        std::vector< LinkEnds > linkEnds;
        for( auto observableIt: observationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                if( std::count( linkEnds.begin( ), linkEnds.end( ), linkEndsIt.first ) == 0 )
                {
                    linkEnds.push_back( linkEndsIt.first );
                }
            }
        }
        return linkEnds;
    }

    std::vector< std::string > getBodiesInLinkEnds(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        std::vector< std::string > bodiesInLinkEnds;
        for( auto observableIt: observationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto it: linkEndsIt.first )
                {
                    if( std::count( bodiesInLinkEnds.begin( ),
                                    bodiesInLinkEnds.end( ),
                                    it.second.bodyName_ ) == 0 )
                    {
                        bodiesInLinkEnds.push_back( it.second.bodyName_ );
                    }
                }
            }
        }
        return bodiesInLinkEnds;
    }

    std::vector< std::string > getReferencePointsInLinkEnds(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        std::vector< std::string > referencePoints;
        for( auto observableIt: observationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto it: linkEndsIt.first )
                {
                    if( ( it.second.stationName_ != "" ) &&
                        ( std::count( referencePoints.begin( ),
                                      referencePoints.end( ),
                                      it.second.stationName_ ) == 0 ) )
                    {
                        referencePoints.push_back( it.second.stationName_ );
                    }
                }
            }
        }
        return referencePoints;
    }

    std::vector< std::pair< TimeType, TimeType > > getTimeBoundsList(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        std::vector< std::pair< TimeType, TimeType > > timeBounds;
        for( auto observableIt: observationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto index: linkEndsIt.second )
                {
                    std::pair< TimeType, TimeType > currentTimeBounds =
                            observationSetList_.at( observableIt.first )
                                    .at( linkEndsIt.first )
                                    .at( index )
                                    ->getTimeBounds( );
                    if( std::count( timeBounds.begin( ), timeBounds.end( ), currentTimeBounds ) ==
                        0 )
                    {
                        timeBounds.push_back( currentTimeBounds );
                    }
                }
            }
        }
        return timeBounds;
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservations(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for( auto observableIt: observationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto index: linkEndsIt.second )
                {
                    observations.push_back( observationSetList_.at( observableIt.first )
                                                    .at( linkEndsIt.first )
                                                    .at( index )
                                                    ->getObservationsVector( ) );
                }
            }
        }
        return observations;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getConcatenatedObservations(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations =
                getObservations( observationParser );

        unsigned int totalObsSize = 0;
        for( auto obs: observations )
        {
            totalObsSize += obs.size( );
        }

        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > concatenatedObservations =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObsSize, 1 );
        unsigned int obsIndex = 0;
        for( unsigned int k = 0; k < observations.size( ); k++ )
        {
            concatenatedObservations.block( obsIndex, 0, observations.at( k ).size( ), 1 ) =
                    observations.at( k );
            obsIndex += observations.at( k ).size( );
        }

        return concatenatedObservations;
    }

    std::vector< std::vector< TimeType > > getObservationTimes(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< std::vector< TimeType > > observationsTimes;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for( auto observableIt: observationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto index: linkEndsIt.second )
                {
                    observationsTimes.push_back( observationSetList_.at( observableIt.first )
                                                         .at( linkEndsIt.first )
                                                         .at( index )
                                                         ->getObservationTimes( ) );
                }
            }
        }
        return observationsTimes;
    }

    std::vector< TimeType > getConcatenatedObservationTimes(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< std::vector< TimeType > > observationTimes =
                getObservationTimes( observationParser );

        std::vector< TimeType > concatenatedObservationsTimes;
        for( auto times: observationTimes )
        {
            for( unsigned int k = 0; k < times.size( ); k++ )
            {
                concatenatedObservationsTimes.push_back( times.at( k ) );
            }
        }

        return concatenatedObservationsTimes;
    }

    std::vector< double > getConcatenatedDoubleObservationTimes(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        return utilities::staticCastVector< double, TimeType >(
                getConcatenatedObservationTimes( observationParser ) );
    }

    std::pair< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >,
               std::vector< std::vector< TimeType > > >
    getObservationsAndTimes( std::shared_ptr< ObservationCollectionParser > observationParser =
                                     std::make_shared< ObservationCollectionParser >( ) )
    {
        return std::make_pair( getObservations( observationParser ),
                               getObservationTimes( observationParser ) );
    }

    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > >
    getConcatenatedObservationsAndTimes(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        return std::make_pair( getConcatenatedObservations( observationParser ),
                               getConcatenatedObservationTimes( observationParser ) );
    }

    std::vector< Eigen::VectorXd > getWeights(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::VectorXd > weights;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for( auto observableIt: observationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto index: linkEndsIt.second )
                {
                    weights.push_back( observationSetList_.at( observableIt.first )
                                               .at( linkEndsIt.first )
                                               .at( index )
                                               ->getWeightsVector( ) );
                }
            }
        }
        return weights;
    }

    Eigen::VectorXd getConcatenatedWeights(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::VectorXd > weights = getWeights( observationParser );

        unsigned int totalObsSize = 0;
        for( auto weight: weights )
        {
            totalObsSize += weight.size( );
        }

        Eigen::VectorXd concatenatedWeights = Eigen::VectorXd::Zero( totalObsSize );
        unsigned int obsIndex = 0;
        for( unsigned int k = 0; k < weights.size( ); k++ )
        {
            concatenatedWeights.block( obsIndex, 0, weights.at( k ).size( ), 1 ) = weights.at( k );
            obsIndex += weights.at( k ).size( );
        }

        return concatenatedWeights;
    }

    Eigen::VectorXd getUnparsedConcatenatedWeights( ) const
    {
        return getConcatenatedWeights( );
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getResiduals(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for( auto observableIt: observationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto index: linkEndsIt.second )
                {
                    residuals.push_back( observationSetList_.at( observableIt.first )
                                                 .at( linkEndsIt.first )
                                                 .at( index )
                                                 ->getResidualsVector( ) );
                }
            }
        }
        return residuals;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getConcatenatedResiduals(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals =
                getResiduals( observationParser );

        unsigned int totalObsSize = 0;
        for( auto residual: residuals )
        {
            totalObsSize += residual.size( );
        }

        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > concatenatedResiduals =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObsSize, 1 );
        unsigned int obsIndex = 0;
        for( unsigned int k = 0; k < residuals.size( ); k++ )
        {
            concatenatedResiduals.block( obsIndex, 0, residuals.at( k ).size( ), 1 ) =
                    residuals.at( k );
            obsIndex += residuals.at( k ).size( );
        }

        return concatenatedResiduals;
    }

    std::vector< Eigen::VectorXd > getRmsResiduals(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::VectorXd > rmsResiduals;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for( auto observableIt: observationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto index: linkEndsIt.second )
                {
                    rmsResiduals.push_back( observationSetList_.at( observableIt.first )
                                                    .at( linkEndsIt.first )
                                                    .at( index )
                                                    ->getRmsResiduals( ) );
                }
            }
        }
        return rmsResiduals;
    }

    std::vector< Eigen::VectorXd > getMeanResiduals(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::VectorXd > meanResiduals;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for( auto observableIt: observationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto index: linkEndsIt.second )
                {
                    meanResiduals.push_back( observationSetList_.at( observableIt.first )
                                                     .at( linkEndsIt.first )
                                                     .at( index )
                                                     ->getMeanResiduals( ) );
                }
            }
        }
        return meanResiduals;
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
    getComputedObservations( std::shared_ptr< ObservationCollectionParser > observationParser =
                                     std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                computedObservations;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for( auto observableIt: observationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto index: linkEndsIt.second )
                {
                    computedObservations.push_back( observationSetList_.at( observableIt.first )
                                                            .at( linkEndsIt.first )
                                                            .at( index )
                                                            ->getComputedObservationsVector( ) );
                }
            }
        }
        return computedObservations;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getConcatenatedComputedObservations(
            std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                computedObservations = getComputedObservations( observationParser );

        unsigned int totalObsSize = 0;
        for( auto obs: computedObservations )
        {
            totalObsSize += obs.size( );
        }

        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > concatenatedComputedObservations =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObsSize );
        unsigned int obsIndex = 0;
        for( unsigned int k = 0; k < computedObservations.size( ); k++ )
        {
            concatenatedComputedObservations.block(
                    obsIndex, 0, computedObservations.at( k ).size( ), 1 ) =
                    computedObservations.at( k );
            obsIndex += computedObservations.at( k ).size( );
        }

        return concatenatedComputedObservations;
    }

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > >
    getObservationSetStartAndSize(
            std::shared_ptr< ObservationCollectionParser > observationParser )
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > >
                observationSetStartAndSize;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for( auto observableIt: observationSetsIndices )
        {
            std::map< LinkEnds, std::vector< std::pair< int, int > > > currentObservableTypeIndices;
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto index: linkEndsIt.second )
                {
                    currentObservableTypeIndices[ linkEndsIt.first ].push_back(
                            observationSetStartAndSize_.at( observableIt.first )
                                    .at( linkEndsIt.first )
                                    .at( index ) );
                }
            }
            observationSetStartAndSize[ observableIt.first ] = currentObservableTypeIndices;
        }
        return observationSetStartAndSize;
    }

    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
    getSingleLinkAndTypeObservationSets( const ObservableType observableType,
                                         const LinkDefinition linkEndsDefinition )
    {
        std::vector< std::shared_ptr< ObservationCollectionParser > > multiTypeParserList;
        multiTypeParserList.push_back( observationParser( observableType ) );
        multiTypeParserList.push_back( observationParser( linkEndsDefinition.linkEnds_ ) );
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                observationSets =
                        getSingleObservationSets( observationParser( multiTypeParserList, true ) );
        if( observationSets.empty( ) )
        {
            throw std::runtime_error(
                    " Error when getting single link observation sets, not observations of type " +
                    std::to_string( observableType ) + " for given link ends." );
        }
        return observationSets;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getSingleLinkObservations(
            const ObservableType observableType,
            const LinkDefinition& linkEndsDefinition )
    {
        std::vector< std::shared_ptr< ObservationCollectionParser > > multiTypeParserList;
        multiTypeParserList.push_back( observationParser( observableType ) );
        multiTypeParserList.push_back( observationParser( linkEndsDefinition.linkEnds_ ) );
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observations =
                getConcatenatedObservations( observationParser( multiTypeParserList, true ) );
        if( observations.size( ) == 0 )
        {
            throw std::runtime_error(
                    " Error when getting single link observations, not observations of type " +
                    std::to_string( observableType ) + " for given link ends." );
        }
        return observations;
    }

    std::vector< TimeType > getSingleLinkTimes( const ObservableType observableType,
                                                const LinkDefinition& linkEndsDefinition )
    {
        std::vector< std::shared_ptr< ObservationCollectionParser > > multiTypeParserList;
        multiTypeParserList.push_back( observationParser( observableType ) );
        multiTypeParserList.push_back( observationParser( linkEndsDefinition.linkEnds_ ) );
        std::vector< TimeType > observationTimes =
                getConcatenatedObservationTimes( observationParser( multiTypeParserList, true ) );
        if( observationTimes.empty( ) )
        {
            throw std::runtime_error(
                    " Error when getting single link observation times, not observations of type " +
                    std::to_string( observableType ) + " for given link ends." );
        }
        return observationTimes;
    }

    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > >
    getSingleLinkObservationsAndTimes( const ObservableType observableType,
                                       const LinkDefinition& linkEndsDefinition )
    {
        std::vector< std::shared_ptr< ObservationCollectionParser > > multiTypeParserList;
        multiTypeParserList.push_back( observationParser( observableType ) );
        multiTypeParserList.push_back( observationParser( linkEndsDefinition.linkEnds_ ) );
        std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
                   std::vector< TimeType > >
                observationsAndTimes = getConcatenatedObservationsAndTimes(
                        observationParser( multiTypeParserList, true ) );
        if( observationsAndTimes.second.empty( ) )
        {
            throw std::runtime_error(
                    " Error when getting single link observations and times, not observation set "
                    "of type " +
                    std::to_string( observableType ) + " for given link ends." );
        }
        return observationsAndTimes;
    }

    void setConstantWeight( const double weight = 1.0,
                            const std::shared_ptr< ObservationCollectionParser > observationParser =
                                    std::make_shared< ObservationCollectionParser >( ) )
    {
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObsSets = getSingleObservationSets( observationParser );
        if( singleObsSets.empty( ) )
        {
            std::cerr << "Warning when setting constant weights, no single observation set found "
                         "for specified observation parser. Weights not set";
        }
        for( auto singleObsSet: singleObsSets )
        {
            singleObsSet->setConstantWeight( weight );
        }
    }

    void setConstantWeight( const Eigen::VectorXd weight,
                            const std::shared_ptr< ObservationCollectionParser > observationParser =
                                    std::make_shared< ObservationCollectionParser >( ) )
    {
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObsSets = getSingleObservationSets( observationParser );
        if( singleObsSets.empty( ) )
        {
            std::cerr << "Warning when setting constant weights, no single observation set found "
                         "for specified observation parser. Weights not set";
        }
        for( auto singleObsSet: singleObsSets )
        {
            singleObsSet->setConstantWeight( weight );
        }
    }

    void setConstantWeightPerObservable(
            const std::map< std::shared_ptr< ObservationCollectionParser >, double >
                    weightsPerObservationParser )
    {
        for( auto parserIt: weightsPerObservationParser )
        {
            setConstantWeight( parserIt.second, parserIt.first );
        }
    }

    void setConstantWeightPerObservable(
            const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd >
                    weightsPerObservationParser )
    {
        for( auto parserIt: weightsPerObservationParser )
        {
            setConstantWeight( parserIt.second, parserIt.first );
        }
    }

    void setTabulatedWeights(
            const Eigen::VectorXd tabulatedWeights,
            const std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObsSets = getSingleObservationSets( observationParser );
        if( singleObsSets.empty( ) )
        {
            std::cerr << "Warning when setting tabulated weights, no single observation set found "
                         "for specified observation parser. Weights not set";
        }

        bool areObsSetsSameSize = true;
        int totalSizeAllObsSets = singleObsSets.at( 0 )->getTotalObservationSetSize( );
        for( unsigned int k = 1; k < singleObsSets.size( ); k++ )
        {
            unsigned int currentObsSetSize = singleObsSets.at( k )->getTotalObservationSetSize( );
            totalSizeAllObsSets += currentObsSetSize;
            if( currentObsSetSize != singleObsSets.at( 0 )->getTotalObservationSetSize( ) )
            {
                areObsSetsSameSize = false;
            }
        }

        unsigned int startObsSet = 0;
        for( auto obsSet: singleObsSets )
        {
            if( tabulatedWeights.size( ) == totalSizeAllObsSets )
            {
                Eigen::VectorXd singleSetWeights = tabulatedWeights.segment(
                        startObsSet, obsSet->getTotalObservationSetSize( ) );
                startObsSet += obsSet->getTotalObservationSetSize( );
                obsSet->setTabulatedWeights( singleSetWeights );
            }
            else if( areObsSetsSameSize &&
                     ( tabulatedWeights.size( ) ==
                       singleObsSets.at( 0 )->getTotalObservationSetSize( ) ) )
            {
                obsSet->setTabulatedWeights( tabulatedWeights );
            }
            else
            {
                throw std::runtime_error(
                        "Error when setting tabulated weights, the size of the input weight vector "
                        "should be consistent with "
                        "either the size of each individual observation set, or the combined size "
                        "of all required observation sets." );
            }
        }
    }

    void setTabulatedWeights( const std::map< std::shared_ptr< ObservationCollectionParser >,
                                              Eigen::VectorXd > weightsPerObservationParser )
    {
        for( auto parserIt: weightsPerObservationParser )
        {
            setTabulatedWeights( parserIt.second, parserIt.first );
        }
    }

    void appendObservationCollection(
        std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollectionToAppend )
    {
        const SortedObservationSets& setsToAppend = observationCollectionToAppend->getObservationsReference( );

        for( auto obs_it : setsToAppend )
        {
            if( observationSetList_.count( obs_it.first ) == 0 )
            {
                observationSetList_[ obs_it.first ] == obs_it.second;
            }
            else
            {
                for ( auto link_end_it : obs_it.second )
                {
                    if( observationSetList_.at( obs_it.first ).count( link_end_it.first ) == 0 )
                    {
                        observationSetList_[ obs_it.first ][ link_end_it.first ] == link_end_it.second;
                    }
                    else
                    {
                        auto listToAdd = link_end_it.second;
                        observationSetList_[ obs_it.first ][ link_end_it.first ].insert(
                            observationSetList_[ obs_it.first ][ link_end_it.first ].end( ),
                            listToAdd.begin( ), listToAdd.end( ));
                    }
                }
            }
        }
    }

    void filterObservations(
            const std::map< std::shared_ptr< ObservationCollectionParser >,
                            std::shared_ptr< ObservationFilterBase > >& observationFilters,
            const bool saveFilteredObservations = true )
    {
        // Parse all observation filters
        for( auto filterIt: observationFilters )
        {
            std::shared_ptr< ObservationCollectionParser > observationParser = filterIt.first;
            std::shared_ptr< ObservationFilterBase > filter = filterIt.second;

            // Retrieve single observation sets based on observation parser
            std::vector<
                    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                    singleObsSets = getSingleObservationSets( observationParser );

            // Filter observations for all single observation sets
            for( auto obsSet: singleObsSets )
            {
                obsSet->filterObservations( filter, saveFilteredObservations );
            }
        }

        // Reset observation set indices and concatenated observations and times
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    void filterObservations( std::shared_ptr< ObservationFilterBase > observationFilter,
                             std::shared_ptr< ObservationCollectionParser > observationParser =
                                     std::make_shared< ObservationCollectionParser >( ),
                             const bool saveFilteredObservations = true )
    {
        std::map< std::shared_ptr< ObservationCollectionParser >,
                  std::shared_ptr< ObservationFilterBase > >
                observationFilters{ { observationParser, observationFilter } };

        filterObservations( observationFilters, saveFilteredObservations );
    }

    void splitObservationSets( std::shared_ptr< ObservationSetSplitterBase > observationSetSplitter,
                               std::shared_ptr< ObservationCollectionParser > observationParser =
                                       std::make_shared< ObservationCollectionParser >( ) )
    {
        // Retrieve single observation sets based on observation parser
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                singleObsSetsIndices = getSingleObservationSetsIndices( observationParser );

        // Boolean denoting whether the warning about splitting observation sets which contain
        // filtered observations has already been printed or not (this warning should only been
        // thrown out once).
        bool warningPrinted = false;

        for( auto observableIt: singleObsSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                std::vector<
                        std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                        singleObsSets =
                                observationSetList_.at( observableIt.first ).at( linkEndsIt.first );

                unsigned int obsSetCounter = 0;
                for( auto indexSetToSplit: linkEndsIt.second )
                {
                    // Get new observation sets after splitting
                    std::vector< std::shared_ptr<
                            SingleObservationSet< ObservationScalarType, TimeType > > >
                            newObsSets = splitObservationSet(
                                    singleObsSets.at( indexSetToSplit + obsSetCounter ),
                                    observationSetSplitter,
                                    !warningPrinted );
                    warningPrinted = true;

                    // Remove original observation set
                    singleObsSets.erase( singleObsSets.begin( ) +
                                         ( indexSetToSplit + obsSetCounter ) );

                    // Replace by new sets
                    for( auto newSet: newObsSets )
                    {
                        singleObsSets.insert(
                                singleObsSets.begin( ) + ( indexSetToSplit + obsSetCounter ),
                                newSet );
                        obsSetCounter++;
                    }

                    obsSetCounter -= 1;
                }

                observationSetList_.at( observableIt.first ).at( linkEndsIt.first ) = singleObsSets;
            }
        }

        // Reset observation set indices and concatenated observations and times
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    void replaceSingleObservationSet(
            const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >&
                    newSet,
            const unsigned int setIndex )
    {
        if( observationSetList_.count( newSet->getObservableType( ) ) == 0 )
        {
            throw std::runtime_error(
                    "Error when replacing single observation set in observation collection, "
                    "observable type not found." );
        }
        if( observationSetList_.at( newSet->getObservableType( ) )
                    .count( newSet->getLinkEnds( ).linkEnds_ ) == 0 )
        {
            throw std::runtime_error(
                    "Error when replacing single observation set in observation collection, link "
                    "ends not found for given observable type." );
        }
        if( setIndex > observationSetList_.at( newSet->getObservableType( ) )
                               .at( newSet->getLinkEnds( ).linkEnds_ )
                               .size( ) -
                    1 )
        {
            throw std::runtime_error(
                    "Error when replacing single observation set in observation collection, "
                    "set index exceeds number of sets found for given observable and link ends." );
        }
        observationSetList_.at( newSet->getObservableType( ) )
                .at( newSet->getLinkEnds( ).linkEnds_ )
                .at( setIndex ) = newSet;

        // Reset observation set indices and concatenated observations and times
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    void removeSingleObservationSets(
            const std::shared_ptr< ObservationCollectionParser > observationParser )
    {
        // Retrieve single observation sets based on observation parser
        removeSingleObservationSets( getSingleObservationSetsIndices( observationParser ) );
    }

    void removeEmptySingleObservationSets( )
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                indicesSetsToRemove;

        for( auto observableIt: observationSetList_ )
        {
            std::map< LinkEnds, std::vector< unsigned int > > indicesToRemovePerLinkEnds;
            for( auto linkEndsIt: observableIt.second )
            {
                // Identify empty observation sets
                for( unsigned int k = 0; k < linkEndsIt.second.size( ); k++ )
                {
                    if( linkEndsIt.second.at( k )->getNumberOfObservables( ) == 0 )
                    {
                        indicesToRemovePerLinkEnds[ linkEndsIt.first ].push_back( k );
                    }
                }
            }
            indicesSetsToRemove[ observableIt.first ] = indicesToRemovePerLinkEnds;
        }

        // Remove empty observation sets
        removeSingleObservationSets( indicesSetsToRemove );
    }

    void removeSingleObservationSets(
            const std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >&
                    indicesSetsToRemove )
    {
        // Parse observation set list and remove selected sets
        for( auto observableIt: indicesSetsToRemove )
        {
            if( observationSetList_.count( observableIt.first ) == 0 )
            {
                throw std::runtime_error(
                        "Error when removing single observation sets, given observable type not "
                        "found." );
            }
            for( auto linkEndsIt: observableIt.second )
            {
                if( observationSetList_.at( observableIt.first ).count( linkEndsIt.first ) == 0 )
                {
                    throw std::runtime_error(
                            "Error when removing single observation sets, link ends not found for "
                            "given observable type." );
                }
                unsigned int counterRemovedSets = 0;
                for( auto indexToRemove: linkEndsIt.second )
                {
                    if( indexToRemove - counterRemovedSets >
                        observationSetList_.at( observableIt.first )
                                        .at( linkEndsIt.first )
                                        .size( ) -
                                1 )
                    {
                        throw std::runtime_error(
                                "Error when removing single observation sets, set index exceeds "
                                "number of observation sets"
                                " for given observation type and link ends." );
                    }
                    observationSetList_.at( observableIt.first )
                            .at( linkEndsIt.first )
                            .erase( observationSetList_.at( observableIt.first )
                                            .at( linkEndsIt.first )
                                            .begin( ) +
                                    indexToRemove - counterRemovedSets );
                    counterRemovedSets += 1;
                }
            }
        }

        // Reset observation set indices and concatenated observations and times
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
    getSingleObservationSetsIndices(
            const std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                observationSetsIndices;

        ObservationParserType parserType = observationParser->getObservationParserType( );
        switch( parserType )
        {
            case empty_parser: {
                for( auto observableIt: observationSetList_ )
                {
                    std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                    for( auto linkEndsIt: observableIt.second )
                    {
                        for( unsigned int k = 0; k < linkEndsIt.second.size( ); k++ )
                        {
                            indicesPerObservable[ linkEndsIt.first ].push_back( k );
                        }
                    }
                    observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                }
                break;
            }
            case observable_type_parser: {
                std::vector< ObservableType > observableTypes =
                        std::dynamic_pointer_cast< ObservationCollectionObservableTypeParser >(
                                observationParser )
                                ->getObservableTypes( );
                for( auto observableIt: observationSetList_ )
                {
                    if( ( ( std::count( observableTypes.begin( ),
                                        observableTypes.end( ),
                                        observableIt.first ) ) &&
                          ( !observationParser->useOppositeCondition( ) ) ) ||
                        ( ( std::count( observableTypes.begin( ),
                                        observableTypes.end( ),
                                        observableIt.first ) == 0 ) &&
                          ( observationParser->useOppositeCondition( ) ) ) )
                    {
                        std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                        for( auto linkEndsIt: observableIt.second )
                        {
                            for( unsigned int k = 0; k < linkEndsIt.second.size( ); k++ )
                            {
                                indicesPerObservable[ linkEndsIt.first ].push_back( k );
                            }
                        }
                        observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                    }
                }

                break;
            }
            case link_ends_parser: {
                std::vector< LinkEnds > linkEndsVector =
                        std::dynamic_pointer_cast< ObservationCollectionLinkEndsParser >(
                                observationParser )
                                ->getLinkEndsVector( );
                std::vector< LinkEnds > allLinkEnds = getLinkEnds( );
                for( auto linkEnds: allLinkEnds )
                {
                    if( ( ( std::count(
                                  linkEndsVector.begin( ), linkEndsVector.end( ), linkEnds ) ) &&
                          ( !observationParser->useOppositeCondition( ) ) ) ||
                        ( ( std::count( linkEndsVector.begin( ),
                                        linkEndsVector.end( ),
                                        linkEnds ) == 0 ) &&
                          ( observationParser->useOppositeCondition( ) ) ) )
                    {
                        for( auto observableIt: observationSetList_ )
                        {
                            std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                            if( observationSetsIndices.count( observableIt.first ) )
                            {
                                indicesPerObservable =
                                        observationSetsIndices.at( observableIt.first );
                            }

                            if( observableIt.second.count( linkEnds ) )
                            {
                                for( unsigned int k = 0;
                                     k < observableIt.second.at( linkEnds ).size( );
                                     k++ )
                                {
                                    indicesPerObservable[ linkEnds ].push_back( k );
                                }
                            }
                            if( indicesPerObservable.size( ) > 0 )
                            {
                                observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                            }
                        }
                    }
                }
                break;
            }
            case link_end_string_parser: {
                std::shared_ptr< ObservationCollectionLinkEndStringParser >
                        linkEndStringObservationParser = std::dynamic_pointer_cast<
                                ObservationCollectionLinkEndStringParser >( observationParser );
                std::vector< std::string > linkEndsNames =
                        linkEndStringObservationParser->getLinkEndNames( );
                bool isReferencePoint = linkEndStringObservationParser->isReferencePoint( );

                // Retrieve all body names
                std::vector< std::string > allLinkEndsNames;
                if( !isReferencePoint )
                {
                    allLinkEndsNames = getBodiesInLinkEnds( );
                }
                else
                {
                    allLinkEndsNames = getReferencePointsInLinkEnds( );
                }

                for( auto name: allLinkEndsNames )
                {
                    if( ( ( std::count( linkEndsNames.begin( ), linkEndsNames.end( ), name ) ) &&
                          ( !observationParser->useOppositeCondition( ) ) ) ||
                        ( ( std::count( linkEndsNames.begin( ), linkEndsNames.end( ), name ) ==
                            0 ) &&
                          ( observationParser->useOppositeCondition( ) ) ) )
                    {
                        for( auto observableIt: observationSetList_ )
                        {
                            std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                            if( observationSetsIndices.count( observableIt.first ) )
                            {
                                indicesPerObservable =
                                        observationSetsIndices.at( observableIt.first );
                            }

                            for( auto linkEndsIt: observableIt.second )
                            {
                                // FIX!
                                bool isBodyInLinkEnds = false;
                                bool isGroundStationInLinkEnds = false;
                                for( auto it: linkEndsIt.first )
                                {
                                    if( it.second.bodyName_ == name )
                                    {
                                        isBodyInLinkEnds = true;
                                    }
                                    if( it.second.stationName_ == name )
                                    {
                                        isGroundStationInLinkEnds = true;
                                    }
                                }

                                if( ( !isReferencePoint &&
                                      isBodyInLinkEnds /*( linkEndsIt.first, name )*/ ) ||
                                    ( isReferencePoint &&
                                      isGroundStationInLinkEnds /*( linkEndsIt.first, name )*/ ) )
                                {
                                    for( unsigned int k = 0; k < linkEndsIt.second.size( ); k++ )
                                    {
                                        indicesPerObservable[ linkEndsIt.first ].push_back( k );
                                    }
                                }
                            }
                            if( indicesPerObservable.size( ) > 0 )
                            {
                                observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                            }
                        }
                    }
                }
                break;
            }
            case link_end_id_parser: {
                std::shared_ptr< ObservationCollectionLinkEndIdParser > linkEndIdObservationParser =
                        std::dynamic_pointer_cast< ObservationCollectionLinkEndIdParser >(
                                observationParser );
                std::vector< LinkEndId > linkEndIds = linkEndIdObservationParser->getLinkEndIds( );

                for( auto linkEndId: linkEndIds )
                {
                    for( auto observableIt: observationSetList_ )
                    {
                        std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                        if( observationSetsIndices.count( observableIt.first ) )
                        {
                            indicesPerObservable = observationSetsIndices.at( observableIt.first );
                        }

                        for( auto linkEndsIt: observableIt.second )
                        {
                            LinkEnds linkEnds = linkEndsIt.first;

                            // Check whether relevant link end id is present in linkEnds
                            bool isLinkEndIdInLinkEnds = false;
                            for( auto it: linkEnds )
                            {
                                if( it.second == linkEndId )
                                {
                                    isLinkEndIdInLinkEnds = true;
                                }
                            }

                            if( ( isLinkEndIdInLinkEnds &&
                                  ( !observationParser->useOppositeCondition( ) ) ) ||
                                ( !isLinkEndIdInLinkEnds &&
                                  ( observationParser->useOppositeCondition( ) ) ) )
                            {
                                for( unsigned int k = 0; k < linkEndsIt.second.size( ); k++ )
                                {
                                    indicesPerObservable[ linkEndsIt.first ].push_back( k );
                                }
                            }
                        }
                        if( indicesPerObservable.size( ) > 0 )
                        {
                            observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                        }
                    }
                }
                break;
            }
            case link_end_type_parser: {
                std::shared_ptr< ObservationCollectionLinkEndTypeParser >
                        linkEndTypeObservationParser =
                                std::dynamic_pointer_cast< ObservationCollectionLinkEndTypeParser >(
                                        observationParser );
                std::vector< LinkEndType > linkEndTypes =
                        linkEndTypeObservationParser->getLinkEndTypes( );

                for( auto linkEndType: linkEndTypes )
                {
                    for( auto observableIt: observationSetList_ )
                    {
                        std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                        if( observationSetsIndices.count( observableIt.first ) )
                        {
                            indicesPerObservable = observationSetsIndices.at( observableIt.first );
                        }

                        for( auto linkEndsIt: observableIt.second )
                        {
                            LinkEnds linkEnds = linkEndsIt.first;

                            // Check whether relevant link end id is present in linkEnds
                            bool isLinkEndTypeInLinkEnds = false;
                            for( auto it: linkEnds )
                            {
                                if( it.first == linkEndType )
                                {
                                    isLinkEndTypeInLinkEnds = true;
                                }
                            }

                            if( ( isLinkEndTypeInLinkEnds &&
                                  ( !observationParser->useOppositeCondition( ) ) ) ||
                                ( !isLinkEndTypeInLinkEnds &&
                                  ( observationParser->useOppositeCondition( ) ) ) )
                            {
                                for( unsigned int k = 0; k < linkEndsIt.second.size( ); k++ )
                                {
                                    indicesPerObservable[ linkEndsIt.first ].push_back( k );
                                }
                            }
                        }
                        if( indicesPerObservable.size( ) > 0 )
                        {
                            observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                        }
                    }
                }
                break;
            }
            case single_link_end_parser: {
                std::shared_ptr< ObservationCollectionSingleLinkEndParser >
                        singleLinkEndObservationParser = std::dynamic_pointer_cast<
                                ObservationCollectionSingleLinkEndParser >( observationParser );
                std::vector< std::pair< LinkEndType, LinkEndId > > singleLinkEnds =
                        singleLinkEndObservationParser->getSingleLinkEnds( );

                for( auto singleLinkEnd: singleLinkEnds )
                {
                    for( auto observableIt: observationSetList_ )
                    {
                        std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                        if( observationSetsIndices.count( observableIt.first ) )
                        {
                            indicesPerObservable = observationSetsIndices.at( observableIt.first );
                        }

                        for( auto linkEndsIt: observableIt.second )
                        {
                            LinkEnds linkEnds = linkEndsIt.first;

                            // Check whether relevant link end type is present in linkEnds
                            bool isInLinkEnds = false;
                            for( auto it: linkEnds )
                            {
                                if( ( it.first == singleLinkEnd.first ) &&
                                    ( it.second == singleLinkEnd.second ) )
                                {
                                    isInLinkEnds = true;
                                }
                            }

                            if( ( isInLinkEnds &&
                                  ( !observationParser->useOppositeCondition( ) ) ) ||
                                ( !isInLinkEnds &&
                                  ( observationParser->useOppositeCondition( ) ) ) )
                            {
                                for( unsigned int k = 0; k < linkEndsIt.second.size( ); k++ )
                                {
                                    indicesPerObservable[ linkEndsIt.first ].push_back( k );
                                }
                            }
                        }
                        if( indicesPerObservable.size( ) > 0 )
                        {
                            observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                        }
                    }
                }
                break;
            }
            case time_bounds_parser: {
                std::vector< std::pair< double, double > > timeBoundsVector =
                        std::dynamic_pointer_cast< ObservationCollectionTimeBoundsParser >(
                                observationParser )
                                ->getTimeBoundsVector( );

                for( auto timeBounds: timeBoundsVector )
                {
                    for( auto observableIt: observationSetList_ )
                    {
                        std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                        if( observationSetsIndices.count( observableIt.first ) )
                        {
                            indicesPerObservable = observationSetsIndices.at( observableIt.first );
                        }

                        for( auto linkEndsIt: observableIt.second )
                        {
                            for( unsigned int k = 0; k < linkEndsIt.second.size( ); k++ )
                            {
                                bool isInTimeBounds = false;
                                if( ( linkEndsIt.second.at( k )->getTimeBounds( ).first >=
                                      timeBounds.first ) &&
                                    ( linkEndsIt.second.at( k )->getTimeBounds( ).second <=
                                      timeBounds.second ) )
                                {
                                    isInTimeBounds = true;
                                }
                                if( ( isInTimeBounds &&
                                      ( !observationParser->useOppositeCondition( ) ) ) ||
                                    ( !isInTimeBounds &&
                                      ( observationParser->useOppositeCondition( ) ) ) )
                                {
                                    indicesPerObservable[ linkEndsIt.first ].push_back( k );
                                }
                            }
                        }
                        if( indicesPerObservable.size( ) > 0 )
                        {
                            observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                        }
                    }
                }
                break;
            }
            case ancillary_settings_parser: {
                std::vector< std::shared_ptr< ObservationAncilliarySimulationSettings > >
                        ancillarySettings = std::dynamic_pointer_cast<
                                                    ObservationCollectionAncillarySettingsParser >(
                                                    observationParser )
                                                    ->getAncillarySettings( );

                for( auto setting: ancillarySettings )
                {
                    for( auto observableIt: observationSetList_ )
                    {
                        std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                        if( observationSetsIndices.count( observableIt.first ) )
                        {
                            indicesPerObservable = observationSetsIndices.at( observableIt.first );
                        }

                        for( auto linkEndsIt: observableIt.second )
                        {
                            for( unsigned int k = 0; k < linkEndsIt.second.size( ); k++ )
                            {
                                bool identicalAncillarySettings = false;
                                if( ( linkEndsIt.second.at( k )
                                              ->getAncilliarySettings( )
                                              ->getDoubleData( ) == setting->getDoubleData( ) ) &&
                                    ( linkEndsIt.second.at( k )
                                              ->getAncilliarySettings( )
                                              ->getDoubleVectorData( ) ==
                                      setting->getDoubleVectorData( ) ) )
                                {
                                    identicalAncillarySettings = true;
                                }
                                if( ( identicalAncillarySettings &&
                                      !observationParser->useOppositeCondition( ) ) ||
                                    ( !identicalAncillarySettings &&
                                      observationParser->useOppositeCondition( ) ) )
                                {
                                    indicesPerObservable[ linkEndsIt.first ].push_back( k );
                                }
                            }
                        }
                        if( indicesPerObservable.size( ) > 0 )
                        {
                            observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                        }
                    }
                }
                break;
            }
            case multi_type_parser: {
                std::shared_ptr< ObservationCollectionMultiTypeParser > multiTypeParser =
                        std::dynamic_pointer_cast< ObservationCollectionMultiTypeParser >(
                                observationParser );
                std::vector< std::shared_ptr< ObservationCollectionParser > > observationParsers =
                        multiTypeParser->getObservationParsers_( );

                bool areConditionsCombined = multiTypeParser->areConditionsCombined( );

                if( !areConditionsCombined )
                {
                    for( auto parser: observationParsers )
                    {
                        std::map< ObservableType,
                                  std::map< LinkEnds, std::vector< unsigned int > > >
                                currentObservationSetsIndices =
                                        getSingleObservationSetsIndices( parser );

                        for( auto observableIt: currentObservationSetsIndices )
                        {
                            if( observationSetsIndices.count( observableIt.first ) == 0 )
                            {
                                observationSetsIndices[ observableIt.first ] = observableIt.second;
                            }
                            else
                            {
                                for( auto linkEndsIt: observableIt.second )
                                {
                                    if( observationSetsIndices.at( observableIt.first )
                                                .count( linkEndsIt.first ) == 0 )
                                    {
                                        observationSetsIndices.at(
                                                observableIt.first )[ linkEndsIt.first ] =
                                                linkEndsIt.second;
                                    }
                                    else
                                    {
                                        std::vector< unsigned int > indices =
                                                observationSetsIndices.at( observableIt.first )
                                                        .at( linkEndsIt.first );
                                        for( auto index: linkEndsIt.second )
                                        {
                                            if( std::count( indices.begin( ),
                                                            indices.end( ),
                                                            index ) == 0 )
                                            {
                                                observationSetsIndices.at( observableIt.first )
                                                        .at( linkEndsIt.first )
                                                        .push_back( index );
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                else
                {
                    // First retrieve all observation sets
                    std::vector< std::map< ObservableType,
                                           std::map< LinkEnds, std::vector< unsigned int > > > >
                            allObservationSetsIndices;
                    for( auto parser: observationParsers )
                    {
                        allObservationSetsIndices.push_back(
                                getSingleObservationSetsIndices( parser ) );
                    }

                    if( allObservationSetsIndices.size( ) > 0 )
                    {
                        std::map< ObservableType,
                                  std::map< LinkEnds, std::vector< unsigned int > > >
                                originalSetsIndices = allObservationSetsIndices.at( 0 );

                        // Retrieve common observable types
                        std::vector< ObservableType > commonObsTypes =
                                utilities::createVectorFromMapKeys( originalSetsIndices );
                        for( unsigned int k = 1; k < allObservationSetsIndices.size( ); k++ )
                        {
                            std::vector< ObservableType > currentObsTypes =
                                    utilities::createVectorFromMapKeys(
                                            allObservationSetsIndices.at( k ) );
                            std::vector< ObservableType > newCommonTypes;
                            std::set_intersection( commonObsTypes.begin( ),
                                                   commonObsTypes.end( ),
                                                   currentObsTypes.begin( ),
                                                   currentObsTypes.end( ),
                                                   std::back_inserter( newCommonTypes ) );
                            commonObsTypes = newCommonTypes;
                        }

                        for( auto obsType: commonObsTypes )
                        {
                            std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;

                            // For given observable type, retrieve common link ends
                            std::vector< LinkEnds > commonLinkEndsList =
                                    utilities::createVectorFromMapKeys(
                                            originalSetsIndices.at( obsType ) );
                            for( unsigned int k = 1; k < allObservationSetsIndices.size( ); k++ )
                            {
                                std::vector< LinkEnds > currentLinkEndsList =
                                        utilities::createVectorFromMapKeys(
                                                allObservationSetsIndices.at( k ).at( obsType ) );
                                std::vector< LinkEnds > newCommonLinkEndsList;
                                std::set_intersection(
                                        commonLinkEndsList.begin( ),
                                        commonLinkEndsList.end( ),
                                        currentLinkEndsList.begin( ),
                                        currentLinkEndsList.end( ),
                                        std::back_inserter( newCommonLinkEndsList ) );
                                commonLinkEndsList = newCommonLinkEndsList;
                            }

                            for( auto linkEnds: commonLinkEndsList )
                            {
                                // For given observable type and link ends, retrieve common
                                // observation set indices
                                std::vector< unsigned int > commonIndices =
                                        originalSetsIndices.at( obsType ).at( linkEnds );
                                for( unsigned int k = 1; k < allObservationSetsIndices.size( );
                                     k++ )
                                {
                                    std::vector< unsigned int > currentIndices =
                                            allObservationSetsIndices.at( k ).at( obsType ).at(
                                                    linkEnds );
                                    std::vector< unsigned int > newCommonIndices;
                                    std::set_intersection( commonIndices.begin( ),
                                                           commonIndices.end( ),
                                                           currentIndices.begin( ),
                                                           currentIndices.end( ),
                                                           std::back_inserter( newCommonIndices ) );
                                    commonIndices = newCommonIndices;
                                }
                                indicesPerObservable[ linkEnds ] = commonIndices;
                            }
                            observationSetsIndices[ obsType ] = indicesPerObservable;
                        }
                    }
                }

                break;
            }
            default:
                throw std::runtime_error( "Observation parser type not recognised." );
        }

        return observationSetsIndices;
    }

    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
    getSingleObservationSets(
            const std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >
                singleObservationSetsIndices = getSingleObservationSetsIndices( observationParser );

        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObservationSets;
        for( auto observableIt: singleObservationSetsIndices )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto index: linkEndsIt.second )
                {
                    singleObservationSets.push_back( observationSetList_.at( observableIt.first )
                                                             .at( linkEndsIt.first )
                                                             .at( index ) );
                }
            }
        }
        return singleObservationSets;
    }

    void printObservationSetsStartAndSize( ) const
    {
        std::cout
                << "Observation collection structure: start and size of individual observation sets"
                << std::endl;
        for( auto observableIt: observationSetStartAndSize_ )
        {
            for( auto linkEndsIt: observableIt.second )
            {
                for( auto indices: linkEndsIt.second )
                {
                    std::cout << "start index " << indices.first << " - size " << indices.second
                              << " : " << "type " << observableIt.first << ", link ends "
                              << getLinkEndsString( linkEndsIt.first ) << std::endl;
                }
            }
        }
        std::cout << "\n\n";
    }

    void setReferencePoint( simulation_setup::SystemOfBodies& bodies,
                            const Eigen::Vector3d& antennaPosition,
                            const std::string& antennaName,
                            const std::string& spacecraftName,
                            const LinkEndType linkEndType,
                            const std::shared_ptr< ObservationCollectionParser > observationParser =
                                    std::make_shared< ObservationCollectionParser >( ) )
    {
        // Retrieve existing reference points with fixed position in the spacecraft-fixed frame
        std::map< std::string, std::shared_ptr< ephemerides::ConstantEphemeris > >
                existingReferencePoints = bodies.at( spacecraftName )
                                                  ->getVehicleSystems( )
                                                  ->getFixedReferencePoints( );

        // Check if reference point exists and create missing antenna if necessary
        bool antennaDetected = false;
        for( auto refPointsIt: existingReferencePoints )
        {
            if( ( refPointsIt.second->getCartesianState( ).segment( 0, 3 ) == antennaPosition ) &&
                ( refPointsIt.first == antennaName ) )
            {
                antennaDetected = true;
            }
            else if( ( refPointsIt.second->getCartesianState( ).segment( 0, 3 ) ==
                       antennaPosition ) &&
                     ( refPointsIt.first != antennaName ) )
            {
                throw std::runtime_error(
                        "Error when setting reference point in observation collection, the "
                        "required antenna position is already defined "
                        "as a reference point, but associated with different antenna name (" +
                        antennaName + ")." );
            }
            else if( ( refPointsIt.second->getCartesianState( ).segment( 0, 3 ) !=
                       antennaPosition ) &&
                     ( refPointsIt.first == antennaName ) )
            {
                throw std::runtime_error(
                        "Error when setting reference point in observation collection, the "
                        "required antenna name is already defined "
                        "as a reference point, but associated with different position." );
            }
        }
        if( !antennaDetected )
        {
            bodies.at( spacecraftName )
                    ->getVehicleSystems( )
                    ->setReferencePointPosition( antennaName, antennaPosition );
        }

        // Set reference point in observation collection
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleSets = getSingleObservationSets( observationParser );
        for( auto set: singleSets )
        {
            std::map< LinkEndType, LinkEndId > newLinkEnds = set->getLinkEnds( ).linkEnds_;
            newLinkEnds[ linkEndType ] = LinkEndId(
                    std::make_pair( newLinkEnds[ linkEndType ].bodyName_, antennaName ) );
            LinkDefinition newLinkDefinition = LinkDefinition( newLinkEnds );
            set->setLinkEnds( newLinkDefinition );
        }

        observationSetList_ =
                createSortedObservationSetList< ObservationScalarType, TimeType >( singleSets );
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    void setReferencePoints(
            simulation_setup::SystemOfBodies& bodies,
            const std::map< double, Eigen::Vector3d >& antennaSwitchHistory,
            const std::string& spacecraftName,
            const LinkEndType linkEndType,
            const std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        // Check if reference points exist and create missing antenna if necessary
        std::map< double, std::string > antennaNames;

        unsigned int counter = 0;
        for( auto antennaIt: antennaSwitchHistory )
        {
            // Retrieve existing reference points with fixed position in the spacecraft-fixed frame
            std::map< std::string, std::shared_ptr< ephemerides::ConstantEphemeris > >
                    existingReferencePoints = bodies.at( spacecraftName )
                                                      ->getVehicleSystems( )
                                                      ->getFixedReferencePoints( );

            bool antennaDetected = false;
            for( auto refPointsIt: existingReferencePoints )
            {
                if( refPointsIt.second->getCartesianState( ).segment( 0, 3 ) == antennaIt.second )
                {
                    antennaDetected = true;
                    antennaNames[ antennaIt.first ] = refPointsIt.first;
                }
            }
            if( !antennaDetected )
            {
                counter++;
                std::string newReferencePointName = "Antenna" + std::to_string( counter );
                bodies.at( spacecraftName )
                        ->getVehicleSystems( )
                        ->setReferencePointPosition( newReferencePointName, antennaIt.second );
                antennaNames[ antennaIt.first ] = newReferencePointName;
            }
        }

        // Retrieve switch times
        std::vector< double > originalSwitchTimes =
                utilities::createVectorFromMapKeys( antennaSwitchHistory );
        std::pair< TimeType, TimeType > timeBounds = getTimeBounds( );
        if( originalSwitchTimes.front( ) > timeBounds.first ||
            originalSwitchTimes.back( ) < timeBounds.second )
        {
            throw std::runtime_error(
                    "Error when setting reference points in ObservationCollection, the antenna "
                    "switch history does not cover the required observation time interval." );
        }

        // Split single observation sets accordingly
        splitObservationSets( observationSetSplitter( time_tags_splitter, originalSwitchTimes ),
                              observationParser );

        // Set reference points
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleSets = getSingleObservationSets( observationParser );

        for( auto set: singleSets )
        {
            TimeType setStartTime = set->getTimeBounds( ).first;
            TimeType setEndTime = set->getTimeBounds( ).second;

            for( unsigned int k = 0; k < originalSwitchTimes.size( ) - 1; k++ )
            {
                if( setStartTime >= originalSwitchTimes[ k ] &&
                    setEndTime <= originalSwitchTimes[ k + 1 ] )
                {
                    std::string currentAntenna = antennaNames[ originalSwitchTimes[ k ] ];
                    std::map< LinkEndType, LinkEndId > newLinkEnds = set->getLinkEnds( ).linkEnds_;
                    newLinkEnds[ linkEndType ] = LinkEndId( std::make_pair(
                            newLinkEnds[ linkEndType ].bodyName_, currentAntenna ) );
                    LinkDefinition newLinkDefinition = LinkDefinition( newLinkEnds );
                    set->setLinkEnds( newLinkDefinition );
                }
            }
        }

        observationSetList_ =
                createSortedObservationSetList< ObservationScalarType, TimeType >( singleSets );
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    void setReferencePoint(
            simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< ephemerides::Ephemeris > antennaBodyFixedEphemeris,
            const std::string& antennaName,
            const std::string& spacecraftName,
            const LinkEndType linkEndType,
            const std::shared_ptr< ObservationCollectionParser > observationParser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        // Retrieve existing reference points
        std::map< std::string, std::shared_ptr< ephemerides::Ephemeris > > existingReferencePoints =
                bodies.at( spacecraftName )->getVehicleSystems( )->getReferencePoints( );

        // Check if ephemeris already defined for given reference point
        if( existingReferencePoints.count( antennaName ) > 0 )
        {
            std::cerr
                    << "Warning when setting reference point (" << spacecraftName << ", "
                    << antennaName
                    << "), the reference point already exists. Its ephemeris will be overwritten.";
        }

        // Set reference point body-fixed ephemeris
        bodies.at( spacecraftName )
                ->getVehicleSystems( )
                ->setReferencePointPosition( antennaName, antennaBodyFixedEphemeris );

        // Set reference point
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleSets = getSingleObservationSets( observationParser );

        for( auto set: singleSets )
        {
            std::map< LinkEndType, LinkEndId > newLinkEnds = set->getLinkEnds( ).linkEnds_;
            newLinkEnds[ linkEndType ] = LinkEndId(
                    std::make_pair( newLinkEnds[ linkEndType ].bodyName_, antennaName ) );
            LinkDefinition newLinkDefinition = LinkDefinition( newLinkEnds );
            set->setLinkEnds( newLinkDefinition );
        }

        observationSetList_ =
                createSortedObservationSetList< ObservationScalarType, TimeType >( singleSets );
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    void setTransponderDelay(
            const std::string& spacecraftName,
            const double transponderDelay,
            const std::shared_ptr< ObservationCollectionParser > inputObservationParser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        // Create observation parser with the spacecraft name
        std::shared_ptr< ObservationCollectionParser > spacecraftParser =
                observationParser( spacecraftName );

        // Create combined observation parser
        std::shared_ptr< ObservationCollectionMultiTypeParser > jointParser =
                std::make_shared< ObservationCollectionMultiTypeParser >(
                        std::vector< std::shared_ptr< ObservationCollectionParser > >(
                                { spacecraftParser, inputObservationParser } ),
                        true );

        // Set transponder delay
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleSets = getSingleObservationSets( jointParser );

        for( auto set: singleSets )
        {
            std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >
                    ancilliarySettings = set->getAncilliarySettings( );

            if( ancilliarySettings != nullptr )
            {
                std::vector< double > linkEndsDelays_ =
                        ancilliarySettings->getAncilliaryDoubleVectorData( link_ends_delays,
                                                                           false );
                linkEndsDelays_[ 1 ] = transponderDelay;
                ancilliarySettings->setAncilliaryDoubleVectorData( link_ends_delays,
                                                                   linkEndsDelays_ );
            }
        }
    }

    // Function to add an observation dependent variable to (a subset of) the single observation
    // sets
    std::shared_ptr< ObservationCollectionParser > addDependentVariable(
            std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > settings,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< ObservationCollectionParser > parser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        // Create observation collection parser corresponding to the required dependent variable
        // settings
        std::shared_ptr< ObservationCollectionParser > dependentVariablesParser =
                getObservationParserFromDependentVariableSettings( settings );

        // Create combined observation parser
        std::shared_ptr< ObservationCollectionMultiTypeParser > jointParser =
                std::make_shared< ObservationCollectionMultiTypeParser >(
                        std::vector< std::shared_ptr< ObservationCollectionParser > >(
                                { dependentVariablesParser, parser } ),
                        true );

        // Retrieve single observation sets from joint parser
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObsSets = getSingleObservationSets( jointParser );

        // Parse single observation sets
        for( auto set: singleObsSets )
        {
            // Retrieve single observation set observable type and link ends
            ObservableType observableType = set->getObservableType( );
            LinkEnds linkEnds = set->getLinkEnds( ).linkEnds_;

            // Retrieve complete list of all dependent variables settings compatible with the input
            // settings (which might not be fully defined, i.e. with missing link ends information,
            // etc.)
            std::vector< std::shared_ptr< ObservationDependentVariableSettings > >
                    allSettingsToCreate = createAllCompatibleDependentVariableSettings(
                            observableType, linkEnds, settings );

            // Add dependent variables for all compatible settings
            set->addDependentVariables( allSettingsToCreate, bodies );
        }

        // Return joint observation collection parser
        return jointParser;
    }

    // Function to retrieve the values of a given dependent variable for (a subset of) the single
    // observation sets, sorted per single observation set
    std::pair< std::vector< Eigen::MatrixXd >, std::shared_ptr< ObservationCollectionParser > >
    getDependentVariables( std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >
                                   dependentVariableSettings,
                           const bool returnFirstCompatibleSettings = false,
                           const std::shared_ptr< ObservationCollectionParser > parser =
                                   std::make_shared< ObservationCollectionParser >( ) )
    {
        // Create observation collection parser corresponding to the required dependent variable
        // settings
        std::shared_ptr< ObservationCollectionParser > dependentVariablesParser =
                getObservationParserFromDependentVariableSettings( dependentVariableSettings );

        // Create combined observation parser (for which the relevant dependent variable is also
        // defined)
        std::shared_ptr< ObservationCollectionMultiTypeParser > jointParser =
                std::make_shared< ObservationCollectionMultiTypeParser >(
                        std::vector< std::shared_ptr< ObservationCollectionParser > >(
                                { dependentVariablesParser, parser } ),
                        true );

        // Retrieve single observation sets from joint parser
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObsSets = getSingleObservationSets( jointParser );

        // Retrieve dependent variables from the joint observation parser
        std::vector< Eigen::MatrixXd > dependentVariablesValues;
        for( auto set: singleObsSets )
        {
            dependentVariablesValues.push_back( set->getSingleDependentVariable(
                    dependentVariableSettings, returnFirstCompatibleSettings ) );
        }

        // Return dependent variables values and joint parser
        return std::make_pair( dependentVariablesValues, jointParser );
    }

    //! Function that retrieves the time history of a given dependent variable, sorted per observation set. It must be noted that the reported epochs are the times at which the
    //! observations are computed/acquired, which might differ from the times at which the dependent variables are evaluated.
    std::vector< std::map< TimeType, Eigen::VectorXd > >
    getDependentVariableHistoryPerObservationSet(
            std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >
                    dependentVariableSettings,
            const bool returnFirstCompatibleSettings = false,
            const std::shared_ptr< ObservationCollectionParser > parser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        std::pair< std::vector< Eigen::MatrixXd >, std::shared_ptr< ObservationCollectionParser > >
                dependentVariableResultList = getDependentVariables(
                        dependentVariableSettings, returnFirstCompatibleSettings, parser );

        // Retrieve single observation sets for which the dependent variable of interest is computed
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObservationSets =
                        getSingleObservationSets( dependentVariableResultList.second );

        // Parse all relevant single observation sets
        std::vector< std::map< TimeType, Eigen::VectorXd > > dependentVariableResult;
        for( auto set: singleObservationSets )
        {
            dependentVariableResult.push_back( set->getSingleDependentVariableHistory(
                    dependentVariableSettings, returnFirstCompatibleSettings ) );
        }
        return dependentVariableResult;
    }

    //! Function that retrieves the time history of a given dependent variable, concatenated over all relevant observation sets.
    //! It must be noted that the reported epochs are the times at which the observations are computed/acquired, which might
    //! differ from the times at which the dependent variables are evaluated.
    std::map< TimeType, Eigen::VectorXd > getDependentVariableHistory(
            std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >
                    dependentVariableSettings,
            const bool returnFirstCompatibleSettings = false,
            const std::shared_ptr< ObservationCollectionParser > parser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        std::vector< std::map< TimeType, Eigen::VectorXd > >
                dependentVariableResultPerObservationSet =
                        getDependentVariableHistoryPerObservationSet(
                                dependentVariableSettings, returnFirstCompatibleSettings, parser );
        return utilities::concatenateMaps( dependentVariableResultPerObservationSet );
    }

    //! Function returning the list of all dependent variable settings compatible with the settings provided as inputs
    //! (which might not be fully defined, i.e. with missing link ends information, etc.). The inner vector of the first element of the pair
    //! contains the list of compatible settings, for each observation sets (outer vector).
    std::pair<
            std::vector< std::vector< std::shared_ptr< ObservationDependentVariableSettings > > >,
            std::shared_ptr< ObservationCollectionParser > >
    getCompatibleDependentVariablesSettingsList(
            std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings,
            std::shared_ptr< ObservationCollectionParser > parser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        // Create observation collection parser corresponding to the required dependent variable
        // settings
        std::shared_ptr< ObservationCollectionParser > dependentVariablesParser =
                getObservationParserFromDependentVariableSettings( dependentVariableSettings );

        // Define joint observation collection parser
        std::shared_ptr< ObservationCollectionMultiTypeParser > jointParser =
                std::make_shared< ObservationCollectionMultiTypeParser >(
                        std::vector< std::shared_ptr< ObservationCollectionParser > >(
                                { dependentVariablesParser, parser } ),
                        true );

        // Retrieve relevant single observation sets
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObservationSets = getSingleObservationSets( jointParser );

        // Parse all single observation sets
        std::vector< std::vector< std::shared_ptr< ObservationDependentVariableSettings > > >
                dependentVariablesListPerSet;
        for( auto set: singleObservationSets )
        {
            std::vector< std::shared_ptr< ObservationDependentVariableSettings > >
                    currentVariableSettingsList = set->getCompatibleDependentVariablesSettingsList(
                            dependentVariableSettings );
            if( currentVariableSettingsList.size( ) > 0 )
            {
                dependentVariablesListPerSet.push_back( currentVariableSettingsList );
            }
        }

        return std::make_pair( dependentVariablesListPerSet, jointParser );
    }

    //! Function returning the values of all dependent variables compatible with the settings provided as input, sorted per observation sets (outer
    //! vector of the first element of the pair). The order in which the settings are provided in the inner vector matches the settings list output
    //! of the getCompatibleDependentVariablesSettingsList function
    std::pair< std::vector< std::vector< Eigen::MatrixXd > >,
               std::shared_ptr< ObservationCollectionParser > >
    getAllCompatibleDependentVariables(
            std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings,
            std::shared_ptr< ObservationCollectionParser > parser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        // Create observation collection parser corresponding to the required dependent variable
        // settings
        std::shared_ptr< ObservationCollectionParser > dependentVariablesParser =
                getObservationParserFromDependentVariableSettings( dependentVariableSettings );

        // Define joint observation collection parser
        std::shared_ptr< ObservationCollectionMultiTypeParser > jointParser =
                std::make_shared< ObservationCollectionMultiTypeParser >(
                        std::vector< std::shared_ptr< ObservationCollectionParser > >(
                                { dependentVariablesParser, parser } ),
                        true );

        // Retrieve relevant single observation sets
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObservationSets = getSingleObservationSets( jointParser );

        // Parse all single observation sets
        std::vector< std::vector< Eigen::MatrixXd > > dependentVariablesListPerSet;
        for( auto set: singleObservationSets )
        {
            std::vector< Eigen::MatrixXd > currentVariablesList =
                    set->getAllCompatibleDependentVariables( dependentVariableSettings );
            if( currentVariablesList.size( ) > 0 )
            {
                dependentVariablesListPerSet.push_back( currentVariablesList );
            }
        }

        return std::make_pair( dependentVariablesListPerSet, jointParser );
    }

    // Function to retrieve the values of a given dependent variable for (a subset of) the single
    // observation sets, concatenated over all relevant single observation sets
    std::pair< Eigen::MatrixXd, std::shared_ptr< ObservationCollectionParser > >
    getConcatenatedDependentVariables(
            std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >
                    dependentVariableSettings,
            const bool returnFirstCompatibleSettings = false,
            const std::shared_ptr< ObservationCollectionParser > parser =
                    std::make_shared< ObservationCollectionParser >( ) )
    {
        // Retrieve dependent variables stored per single observation sets
        std::pair< std::vector< Eigen::MatrixXd >, std::shared_ptr< ObservationCollectionParser > >
                dependentVariables = getDependentVariables(
                        dependentVariableSettings, returnFirstCompatibleSettings, parser );

        // Compute size concatenated dependent variables matrix
        unsigned int dependentVariableSize = 0;
        if( dependentVariables.first.size( ) > 0 )
        {
            dependentVariableSize = dependentVariables.first[ 0 ].cols( );
        }

        unsigned int dependentVariablesRows = 0;
        for( auto it: dependentVariables.first )
        {
            if( it.cols( ) != dependentVariableSize )
            {
                throw std::runtime_error(
                        "Error when concatenated dependent variable over obs collection, dependent "
                        "variable size is inconsistent between single observation sets." );
            }
            dependentVariablesRows += it.rows( );
        }

        // Concatenate dependent variables over all single observation sets
        Eigen::MatrixXd concatenatedDependentVariables =
                Eigen::MatrixXd::Zero( dependentVariablesRows, dependentVariableSize );
        unsigned int index = 0;
        for( unsigned int k = 0; k < dependentVariables.first.size( ); k++ )
        {
            concatenatedDependentVariables.block(
                    index, 0, dependentVariables.first.at( k ).rows( ), dependentVariableSize ) =
                    dependentVariables.first.at( k );
            index += dependentVariables.first.at( k ).rows( );
        }

        return std::make_pair( concatenatedDependentVariables, dependentVariables.second );
    }

   private:
    void setObservationSetIndices( )
    {
        int currentStartIndex = 0;
        int currentTypeStartIndex = 0;
        totalNumberOfObservables_ = 0;
        totalObservableSize_ = 0;

        observationSetStartAndSize_.clear( );
        concatenatedObservationSetStartAndSize_.clear( );
        observationTypeAndLinkEndStartAndSize_.clear( );
        observationTypeStartAndSize_.clear( );

        for( auto observationIterator: observationSetList_ )
        {
            ObservableType currentObservableType = observationIterator.first;

            currentTypeStartIndex = currentStartIndex;
            int observableSize = getObservableSize( currentObservableType );

            int currentObservableTypeSize = 0;

            for( auto linkEndIterator: observationIterator.second )
            {
                LinkEnds currentLinkEnds = linkEndIterator.first;
                int currentLinkEndStartIndex = currentStartIndex;
                int currentLinkEndSize = 0;
                for( unsigned int i = 0; i < linkEndIterator.second.size( ); i++ )
                {
                    int currentNumberOfObservables =
                            linkEndIterator.second.at( i )->getNumberOfObservables( );
                    int currentObservableVectorSize = currentNumberOfObservables * observableSize;

                    observationSetStartAndSize_[ currentObservableType ][ currentLinkEnds ]
                            .push_back( std::make_pair( currentStartIndex,
                                                        currentObservableVectorSize ) );
                    concatenatedObservationSetStartAndSize_.push_back(
                            std::make_pair( currentStartIndex, currentObservableVectorSize ) );
                    currentStartIndex += currentObservableVectorSize;
                    currentObservableTypeSize += currentObservableVectorSize;
                    currentLinkEndSize += currentObservableVectorSize;

                    totalObservableSize_ += currentObservableVectorSize;
                    totalNumberOfObservables_ += currentNumberOfObservables;
                }
                observationTypeAndLinkEndStartAndSize_[ currentObservableType ][ currentLinkEnds ] =
                        std::make_pair( currentLinkEndStartIndex, currentLinkEndSize );
            }
            observationTypeStartAndSize_[ currentObservableType ] =
                    std::make_pair( currentTypeStartIndex, currentObservableTypeSize );
        }
    }

    void setConcatenatedObservationsAndTimes( )
    {
        concatenatedObservations_ = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero(
                totalObservableSize_ );
        concatenatedTimes_.resize( totalObservableSize_ );
        concatenatedLinkEndIds_.resize( totalObservableSize_ );
        concatenatedLinkEndIdNames_.resize( totalObservableSize_ );

        linkDefinitionsPerObservable_.clear( );
        linkEndIds_.clear( );
        inverseLinkEndIds_.clear( );
        observationSetStartAndSizePerLinkEndIndex_.clear( );

        int observationCounter = 0;
        int maximumStationId = 0;
        int currentStationId;

        for( auto observationIterator: observationSetList_ )
        {
            ObservableType currentObservableType = observationIterator.first;
            int observableSize = getObservableSize( currentObservableType );

            for( auto linkEndIterator: observationIterator.second )
            {
                LinkEnds currentLinkEnds = linkEndIterator.first;
                LinkDefinition firstLinkDefinition;

                if( linkEndIterator.second.size( ) > 0 )
                {
                    firstLinkDefinition = linkEndIterator.second.at( 0 )->getLinkEnds( );
                    linkDefinitionsPerObservable_[ currentObservableType ].push_back(
                            firstLinkDefinition );
                }

                if( linkEndIds_.count( currentLinkEnds ) == 0 )
                {
                    linkEndIds_[ currentLinkEnds ] = maximumStationId;
                    inverseLinkEndIds_[ maximumStationId ] = currentLinkEnds;
                    currentStationId = maximumStationId;
                    maximumStationId++;
                }
                else
                {
                    currentStationId = linkEndIds_[ currentLinkEnds ];
                }

                for( unsigned int i = 0; i < linkEndIterator.second.size( ); i++ )
                {
                    LinkDefinition currentLinkDefinition =
                            linkEndIterator.second.at( i )->getLinkEnds( );
                    if( !( currentLinkDefinition == firstLinkDefinition ) )
                    {
                        throw std::runtime_error(
                                "Error when creating ObservationCollection, link definitions of "
                                "same link ends are not equal " );
                    }

                    std::pair< int, int > startAndSize =
                            observationSetStartAndSize_.at( currentObservableType )
                                    .at( currentLinkEnds )
                                    .at( i );
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > currentObservables =
                            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero(
                                    startAndSize.second );

                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                            currentObservationSet =
                                    linkEndIterator.second.at( i )->getObservations( );
                    std::vector< TimeType > currentObservationTimes =
                            linkEndIterator.second.at( i )->getObservationTimes( );
                    for( unsigned int j = 0; j < currentObservationSet.size( ); j++ )
                    {
                        currentObservables.segment( j * observableSize, observableSize ) =
                                currentObservationSet.at( j );
                        for( int k = 0; k < observableSize; k++ )
                        {
                            concatenatedTimes_[ observationCounter ] =
                                    currentObservationTimes.at( j );
                            concatenatedLinkEndIds_[ observationCounter ] = currentStationId;
                            concatenatedLinkEndIdNames_[ observationCounter ] = currentLinkEnds;
                            observationCounter++;
                        }
                    }
                    concatenatedObservations_.segment( startAndSize.first, startAndSize.second ) =
                            currentObservables;
                }
            }
        }

        for( auto it1: observationSetStartAndSize_ )
        {
            for( auto it2: it1.second )
            {
                observationSetStartAndSizePerLinkEndIndex_[ it1.first ]
                                                          [ linkEndIds_[ it2.first ] ] = it2.second;
            }
        }
    }

    std::shared_ptr< ObservationCollectionParser >
    getObservationParserFromDependentVariableSettings(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >
                    dependentVariableSettings )
    {
        std::shared_ptr< ObservationCollectionParser > observationParser;

        if( !isObservationDependentVariableAncilliarySetting(
                    dependentVariableSettings->variableType_ ) )
        {
            std::vector< std::shared_ptr< ObservationCollectionParser > > parserList;

            // Check if relevant link end id and type are both specified
            if( ( dependentVariableSettings->linkEndId_ != LinkEndId( "", "" ) ) &&
                ( dependentVariableSettings->linkEndType_ != unidentified_link_end ) )
            {
                parserList.push_back( std::make_shared< ObservationCollectionSingleLinkEndParser >(
                        std::make_pair( dependentVariableSettings->linkEndType_,
                                        dependentVariableSettings->linkEndId_ ) ) );
            }
            // if only relevant link end id is specified
            else if( dependentVariableSettings->linkEndId_ != LinkEndId( "", "" ) )
            {
                parserList.push_back( std::make_shared< ObservationCollectionLinkEndIdParser >(
                        dependentVariableSettings->linkEndId_ ) );
            }
            // if only relevant link end type is specified
            else if( dependentVariableSettings->linkEndType_ != unidentified_link_end )
            {
                parserList.push_back( std::make_shared< ObservationCollectionLinkEndTypeParser >(
                        dependentVariableSettings->linkEndType_ ) );
            }

            // Check if originating link end id and type are both specified
            if( ( dependentVariableSettings->originatingLinkEndId_ != LinkEndId( "", "" ) ) &&
                ( dependentVariableSettings->originatingLinkEndType_ != unidentified_link_end ) )
            {
                parserList.push_back( std::make_shared< ObservationCollectionSingleLinkEndParser >(
                        std::make_pair( dependentVariableSettings->originatingLinkEndType_,
                                        dependentVariableSettings->originatingLinkEndId_ ) ) );
            }
            // if only originating link end id is specified
            else if( dependentVariableSettings->originatingLinkEndId_ != LinkEndId( "", "" ) )
            {
                parserList.push_back( std::make_shared< ObservationCollectionLinkEndIdParser >(
                        dependentVariableSettings->originatingLinkEndId_ ) );
            }
            // if only originating link end type is specified
            else if( dependentVariableSettings->originatingLinkEndType_ != unidentified_link_end )
            {
                parserList.push_back( std::make_shared< ObservationCollectionLinkEndTypeParser >(
                        dependentVariableSettings->originatingLinkEndType_ ) );
            }

            // Create multi-type observation collection parser
            if( parserList.size( ) > 0 )
            {
                observationParser = std::make_shared< ObservationCollectionMultiTypeParser >(
                        parserList, true );
            }
            else
            {
                observationParser = std::make_shared< ObservationCollectionParser >( );
            }
        }
        else
        {
            std::vector< std::shared_ptr< ObservationCollectionParser > > parserList;

            // Check if indirect condition on observable type via ancillary settings
            if( std::dynamic_pointer_cast<
                        simulation_setup::AncillaryObservationDependentVariableSettings >(
                        dependentVariableSettings ) != nullptr )
            {
                std::shared_ptr< simulation_setup::AncillaryObservationDependentVariableSettings >
                        ancillaryDependentVariables = std::dynamic_pointer_cast<
                                simulation_setup::AncillaryObservationDependentVariableSettings >(
                                dependentVariableSettings );
                if( ancillaryDependentVariables->observableType_ != undefined_observation_model )
                {
                    parserList.push_back(
                            std::make_shared< ObservationCollectionObservableTypeParser >(
                                    ancillaryDependentVariables->observableType_ ) );
                }
                else
                {
                    std::vector< ObservableType > allObservableTypes =
                            utilities::createVectorFromMapKeys( observationSetList_ );
                    for( auto observableIt: allObservableTypes )
                    {
                        if( ancillaryDependentVariables->isObservableTypeCompatible_(
                                    observableIt ) )
                        {
                            parserList.push_back(
                                    std::make_shared< ObservationCollectionObservableTypeParser >(
                                            observableIt ) );
                        }
                    }
                }
            }

            // Create multi-type observation collection parser
            if( parserList.size( ) > 0 )
            {
                observationParser = std::make_shared< ObservationCollectionMultiTypeParser >(
                        parserList, false );
            }
            else
            {
                observationParser = std::make_shared< ObservationCollectionParser >( );
            }
        }

        return observationParser;
    }

    SortedObservationSets observationSetList_;

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > concatenatedObservations_;

    std::vector< TimeType > concatenatedTimes_;

    Eigen::Matrix< double, Eigen::Dynamic, 1 > concatenatedWeights_;

    std::vector< int > concatenatedLinkEndIds_;

    std::vector< LinkEnds > concatenatedLinkEndIdNames_;

    std::map< ObservableType, std::vector< LinkDefinition > > linkDefinitionsPerObservable_;

    std::map< observation_models::LinkEnds, int > linkEndIds_;

    std::map< int, observation_models::LinkEnds > inverseLinkEndIds_;

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > >
            observationSetStartAndSize_;

    std::vector< std::pair< int, int > > concatenatedObservationSetStartAndSize_;

    std::map< ObservableType, std::map< int, std::vector< std::pair< int, int > > > >
            observationSetStartAndSizePerLinkEndIndex_;

    std::map< ObservableType, std::map< LinkEnds, std::pair< int, int > > >
            observationTypeAndLinkEndStartAndSize_;

    std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize_;

    int totalObservableSize_;

    int totalNumberOfObservables_;
};

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if<
                  is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value,
                  int >::type = 0 >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
createNewObservationCollection(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
                observationCollection,
        const std::shared_ptr< ObservationCollectionParser > observationParser =
                std::make_shared< ObservationCollectionParser >( ) )
{
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
            oldObservationSets =
                    observationCollection->getSingleObservationSets( observationParser );

    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
            newSingleObservationSets;
    for( auto oldObsSet: oldObservationSets )
    {
        std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > newObsSet =
                std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                        oldObsSet->getObservableType( ),
                        oldObsSet->getLinkEnds( ),
                        oldObsSet->getObservations( ),
                        oldObsSet->getObservationTimes( ),
                        oldObsSet->getReferenceLinkEnd( ),
                        oldObsSet->getObservationsDependentVariablesReference( ),
                        oldObsSet->getDependentVariableCalculator( ),
                        oldObsSet->getAncilliarySettings( ) );
        newObsSet->setTabulatedWeights( oldObsSet->getWeightsVector( ) );
        newObsSet->setResiduals( oldObsSet->getResiduals( ) );

        newSingleObservationSets.push_back( newObsSet );
    }
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >(
            newSingleObservationSets );
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if<
                  is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value,
                  int >::type = 0 >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > filterObservations(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
                observationCollection,
        const std::map< std::shared_ptr< ObservationCollectionParser >,
                        std::shared_ptr< ObservationFilterBase > >& observationFilters )
{
    // Create new observation collection
    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
            newObservationCollection = createNewObservationCollection( observationCollection );

    // Apply filtering to observation collection
    newObservationCollection->filterObservations( observationFilters, false );

    return newObservationCollection;
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if<
                  is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value,
                  int >::type = 0 >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > filterObservations(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
                observationCollection,
        const std::shared_ptr< ObservationFilterBase > observationFilter,
        const std::shared_ptr< ObservationCollectionParser > observationParser =
                std::make_shared< ObservationCollectionParser >( ) )
{
    // Create new observation collection
    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
            newObservationCollection = createNewObservationCollection( observationCollection );

    // Apply filtering to observation collection
    newObservationCollection->filterObservations( observationFilter, observationParser, false );

    return newObservationCollection;
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if<
                  is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value,
                  int >::type = 0 >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > splitObservationSets(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
                observationCollection,
        const std::shared_ptr< ObservationSetSplitterBase > observationSetSplitter,
        const std::shared_ptr< ObservationCollectionParser > observationParser =
                std::make_shared< ObservationCollectionParser >( ) )
{
    // Create new observation collection
    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
            newObservationCollection = createNewObservationCollection( observationCollection );

    // Split observation sets within new observation collection
    newObservationCollection->splitObservationSets( observationSetSplitter, observationParser );

    return newObservationCollection;
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if<
                  is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value,
                  int >::type = 0 >
inline std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >
createSingleObservationSet(
        const ObservableType observableType,
        const LinkEnds& linkEnds,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
                observations,
        const std::vector< TimeType > observationTimes,
        const LinkEndType referenceLinkEnd,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >
                ancilliarySettings )
{
    return std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
            observableType,
            linkEnds,
            observations,
            observationTimes,
            referenceLinkEnd,
            std::vector< Eigen::VectorXd >( ),
            nullptr,
            ancilliarySettings );
}

// template< typename ObservationScalarType = double, typename TimeType = double,
//           typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType
//           >::value, int >::type = 0 >
// inline std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
// createManualObservationCollection(
//         const std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr<
//         SingleObservationSet< ObservationScalarType, TimeType > > > > >& observationSetList )
//{
//     return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >(
//     observationSetList );
// }

// template< typename ObservationScalarType = double, typename TimeType = double,
//           typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType
//           >::value, int >::type = 0 >
// inline std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
// createManualObservationCollection(
//         const ObservableType observableType,
//         const LinkEnds& linkEnds,
//         const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >&
//         observationSetList )
//{
//     return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >(
//     observationSetList );
// }

template< typename ObservationScalarType = double, typename TimeType = double >
inline std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
createManualObservationCollection(
        const ObservableType observableType,
        const LinkDefinition& linkEnds,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
                observations,
        const std::vector< TimeType > observationTimes,
        const LinkEndType referenceLinkEnd,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >
                ancilliarySettings = nullptr )
{
    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >
            singleObservationSet = createSingleObservationSet( observableType,
                                                               linkEnds.linkEnds_,
                                                               observations,
                                                               observationTimes,
                                                               referenceLinkEnd,
                                                               ancilliarySettings );

    std::map< ObservableType,
              std::map< LinkEnds,
                        std::vector< std::shared_ptr<
                                SingleObservationSet< ObservationScalarType, TimeType > > > > >
            observationSetList;
    observationSetList[ observableType ][ linkEnds.linkEnds_ ].push_back( singleObservationSet );
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >(
            observationSetList );
}

template< typename ObservationScalarType = double, typename TimeType = double >
inline std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
createManualObservationCollection(
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
                singleObservationSets )
{
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >(
            singleObservationSets );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >
mergeObservationCollections(
    std::vector< std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > > observationCollectionList )
{
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > combinedObservationSets;
    for( unsigned int i = 0; i < observationCollectionList.size( ); i++ )
    {
        auto currentObservationSets = observationCollectionList.at( i )->getObservationsReference( );
        {
            for( auto obs_it : currentObservationSets )
            {
                for( auto link_end_it : obs_it.second )
                {
                    auto listToAdd = link_end_it.second;
                    combinedObservationSets.insert( combinedObservationSets.end( ), listToAdd.begin( ), listToAdd.end( ) );
                }
            }
        }
    }
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( combinedObservationSets );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_COLLECTION_H
