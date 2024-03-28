/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONS_H
#define TUDAT_OBSERVATIONS_H

#include <vector>

#include <memory>
#include <functional>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class SingleObservationSet
{
public:
    SingleObservationSet(
            const ObservableType observableType,
            const LinkDefinition& linkEnds,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType > observationTimes,
            const LinkEndType referenceLinkEnd,
            const std::vector< Eigen::VectorXd >& observationsDependentVariables = std::vector< Eigen::VectorXd >( ),
            const std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr ):
        observableType_( observableType ),
        linkEnds_( linkEnds ),
        observations_( observations ),
        observationTimes_( observationTimes ),
        referenceLinkEnd_( referenceLinkEnd ),
        observationsDependentVariables_( observationsDependentVariables ),
        dependentVariableCalculator_( dependentVariableCalculator ),
        ancilliarySettings_( ancilliarySettings ),
        numberOfObservations_( observations_.size( ) )
    {
        if( dependentVariableCalculator_ != nullptr )
        {
            if( dependentVariableCalculator_->getObservableType( ) != observableType_ )
            {
                throw std::runtime_error( "Error when creating SingleObservationSet, ObservationDependentVariableCalculator has incompatible type " );
            }

            if( !( dependentVariableCalculator_->getLinkEnds( ) == linkEnds ) )
            {
                throw std::runtime_error( "Error when creating SingleObservationSet, ObservationDependentVariableCalculator has incompatible link ends " );
            }
        }

        if( observations_.size( ) != observationTimes_.size( ) )
        {
            throw std::runtime_error( "Error when making SingleObservationSet, input sizes are inconsistent." +
                std::to_string( observations_.size( ) ) + ", " + std::to_string( observationTimes_.size( ) ) );
        }

        for( unsigned int i = 1; i < observations.size( ); i++ )
        {
            if( observations.at( i ).rows( ) != observations.at( i - 1 ).rows( ) )
            {
                throw std::runtime_error( "Error when making SingleObservationSet, input observables not of consistent size." );
            }
        }

        if( !std::is_sorted( observationTimes_.begin( ), observationTimes_.end( ) ) )
        {
            std::multimap< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observationsMap;
            for( unsigned int i = 0; i < observations_.size( ); i++ )
            {
                observationsMap.insert( { observationTimes_.at( i ), observations_.at( i ) } );
            }
            observationTimes_ = utilities::createVectorFromMultiMapKeys( observationsMap );
            observations_ = utilities::createVectorFromMultiMapValues( observationsMap );
            if( observationsDependentVariables_.size( ) > 0 )
            {
                std::map< TimeType, Eigen::VectorXd > observationsDependentVariablesMap;
                for( unsigned int i = 0; i < observationsDependentVariables_.size( ); i++ )
                {
                    observationsDependentVariablesMap[ observationTimes_.at( i ) ] = observationsDependentVariables_.at( i );
                }
                observationsDependentVariables_ = utilities::createVectorFromMapValues( observationsDependentVariablesMap );
            }
            if( static_cast< int >( observations_.size( ) ) != numberOfObservations_ )
            {
                throw std::runtime_error( "Error when making SingleObservationSet, number of observations is incompatible after time ordering" );
            }
        }
    }



    ObservableType getObservableType( )
    {
        return observableType_;
    }

    LinkDefinition getLinkEnds( )
    {
        return linkEnds_;
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservations( )
    {
        return observations_;
    }

    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& getObservationsReference( )
    {
        return observations_;
    }


    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservation( const int index )
    {
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when retrieving single observation, index is out of bounds" );
        }
        return observations_.at( index );
    }

    std::vector< TimeType > getObservationTimes( )
    {
        return observationTimes_;
    }

    const std::vector< TimeType >& getObservationTimesReference( )
    {
        return observationTimes_;
    }

    LinkEndType getReferenceLinkEnd( )
    {
        return referenceLinkEnd_;
    }

    int getNumberOfObservables( )
    {
        return numberOfObservations_;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationsVector( )
    {
        int singleObservableSize = 0;
        if( numberOfObservations_ != 0 )
        {
            singleObservableSize = observations_.at( 0 ).rows( );
        }

        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero(
                    singleObservableSize * numberOfObservations_ );
        for( unsigned int i = 0; i < observations_.size( ); i++ )
        {
            observationsVector.segment( i * singleObservableSize, singleObservableSize ) =
                    observations_.at( i );
        }
        return observationsVector;
    }

    std::map< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservationsHistory( )
    {
        return utilities::createMapFromVectors< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >(
                    observationTimes_, observations_ );
    }


    std::vector< Eigen::VectorXd > getObservationsDependentVariables( )
    {
        return observationsDependentVariables_;
    }

    std::vector< Eigen::VectorXd >& getObservationsDependentVariablesReference( )
    {
        return observationsDependentVariables_;
    }


    std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator > getDependentVariableCalculator( )
    {
        return dependentVariableCalculator_;
    }

    std::map< TimeType, Eigen::VectorXd > getDependentVariableHistory( )
    {
        return utilities::createMapFromVectors< TimeType, Eigen::VectorXd >(
                    observationTimes_, observationsDependentVariables_ );
    }

    std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > getAncilliarySettings( )
    {
        return ancilliarySettings_;
    }

    Eigen::VectorXd getWeightsVector( )
    {
        return weightsVector_;
    }

    Eigen::VectorXd& getWeightsVectorReference( )
    {
        return weightsVector_;
    }

    void setWeightsVector( const Eigen::VectorXd& weightsVector )
    {
        int singleObservableSize = 0;
        if( numberOfObservations_ != 0 )
        {
            singleObservableSize = observations_.at( 0 ).rows( );
            if( weightsVector.rows( ) != singleObservableSize * observations_.size( ) )
            {
                throw std::runtime_error( "Error when setting weights in single observation set, sizes are incompatible." );
            }
        }
        else if( weightsVector.rows( ) > 0 )
        {
            throw std::runtime_error( "Error when setting weights in single observation set, observation set has no data." );
        }
        weightsVector_ = weightsVector;
    }

    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > createFilteredObservationSet(
        std::vector< int > indices )
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > filteredObservations = observations_;
        std::vector< TimeType > filteredObservationTimes = observationTimes_;
        std::vector< Eigen::VectorXd > filteredObservationsDependentVariables = observationsDependentVariables_;

        sort(indices.begin(), indices.end(), std::greater<int>());
        for( unsigned int i = 0; i < indices.size( ); i++ )
        {
            filteredObservationTimes.erase( filteredObservationTimes.begin( ) + indices.at( i ) );
            filteredObservations.erase( filteredObservations.begin( ) + indices.at( i ) );
            if( filteredObservationsDependentVariables.size( ) >  0 )
            {
                filteredObservationsDependentVariables.erase( filteredObservationsDependentVariables.begin( ) + indices.at( i ) );
            }
        }
        return std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
            observableType_, linkEnds_, filteredObservations, filteredObservationTimes,
            referenceLinkEnd_, filteredObservationsDependentVariables, dependentVariableCalculator_, ancilliarySettings_ );
    }


private:

    const ObservableType observableType_;

    const LinkDefinition linkEnds_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations_;

    std::vector< TimeType > observationTimes_;

    const LinkEndType referenceLinkEnd_;

    std::vector< Eigen::VectorXd > observationsDependentVariables_;

    const std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator > dependentVariableCalculator_;

    const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings_;

    const int numberOfObservations_;

    Eigen::VectorXd weightsVector_;

};


template< typename ObservationScalarType = double, typename TimeType = double,
    typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >
createResidualObservationSet(
    const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > observedObservationSet,
    const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > computedObservationSet )
{
    if( observedObservationSet->getObservableType( ) != computedObservationSet->getObservableType( ) )
    {
        throw std::runtime_error( "Error when computing residual observation set, observable is not equal" );
    }

    if( observedObservationSet->getReferenceLinkEnd( ) != computedObservationSet->getReferenceLinkEnd( ) )
    {
        throw std::runtime_error( "Error when computing residual observation set, reference link end is not equal" );
    }

    if( observedObservationSet->getLinkEnds( ).linkEnds_ != computedObservationSet->getLinkEnds( ).linkEnds_ )
    {
        throw std::runtime_error( "Error when computing residual observation set, link ends are not equal" );
    }

    if( observedObservationSet->getNumberOfObservables( ) != computedObservationSet->getNumberOfObservables( ) )
    {
        throw std::runtime_error( "Error when computing residual observation set, number of observable are not equal" );
    }

    //ESTIMATION-TODO: Add comparison of ancilliary settings

    std::vector< TimeType > observedTimes = observedObservationSet->getObservationTimes( );
    std::vector< TimeType > computedTimes = computedObservationSet->getObservationTimes( );

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observedData = observedObservationSet->getObservations( );
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > computedData = computedObservationSet->getObservations( );


    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;
    for( unsigned int i = 0; i < observedTimes.size( ); i++ )
    {
        if( observedTimes.at( i ) != computedTimes.at( i ) )
        {
            throw std::runtime_error( "Error when computing residual observation set, observation time of index " + std::to_string( i ) +
            " is not equal: " + std::to_string( static_cast< double >( observedTimes.at( i ) ) ) + ", " + std::to_string( static_cast< double >( computedTimes.at( i ) ) ) );
        }
        residuals.push_back( observedData.at( i ) - computedData.at( i ) );
    }

    return std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
        observedObservationSet->getObservableType( ),
        observedObservationSet->getLinkEnds( ),
        residuals,
        observedObservationSet->getObservationTimes( ),
        observedObservationSet->getReferenceLinkEnd( ),
        std::vector< Eigen::VectorXd >( ),
        nullptr,
        observedObservationSet->getAncilliarySettings( ) );
}

template< typename ObservationScalarType = double, typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > >
createSortedObservationSetList( const std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > observationSetList )
{
   std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > > sortedObservations;
   for( unsigned int i = 0; i < observationSetList.size( ); i++ )
   {

       sortedObservations[ observationSetList.at( i )->getObservableType( ) ][ observationSetList.at( i )->getLinkEnds( ).linkEnds_ ].push_back(
               observationSetList.at( i ) );
   }
   return sortedObservations;
}


template< typename ObservationScalarType = double, typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class ObservationCollection
{
public:

    typedef std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > > SortedObservationSets;

    ObservationCollection(
            const SortedObservationSets& observationSetList = SortedObservationSets( ) ):
        observationSetList_( observationSetList )
    {
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    ObservationCollection(
            const std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > observationSetList ):
        observationSetList_( createSortedObservationSetList< ObservationScalarType, TimeType >( observationSetList ) )
    {
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationVector( )
    {
        return concatenatedObservations_;
    }

    const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& getObservationVectorReference( )
    {
        return concatenatedObservations_;
    }

    std::vector< TimeType > getConcatenatedTimeVector( )
    {
        return concatenatedTimes_;
    }

    std::pair< TimeType, TimeType > getTimeBounds( )
    {
        return std::make_pair ( *std::min_element( concatenatedTimes_.begin( ), concatenatedTimes_.end( ) ),
                                *std::max_element( concatenatedTimes_.begin( ), concatenatedTimes_.end( ) ) );
    }

    std::vector< int > getConcatenatedLinkEndIds( )
    {
        return concatenatedLinkEndIds_;
    }

    std::map< observation_models::LinkEnds, int > getLinkEndIdentifierMap( )
    {
        return linkEndIds_;
    }

    std::map< int, observation_models::LinkEnds > getInverseLinkEndIdentifierMap( )
    {
        return inverseLinkEndIds_;
    }





    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > getObservationSetStartAndSize( )
    {
        return observationSetStartAndSize_;
    }

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > >& getObservationSetStartAndSizeReference( )
    {
        return observationSetStartAndSize_;
    }

    std::vector< std::pair< int, int > > getConcatenatedObservationSetStartAndSize( )
    {
        return concatenatedObservationSetStartAndSize_;
    }

    std::map< ObservableType, std::map< LinkEnds, std::pair< int, int > > > getObservationTypeAndLinkEndStartAndSize( )
    {
        return observationTypeAndLinkEndStartAndSize_;
    }

    std::map< ObservableType, std::map< int, std::vector< std::pair< int, int > > > > getObservationSetStartAndSizePerLinkEndIndex( )
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

    SortedObservationSets getObservations( )
    {
        return observationSetList_;
    }

    const SortedObservationSets& getObservationsReference ( ) const
    {
        return observationSetList_;
    }

    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > getSingleLinkAndTypeObservationSets(
            const ObservableType observableType,
            const LinkDefinition linkEnds )
    {
        if( observationSetList_.count( observableType ) == 0 )
        {
            throw std::runtime_error( "Error when retrieving observable of type " + observation_models::getObservableName( observableType ) +
                                      " from observation collection, no such observable exists" );
        }
        else if( observationSetList_.at( observableType ).count( linkEnds.linkEnds_ ) == 0 )
        {
            throw std::runtime_error( "Error when retrieving observable of type " + observation_models::getObservableName( observableType ) +
                                      " and link ends " + observation_models::getLinkEndsString( linkEnds.linkEnds_ ) +
                                      " from observation collection, no such link ends found for observable" );
        }

        return observationSetList_.at( observableType ).at( linkEnds.linkEnds_ );
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getSingleLinkObservations(
            const ObservableType observableType,
            const LinkDefinition& linkEnds )
    {
        if( observationSetStartAndSize_.count( observableType ) == 0 )
        {
            throw std::runtime_error( " Error when getting single link observations, not observations of type "
                                      + std::to_string( observableType ) );
        }
        else
        {
            if( observationSetStartAndSize_.at( observableType ).count( linkEnds.linkEnds_ ) == 0 )
            {
                throw std::runtime_error( " Error when getting single link observations, not observations of type "
                                          + std::to_string( observableType ) + " for given link ends." );
            }
            else
            {
                std::vector< std::pair< int, int > > combinedIndices =
                        observationSetStartAndSize_.at( observableType ).at( linkEnds.linkEnds_ );
                int startIndex = combinedIndices.at( 0 ).first;
                int finalEntry = combinedIndices.size( ) - 1;

                int numberOfObservables = ( combinedIndices.at( finalEntry ).first - startIndex ) +
                        combinedIndices.at( finalEntry ).second;
                return concatenatedObservations_.segment( startIndex, numberOfObservables );
            }
        }
    }

    std::vector< TimeType > getSingleLinkTimes(
            const ObservableType observableType,
            const LinkDefinition& linkEnds )
    {
        if( observationSetStartAndSize_.count( observableType ) == 0 )
        {
            throw std::runtime_error( " Error when getting single link observations, not observations of type "
                                      + std::to_string( observableType ) );
        }
        else
        {
            if( observationSetStartAndSize_.at( observableType ).count( linkEnds.linkEnds_ ) == 0 )
            {
                throw std::runtime_error( " Error when getting single link observations, not observations of type "
                                          + std::to_string( observableType ) + " for given link ends." );
            }
            else
            {
                std::vector< std::pair< int, int > > combinedIndices =
                        observationSetStartAndSize_.at( observableType ).at( linkEnds.linkEnds_ );
                int startIndex = combinedIndices.at( 0 ).first;
                int finalEntry = combinedIndices.size( ) - 1;

                int numberOfObservables = ( combinedIndices.at( finalEntry ).first - startIndex ) +
                        combinedIndices.at( finalEntry ).second;
                return std::vector< TimeType >( concatenatedTimes_.begin( ) + startIndex,
                                                concatenatedTimes_.begin( ) + ( startIndex + numberOfObservables ) );
            }
        }
    }

    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > getSingleLinkObservationsAndTimes(
            const ObservableType observableType,
            const LinkDefinition& linkEnds )
    {
        return std::make_pair( getSingleLinkObservations( observableType, linkEnds ),
                               getSingleLinkTimes( observableType, linkEnds ) );
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


    std::map< ObservableType, std::map< int, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > > getSortedObservationSets( )
    {
        std::map< ObservableType, std::map< int, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > > observationSetListIndexSorted;
        for( auto it1 : observationSetList_ )
        {
            for( auto it2 : it1.second )
            {
                observationSetListIndexSorted[ it1.first ][ linkEndIds_[ it2.first ] ] = it2.second;
            }
        }
        return observationSetListIndexSorted;
    }

    std::map < ObservableType, std::vector< LinkEnds > > getLinkEndsPerObservableType( )
    {
        std::map < ObservableType, std::vector< LinkEnds > > linkEndsPerObservableType;

        for ( auto observableTypeIt = observationSetList_.begin( ); observableTypeIt != observationSetList_.end( );
            ++observableTypeIt )
        {
            std::vector< LinkEnds > linkEndsVector;

            for ( auto linkEndsIt = observableTypeIt->second.begin( ); linkEndsIt != observableTypeIt->second.end( );
                ++linkEndsIt )
            {
                linkEndsVector.push_back( linkEndsIt->first );
            }

            linkEndsPerObservableType[ observableTypeIt->first ] = linkEndsVector;
        }

        return linkEndsPerObservableType;
    }

    Eigen::VectorXd getWeightsFromSingleObservationSets( )
    {
        Eigen::VectorXd weightsVector = Eigen::VectorXd::Zero( totalObservableSize_ );

        for( auto observationIterator : observationSetList_ )
        {
            ObservableType currentObservableType = observationIterator.first;

            for( auto linkEndIterator : observationIterator.second )
            {
                LinkEnds currentLinkEnds = linkEndIterator.first;
                for( unsigned int i = 0; i < linkEndIterator.second.size( ); i++ )
                {
                    std::pair< int, int > startAndSize = observationSetStartAndSize_.at( currentObservableType ).at( currentLinkEnds ).at( i );
                    if( observationSetList_.at( currentObservableType ).at( currentLinkEnds ).at( i )->getWeightsVectorReference( ).rows( ) ==
                        startAndSize.second )
                    {
                        weightsVector.segment( startAndSize.first, startAndSize.second ) =
                            observationSetList_.at( currentObservableType ).at( currentLinkEnds ).at(
                                i )->getWeightsVector( );
                    }
                    else
                    {
                        throw std::runtime_error( "Error when compiling full weights vector from single observation set, sizes are inconsistent" );
                    }
                }
            }
        }
        return weightsVector;
    }


private:

    void setObservationSetIndices( )
    {
        int currentStartIndex = 0;
        int currentTypeStartIndex = 0;
        totalNumberOfObservables_ = 0;
        totalObservableSize_ = 0;

        for( auto observationIterator : observationSetList_ )
        {
            ObservableType currentObservableType = observationIterator.first;

            currentTypeStartIndex = currentStartIndex;
            int observableSize = getObservableSize(  currentObservableType );

            int currentObservableTypeSize = 0;

            for( auto linkEndIterator : observationIterator.second )
            {
                LinkEnds currentLinkEnds = linkEndIterator.first;
                int currentLinkEndStartIndex = currentStartIndex;
                int currentLinkEndSize = 0;
                for( unsigned int i = 0; i < linkEndIterator.second.size( ); i++ )
                {
                    int currentNumberOfObservables = linkEndIterator.second.at( i )->getNumberOfObservables( );
                    int currentObservableVectorSize = currentNumberOfObservables * observableSize;

                    observationSetStartAndSize_[ currentObservableType ][ currentLinkEnds ].push_back(
                                std::make_pair( currentStartIndex, currentObservableVectorSize ) );
                    concatenatedObservationSetStartAndSize_.push_back( std::make_pair( currentStartIndex, currentObservableVectorSize ) );
                    currentStartIndex += currentObservableVectorSize;
                    currentObservableTypeSize += currentObservableVectorSize;
                    currentLinkEndSize += currentObservableVectorSize;

                    totalObservableSize_ += currentObservableVectorSize;
                    totalNumberOfObservables_ += currentNumberOfObservables;
                }
                observationTypeAndLinkEndStartAndSize_[ currentObservableType ][ currentLinkEnds ] = std::make_pair(
                    currentLinkEndStartIndex, currentLinkEndSize );
            }
            observationTypeStartAndSize_[ currentObservableType ] = std::make_pair(
                        currentTypeStartIndex, currentObservableTypeSize );
        }
    }

    void setConcatenatedObservationsAndTimes( )
    {
        concatenatedObservations_ = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObservableSize_ );
        concatenatedTimes_.resize( totalObservableSize_ );
        concatenatedLinkEndIds_.resize( totalObservableSize_ );
        concatenatedLinkEndIdNames_.resize( totalObservableSize_ );


        int observationCounter = 0;
        int maximumStationId = 0;
        int currentStationId;

        for( auto observationIterator : observationSetList_ )
        {
            ObservableType currentObservableType = observationIterator.first;
            int observableSize = getObservableSize(  currentObservableType );

            for( auto linkEndIterator : observationIterator.second )
            {
                LinkEnds currentLinkEnds = linkEndIterator.first;
                LinkDefinition firstLinkDefinition;

                if( linkEndIterator.second.size( ) > 0 )
                {
                    firstLinkDefinition = linkEndIterator.second.at( 0 )->getLinkEnds( );
                    linkDefinitionsPerObservable_[ currentObservableType].push_back( firstLinkDefinition );
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
                    LinkDefinition currentLinkDefinition = linkEndIterator.second.at( i )->getLinkEnds( );
                    if( !( currentLinkDefinition == firstLinkDefinition ) )
                    {
                        throw std::runtime_error( "Error when creating ObservationCollection, link definitions of same link ends are not equal " );
                    }

                    std::pair< int, int > startAndSize =
                            observationSetStartAndSize_.at( currentObservableType ).at( currentLinkEnds ).at( i );
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > currentObservables =
                            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( startAndSize.second );

                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > currentObservationSet =
                            linkEndIterator.second.at( i )->getObservations( );
                    std::vector< TimeType > currentObservationTimes =
                            linkEndIterator.second.at( i )->getObservationTimes( );
                    for( unsigned int j = 0; j < currentObservationSet.size( ); j++ )
                    {
                        currentObservables.segment( j * observableSize, observableSize ) = currentObservationSet.at( j );
                        for( int k = 0; k < observableSize; k++ )
                        {
                            concatenatedTimes_[ observationCounter ] = currentObservationTimes.at( j );
                            concatenatedLinkEndIds_[ observationCounter ] = currentStationId;
                            concatenatedLinkEndIdNames_[ observationCounter ] = currentLinkEnds;
                            observationCounter++;
                        }
                    }
                    concatenatedObservations_.segment( startAndSize.first, startAndSize.second ) = currentObservables;
                }
            }
        }

        for( auto it1 : observationSetStartAndSize_ )
        {
            for( auto it2 : it1.second )
            {
                observationSetStartAndSizePerLinkEndIndex_[ it1.first ][ linkEndIds_[ it2.first ] ] = it2.second;
            }
        }
    }

    const SortedObservationSets observationSetList_;

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > concatenatedObservations_;

    std::vector< TimeType > concatenatedTimes_;

    std::vector< int > concatenatedLinkEndIds_;

    std::vector< LinkEnds > concatenatedLinkEndIdNames_;

    std::map< ObservableType, std::vector< LinkDefinition > > linkDefinitionsPerObservable_;

    std::map< observation_models::LinkEnds, int > linkEndIds_;

    std::map< int, observation_models::LinkEnds > inverseLinkEndIds_;

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > observationSetStartAndSize_;

    std::vector< std::pair< int, int > > concatenatedObservationSetStartAndSize_;

    std::map< ObservableType, std::map< int, std::vector< std::pair< int, int > > > > observationSetStartAndSizePerLinkEndIndex_;

    std::map< ObservableType, std::map< LinkEnds, std::pair< int, int > > > observationTypeAndLinkEndStartAndSize_;

    std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize_;

    int totalObservableSize_;

    int totalNumberOfObservables_;
};


template< typename ObservationScalarType = double, typename TimeType = double,
    typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > splitSingleObservationSetIntoArcs(
    const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > originalObservationSet,
    const double arcSplitInterval,
    const int minimumNumberOfObservations )
{
    std::vector< int > rawArcStartIndices = { 0 };
    const std::vector< TimeType > originalObservationTimes = originalObservationSet->getObservationTimesReference( );
    for( unsigned int i = 1; i < originalObservationTimes.size( ); i++ )
    {
        if( ( originalObservationTimes.at( i ) - originalObservationTimes.at( i - 1 ) ) > arcSplitInterval )
        {
            rawArcStartIndices.push_back( i );
        }
    }
    rawArcStartIndices.push_back( originalObservationTimes.size( ) );

    std::vector< std::pair< int, int > > arcSplitIndices;

    for( unsigned int j = 1; j < rawArcStartIndices.size( ); j++ )
    {
        if( ( rawArcStartIndices.at( j ) - rawArcStartIndices.at( j - 1 ) ) > minimumNumberOfObservations )
        {
            arcSplitIndices.push_back( std::make_pair( rawArcStartIndices.at( j - 1 ), rawArcStartIndices.at( j ) - rawArcStartIndices.at( j - 1 ) ) );
        }
    }

    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > splitSingleObervationSet;
    for( unsigned int i = 0; i < arcSplitIndices.size( ); i++ )
    {
        std::vector< Eigen::VectorXd > currentSplitDependentVariables;
        if( originalObservationSet->getObservationsDependentVariablesReference( ).size( ) > 0 )
        {
            currentSplitDependentVariables =
                utilities::getStlVectorSegment( originalObservationSet->getObservationsDependentVariablesReference( ),
                                  arcSplitIndices.at( i ).first, arcSplitIndices.at( i ).second );
        }

        splitSingleObervationSet.push_back(
            std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                originalObservationSet->getObservableType( ),
                originalObservationSet->getLinkEnds( ),
                utilities::getStlVectorSegment( originalObservationSet->getObservationsReference( ), arcSplitIndices.at( i ).first, arcSplitIndices.at( i ).second ),
                utilities::getStlVectorSegment( originalObservationSet->getObservationTimesReference( ), arcSplitIndices.at( i ).first, arcSplitIndices.at( i ).second ),
                originalObservationSet->getReferenceLinkEnd( ),
                currentSplitDependentVariables,
                originalObservationSet->getDependentVariableCalculator( ),
                originalObservationSet->getAncilliarySettings( ) ) );
    }
    return splitSingleObervationSet;
}


template< typename ObservationScalarType = double, typename TimeType = double,
    typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > splitObservationSetsIntoArcs(
    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > originalObservationCollection,
    const double arcSplitInterval,
    const int minimumNumberOfObservations )
{
    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets originalObservationSets =
        originalObservationCollection->getObservations( );
    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets splitObservationSets;
    for( auto observationIt : originalObservationSets )
    {
        for( auto linkEndIt : observationIt.second )
        {
            std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > splitSingleObervationSet;
            for( unsigned int i = 0; i < linkEndIt.second.size( ); i++ )
            {
                auto singleSplitObservationSet = splitSingleObservationSetIntoArcs( linkEndIt.second.at( i ), arcSplitInterval, minimumNumberOfObservations );
                splitSingleObervationSet.insert( splitSingleObervationSet.end( ), singleSplitObservationSet.begin( ), singleSplitObservationSet.end( ) );
            }
            splitObservationSets[ observationIt.first ][ linkEndIt.first ] = splitSingleObervationSet;
        }
    }
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( splitObservationSets );
}


template< typename ObservationScalarType = double, typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > getObservationListWithDependentVariables(
        const std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > fullObservationList,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > dependentVariableToRetrieve )
{
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > observationList;
    for( unsigned int i = 0; i < fullObservationList.size( ); i++ )
    {
        std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator > dependentVariableCalculator =
                fullObservationList.at( i )->getDependentVariableCalculator( );
        if( dependentVariableCalculator != nullptr )
        {
            std::pair< int, int > variableIndices = dependentVariableCalculator->getDependentVariableIndices(
                        dependentVariableToRetrieve );

            if( variableIndices.second != 0 )
            {
                observationList.push_back( fullObservationList.at( i ) );
            }
        }
    }
    return observationList;
}

template< typename ObservationScalarType = double, typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > getObservationListWithDependentVariables(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > dependentVariableToRetrieve,
        const ObservableType observableType )
{
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > observationList;
    if( observationCollection->getObservations( ).count( observableType ) != 0 )
    {
        std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > >
                observationsOfGivenType = observationCollection->getObservations( ).at( observableType );

        for( auto linkEndIterator : observationsOfGivenType )
        {
            std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > fullObservationList =
                    linkEndIterator.second;
            std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > partialObservationList =
                    getObservationListWithDependentVariables( fullObservationList, dependentVariableToRetrieve );
            observationList.insert( observationList.end( ), partialObservationList.begin( ), partialObservationList.end( ) );
        }
    }
    return observationList;
}

template< typename ObservationScalarType = double, typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > getObservationListWithDependentVariables(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > dependentVariableToRetrieve,
        const ObservableType observableType ,
        const LinkEnds& linkEnds )
{
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > observationList;
    if( observationCollection->getObservations( ).count( observableType ) != 0 )
    {
        if( observationCollection->getObservations( ).at( observableType ).count( linkEnds ) != 0 )
        {
            std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > fullObservationList =
                    observationCollection->getObservations( ).at( observableType ).at( linkEnds );
            std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > partialObservationList =
                    getObservationListWithDependentVariables( fullObservationList, dependentVariableToRetrieve );
            observationList.insert( observationList.end( ), partialObservationList.begin( ), partialObservationList.end( ) );
        }
    }
    return observationList;
}

template< typename ObservationScalarType = double, typename TimeType = double, typename... ArgTypes,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::vector< std::map< double, Eigen::VectorXd > > getDependentVariableResultPerObservationSet(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > dependentVariableToRetrieve,
        ArgTypes... args )
{
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > observationsWithVariable =
            getObservationListWithDependentVariables(
                observationCollection, dependentVariableToRetrieve, args ... );

    std::vector< std::map< double, Eigen::VectorXd > > dependentVariableList;
    for( unsigned int i = 0; i < observationsWithVariable.size( ); i++ )
    {
        std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator > dependentVariableCalculator =
                observationsWithVariable.at( i )->getDependentVariableCalculator( );

        std::pair< int, int > variableIndices = dependentVariableCalculator->getDependentVariableIndices(
                    dependentVariableToRetrieve );
        std::map< double, Eigen::VectorXd > slicedHistory =
                utilities::sliceMatrixHistory(
                    observationsWithVariable.at( i )->getDependentVariableHistory( ), variableIndices );

        dependentVariableList.push_back( slicedHistory );
    }

    return dependentVariableList;
}

template< typename ObservationScalarType = double, typename TimeType = double, typename... ArgTypes,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::map< double, Eigen::VectorXd > getDependentVariableResultList(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > dependentVariableToRetrieve,
        ArgTypes... args )
{
    std::vector< std::map< double, Eigen::VectorXd > > dependentVariableResultPerObservationSet =
            getDependentVariableResultPerObservationSet(
                observationCollection, dependentVariableToRetrieve, args ... );
    return utilities::concatenateMaps( dependentVariableResultPerObservationSet );


}

template< typename ObservationScalarType = double, typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
inline std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > createSingleObservationSet(
        const ObservableType observableType,
        const LinkEnds& linkEnds,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&  observations,
        const std::vector< TimeType > observationTimes,
        const LinkEndType referenceLinkEnd,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings )
{
    return std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                observableType, linkEnds, observations, observationTimes, referenceLinkEnd,
                std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ), nullptr, ancilliarySettings );
}

//template< typename ObservationScalarType = double, typename TimeType = double,
//          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
//inline std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >  createManualObservationCollection(
//        const std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > >& observationSetList )
//{
//    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( observationSetList );
//}

//template< typename ObservationScalarType = double, typename TimeType = double,
//          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
//inline std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >  createManualObservationCollection(
//        const ObservableType observableType,
//        const LinkEnds& linkEnds,
//        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >& observationSetList )
//{
//    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( observationSetList );
//}

template< typename ObservationScalarType = double, typename TimeType = double >
inline std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >  createManualObservationCollection(
        const ObservableType observableType,
        const LinkDefinition& linkEnds,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType > observationTimes,
        const LinkEndType referenceLinkEnd,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr )
{
    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > singleObservationSet =
            createSingleObservationSet( observableType, linkEnds.linkEnds_, observations, observationTimes, referenceLinkEnd,
                                        ancilliarySettings );

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > > observationSetList;
    observationSetList[ observableType ][ linkEnds.linkEnds_ ].push_back( singleObservationSet );
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( observationSetList );
}

template< typename ObservationScalarType = double, typename TimeType = double >
inline std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >  createManualObservationCollection(
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObservationSets )
{
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( singleObservationSets );
}


template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >  createResidualCollection(
    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observedData,
    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > computedData )
{
//    std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > >
    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets observedObservationSets = observedData->getObservations( );
    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets computedObservationSets = computedData->getObservations( );
    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets residualObservationSets;

    for( auto observationIt : observedObservationSets )
    {
        for( auto linkEndIt : observationIt.second )
        {
            for( unsigned int i = 0; i < linkEndIt.second.size( ); i++ )
            {
                residualObservationSets[ observationIt.first ][ linkEndIt.first ].push_back(
                    createResidualObservationSet( linkEndIt.second.at( i ), computedObservationSets.at( observationIt.first ).at( linkEndIt.first ).at( i ) ) );
            }
        }
    }
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( residualObservationSets );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::map< observation_models::ObservableType, std::vector< std::pair< LinkEnds, std::vector< std::vector< int > > > > >
    getObservationCollectionEntriesToFiler(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > dataToFiler,
        const Eigen::VectorXd& residualVector,
        const std::map< ObservableType, double > residualCutoffValuePerObservable )
{
    // Check if input data is compatible
    if( residualVector.rows( ) != dataToFiler->getTotalObservableSize( ) )
    {
        throw std::runtime_error( "Error when filtering observations, input size is incompatible" );
    }

    // Retrieve observations to filter
    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets observationSetsToFilter = dataToFiler->getObservations( );

    // Create data structure with filtered results
    std::map< observation_models::ObservableType, std::vector< std::pair< LinkEnds, std::vector< std::vector< int > > > > > filterEntries;

    // Iterate over all observable types
    for( auto observationIt : observationSetsToFilter )
    {
        observation_models::ObservableType observableType = observationIt.first;

        // Get filter value for current observable
        double filterValue = residualCutoffValuePerObservable.at( observationIt.first );

        // Iterate over all link ends
        std::vector< std::pair< LinkEnds, std::vector< std::vector< int > > > > currentObservableEntriesToFilter;
        for ( auto linkEndIt: observationIt.second )
        {
            observation_models::LinkEnds linkEnds = linkEndIt.first;

            // Iterate over all observations with current link ends and type
            std::vector< std::vector< int > > currentLinkEndsEntriesToFiler;
            for ( unsigned int i = 0; i < linkEndIt.second.size( ); i++ )
            {
                std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > currentObservationSet =
                    observationSetsToFilter.at( observationIt.first ).at( linkEndIt.first ).at( i );

                // Get residuals for current set
                std::pair< int, int > fullVectorStartAndSize =
                    dataToFiler->getObservationSetStartAndSizeReference( ).at( observableType ).at( linkEnds ).at( i );
                Eigen::VectorXd currentSetResiduals = residualVector.segment( fullVectorStartAndSize.first, fullVectorStartAndSize.second );

                int currentObservableSize = getObservableSize( observableType );

                // Check if data is compatible
                if( currentSetResiduals.rows( ) != currentObservableSize * currentObservationSet->getNumberOfObservables( ) )
                {
                    throw std::runtime_error( "Error when filtering observations, input size of single observation set for " +
                        getObservableName( observableType ) +", " + getLinkEndsString( linkEnds ) + ", set " +
                        std::to_string( i ) + " is incompatible" );
                }

                std::vector< int > indicesToRemove;
                for( unsigned int j = 0; j < currentObservationSet->getObservationTimes( ).size( ); j++ )
                {
                    bool removeObservation = false;
                    for( int k = 0; k < currentObservableSize; k++ )
                    {
                        if( currentSetResiduals( j * currentObservableSize + k ) > filterValue )
                        {
                            removeObservation = true;
                        }
                    }
                    if( removeObservation )
                    {
                        indicesToRemove.push_back( j );
                    }
                }
                currentLinkEndsEntriesToFiler.push_back( indicesToRemove );
            }
            currentObservableEntriesToFilter.push_back( std::make_pair( linkEnds, currentLinkEndsEntriesToFiler ) );
        }
        filterEntries[ observableType ] = currentObservableEntriesToFilter;
    }
    return filterEntries;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > filterData(
    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
    const std::map< observation_models::ObservableType, std::vector< std::pair< LinkEnds, std::vector< std::vector< int > > > > >& filterEntries )
{
    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets filteredObservedObservationSets;
    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > filteredObservationCollection;

    for( auto observableIt : filterEntries )
    {
        ObservableType observableType = observableIt.first;
        for( unsigned int i = 0; i < observableIt.second.size( ); i++ )
        {
            LinkEnds currentLinkEnds = observableIt.second.at( i ).first;
            std::vector< std::vector< int > > linkEndListEntriesToRemove = observableIt.second.at( i ).second;

            std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > observedSets =
                observationCollection->getObservationsReference( ).at( observableType ).at( currentLinkEnds );
            if( observedSets.size( ) != linkEndListEntriesToRemove.size( ) )
            {
                throw std::runtime_error( "Error when filtering observations, number of observation sets and filter list for " +
                                          getObservableName( observableType ) +", " + getLinkEndsString( currentLinkEnds ) + " is incompatible" );
            }
            for( unsigned int j = 0; j < observedSets.size( ); j++ )
            {
                filteredObservedObservationSets[ observableType ][ currentLinkEnds ].push_back(
                    observedSets.at( j )->createFilteredObservationSet(
                        linkEndListEntriesToRemove.at( j ) ) );
            }
        }
    }
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( filteredObservedObservationSets );
}

template< typename ObservationScalarType = double, typename TimeType = double >
void filterObservedAndComputedData(
    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observedDataCollection,
    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > computedDataCollection,
    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >& observedFilteredDataCollection,
    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >& computedFilteredDataCollection,
    const std::map< ObservableType, double > residualCutoffValuePerObservable )
{
    Eigen::VectorXd residualVector =
        ( observedDataCollection->getObservationVector( ) - computedDataCollection->getObservationVector( ) ).template cast< double >( );

    std::map< observation_models::ObservableType, std::vector< std::pair< LinkEnds, std::vector< std::vector< int > > > > > filterEntries =
        getObservationCollectionEntriesToFiler( observedDataCollection, residualVector, residualCutoffValuePerObservable );

    observedFilteredDataCollection = filterData( observedDataCollection, filterEntries );
    computedFilteredDataCollection = filterData( computedDataCollection, filterEntries );
}

//////////////////// DEPRECATEC

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > filterResidualOutliers(
    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observedData,
    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > residualData,
    const std::map< ObservableType, double > residualCutoffValuePerObservable )
{
    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets observedObservationSets = observedData->getObservations( );
    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets filteredObservedObservationSets;
    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets residualObservationSets = residualData->getObservations( );

    for( auto observationIt : residualObservationSets )
    {
        if( residualObservationSets.count( observationIt.first ) )
        {
            double filterValue = residualCutoffValuePerObservable.at( observationIt.first );
            for ( auto linkEndIt: observationIt.second )
            {
                for ( unsigned int i = 0; i < linkEndIt.second.size( ); i++ )
                {
                    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > currentObservationSet =
                        residualObservationSets.at( observationIt.first ).at( linkEndIt.first ).at( i );
                    std::vector<int> indicesToRemove;
                    std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> >
                        residualObservationSet = linkEndIt.second.at( i );
                    for ( unsigned int j = 0; j < currentObservationSet->getObservationTimes( ).size( ); j++ )
                    {
                        if ( currentObservationSet->getObservation( j ).array( ).abs( ).maxCoeff( ) > filterValue )
                        {
                            indicesToRemove.push_back( j );
                        }
                    }

                    filteredObservedObservationSets[ observationIt.first ][ linkEndIt.first ].push_back(
                        observedObservationSets[ observationIt.first ][ linkEndIt.first ].at( i )->createFilteredObservationSet(
                            indicesToRemove ) );
                }
            }
        }
    }
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( filteredObservedObservationSets );
}

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_OBSERVATIONS_H
