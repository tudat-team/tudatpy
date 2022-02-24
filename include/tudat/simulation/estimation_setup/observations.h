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

#include <boost/bind.hpp>
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
            const LinkEnds& linkEnds,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType > observationTimes,
            const LinkEndType referenceLinkEnd,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observationsDependentVariables =
            std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
            const std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator > dependentVariableCalculator = nullptr ):
        observableType_( observableType ),
        linkEnds_( linkEnds ),
        observations_( observations ),
        observationTimes_( observationTimes ),
        referenceLinkEnd_( referenceLinkEnd ),
        observationsDependentVariables_( observationsDependentVariables ),
        dependentVariableCalculator_( dependentVariableCalculator ),
        numberOfObservations_( observations_.size( ) )
    {
        if( observations_.size( ) != observationTimes_.size( ) )
        {
            throw std::runtime_error( "Error when making SingleObservationSet, input sizes are inconsistent." );
        }

        for( unsigned int i = 1; i < observations.size( ); i++ )
        {
            if( observations.at( i ).rows( ) != observations.at( i - 1 ).rows( ) )
            {
                throw std::runtime_error( "Error when making SingleObservationSet, input observables not of consistent size." );
            }
        }
    }

    ObservableType getObservableType( )
    {
        return observableType_;
    }

    LinkEnds getLinkEnds( )
    {
        return linkEnds_;
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservations( )
    {
        return observations_;
    }

    std::vector< TimeType > getObservationTimes( )
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


    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservationsDependentVariables( )
    {
        return observationsDependentVariables_;
    }


    std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator > getDependentVariableCalculator( )
    {
        return dependentVariableCalculator_;
    }

    std::map< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getDependentVariableHistory( )
    {
        return utilities::createMapFromVectors< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >(
                    observationTimes_, observationsDependentVariables_ );
    }



private:

    const ObservableType observableType_;

    const LinkEnds linkEnds_;

    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations_;

    const std::vector< TimeType > observationTimes_;

    const LinkEndType referenceLinkEnd_;

    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observationsDependentVariables_;

    const std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator > dependentVariableCalculator_;

    const int numberOfObservations_;

};


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

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationVector( )
    {
        return concatenatedObservations_;
    }

    std::vector< TimeType > getConcatenatedTimeVector( )
    {
        return concatenatedTimes_;
    }

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > getObservationSetStartAndSize( )
    {
        return observationSetStartAndSize_;
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

    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > getSingleLinkAndTypeObservationSets(
            const ObservableType observableType,
            const LinkEnds linkEnds )
    {
        return observationSetList_.at( observableType ).at( linkEnds );
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getSingleLinkObservations(
            const ObservableType observableType,
            const LinkEnds& linkEnds )
    {
        if( observationSetStartAndSize_.count( observableType ) == 0 )
        {
            throw std::runtime_error( " Error when getting single link observations, not observations of type "
                                      + std::to_string( observableType ) );
        }
        else
        {
            if( observationSetStartAndSize_.at( observableType ).count( linkEnds ) == 0 )
            {
                throw std::runtime_error( " Error when getting single link observations, not observations of type "
                                          + std::to_string( observableType ) + " for given link ends." );
            }
            else
            {
                std::vector< std::pair< int, int > > combinedIndices =
                        observationSetStartAndSize_.at( observableType ).at( linkEnds );
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
            const LinkEnds& linkEnds )
    {
        if( observationSetStartAndSize_.count( observableType ) == 0 )
        {
            throw std::runtime_error( " Error when getting single link observations, not observations of type "
                                      + std::to_string( observableType ) );
        }
        else
        {
            if( observationSetStartAndSize_.at( observableType ).count( linkEnds ) == 0 )
            {
                throw std::runtime_error( " Error when getting single link observations, not observations of type "
                                          + std::to_string( observableType ) + " for given link ends." );
            }
            else
            {
                std::vector< std::pair< int, int > > combinedIndices =
                        observationSetStartAndSize_.at( observableType ).at( linkEnds );
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
            const LinkEnds& linkEnds )
    {
        return std::make_pair( getSingleLinkObservations( observableType, linkEnds ),
                               getSingleLinkTimes( observableType, linkEnds ) );
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
                for( unsigned int i = 0; i < linkEndIterator.second.size( ); i++ )
                {
                    int currentNumberOfObservables = linkEndIterator.second.at( i )->getNumberOfObservables( );
                    int currentObservableVectorSize = currentNumberOfObservables * observableSize;

                    observationSetStartAndSize_[ currentObservableType ][ currentLinkEnds ].push_back(
                                std::make_pair( currentStartIndex, currentObservableVectorSize ) );
                    currentStartIndex += currentObservableVectorSize;
                    currentObservableTypeSize += currentObservableVectorSize;

                    totalObservableSize_ += currentObservableVectorSize;
                    totalNumberOfObservables_ += currentNumberOfObservables;
                }
            }
            observationTypeStartAndSize_[ currentObservableType ] = std::make_pair(
                        currentTypeStartIndex, currentObservableTypeSize );
        }
    }

    void setConcatenatedObservationsAndTimes( )
    {
        concatenatedObservations_ = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObservableSize_ );
        concatenatedTimes_.resize( totalObservableSize_ );

        int observationCounter = 0;

        for( auto observationIterator : observationSetList_ )
        {
            ObservableType currentObservableType = observationIterator.first;
            int observableSize = getObservableSize(  currentObservableType );

            for( auto linkEndIterator : observationIterator.second )
            {
                LinkEnds currentLinkEnds = linkEndIterator.first;
                for( unsigned int i = 0; i < linkEndIterator.second.size( ); i++ )
                {
                    std::pair< int, int > startAndSize =
                            observationSetStartAndSize_.at( currentObservableType ).at( currentLinkEnds ).at( i );
                    Eigen::VectorXd currentObservables = Eigen::VectorXd::Zero( startAndSize.second );

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
                            observationCounter++;
                        }
                    }
                    concatenatedObservations_.segment( startAndSize.first, startAndSize.second ) = currentObservables;
                }
            }
        }
    }

    const SortedObservationSets observationSetList_;

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > concatenatedObservations_;

    std::vector< TimeType > concatenatedTimes_;


    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > observationSetStartAndSize_;

    std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize_;

    int totalObservableSize_;

    int totalNumberOfObservables_;

};

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
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType > observationTimes,
        const LinkEndType referenceLinkEnd )
{
    return std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                observableType, linkEnds, observations, observationTimes, referenceLinkEnd );
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
        const LinkEnds& linkEnds,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType > observationTimes,
        const LinkEndType referenceLinkEnd )
{
    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > singleObservationSet =
            createSingleObservationSet( observableType, linkEnds, observations, observationTimes, referenceLinkEnd );

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > > observationSetList;
    observationSetList[ observableType ][ linkEnds ].push_back( singleObservationSet );
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( observationSetList );
}



} // namespace observation_models

} // namespace tudat

#endif // TUDAT_OBSERVATIONS_H
