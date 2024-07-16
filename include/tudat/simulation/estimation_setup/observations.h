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

        singleObservationSize_ = 0;
        if ( numberOfObservations_ > 0 )
        {
            singleObservationSize_ = observations_.at( 0 ).rows( );
        }

        // Initialise weights
        for ( unsigned int k = 0 ; k < numberOfObservations_ ; k++ )
        {
            weights_.push_back( Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( singleObservationSize_, 1 ) );
        }

        // Initialise time bounds
       updateTimeBounds( );

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

    void setObservations( std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations )
    {
        if ( observations.size( ) != observations_.size( ) )
        {
            throw std::runtime_error( "Error when resetting observations, number of observations is incompatible." );
        }
        observations_ = observations;
        if ( !observations_.empty( ) )
        {
            singleObservationSize_ = observations_.at( 0 ).size( );
        }
    }

    void setResiduals( std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals )
    {
        if ( residuals.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Error when setting residuals, number of observations is inconsistent." );
        }
        residuals_ = residuals;
    }

    void setResiduals( Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualsVector )
    {
        if ( residualsVector.size( ) != numberOfObservations_ * singleObservationSize_ )
        {
            throw std::runtime_error( "Error when setting residuals, number of observations is inconsistent." );
        }

        residuals_.clear( );
        for ( unsigned int k = 0 ; k < numberOfObservations_ ; k++ )
        {
            residuals_.push_back( residualsVector.segment( k * singleObservationSize_, singleObservationSize_ ) );
        }
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

    void setObservation( const int index,
                         Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observation )
    {
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when setting single observation value, index is out of bounds" );
        }
        if ( observation.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error( "Error when setting single observation value, the observation size is inconsistent." );
        }
        observations_.at( index ) = observation;
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

    unsigned int getSingleObservableSize( ) const
    {
        return singleObservationSize_;
    }

    unsigned int getTotalObservationSetSize( ) const
    {
        return numberOfObservations_ * singleObservationSize_;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationsVector( )
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_ * numberOfObservations_ );
        for( unsigned int i = 0; i < observations_.size( ); i++ )
        {
            observationsVector.segment( i * singleObservationSize_, singleObservationSize_ ) = observations_.at( i );
        }
        return observationsVector;
    }


    std::pair< TimeType, TimeType > getTimeBounds( )
    {
        if( observationTimes_.size( ) == 0 )
        {
            throw std::runtime_error( "Error when getting time bounds of observation set, no observations found" );
        }
        return timeBounds_; /*std::make_pair ( *std::min_element( observationTimes_.begin( ), observationTimes_.end( ) ),
                                *std::max_element( observationTimes_.begin( ), observationTimes_.end( ) ) );*/
    }

    void updateTimeBounds( )
    {
        if( observationTimes_.size( ) == 0 )
        {
            timeBounds_ = std::make_pair( TUDAT_NAN, TUDAT_NAN );
        }
        else
        {
            timeBounds_ = std::make_pair ( *std::min_element( observationTimes_.begin( ), observationTimes_.end( ) ),
                                           *std::max_element( observationTimes_.begin( ), observationTimes_.end( ) ) );
        }
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

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getWeights( ) const
    {
        return weights_;
    }

    Eigen::VectorXd getWeightsVector( ) const
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > weightsVector =
                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( singleObservationSize_ * numberOfObservations_, 1 );
        for( unsigned int i = 0; i < numberOfObservations_ ; i++ )
        {
            weightsVector.block( i * singleObservationSize_, 0, singleObservationSize_, 1 ) = weights_.at( i );
        }
        return weightsVector;
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getResiduals( ) const
    {
        return residuals_;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidualsVector( ) const
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualsVector =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_ * numberOfObservations_, 1 );
        for( unsigned int i = 0; i < numberOfObservations_ ; i++ )
        {
            residualsVector.block( i * singleObservationSize_, 0, singleObservationSize_, 1 ) = residuals_.at( i );
        }
        return residualsVector;
    }

//    Eigen::VectorXd& getWeightsVectorReference( )
//    {
//        return weightsVector_;
//    }

//    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& getWeightsReference( )
//    {
//        return weights_;
//    }

    void setConstantWeight( const double weight )
    {
        for ( unsigned int k = 0 ; k < numberOfObservations_ ; k++ )
        {
            weights_.at( k ) = weight * Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( singleObservationSize_, 1 );
        }
    }

    void setConstantWeight( const Eigen::Matrix< double, Eigen::Dynamic, 1 >& weight )
    {
        if ( weight.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error( "Error when setting constant weight in single observation set, weight size is inconsistent with single observation size." );
        }
        for ( unsigned int k = 0 ; k < weights_.size( ) ; k++ )
        {
            weights_.at( k ) = weight;
        }
    }

    void setTabulatedWeights( const Eigen::VectorXd& weightsVector )
    {
        if( weightsVector.rows( ) != singleObservationSize_ * observations_.size( ) )
        {
            throw std::runtime_error( "Error when setting weights in single observation set, sizes are incompatible." );
        }
//        else if( weightsVector.rows( ) > 0 )
//        {
//            throw std::runtime_error( "Error when setting weights in single observation set, observation set has no data." );
//        }
        for ( unsigned int k = 0 ; k < numberOfObservations_ ; k++ )
        {
            for ( unsigned int i = 0 ; i < singleObservationSize_ ; i++ )
            {
                weights_.at( k )[ i ] = weightsVector[ k * singleObservationSize_ + i ];
            }
        }
//        weightsVector_ = weightsVector;
    }


//    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > createFilteredObservationSet(
//        const std::vector< int > indices )
//    {
//        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > filteredObservations = observations_;
//        std::vector< TimeType > filteredObservationTimes = observationTimes_;
//        std::vector< Eigen::VectorXd > filteredObservationsDependentVariables = observationsDependentVariables_;
//
//        sort(indices.begin(), indices.end(), std::greater<int>());
//        for( unsigned int i = 0; i < indices.size( ); i++ )
//        {
//            filteredObservationTimes.erase( filteredObservationTimes.begin( ) + indices.at( i ) );
//            filteredObservations.erase( filteredObservations.begin( ) + indices.at( i ) );
//            if( filteredObservationsDependentVariables.size( ) >  0 )
//            {
//                filteredObservationsDependentVariables.erase( filteredObservationsDependentVariables.begin( ) + indices.at( i ) );
//            }
//        }
//        return std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
//            observableType_, linkEnds_, filteredObservations, filteredObservationTimes,
//            referenceLinkEnd_, filteredObservationsDependentVariables, dependentVariableCalculator_, ancilliarySettings_ );
//    }


    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getComputedObservations( ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > computedObservations;
        for ( unsigned int k = 0 ; k < observations_.size( ) ; k++ )
        {
            computedObservations.push_back( observations_.at( k ) - residuals_.at( k ) );
        }
        return computedObservations;
    }


    void removeSingleObservation( unsigned int indexToRemove )
    {
        observations_.erase( observations_.begin( ) + indexToRemove );
        residuals_.erase( residuals_.begin( ) + indexToRemove );
        weights_.erase( weights_.begin( ) + indexToRemove );

        // Update number of observations and time bounds
        numberOfObservations_ = observations_.size( );
        updateTimeBounds( );
    }

    void removeObservations( std::vector< unsigned int >& indicesToRemove )
    {
        unsigned int counter = 0;
        for ( auto ind : indicesToRemove )
        {
            removeSingleObservation( ind - counter );
            counter += 1;
        }
    }

    void filterObservations( double residualCutOffValue )
    {
        std::vector< unsigned int > indicesToRemove;
        for( unsigned int j = 0 ; j < numberOfObservations_ ; j++ )
        {
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservationResidual = residuals_.at( j );
            bool removeObservation = false;
            for( int k = 0 ; k < singleObservationSize_ ; k++ )
            {
                if( singleObservationResidual[ k ] > residualCutOffValue )
                {
                    removeObservation = true;
                }
            }
            if( removeObservation )
            {
                indicesToRemove.push_back( j );
            }
        }
        removeObservations( indicesToRemove );
    }

//    auto getComputedObservable( )
//    {
//        return residuals_ + observations_;
//    }

private:

    const ObservableType observableType_;

    const LinkDefinition linkEnds_;

    std::pair< TimeType, TimeType > timeBounds_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations_;

    std::vector< TimeType > observationTimes_;

    const LinkEndType referenceLinkEnd_;

    std::vector< Eigen::VectorXd > observationsDependentVariables_;

    const std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator > dependentVariableCalculator_;

    const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings_;

    int numberOfObservations_;

    int singleObservationSize_;

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights_;

//    Eigen::VectorXd weightsVector_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals_;

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

    Eigen::Matrix< double, Eigen::Dynamic, 1 > getConcatenatedWeightVector( )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > weightVector =
                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( totalObservableSize_, 1 );
        unsigned int indexObs = 0;
        for ( auto observableIt : observationSetList_ )
        {
            for ( auto linkEndIt : observableIt.second )
            {
                for ( auto singleObsSet : linkEndIt.second )
                {
                    unsigned int singleObsSetSize = singleObsSet->getTotalObservationSetSize( );
                    weightVector.segment( indexObs, singleObsSetSize ) = singleObsSet->getWeightsVector( );
                    indexObs += singleObsSetSize;
                }
            }
        }
        return weightVector;
    }

//    std::vector< TimeType > getConcatenatedWeightVector( )
//    {
//        // for now, this only takes the weights set from the single observation sets.
//        // TODO may require a change later to accomodate other sources.
//
//        // lazy evaluation, check if the weights has been set, otherwise set it then return
//        if ( concatenatedWeights_.size( ) != totalObservableSize_ )
//        {
//            concatenatedWeights_.resize( totalObservableSize_ );
//            Eigen::VectorXd temp = getWeightsFromSingleObservationSets();
//            if (temp.size() > 0) {
//                Eigen::Map<Eigen::VectorXd> tempMap(temp.data(), temp.size());
//                concatenatedWeights_ = std::vector<ObservationScalarType>(
//                    tempMap.data(), tempMap.data() + tempMap.size());
//            }
//        }
//
//        return concatenatedWeights_;
//    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getConcatenatedResidualsVector( )
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualsVector =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObservableSize_, 1 );
        unsigned int indexObs = 0;
        for ( auto observableIt : observationSetList_ )
        {
            for ( auto linkEndIt : observableIt.second )
            {
                for ( auto singleObsSet : linkEndIt.second )
                {
                    unsigned int singleObsSetSize = singleObsSet->getTotalObservationSetSize( );
                    residualsVector.segment( indexObs, singleObsSetSize ) = singleObsSet->getResidualsVector( );
                    indexObs += singleObsSetSize;
                }
            }
        }
        return residualsVector;
    }

    std::vector< double > getConcatenatedDoubleTimeVector( )
    {
        return utilities::staticCastVector<double, TimeType>( concatenatedTimes_ );
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

    std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > getSingleTypeObservationSets(
            const ObservableType observableType ) const
    {
        if( observationSetList_.count( observableType ) == 0 )
        {
            throw std::runtime_error( "Error when retrieving observable of type " + observation_models::getObservableName( observableType ) +
                                      " from observation collection, no such observable exists" );
        }
        return observationSetList_.at( observableType );
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

    std::map< ObservableType, std::map< int, std::vector< std::pair< double, double > > > > getSortedObservationSetsTimeBounds( )
    {
        std::map< ObservableType, std::map< int, std::vector< std::pair< double, double > > > > observationSetTimeBounds;
        for( auto it1 : observationSetList_ )
        {
            for( auto it2 : it1.second )
            {
                int currentLinkEndId = linkEndIds_[ it2.first ];
                for( unsigned int i = 0; i < it2.second.size( ); i++ )
                {
                    std::pair< TimeType, TimeType > timeBounds = it2.second.at( i )->getTimeBounds( );
                    observationSetTimeBounds[ it1.first ][ currentLinkEndId  ].push_back(
                        { static_cast< double >( timeBounds.first ), static_cast< double >( timeBounds.second ) } );
                }
            }
        }
        return observationSetTimeBounds;
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

    void addSingleObservationSet(
            const ObservableType observableType,
            const LinkEnds& linkEnds,
            std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > observationSet )
    {
        observationSetList_[ observableType ][ linkEnds ].push_back( observationSet );
    }

    void removeSingleObservationSet(
            const unsigned int indexSetToRemove )
    {
//        observationSetList_.
    }

    void removeSingleObservationSet(
            const ObservableType observableType,
            const LinkDefinition& linkEnds,
            const std::pair< TimeType, TimeType > timeBounds )
    {

    }

//    Eigen::VectorXd getWeightsFromSingleObservationSets( )
//    {
//        Eigen::VectorXd weightsVector = Eigen::VectorXd::Zero( totalObservableSize_ );
//
//        for( auto observationIterator : observationSetList_ )
//        {
//            ObservableType currentObservableType = observationIterator.first;
//
//            for( auto linkEndIterator : observationIterator.second )
//            {
//                LinkEnds currentLinkEnds = linkEndIterator.first;
//                for( unsigned int i = 0; i < linkEndIterator.second.size( ); i++ )
//                {
//                    std::pair< int, int > startAndSize = observationSetStartAndSize_.at( currentObservableType ).at( currentLinkEnds ).at( i );
////                    if( observationSetList_.at( currentObservableType ).at( currentLinkEnds ).at( i )->getWeightsVectorReference( ).rows( ) ==
////                        startAndSize.second )
////                    {
//                        weightsVector.segment( startAndSize.first, startAndSize.second ) =
//                            observationSetList_.at( currentObservableType ).at( currentLinkEnds ).at(
//                                i )->getWeightsVector( );
////                    }
////                    else
////                    {
////                        throw std::runtime_error( "Error when compiling full weights vector from single observation set, sizes are inconsistent" );
////                    }
//                }
//            }
//        }
//        return weightsVector;
//    }

    //! Function to set a constant values for all observation weights
    /*!
     * Function to set a constant values for all observation weights
     * \param constantWeight Constant weight that is to be set for all observations
     */
    void setConstantWeightsMatrix( const double weight = 1.0 )
    {
        for ( auto observableIt : observationSetList_ )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto singleObsSet : linkEndsIt.second )
                {
                    singleObsSet->setConstantWeight( weight );
                }
            }
        }
    }

    //! Set constant scalar weight for all observables of given type
    void setConstantSingleObservableWeights(
            const observation_models::ObservableType currentObservable,
            const double weight )
    {
        if ( observationSetList_.count( currentObservable ) )
        {
            for ( auto linkEndsIt : observationSetList_.at( currentObservable ) )
            {
                for ( auto singleObsSet : linkEndsIt.second )
                {
                    singleObsSet->setConstantWeight( weight );
                }
            }
        }
        else
        {
            std::cerr<< "Warning when setting weights for data type "<< std::to_string( currentObservable) <<  ". " <<
                     " No data of given type found." <<std::endl;
        }
    }

    //! Set constant vector weight for all observables of given type
    void setConstantSingleObservableVectorWeights(
            const observation_models::ObservableType currentObservable,
            const Eigen::VectorXd& weight )
    {
        if ( observationSetList_.count( currentObservable ) )
        {
            for ( auto linkEndsIt : observationSetList_.at( currentObservable ) )
            {
                for ( auto singleObsSet : linkEndsIt.second )
                {
                    if ( singleObsSet->getSingleObservableSize( ) != weight.size( ) )
                    {
                        throw std::runtime_error( "Error when setting weight vector for observable, size should be " +
                                                  std::to_string( singleObsSet->getSingleObservableSize( ) ) + " but is " + std::to_string( weight.rows( ) ) );
                    }
                    singleObsSet->setConstantWeight( weight );
                }
            }
        }
        else
        {
            std::cerr<< "Warning when setting weights for data type "<< std::to_string( currentObservable) <<  ". " <<
                     " No data of given type found." <<std::endl;
        }
    }

    //! Set constant scalar weight for all observables of given type and link ends
    void setConstantSingleObservableAndLinkEndsWeights(
            const observation_models::ObservableType currentObservable,
            const observation_models::LinkEnds currentLinkEnds,
            const double weight )
    {
        if( observationSetList_.count( currentObservable) == 0 )
        {
            std::cerr<< "Warning when setting weights for data type "<< std::to_string( currentObservable) <<  ". " <<
                     " No data of given type found." <<std::endl;
        }
        else if( observationSetList_.at( currentObservable ).count( currentLinkEnds ) == 0 )
        {

            std::cerr<< "Warning when setting weights for data type " << std::to_string( currentObservable)<< " and link ends " <<
                     ". No data of given type and link ends found." <<std::endl;
        }
        else
        {
            for ( auto singleObsSet : observationSetList_.at( currentObservable ).at( currentLinkEnds ) )
            {
                singleObsSet->setConstantWeight( weight );
            }
        }
    }

    //! Set constant vector weight for all observables of given type and link ends
    void setConstantSingleObservableAndLinkEndsVectorWeights(
            const observation_models::ObservableType currentObservable,
            const observation_models::LinkEnds currentLinkEnds,
            const Eigen::VectorXd& weight )
    {
        if( observationSetList_.count( currentObservable) == 0 )
        {
            std::cerr<< "Warning when setting weights for data type "<< std::to_string( currentObservable) <<  ". " <<
                     " No data of given type found." <<std::endl;
        }
        else if( observationSetList_.at( currentObservable ).count( currentLinkEnds ) == 0 )
        {

            std::cerr<< "Warning when setting weights for data type " << std::to_string( currentObservable)<< " and link ends " <<
                     ". No data of given type and link ends found." <<std::endl;
        }
        else
        {
            for ( auto singleObsSet : observationSetList_.at( currentObservable ).at( currentLinkEnds ) )
            {
                if ( singleObsSet->getSingleObservableSize( ) != weight.size( ) )
                {
                    throw std::runtime_error( "Error when setting weight vector for observable, size should be " +
                                              std::to_string( singleObsSet->getSingleObservableSize( ) ) + " but is " + std::to_string( weight.rows( ) ) );
                }
                singleObsSet->setConstantWeight( weight );
            }
        }
    }

    //! Set constant vector weight for all observables of given type and link ends
    void setTabulatedSingleObservableAndLinkEndsWeights(
            const observation_models::ObservableType currentObservable,
            const observation_models::LinkEnds currentLinkEnds,
            const Eigen::VectorXd& weight )
    {
        if( observationSetList_.count( currentObservable) == 0 )
        {
            std::cerr<< "Warning when setting weights for data type "<< std::to_string( currentObservable) <<  ". " <<
                     " No data of given type found." <<std::endl;
        }
        else if( observationSetList_.at( currentObservable ).count( currentLinkEnds ) == 0 )
        {

            std::cerr<< "Warning when setting weights for data type " << std::to_string( currentObservable)<< " and link ends " <<
                     ". No data of given type and link ends found." <<std::endl;
        }
        else
        {
            std::pair< int, int > indicesToUse = observationSetList_.at( currentObservable ).at( currentLinkEnds );
            if( indicesToUse.second != weight.rows( ) )
            {
                throw std::runtime_error( "Error when setting total weight vector for observable and link ends, size should be " +
                                          std::to_string( indicesToUse.second ) + " but is " + std::to_string( weight.rows( ) ) );
            }

            unsigned int weightStartIndex = 0;
            for ( auto singleObsSet : observationSetList_.at( currentObservable ).at( currentLinkEnds ) )
            {
                Eigen::VectorXd currentObsSetWeights = weight.segment( weightStartIndex, singleObsSet->getTotalObservationSetSize( ) );
                singleObsSet->setTabulatedWeights( currentObsSetWeights );
                weightStartIndex += singleObsSet->getTotalObservationSetSize( );
            }
        }
    }


    //! Function to set a values for observation weights, constant per observable type
    /*!
     * Function to set a values for observation weights, constant per observable type
     * \param weightPerObservable Values for observation weights, constant per observable type
     */
    void setConstantPerObservableWeightsMatrix(
            const std::map< observation_models::ObservableType, double >& weightPerObservable )
    {
        for( auto observableIt : weightPerObservable )
        {
            setConstantSingleObservableWeights( observableIt.first, observableIt.second );
        }
    }

    void setConstantPerObservableVectorWeightsMatrix(
            const std::map< observation_models::ObservableType, Eigen::VectorXd >& weightPerObservable )
    {
        for( auto observableIt : weightPerObservable )
        {
            setConstantSingleObservableVectorWeights( observableIt.first, observableIt.second );
        }
    }

    //! Function to set a values for observation weights, constant per observable type and link ends type
    /*!
     * Function to set a values for observation weights, constant per observable type and link ends type
     * \param weightPerObservableAndLinkEnds Values for observation weights, constant per observable type and link ends type
     */
    void setConstantPerObservableAndLinkEndsWeights(
            const std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, double > >& weightPerObservableAndLinkEnds )
    {
        for( auto observableIt : weightPerObservableAndLinkEnds )
        {
            for( auto linkEndIt : observableIt.second )
            {
                setConstantSingleObservableAndLinkEndsWeights( observableIt.first, linkEndIt.first, linkEndIt.second );
            }
        }
    }

    void setConstantPerObservableAndLinkEndsVectorWeights(
            const std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > >& weightPerObservableAndLinkEnds )
    {
        for( auto observableIt : weightPerObservableAndLinkEnds )
        {
            for( auto linkEndIt : observableIt.second )
            {
                setConstantSingleObservableAndLinkEndsVectorWeights( observableIt.first, linkEndIt.first, linkEndIt.second );
            }
        }
    }

    void setConstantPerObservableAndLinkEndsWeights(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const double weight )
    {
        std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, double > > weightPerObservableAndLinkEnds;
        for( unsigned int i = 0; i < linkEnds.size( ); i++ )
        {
            weightPerObservableAndLinkEnds[ observableType ][ linkEnds.at( i ) ] =  weight;
        }
        setConstantPerObservableAndLinkEndsWeights( weightPerObservableAndLinkEnds );
    }

    void setConstantPerObservableAndLinkEndsVectorWeights(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const Eigen::VectorXd& weight )
    {
        std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > > weightPerObservableAndLinkEnds;
        for( unsigned int i = 0; i < linkEnds.size( ); i++ )
        {
            weightPerObservableAndLinkEnds[ observableType ][ linkEnds.at( i ) ] =  weight;
        }
        setConstantPerObservableAndLinkEndsVectorWeights( weightPerObservableAndLinkEnds );
    }

    void setTabulatedPerObservableAndLinkEndsWeights(
            const std::map< observation_models::ObservableType,
            const std::map< observation_models::LinkEnds, Eigen::VectorXd > >& weightsPerObservableAndLinkEnds )
    {
        for( auto observableIt : weightsPerObservableAndLinkEnds )
        {
            for( auto linkEndIt : observableIt.second )
            {
               setTabulatedSingleObservableAndLinkEndsWeights( observableIt.first, linkEndIt.first, linkEndIt.second );
            }
        }
    }

    void setTabulatedPerObservableAndLinkEndsWeights(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const Eigen::VectorXd& weights )
    {
        std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > > weightsPerObservableAndLinkEnds;
        for ( unsigned int i = 0; i < linkEnds.size( ); i++ )
        {
            weightsPerObservableAndLinkEnds[ observableType ][ linkEnds.at( i ) ] =  weights;
        }
        setTabulatedPerObservableAndLinkEndsWeights( weightsPerObservableAndLinkEnds );
    }


    std::tuple< Eigen::VectorXd, Eigen::VectorXd > getFullDependentVariableVector( const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > dependentVariable )
    {
        Eigen::VectorXd dependentVariablesVector = Eigen::VectorXd::Zero( concatenatedTimes_.size( ) );
        Eigen::VectorXd dependentVariableTimes = Eigen::VectorXd::Zero( concatenatedTimes_.size( ) );

        for( auto it : observationSetList_ )
        {
            if( observation_models::getObservableSize( it.first ) != 1 )
            {
                throw std::runtime_error( "Error, getting full dependent variable vector not yet supported for observables of size > 1" );
            }
            for( auto it2 : it.second )
            {
                for( unsigned int i = 0; i < it2.second.size( ); i++ )
                {
                    std::shared_ptr<simulation_setup::ObservationDependentVariableCalculator>
                        dependentVariableCalculator =
                        it2.second.at( i )->getDependentVariableCalculator( );
                    std::pair< int, int > currentVectorIndices = observationSetStartAndSize_.at( it.first ).at( it2.first ).at( i );

                    bool addDependentVariables = false;
                    if ( dependentVariableCalculator != nullptr )
                    {
                        std::pair<int, int> variableIndices = dependentVariableCalculator->getDependentVariableIndices(
                            dependentVariable );

                        if( variableIndices.second > 0 )
                        {
                            addDependentVariables = true;
                            std::map< double, Eigen::VectorXd > currentDependentVariableHistory =
                                utilities::sliceMatrixHistory< TimeType, double, double >( it2.second.at( i )->getDependentVariableHistory( ), variableIndices );

                            dependentVariableTimes.segment( currentVectorIndices.first, currentVectorIndices.second ) =
                                utilities::convertStlVectorToEigenVector( utilities::createVectorFromMapKeys( currentDependentVariableHistory ) );
                            Eigen::MatrixXd currentDependentVariableBlock = utilities::convertVectorHistoryToMatrix( currentDependentVariableHistory );
                            dependentVariablesVector.segment( currentVectorIndices.first, currentVectorIndices.second ) =
                                currentDependentVariableBlock.block( 0, 0, currentDependentVariableBlock.rows( ), 1 );
                        }
                    }
                }
            }
        }
        return std::make_pair( dependentVariableTimes, dependentVariablesVector );
    }


    void filterObservations(
            const std::map< ObservableType, double >& residualCutoffValuePerObservable )
    {
        for ( auto observationIt : residualCutoffValuePerObservable )
        {
            if ( observationSetList_.count( observationIt.first ) )
            {
                double filterValue = residualCutoffValuePerObservable.at( observationIt.first );
                for ( auto linkEndIt: observationSetList_.at( observationIt.first ) )
                {
                    for ( auto singleObsSet : linkEndIt.second )
                    {
                        singleObsSet->filterObservations( filterValue );
                    }
                }
            }
        }
    }

    void filterObservations(
            const std::map< ObservableType, std::map< LinkEnds, double > >& residualCutoffValuePerObservablePerLinkEnd )
    {
        for ( auto observationIt : residualCutoffValuePerObservablePerLinkEnd )
        {
            if ( observationSetList_.count( observationIt.first ) )
            {
                for ( auto linkEndIt : observationIt.second )
                {
                    if ( observationSetList_.at( observationIt ).count( linkEndIt.first ) )
                    {
                        double filterValue = residualCutoffValuePerObservablePerLinkEnd.at( observationIt.first ).at( linkEndIt.first );
                        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObservationSets =
                                observationSetList_.at( observationIt.first ).at( linkEndIt.first );
                        for ( auto singleObsSet : singleObservationSets )
                        {
                            singleObsSet->filterObservations( filterValue );
                        }
                    }
                }
            }
        }
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

    Eigen::Matrix< double, Eigen::Dynamic, 1 > concatenatedWeights_;

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

        std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > newObservationSet =
            std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                originalObservationSet->getObservableType( ),
                originalObservationSet->getLinkEnds( ),
                utilities::getStlVectorSegment( originalObservationSet->getObservationsReference( ), arcSplitIndices.at( i ).first, arcSplitIndices.at( i ).second ),
                utilities::getStlVectorSegment( originalObservationSet->getObservationTimesReference( ), arcSplitIndices.at( i ).first, arcSplitIndices.at( i ).second ),
                originalObservationSet->getReferenceLinkEnd( ),
                currentSplitDependentVariables,
                originalObservationSet->getDependentVariableCalculator( ),
                originalObservationSet->getAncilliarySettings( ) );
        if( newObservationSet->getObservationTimes( ).size( ) > 0 )
        {
            splitSingleObervationSet.push_back( newObservationSet );
        }
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

//template< typename ObservationScalarType = double, typename TimeType = double >
//std::map< observation_models::ObservableType, std::vector< std::pair< LinkEnds, std::vector< std::vector< int > > > > >
//    getObservationCollectionEntriesToFiler(
//        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > dataToFiler,
//        const Eigen::VectorXd& residualVector,
//        const std::map< ObservableType, double > residualCutoffValuePerObservable )
//{
//    // Check if input data is compatible
//    if( residualVector.rows( ) != dataToFiler->getTotalObservableSize( ) )
//    {
//        throw std::runtime_error( "Error when filtering observations, input size is incompatible" );
//    }
//
//    // Retrieve observations to filter
//    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets observationSetsToFilter = dataToFiler->getObservations( );
//
//    // Create data structure with filtered results
//    std::map< observation_models::ObservableType, std::vector< std::pair< LinkEnds, std::vector< std::vector< int > > > > > filterEntries;
//
//    // Iterate over all observable types
//    for( auto observationIt : observationSetsToFilter )
//    {
//        observation_models::ObservableType observableType = observationIt.first;
//
//        // Get filter value for current observable
//        double filterValue = residualCutoffValuePerObservable.at( observationIt.first );
//
//        // Iterate over all link ends
//        std::vector< std::pair< LinkEnds, std::vector< std::vector< int > > > > currentObservableEntriesToFilter;
//        for ( auto linkEndIt: observationIt.second )
//        {
//            observation_models::LinkEnds linkEnds = linkEndIt.first;
//
//            // Iterate over all observations with current link ends and type
//            std::vector< std::vector< int > > currentLinkEndsEntriesToFiler;
//            for ( unsigned int i = 0; i < linkEndIt.second.size( ); i++ )
//            {
//                std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > currentObservationSet =
//                    observationSetsToFilter.at( observationIt.first ).at( linkEndIt.first ).at( i );
//
//                // Get residuals for current set
//                std::pair< int, int > fullVectorStartAndSize =
//                    dataToFiler->getObservationSetStartAndSizeReference( ).at( observableType ).at( linkEnds ).at( i );
//                Eigen::VectorXd currentSetResiduals = residualVector.segment( fullVectorStartAndSize.first, fullVectorStartAndSize.second );
//
//                int currentObservableSize = getObservableSize( observableType );
//
//                // Check if data is compatible
//                if( currentSetResiduals.rows( ) != currentObservableSize * currentObservationSet->getNumberOfObservables( ) )
//                {
//                    throw std::runtime_error( "Error when filtering observations, input size of single observation set for " +
//                        getObservableName( observableType ) +", " + getLinkEndsString( linkEnds ) + ", set " +
//                        std::to_string( i ) + " is incompatible" );
//                }
//
//                std::vector< int > indicesToRemove;
//                for( unsigned int j = 0; j < currentObservationSet->getObservationTimes( ).size( ); j++ )
//                {
//                    bool removeObservation = false;
//                    for( int k = 0; k < currentObservableSize; k++ )
//                    {
//                        if( currentSetResiduals( j * currentObservableSize + k ) > filterValue )
//                        {
//                            removeObservation = true;
//                        }
//                    }
//                    if( removeObservation )
//                    {
//                        indicesToRemove.push_back( j );
//                    }
//                }
//                currentLinkEndsEntriesToFiler.push_back( indicesToRemove );
//            }
//            currentObservableEntriesToFilter.push_back( std::make_pair( linkEnds, currentLinkEndsEntriesToFiler ) );
//        }
//        filterEntries[ observableType ] = currentObservableEntriesToFilter;
//    }
//    return filterEntries;
//}

//template< typename ObservationScalarType = double, typename TimeType = double >
//std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > filterData(
//    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
//    const std::map< observation_models::ObservableType, std::vector< std::pair< LinkEnds, std::vector< std::vector< int > > > > >& filterEntries )
//{
//    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets filteredObservedObservationSets;
//    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > filteredObservationCollection;
//
//    for( auto observableIt : filterEntries )
//    {
//        ObservableType observableType = observableIt.first;
//        for( unsigned int i = 0; i < observableIt.second.size( ); i++ )
//        {
//            LinkEnds currentLinkEnds = observableIt.second.at( i ).first;
//            std::vector< std::vector< int > > linkEndListEntriesToRemove = observableIt.second.at( i ).second;
//
//            std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > observedSets =
//                observationCollection->getObservationsReference( ).at( observableType ).at( currentLinkEnds );
//            if( observedSets.size( ) != linkEndListEntriesToRemove.size( ) )
//            {
//                throw std::runtime_error( "Error when filtering observations, number of observation sets and filter list for " +
//                                          getObservableName( observableType ) +", " + getLinkEndsString( currentLinkEnds ) + " is incompatible" );
//            }
//            for( unsigned int j = 0; j < observedSets.size( ); j++ )
//            {
//                filteredObservedObservationSets[ observableType ][ currentLinkEnds ].push_back(
//                    observedSets.at( j )->createFilteredObservationSet(
//                        linkEndListEntriesToRemove.at( j ) ) );
//            }
//        }
//    }
//    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( filteredObservedObservationSets );
//}

//template< typename ObservationScalarType = double, typename TimeType = double >
//void filterObservedAndComputedData(
//    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observedDataCollection,
//    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > computedDataCollection,
//    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >& observedFilteredDataCollection,
//    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >& computedFilteredDataCollection,
//    const std::map< ObservableType, double > residualCutoffValuePerObservable )
//{
//    Eigen::VectorXd residualVector =
//        ( observedDataCollection->getObservationVector( ) - computedDataCollection->getObservationVector( ) ).template cast< double >( );
//
//    std::map< observation_models::ObservableType, std::vector< std::pair< LinkEnds, std::vector< std::vector< int > > > > > filterEntries =
//        getObservationCollectionEntriesToFiler( observedDataCollection, residualVector, residualCutoffValuePerObservable );
//
//    observedFilteredDataCollection = filterData( observedDataCollection, filterEntries );
//    computedFilteredDataCollection = filterData( computedDataCollection, filterEntries );
//}

////////////////////// DEPRECATEC
//
//template< typename ObservationScalarType = double, typename TimeType = double >
//std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > filterResidualOutliers(
//    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observedData,
//    const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > residualData,
//    const std::map< ObservableType, double > residualCutoffValuePerObservable )
//{
//    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets observedObservationSets = observedData->getObservations( );
//    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets filteredObservedObservationSets;
//    typename ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets residualObservationSets = residualData->getObservations( );
//
//    for( auto observationIt : residualObservationSets )
//    {
//        if( residualObservationSets.count( observationIt.first ) )
//        {
//            double filterValue = residualCutoffValuePerObservable.at( observationIt.first );
//            for ( auto linkEndIt: observationIt.second )
//            {
//                for ( unsigned int i = 0; i < linkEndIt.second.size( ); i++ )
//                {
//                    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > currentObservationSet =
//                        residualObservationSets.at( observationIt.first ).at( linkEndIt.first ).at( i );
//                    std::vector<int> indicesToRemove;
//                    std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> >
//                        residualObservationSet = linkEndIt.second.at( i );
//                    for ( unsigned int j = 0; j < currentObservationSet->getObservationTimes( ).size( ); j++ )
//                    {
//                        if ( currentObservationSet->getObservation( j ).array( ).abs( ).maxCoeff( ) > filterValue )
//                        {
//                            indicesToRemove.push_back( j );
//                        }
//                    }
//
//                    filteredObservedObservationSets[ observationIt.first ][ linkEndIt.first ].push_back(
//                        observedObservationSets[ observationIt.first ][ linkEndIt.first ].at( i )->createFilteredObservationSet(
//                            indicesToRemove ) );
//                }
//            }
//        }
//    }
//    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( filteredObservedObservationSets );
//}

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_OBSERVATIONS_H
