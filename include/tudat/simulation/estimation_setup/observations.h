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

enum ObservationFilterType
{
    residual_filtering,
    absolute_value_filtering,
//    epochs,
//    time_interval,
    dependent_variable_filtering
};

struct ObservationFilter
{
public:
    ObservationFilter( const ObservationFilterType filterType, const double filterValue, const bool filterOut = true ) :
            filterType_( filterType ), doubleFilterValue_( filterValue ), vectorFilterValue_( Eigen::VectorXd::Zero( 0 ) ), filterOut_( filterOut )
    {
    //        if ( !( isFilterDoubleValue( filterType_ ) ) )
    //        {
    //            throw std::runtime_error( "Error when creating ObservationFilter of type " + std::to_string( filterType_ ) + ", the input filter value should be a vector and not a double." );
    //        }
    }

    ObservationFilter( const ObservationFilterType filterType, const Eigen::VectorXd& filterValue, const bool filterOut = true ) :
            filterType_( filterType ), doubleFilterValue_( TUDAT_NAN ), vectorFilterValue_( filterValue ), filterOut_( filterOut )
    {
    //        if  ( isFilterDoubleValue( filterType_ ) )
    //        {
    //            throw std::runtime_error( "Error when creating ObservationFilter of type " + std::to_string( filterType_ ) + ", the input filter value should be a double and not a vector." );
    //        }
    }

    virtual ~ObservationFilter( ){ }

    ObservationFilterType getFilterType( ) const
    {
        return filterType_;
    }

    double getDoubleFilterValue( ) const
    {
        return doubleFilterValue_;
    }

    Eigen::VectorXd getVectorFilterValue( ) const
    {
        return vectorFilterValue_;
    }

    const bool filterOut( ) const
    {
        return filterOut_;
    }

    protected:
    ObservationFilterType filterType_;

    double doubleFilterValue_;

    Eigen::VectorXd vectorFilterValue_;

    bool filterOut_;

};

inline std::shared_ptr< ObservationFilter > observationFilter( const ObservationFilterType filterType, const double filterValue, const bool filterOut = true )
{
    return std::make_shared< ObservationFilter >( filterType, filterValue, filterOut );
}

inline std::shared_ptr< ObservationFilter > observationFilter( const ObservationFilterType filterType, const Eigen::VectorXd& filterValue, const bool filterOut = true )
{
    return std::make_shared< ObservationFilter >( filterType, filterValue, filterOut );
}


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

        // Sort observations and metadata per observation time
        orderObservationsAndMetadata( );

        singleObservationSize_ = getObservableSize( observableType );

        // Initialise weights
        for ( unsigned int k = 0 ; k < numberOfObservations_ ; k++ )
        {
            weights_.push_back( Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( singleObservationSize_, 1 ) );
        }

        // Initialise residuals
        for ( unsigned int k = 0 ; k < numberOfObservations_ ; k++ )
        {
            residuals_.push_back( Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_, 1 ) );
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

    TimeType getObservationTime( unsigned int index ) const
    {
        if ( index >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when retrieving single observation time, required index incompatible with number of observations." );
        }
        return observationTimes_.at( index );
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
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_ * numberOfObservations_, 1 );
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
        return timeBounds_;
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

    Eigen::VectorXd getSingleObservationDependentVariable( unsigned int index ) const
    {
        if ( index >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when retrieving single observation dependent variables, required index incompatible with number of observations." );
        }
        return observationsDependentVariables_.at( index );
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

    Eigen::Matrix< double, Eigen::Dynamic, 1 > getWeight( unsigned int index ) const
    {
        if ( index >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when retrieving single observation weight, required index incompatible with number of observations." );
        }
        return weights_.at( index );
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

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidual( unsigned int index ) const
    {
        if ( index >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when retrieving single observation residual, required index incompatible with number of observations." );
        }
        return residuals_.at( index );
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


    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getComputedObservations( ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > computedObservations;
        for ( unsigned int k = 0 ; k < observations_.size( ) ; k++ )
        {
            computedObservations.push_back( observations_.at( k ) - residuals_.at( k ) );
        }
        return computedObservations;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getComputedObservationsVector( ) const
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > computedObservationsVector =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_ * numberOfObservations_, 1 );

        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > computedObservations = getComputedObservations( );
        for( unsigned int i = 0; i < numberOfObservations_ ; i++ )
        {
            computedObservationsVector.segment( i * singleObservationSize_, singleObservationSize_ ) = computedObservations.at( i );
        }
        return computedObservationsVector;
    }

    int getNumberOfFilteredObservations( ) const
    {
        int numberFilteredObservations = 0;
        if ( filteredObservationSet_ != nullptr )
        {
            numberFilteredObservations = filteredObservationSet_->getNumberOfObservables( );
        }
        return numberFilteredObservations;
    }


    void removeSingleObservation( unsigned int indexToRemove )
    {
        if ( indexToRemove >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when removing single observation from SingleObservationSet, index incompatible with number of observations." );
        }

        // Update observations
        observations_.erase( observations_.begin( ) + indexToRemove );
        observationTimes_.erase( observationTimes_.begin( ) + indexToRemove );

        if ( residuals_.size( ) > 0 )
        {
            residuals_.erase( residuals_.begin( ) + indexToRemove );
        }

        if ( weights_.size( ) > 0 )
        {
            weights_.erase( weights_.begin( ) + indexToRemove );
        }

        if ( observationsDependentVariables_.size( ) > 0 )
        {
            observationsDependentVariables_.erase( observationsDependentVariables_.begin( ) + indexToRemove );
        }

        // Update number of observations and time bounds
        numberOfObservations_ = observations_.size( );
        updateTimeBounds( );
    }

    void removeObservations( const std::vector< unsigned int >& indicesToRemove )
    {
        unsigned int counter = 0;
        for ( auto ind : indicesToRemove )
        {
            removeSingleObservation( ind - counter ); // observations are already filtered and sorted
            counter += 1;
        }
    }

    void filterObservations( const std::shared_ptr< ObservationFilter > observationFilter, const bool saveFilteredObservations = true )
    {
        if ( observationFilter->filterOut( ) && filteredObservationSet_ == nullptr )
        {
            // Initialise empty filtered observation set
            filteredObservationSet_ = std::make_shared<SingleObservationSet< ObservationScalarType, TimeType > >(
                    observableType_, linkEnds_, std::vector<Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                    std::vector< TimeType >( ), referenceLinkEnd_, std::vector< Eigen::VectorXd >( ), dependentVariableCalculator_, ancilliarySettings_);
        }
        if ( !observationFilter->filterOut( ) && filteredObservationSet_ == nullptr )
        {
            throw std::runtime_error( "Error when attempting to un-filter observations, filtered observation set is empty." );
        }

        std::vector< unsigned int > indicesToRemove;
        switch ( observationFilter->getFilterType( ) )
        {
            case residual_filtering:
            {
                double residualCutOffValue = observationFilter->getDoubleFilterValue( );

                if ( observationFilter->filterOut( ) )
                {
                    for( unsigned int j = 0 ; j < numberOfObservations_ ; j++ )
                    {
                        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservationResidual = residuals_.at( j );
                        bool removeObservation = false;
                        for( int k = 0 ; k < singleObservationSize_ ; k++ )
                        {
                            if(  singleObservationResidual[ k ] > residualCutOffValue )
                            {
                                removeObservation = true;
                            }
                        }
                        if( removeObservation )
                        {
                            indicesToRemove.push_back( j );
                        }
                    }
                }
                else
                {
                    for( unsigned int j = 0 ; j < filteredObservationSet_->getNumberOfObservables( ) ; j++ )
                    {
                        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservationResidual = filteredObservationSet_->getResidual( j );
                        bool removeObservation = false;

                        for( int k = 0 ; k < singleObservationSize_ ; k++ )
                        {
                            if( singleObservationResidual[ k ] <= residualCutOffValue )
                            {
                                removeObservation = true;
                            }
                        }
                        if( removeObservation )
                        {
                            indicesToRemove.push_back( j );
                        }
                    }
                }
                break;
            }
            case absolute_value_filtering:
            {
                double absoluteValueCutOff = observationFilter->getDoubleFilterValue( );

                if ( observationFilter->filterOut( ) )
                {
                    for( unsigned int j = 0 ; j < numberOfObservations_ ; j++ )
                    {
                        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservation = observations_.at( j );
                        bool removeObservation = false;
                        for( int k = 0 ; k < singleObservationSize_ ; k++ )
                        {
                            if ( singleObservation[ k ] > absoluteValueCutOff )
                            {
                                removeObservation = true;
                            }
                        }
                        if( removeObservation )
                        {
                            indicesToRemove.push_back( j );
                        }
                    }
                }
                else
                {
                    for( unsigned int j = 0 ; j < filteredObservationSet_->getNumberOfObservables( ) ; j++ )
                    {
                        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservation = filteredObservationSet_->getObservation( j );
                        bool removeObservation = false;
                        for( int k = 0 ; k < singleObservationSize_ ; k++ )
                        {
                            if ( singleObservation[ k ] <= absoluteValueCutOff )
                            {
                                removeObservation = true;
                            }
                        }
                        if( removeObservation )
                        {
                            indicesToRemove.push_back( j );
                        }
                    }
                }
                break;
            }
            case dependent_variable_filtering:
            {
                throw std::runtime_error( "Error when filtering observations based on dependent variable, filter type not yet implemented "
                                          "(first requires fixing observation dependent variable feature)." );
                break;
            }
            default:
                throw std::runtime_error( "Observation filter type not recognised." );
        }

        if ( observationFilter->filterOut( ) )
        {
            moveObservationsInOutFilteredSet( indicesToRemove, true, saveFilteredObservations );
        }
        else
        {
            moveObservationsInOutFilteredSet( indicesToRemove, false, true );
        }
    }

    void addObservations(
            std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            std::vector< TimeType >& times,
            std::vector< Eigen::VectorXd > dependentVariables = { },
            std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights = { },
            std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals = { },
            const bool sortObservations = true )
    {
        if ( ( observations.size( ) != times.size( ) ) ||
                ( weights.size( ) > 0 && ( observations.size( ) != weights.size( ) ) ) ||
                ( residuals.size( ) > 0 && ( observations.size( ) != residuals.size( ) ) ) ||
                ( dependentVariables.size( ) > 0 && ( observations.size( ) != dependentVariables.size( ) ) ) )
        {
            throw std::runtime_error( "Error when adding observations to single obs set, input sizes are inconsistent." );
        }

        for ( unsigned int k = 0 ; k < observations.size( ) ; k++ )
        {
            observations_.push_back( observations.at( k ) );
            observationTimes_.push_back( times.at( k ) );

            // If residuals are set
            if ( ( residuals_.size( ) > 0 || numberOfObservations_ == 0 ) && residuals.size( ) > 0 )
            {
                residuals_.push_back( residuals.at( k ) );
            }

            // If weights ares set
            if ( ( weights_.size( ) > 0 || numberOfObservations_ == 0 ) && weights.size( ) > 0 )
            {
                weights_.push_back( weights.at( k ) );
            }

            // if dependent variables are set
            if ( ( observationsDependentVariables_.size( ) > 0 || numberOfObservations_ == 0 ) && dependentVariables.size( ) > 0 )
            {
                observationsDependentVariables_.push_back( dependentVariables.at( k ) );
            }

            numberOfObservations_ += 1;
        }

        // Sort observations
        if ( sortObservations )
        {
            orderObservationsAndMetadata( );
        }
    }


private:

    void orderObservationsAndMetadata( )
    {
        if( !std::is_sorted( observationTimes_.begin( ), observationTimes_.end( ) ) )
        {
            if( static_cast< int >( observations_.size( ) ) != numberOfObservations_ )
            {
                throw std::runtime_error( "Error when making SingleObservationSet, number of observations is incompatible after time ordering" );
            }

            std::multimap< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observationsMap;
            for( unsigned int i = 0; i < observations_.size( ); i++ )
            {
                observationsMap.insert( { observationTimes_.at( i ), observations_.at( i ) } );
            }
            observationTimes_ = utilities::createVectorFromMultiMapKeys( observationsMap );
            observations_ = utilities::createVectorFromMultiMapValues( observationsMap );

            if( observationsDependentVariables_.size( ) > 0 )
            {
                if( static_cast< int >( observationsDependentVariables_.size( ) ) != numberOfObservations_ )
                {
                    throw std::runtime_error( "Error when making SingleObservationSet, dependent variables vector size is incompatible after time ordering" );
                }
                std::map< TimeType, Eigen::VectorXd > observationsDependentVariablesMap;
                for( unsigned int i = 0; i < observationsDependentVariables_.size( ); i++ )
                {
                    observationsDependentVariablesMap[ observationTimes_.at( i ) ] = observationsDependentVariables_.at( i );
                }
                observationsDependentVariables_ = utilities::createVectorFromMapValues( observationsDependentVariablesMap );
            }

            if ( weights_.size( ) > 0 )
            {
                if( static_cast< int >( weights_.size( ) ) != numberOfObservations_ )
                {
                    throw std::runtime_error( "Error when making SingleObservationSet, weights size is incompatible after time ordering" );
                }
                std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > weightsMap;
                for( unsigned int i = 0; i < weights_.size( ); i++ )
                {
                    weightsMap[ observationTimes_.at( i ) ] = weights_.at( i );
                }
                weights_ = utilities::createVectorFromMapValues( weightsMap );
            }

            if ( residuals_.size( ) > 0 )
            {
                if( static_cast< int >( residuals_.size( ) ) != numberOfObservations_ )
                {
                    throw std::runtime_error( "Error when making SingleObservationSet, residuals size is incompatible after time ordering" );
                }
                std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > residualsMap;
                for( unsigned int i = 0; i < residuals_.size( ); i++ )
                {
                    residualsMap[ observationTimes_.at( i ) ] = residuals_.at( i );
                }
                residuals_ = utilities::createVectorFromMapValues( residualsMap );
            }
        }
    }

    void moveObservationsInOutFilteredSet( const std::vector< unsigned int >& indices, const bool moveInFilteredSet = true, const bool saveFilteredObservations = true )
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
        std::vector< TimeType > times;
        std::vector< Eigen::VectorXd > dependentVariables;
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights;
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;

        if ( moveInFilteredSet )
        {
            for ( auto index : indices )
            {
                if ( index >= numberOfObservations_ )
                {
                    throw std::runtime_error( "Error when moving observation to filtered observation set, index incompatible with number of observations." );
                }

                observations.push_back( observations_.at( index ) );
                times.push_back( observationTimes_.at( index ) );

                // If weights are set
                if ( weights_.size( ) > 0 )
                {
                    weights.push_back( weights_.at( index ) );
                }

                // If residuals are set
                if ( residuals_.size( ) > 0 )
                {
                    residuals.push_back( residuals_.at( index ) );
                }

                // If dependent variables are set
                if ( observationsDependentVariables_.size( ) > 0 )
                {
                    dependentVariables.push_back( observationsDependentVariables_.at( index ) );
                }
            }

            if ( saveFilteredObservations )
            {
                filteredObservationSet_->addObservations( observations, times, dependentVariables, weights, residuals, true );
            }
            removeObservations( indices );
        }
        else
        {
            for ( auto index : indices )
            {
                if (  getNumberOfFilteredObservations( ) == 0 )
                {
                    throw std::runtime_error( "Error when moving observation back from filtered observation set, filtered observation set is empty." );
                }
                if ( index >= getNumberOfFilteredObservations( ) )
                {
                    throw std::runtime_error( "Error when moving observation back from filtered observation set, index incompatible with number of observations." );
                }

                observations.push_back( filteredObservationSet_->getObservations( ).at( index ) );
                times.push_back( filteredObservationSet_->getObservationTimes( ).at( index ) );

                // If weights are set
                if ( weights_.size( ) > 0 )
                {
                    weights.push_back( filteredObservationSet_->getWeights( ).at( index ) );
                }

                // If residuals are set
                if ( residuals_.size( ) > 0 )
                {
                    residuals.push_back( filteredObservationSet_->getResiduals( ).at( index ) );
                }

                // If dependent variables are set
                if ( observationsDependentVariables_.size( ) > 0 )
                {
                    dependentVariables.push_back( filteredObservationSet_->getObservationsDependentVariables( ).at( index ) );
                }

            }
            addObservations( observations, times, dependentVariables, weights, residuals, true );
            filteredObservationSet_->removeObservations( indices );
        }
    }

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

    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > filteredObservationSet_;

};


enum ObservationParserType
{
    empty_parser,
    observable_type_parser,
    link_ends_parser,
    body_name_parser,
    time_bounds_parser,
    multi_type_parser
};

struct ObservationCollectionParser
{

public:

    ObservationCollectionParser( ) : parserType_( empty_parser ), useOppositeCondition_( false ){ }

    ObservationCollectionParser( const ObservationParserType parserType, const bool useOppositeCondition = false ) :
            parserType_( parserType ), useOppositeCondition_( useOppositeCondition ){ }

    virtual ~ObservationCollectionParser( ){ }

    ObservationParserType getObservationParserType( ) const
    {
        return parserType_;
    }

    bool useOppositeCondition( ) const
    {
        return useOppositeCondition_;
    }

protected:

    const ObservationParserType parserType_;

    const bool useOppositeCondition_;

};

struct ObservationCollectionObservableTypeParser : public ObservationCollectionParser
{
public:

    ObservationCollectionObservableTypeParser( const ObservableType observableType, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( observable_type_parser, useOppositeCondition ), observableTypes_( std::vector< ObservableType >( { observableType } ) ){ }

    ObservationCollectionObservableTypeParser( const std::vector< ObservableType > observableTypes, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( observable_type_parser, useOppositeCondition ), observableTypes_( observableTypes ){ }

    virtual ~ObservationCollectionObservableTypeParser( ){ }

    std::vector< ObservableType > getObservableTypes( ) const
    {
        return observableTypes_;
    }

protected:

    const std::vector< ObservableType > observableTypes_;

};

struct ObservationCollectionLinkEndsParser : public ObservationCollectionParser
{
public:

    ObservationCollectionLinkEndsParser( const LinkEnds linkEnds, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( link_ends_parser, useOppositeCondition ), linkEndsVector_( std::vector< LinkEnds >( { linkEnds } ) ){ }

    ObservationCollectionLinkEndsParser( const std::vector< LinkEnds > linkEndsVector, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( link_ends_parser, useOppositeCondition ), linkEndsVector_( linkEndsVector ){ }

    virtual ~ObservationCollectionLinkEndsParser( ){ }

    std::vector< LinkEnds > getLinkEndsVector( ) const
    {
        return linkEndsVector_;
    }

protected:

    const std::vector< LinkEnds > linkEndsVector_;

};

struct ObservationCollectionBodyNameParser : public ObservationCollectionParser
{
public:

    ObservationCollectionBodyNameParser( const std::string bodyName,
                                         const bool isStationName = false,
                                         const bool useOppositeCondition = false ) :
            ObservationCollectionParser( body_name_parser, useOppositeCondition ), bodyNames_( std::vector< std::string >( { bodyName } ) ), isStationName_( isStationName ){ }

    ObservationCollectionBodyNameParser( const std::vector< std::string > bodyNames,
                                         const bool isStationName = false,
                                         const bool useOppositeCondition = false ) :
            ObservationCollectionParser( body_name_parser, useOppositeCondition ), bodyNames_( bodyNames ), isStationName_( isStationName ){ }

    virtual ~ObservationCollectionBodyNameParser( ){ }

    std::vector< std::string > getBodyNames( ) const
    {
        return bodyNames_;
    }

    bool isCheckForStationName( ) const
    {
        return isStationName_;
    }

protected:

    const std::vector< std::string > bodyNames_;

    const bool isStationName_;

};

struct ObservationCollectionTimeBoundsParser : public ObservationCollectionParser
{
public:

    ObservationCollectionTimeBoundsParser( const std::pair< double, double > timeBounds, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( time_bounds_parser, useOppositeCondition ), timeBoundsVector_( std::vector< std::pair< double, double > >( { timeBounds } ) ){ }

    ObservationCollectionTimeBoundsParser( const std::vector< std::pair< double, double > > timeBoundsVector, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( time_bounds_parser, useOppositeCondition ), timeBoundsVector_( timeBoundsVector ){ }

    virtual ~ObservationCollectionTimeBoundsParser( ){ }

    std::vector< std::pair< double, double > > getTimeBoundsVector( ) const
    {
        return timeBoundsVector_;
    }

protected:

    const std::vector< std::pair< double, double > > timeBoundsVector_;

};



struct ObservationCollectionMultiTypeParser : public ObservationCollectionParser
{
public:

    ObservationCollectionMultiTypeParser( const std::vector< std::shared_ptr< ObservationCollectionParser > >& observationParsers ) :
            ObservationCollectionParser( multi_type_parser, false ), observationParsers_( observationParsers ){ }

    virtual ~ObservationCollectionMultiTypeParser( ){ }

    std::vector< std::shared_ptr< ObservationCollectionParser > > getObservationParsers_(  ) const
    {
        return observationParsers_;
    }

protected:

    const std::vector< std::shared_ptr< ObservationCollectionParser > > observationParsers_;

};

inline std::shared_ptr< ObservationCollectionParser > observationParser( )
{
    return std::make_shared< ObservationCollectionParser >( );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const ObservableType observableType, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionObservableTypeParser >( observableType, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::vector< ObservableType >& observableTypes, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionObservableTypeParser >( observableTypes, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const LinkEnds linkEnds, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionLinkEndsParser >( linkEnds, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::vector< LinkEnds >& linkEndsVector, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionLinkEndsParser >( linkEndsVector, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::string bodyName, const bool isStationName = false, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionBodyNameParser >( bodyName, isStationName, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::vector< std::string >& bodyNames,
                                                                         const bool isStationName = false, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionBodyNameParser >( bodyNames, isStationName, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::pair< double, double >& timeBounds, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionTimeBoundsParser >( timeBounds, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::vector< std::pair< double, double > >& timeBoundsVector, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionTimeBoundsParser >( timeBoundsVector, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::vector< std::shared_ptr< ObservationCollectionParser > >& observationParsers )
{
    return std::make_shared< ObservationCollectionMultiTypeParser >( observationParsers );
}



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

//    Eigen::Matrix< double, Eigen::Dynamic, 1 > getConcatenatedWeightVector( )
//    {
//        Eigen::Matrix< double, Eigen::Dynamic, 1 > weightVector =
//                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( totalObservableSize_, 1 );
//        unsigned int indexObs = 0;
//        for ( auto observableIt : observationSetList_ )
//        {
//            for ( auto linkEndIt : observableIt.second )
//            {
//                for ( auto singleObsSet : linkEndIt.second )
//                {
//                    unsigned int singleObsSetSize = singleObsSet->getTotalObservationSetSize( );
//                    weightVector.segment( indexObs, singleObsSetSize ) = singleObsSet->getWeightsVector( );
//                    indexObs += singleObsSetSize;
//                }
//            }
//        }
//        return weightVector;
//    }

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

//    std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > getSingleTypeObservationSets(
//            const ObservableType observableType ) const
//    {
//        if( observationSetList_.count( observableType ) == 0 )
//        {
//            throw std::runtime_error( "Error when retrieving observable of type " + observation_models::getObservableName( observableType ) +
//                                      " from observation collection, no such observable exists" );
//        }
//        return observationSetList_.at( observableType );
//    }


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

    std::vector< ObservableType > getObservableTypes( std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        std::vector< ObservableType > observableTypes;
        for ( auto observableIt : observationSetsIndices )
        {
            if ( std::count( observableTypes.begin( ), observableTypes.end( ), observableIt.first ) == 0 )
            {
                observableTypes.push_back( observableIt.first );
            }
        }
        return observableTypes;
    }

    std::vector< LinkEnds > getLinkEnds( std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        std::vector< LinkEnds > linkEnds;
        for ( auto observableIt : observationSetsIndices )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                if ( std::count( linkEnds.begin( ), linkEnds.end( ), linkEndsIt.first ) == 0 )
                {
                    linkEnds.push_back( linkEndsIt.first );
                }
            }
        }
        return linkEnds;
    }

    std::vector< std::string > getBodiesInLinkEnds( std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        std::vector< std::string > bodiesInLinkEnds;
        for ( auto observableIt : observationSetsIndices )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto it : linkEndsIt.first )
                {
                    if ( std::count( bodiesInLinkEnds.begin( ), bodiesInLinkEnds.end( ), it.second.bodyName_ ) == 0 )
                    {
                        bodiesInLinkEnds.push_back( it.second.bodyName_ );
                    }
                }
            }
        }
        return bodiesInLinkEnds;
    }

    std::vector< std::string > getReferencePointsInLinkEnds( std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        std::vector< std::string > referencePoints;
        for ( auto observableIt : observationSetsIndices )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto it : linkEndsIt.first )
                {
                    if ( ( it.second.stationName_ != "" ) && ( std::count( referencePoints.begin( ), referencePoints.end( ), it.second.stationName_ ) == 0 ) )
                    {
                        referencePoints.push_back( it.second.stationName_ );
                    }
                }
            }
        }
        return referencePoints;
    }

    std::vector< std::pair< TimeType, TimeType > > getTimeBoundsList(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        std::vector< std::pair< TimeType, TimeType > > timeBounds;
        for ( auto observableIt : observationSetsIndices )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto index : linkEndsIt.second )
                {
                    std::pair< TimeType, TimeType > currentTimeBounds = observationSetList_.at( observableIt.first ).at( linkEndsIt.first ).at( index )->getTimeBounds( );
                    if ( std::count( timeBounds.begin( ), timeBounds.end( ), currentTimeBounds ) == 0 )
                    {
                        timeBounds.push_back( currentTimeBounds );
                    }
                }
            }
        }
        return timeBounds;
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservations(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for ( auto observableIt : observationSetsIndices )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto index : linkEndsIt.second )
                {
                    observations.push_back( observationSetList_.at( observableIt.first ).at( linkEndsIt.first ).at( index )->getObservationsVector( ) );
                }
            }
        }
        return observations;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getConcatenatedObservations(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations = getObservations( observationParser );

        unsigned int totalObsSize = 0;
        for ( auto obs : observations )
        {
            totalObsSize += obs.size( );
        }

        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > concatenatedObservations = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObsSize, 1 );
        unsigned int obsIndex = 0;
        for ( unsigned int k = 0 ; k < observations.size( ) ; k++ )
        {
            concatenatedObservations.block( obsIndex, 0, observations.at( k ).size( ), 1 ) = observations.at( k );
            obsIndex += observations.at( k ).size( );
        }

        return concatenatedObservations;
    }

    std::vector< std::vector< TimeType > > getObservationTimes(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< std::vector< TimeType > > observationsTimes;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for ( auto observableIt : observationSetsIndices )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto index : linkEndsIt.second )
                {
                    observationsTimes.push_back( observationSetList_.at( observableIt.first ).at( linkEndsIt.first ).at( index )->getObservationTimes( ) );
                }
            }
        }
        return observationsTimes;
    }

    std::vector< TimeType > getConcatenatedObservationTimes(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< std::vector< TimeType > > observationTimes = getObservationTimes( observationParser );

        std::vector< TimeType > concatenatedObservationsTimes;
        for ( auto times : observationTimes )
        {
            for ( unsigned int k = 0 ; k < times.size( ) ; k++ )
            {
                concatenatedObservationsTimes.push_back( times.at( k ) );
            }
        }

        return concatenatedObservationsTimes;
    }

    std::pair< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >, std::vector< std::vector< TimeType > > > getObservationsAndTimes(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
    {
        return std::make_pair( getObservations( observationParser ), getObservationTimes( observationParser ) );
    }

    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > getConcatenatedObservationsAndTimes(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
    {
        return std::make_pair( getConcatenatedObservations( observationParser ), getConcatenatedObservationTimes( observationParser ) );
    }

    std::vector< Eigen::VectorXd > getWeights(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::VectorXd > weights;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for ( auto observableIt : observationSetsIndices )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto index : linkEndsIt.second )
                {
                    weights.push_back( observationSetList_.at( observableIt.first ).at( linkEndsIt.first ).at( index )->getWeightsVector( ) );
                }
            }
        }
        return weights;
    }

    Eigen::VectorXd getConcatenatedWeights(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::VectorXd > weights = getWeights( observationParser );

        unsigned int totalObsSize = 0;
        for ( auto weight : weights )
        {
            totalObsSize += weight.size( );
        }

        Eigen::VectorXd concatenatedWeights = Eigen::VectorXd::Zero( totalObsSize );
        unsigned int obsIndex = 0;
        for ( unsigned int k = 0 ; k < weights.size( ) ; k++ )
        {
            concatenatedWeights.block( obsIndex, 0, weights.at( k ).size( ), 1 ) = weights.at( k );
            obsIndex += weights.at( k ).size( );
        }

        return concatenatedWeights;
    }

    std::vector< Eigen::VectorXd > getResiduals(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::VectorXd > residuals;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for ( auto observableIt : observationSetsIndices )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto index : linkEndsIt.second )
                {
                    residuals.push_back( observationSetList_.at( observableIt.first ).at( linkEndsIt.first ).at( index )->getResidualsVector( ) );
                }
            }
        }
        return residuals;
    }

    Eigen::VectorXd getConcatenatedResiduals(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::VectorXd > residuals = getResiduals( observationParser );

        unsigned int totalObsSize = 0;
        for ( auto residual : residuals )
        {
            totalObsSize += residual.size( );
        }

        Eigen::VectorXd concatenatedResiduals = Eigen::VectorXd::Zero( totalObsSize );
        unsigned int obsIndex = 0;
        for ( unsigned int k = 0 ; k < residuals.size( ) ; k++ )
        {
            concatenatedResiduals.block( obsIndex, 0, residuals.at( k ).size( ), 1 ) = residuals.at( k );
            obsIndex += residuals.at( k ).size( );
        }

        return concatenatedResiduals;
    }

    std::vector< Eigen::VectorXd > getComputedObservations(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::VectorXd > computedObservations;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for ( auto observableIt : observationSetsIndices )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto index : linkEndsIt.second )
                {
                    computedObservations.push_back( observationSetList_.at( observableIt.first ).at( linkEndsIt.first ).at( index )->getComputedObservationsVector( ) );
                }
            }
        }
        return computedObservations;
    }

    Eigen::VectorXd getConcatenatedComputedObservations(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::VectorXd > computedObservations = getComputedObservations( observationParser );

        unsigned int totalObsSize = 0;
        for ( auto obs : computedObservations )
        {
            totalObsSize += obs.size( );
        }

        Eigen::VectorXd concatenatedComputedObservations = Eigen::VectorXd::Zero( totalObsSize );
        unsigned int obsIndex = 0;
        for ( unsigned int k = 0 ; k < computedObservations.size( ) ; k++ )
        {
            concatenatedComputedObservations.block( obsIndex, 0, computedObservations.at( k ).size( ), 1 ) = computedObservations.at( k );
            obsIndex += computedObservations.at( k ).size( );
        }

        return concatenatedComputedObservations;
    }

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > getObservationSetStartAndSize(
            std::shared_ptr< ObservationCollectionParser > observationParser )
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > observationSetStartAndSize;

        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices = getSingleObservationSetsIndices( observationParser );
        for ( auto observableIt : observationSetsIndices )
        {
            std::map< LinkEnds, std::vector< std::pair< int, int > > > currentObservableTypeIndices;
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto index : linkEndsIt.second )
                {
                    currentObservableTypeIndices[ linkEndsIt.first ].push_back( observationSetStartAndSize_.at( observableIt.first ).at( linkEndsIt.first ).at( index ) );
                }
            }
            observationSetStartAndSize[ observableIt.first ] = currentObservableTypeIndices;
        }
        return observationSetStartAndSize;
    }


//    void addSingleObservationSet(
//            const ObservableType observableType,
//            const LinkEnds& linkEnds,
//            std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > observationSet )
//    {
//        observationSetList_[ observableType ][ linkEnds ].push_back( observationSet );
//    }

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


//    void filterObservations(
//            const std::map< ObservableType, double >& residualCutoffValuePerObservable )
//    {
//        for ( auto observationIt : residualCutoffValuePerObservable )
//        {
//            if ( observationSetList_.count( observationIt.first ) )
//            {
//                double filterValue = residualCutoffValuePerObservable.at( observationIt.first );
//                for ( auto linkEndIt: observationSetList_.at( observationIt.first ) )
//                {
//                    for ( auto singleObsSet : linkEndIt.second )
//                    {
//                        singleObsSet->filterObservations( filterValue );
//                    }
//                }
//            }
//        }
//    }

//    void filterObservations(
//            const std::map< ObservableType, std::map< LinkEnds, double > >& residualCutoffValuePerObservablePerLinkEnd )
//    {
//        for ( auto observationIt : residualCutoffValuePerObservablePerLinkEnd )
//        {
//            if ( observationSetList_.count( observationIt.first ) )
//            {
//                for ( auto linkEndIt : observationIt.second )
//                {
//                    if ( observationSetList_.at( observationIt ).count( linkEndIt.first ) )
//                    {
//                        double filterValue = residualCutoffValuePerObservablePerLinkEnd.at( observationIt.first ).at( linkEndIt.first );
//                        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObservationSets =
//                                observationSetList_.at( observationIt.first ).at( linkEndIt.first );
//                        for ( auto singleObsSet : singleObservationSets )
//                        {
//                            singleObsSet->filterObservations( filterValue );
//                        }
//                    }
//                }
//            }
//        }
//    }

    void filterObservations( const std::map< std::shared_ptr< ObservationCollectionParser >, std::shared_ptr< ObservationFilter > >& observationFilters,
                             const bool saveFilteredObservations = true )
    {
        // Parse all observation filters
        for ( auto filterIt : observationFilters )
        {
            std::shared_ptr< ObservationCollectionParser > observationParser = filterIt.first;
            std::shared_ptr< ObservationFilter > filter = filterIt.second;

            // Retrieve single observation sets based on observation parser
            std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObsSets =
                    getSingleObservationSets( observationParser );

            // Filter observations for all single observation sets
            for ( auto obsSet : singleObsSets )
            {
               obsSet->filterObservations( filter, saveFilteredObservations );
            }
        }

        // Reset observation set indices and concatenated observations and times
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }


   void filterObservations( std::shared_ptr< ObservationFilter > observationFilter,
                            const bool saveFilteredObservations = true,
                            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
    {
        std::map< std::shared_ptr< ObservationCollectionParser >, std::shared_ptr< ObservationFilter > > observationFilters
                { { observationParser, observationFilter } };

        filterObservations( observationFilters, saveFilteredObservations );
    }

    std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > getSingleObservationSetsIndices(
            const std::shared_ptr< ObservationCollectionParser > observationParser ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > observationSetsIndices;

        ObservationParserType parserType = observationParser->getObservationParserType( );
        switch( parserType )
        {
            case empty_parser:
            {
                for ( auto observableIt : observationSetList_ )
                {
                    std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                    for ( auto linkEndsIt : observableIt.second )
                    {
                        for ( unsigned int k = 0 ; k < linkEndsIt.second.size( ) ; k++ )
                        {
                            indicesPerObservable[ linkEndsIt.first ].push_back( k );
                        }
                    }
                    observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                }
                break;
            }
            case observable_type_parser:
            {
                std::vector< ObservableType > observableTypes = std::dynamic_pointer_cast< ObservationCollectionObservableTypeParser >( observationParser )->getObservableTypes( );
                for ( auto observableIt : observationSetList_ )
                {
                    if ( ( ( std::count( observableTypes.begin( ), observableTypes.end( ), observableIt.first ) ) && ( !observationParser->useOppositeCondition( ) ) ) ||
                         ( ( std::count( observableTypes.begin( ), observableTypes.end( ), observableIt.first ) == 0 ) && ( observationParser->useOppositeCondition( ) ) ) )
                    {
                        std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                        for ( auto linkEndsIt : observableIt.second )
                        {
                            for ( unsigned int k = 0 ; k < linkEndsIt.second.size( ) ; k++ )
                            {
                                indicesPerObservable[ linkEndsIt.first ].push_back( k );
                            }
                        }
                        observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                    }
                }

                break;
            }
            case link_ends_parser:
            {
                std::vector< LinkEnds > linkEndsVector = std::dynamic_pointer_cast< ObservationCollectionLinkEndsParser >( observationParser )->getLinkEndsVector( );
                std::vector< LinkEnds > allLinkEnds = getLinkEnds( );
                for ( auto linkEnds : allLinkEnds )
                {
                    if ( ( ( std::count( linkEndsVector.begin( ), linkEndsVector.end( ), linkEnds ) ) && ( !observationParser->useOppositeCondition( ) ) ) ||
                         ( ( std::count( linkEndsVector.begin( ), linkEndsVector.end( ), linkEnds ) == 0 ) && ( observationParser->useOppositeCondition( ) ) ) )
                    {
                        for ( auto observableIt : observationSetList_ )
                        {
                            std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                            if ( observationSetsIndices.count( observableIt.first ) )
                            {
                                indicesPerObservable = observationSetsIndices.at( observableIt.first );
                            }

                            if ( observableIt.second.count( linkEnds ) )
                            {
                                for ( unsigned int k = 0 ; k < observableIt.second.at( linkEnds ).size( ) ; k++ )
                                {
                                    indicesPerObservable[ linkEnds ].push_back( k );
                                }
                            }
                            observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                        }
                    }
                }
                break;
            }
            case body_name_parser:
            {
                std::shared_ptr< ObservationCollectionBodyNameParser > bodyNameObservationParser =
                        std::dynamic_pointer_cast< ObservationCollectionBodyNameParser >( observationParser );
                std::vector< std::string > bodyNames = bodyNameObservationParser->getBodyNames( );
                bool isStationName = bodyNameObservationParser->isCheckForStationName( );

                // Retrieve all body names
                std::vector< std::string > allBodyNames;
                if ( !isStationName )
                {
                    allBodyNames = getBodiesInLinkEnds( );
                }
                else
                {
                    allBodyNames = getReferencePointsInLinkEnds( );
                }

                for ( auto name : allBodyNames )
                {
                    if ( ( ( std::count( bodyNames.begin( ), bodyNames.end( ), name ) ) && ( !observationParser->useOppositeCondition( ) ) ) ||
                         ( ( std::count( bodyNames.begin( ), bodyNames.end( ), name ) == 0 ) && ( observationParser->useOppositeCondition( ) ) ) )
                    {
                        for ( auto observableIt : observationSetList_ )
                        {
                            std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                            if ( observationSetsIndices.count( observableIt.first ) )
                            {
                                indicesPerObservable = observationSetsIndices.at( observableIt.first );
                            }

                            for ( auto linkEndsIt : observableIt.second )
                            {
                                // FIX!
                                bool isBodyInLinkEnds = false;
                                bool isGroundStationInLinkEnds = false;
                                for ( auto it : linkEndsIt.first )
                                {
                                    if ( it.second.bodyName_ == name )
                                    {
                                        isBodyInLinkEnds = true;
                                    }
                                    if ( it.second.stationName_ == name )
                                    {
                                        isGroundStationInLinkEnds = true;
                                    }
                                }

                                if ( ( !isStationName && isBodyInLinkEnds/*( linkEndsIt.first, name )*/ ) ||
                                     ( isStationName && isGroundStationInLinkEnds/*( linkEndsIt.first, name )*/ ) )
                                {
                                    for ( unsigned int k = 0 ; k < linkEndsIt.second.size( ) ; k++ )
                                    {
                                        indicesPerObservable[ linkEndsIt.first ].push_back( k );
                                    }
                                }
                            }
                            observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                        }
                    }
                }
                break;
            }
            case time_bounds_parser:
            {
                std::vector< std::pair< double, double > > timeBoundsVector = std::dynamic_pointer_cast< ObservationCollectionTimeBoundsParser >( observationParser )->getTimeBoundsVector( );

                for ( auto timeBounds : timeBoundsVector )
                {
                    for ( auto observableIt : observationSetList_ )
                    {
                        std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;
                        if ( observationSetsIndices.count( observableIt.first ) )
                        {
                            indicesPerObservable = observationSetsIndices.at( observableIt.first );
                        }

                        for ( auto linkEndsIt : observableIt.second )
                        {
                            for ( unsigned int k = 0 ; k < linkEndsIt.second.size( ) ; k++ )
                            {
                                bool isInTimeBounds = false;
                                if ( ( linkEndsIt.second.at( k )->getTimeBounds( ).first >= timeBounds.first ) &&
                                     ( linkEndsIt.second.at( k )->getTimeBounds( ).second <= timeBounds.second ) )
                                {
                                    isInTimeBounds = true;
                                }
                                if ( ( isInTimeBounds && ( !observationParser->useOppositeCondition( ) ) ) ||
                                     ( !isInTimeBounds && ( observationParser->useOppositeCondition( ) ) ) )
                                {
                                    indicesPerObservable[ linkEndsIt.first ].push_back( k );
                                }
                            }
                        }
                        observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                    }
                }
                break;
            }
            case multi_type_parser:
            {
                std::vector< std::shared_ptr< ObservationCollectionParser > > observationParsers =
                        std::dynamic_pointer_cast< ObservationCollectionMultiTypeParser >( observationParser )->getObservationParsers_( );

                for ( auto parser : observationParsers )
                {
                    std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > currentObservationSetsIndices = getSingleObservationSetsIndices( parser );

                    for ( auto observableIt : currentObservationSetsIndices )
                    {
                        if ( observationSetsIndices.count( observableIt.first ) == 0 )
                        {
                            observationSetsIndices[ observableIt.first ] = observableIt.second;
                        }
                        else
                        {
                            for ( auto linkEndsIt : observableIt.second )
                            {
                                if ( observationSetsIndices.at( observableIt.first ).count( linkEndsIt.first ) == 0 )
                                {
                                    observationSetsIndices.at( observableIt.first )[ linkEndsIt.first ] = linkEndsIt.second;
                                }
                                else
                                {
                                    std::vector< unsigned int > indices = observationSetsIndices.at( observableIt.first ).at( linkEndsIt.first );
                                    for ( auto index : linkEndsIt.second )
                                    {
                                        if ( std::count( indices.begin( ), indices.end( ), index ) == 0 )
                                        {
                                            observationSetsIndices.at( observableIt.first ).at( linkEndsIt.first ).push_back( index );
                                        }
                                    }
                                }
                            }

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

    const std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > getSingleObservationSets(
            const std::shared_ptr< ObservationCollectionParser > observationParser )
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > singleObservationSetsIndices =
                getSingleObservationSetsIndices( observationParser );

        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObservationSets;
        for ( auto observableIt : singleObservationSetsIndices )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto index : linkEndsIt.second )
                {
                    singleObservationSets.push_back( observationSetList_.at( observableIt.first ).at( linkEndsIt.first ).at( index ) );
                }
            }
        }
        return singleObservationSets;
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

        linkDefinitionsPerObservable_.clear( );
        linkEndIds_.clear( );
        inverseLinkEndIds_.clear( );
        observationSetStartAndSizePerLinkEndIndex_.clear( );

        int observationCounter = 0;
        int maximumStationId = 0;
        int currentStationId;

        for( auto observationIterator : observationSetList_ )
        {
            ObservableType currentObservableType = observationIterator.first;
            int observableSize = getObservableSize( currentObservableType );

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


//bool isBodyInLinkEnds( const LinkEnds& linkEnds, const std::string bodyName )
//{
//    bool inLinkEnds = false;
//    for ( auto it : linkEnds )
//    {
//        if ( it.second.bodyName_ == bodyName )
//        {
//            inLinkEnds = true;
//        }
//    }
//    return inLinkEnds;
//}
//
//bool isGroundStationInLinkEnds( const LinkEnds& linkEnds, const std::string stationName )
//{
//    bool inLinkEnds = false;
//    for ( auto it : linkEnds )
//    {
//        if ( it.second.stationName_ == stationName )
//        {
//            inLinkEnds = true;
//        }
//    }
//    return inLinkEnds;
//}


template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > createNewObservationCollection(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< ObservationCollectionParser > observationParser )
{
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( observationCollection->getSingleObservationSets( observationParser ) );
}


template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void setConstantWeight(
        std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const double weight = 1.0,
        const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
{
    std::vector< std::pair< int, int > > observationSetIndices;
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObsSets = observationCollection->getSingleObservationSets( observationParser );
    for ( auto singleObsSet : singleObsSets )
    {
        singleObsSet->setConstantWeight( weight );
    }
}

template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void setConstantWeight(
        std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< ObservationCollectionParser > observationParser,
        const Eigen::VectorXd weight )
{
    std::vector< std::pair< int, int > > observationSetIndices;
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObsSets = observationCollection->getSingleObservationSets( observationParser );
    for ( auto singleObsSet : singleObsSets )
    {
        singleObsSet->setConstantWeight( weight );
    }
}

template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void setConstantWeightPerObservable(
        std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser )
{
    for ( auto parserIt : weightsPerObservationParser )
    {
        setConstantWeight( observationCollection, parserIt.second, parserIt.first );
    }
}

template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void setConstantWeightPerObservable(
        std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd > weightsPerObservationParser )
{
    for ( auto parserIt : weightsPerObservationParser )
    {
        setConstantWeight( observationCollection, parserIt.second, parserIt.first );
    }
}

template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void setTabulatedWeights(
        std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< ObservationCollectionParser > observationParser,
        const Eigen::VectorXd tabulatedWeights )
{
    std::vector< std::pair< int, int > > observationSetIndices;
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObsSets = observationCollection->getSingleObservationSets( observationParser );
    for ( auto singleObsSet : singleObsSets )
    {
        singleObsSet->setTabulatedWeights( tabulatedWeights );
    }
}

template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void setTabulatedWeights(
        std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd > weightsPerObservationParser )
{
    for ( auto parserIt : weightsPerObservationParser )
    {
        setTabulatedWeights( observationCollection, parserIt.first, parserIt.second );
    }
}

//void splitObservationCollectionTimeBounds( const std::pair< double, double >& timeBounds,
//                                           std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > observationCollection,
//                                                   std::shared_ptr< ObservationCollectionParser > parser = std::make_shared< ObservationCollecterParser<  > >( ) )
//{
//    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObservationSets = getSingleObservationSets( observationCollection, parser );
//
//
//}




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
