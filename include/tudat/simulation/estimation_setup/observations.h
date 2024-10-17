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
#include "tudat/simulation/estimation_setup/observationsProcessing.h"

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


        singleObservationSize_ = getObservableSize( observableType );

        // Initialise weights
        for ( int k = 0 ; k < numberOfObservations_ ; k++ )
        {
            weights_.push_back( Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( singleObservationSize_, 1 ) );
        }

        // Initialise residuals
        for ( int k = 0 ; k < numberOfObservations_ ; k++ )
        {
            residuals_.push_back( Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_, 1 ) );
        }

        // Sort observations and metadata per observation time
        orderObservationsAndMetadata( );

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

    void setLinkEnds( LinkDefinition& linkEnds )
    {
        linkEnds_ = linkEnds;
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservations( )
    {
        return observations_;
    }

    void setObservations( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations )
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

    void setObservations( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationsVector )
    {
        if ( observationsVector.size( ) != numberOfObservations_ * singleObservationSize_ )
        {
            throw std::runtime_error( "Error when resetting observations, number of observations is incompatible." );
        }

        observations_.clear( );
        for ( int k = 0 ; k < numberOfObservations_ ; k++ )
        {
            observations_.push_back( observationsVector.segment( k * singleObservationSize_, singleObservationSize_ ) );
        }
        if ( !observations_.empty( ) )
        {
            singleObservationSize_ = observations_.at( 0 ).size( );
        }
    }

    void setResiduals( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals )
    {
        if ( residuals.size( ) != static_cast< unsigned int >( numberOfObservations_ ) )
        {
            throw std::runtime_error( "Error when setting residuals, number of observations is inconsistent." );
        }
        residuals_ = residuals;
    }

    void setResiduals( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualsVector )
    {
        if ( residualsVector.size( ) != numberOfObservations_ * singleObservationSize_ )
        {
            throw std::runtime_error( "Error when setting residuals, number of observations is inconsistent." );
        }

        residuals_.clear( );
        for ( int k = 0 ; k < numberOfObservations_ ; k++ )
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

    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& getWeightsReference( )
    {
        return weights_;
    }

    Eigen::VectorXd getWeightsVector( ) const
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > weightsVector =
                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( singleObservationSize_ * numberOfObservations_, 1 );
        for( int i = 0; i < numberOfObservations_ ; i++ )
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
        for( int i = 0; i < numberOfObservations_ ; i++ )
        {
            residualsVector.block( i * singleObservationSize_, 0, singleObservationSize_, 1 ) = residuals_.at( i );
        }
        return residualsVector;
    }

    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& getResidualsReference( )
    {
        return residuals_;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidual( unsigned int index ) const
    {
        if ( index >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when retrieving single observation residual, required index incompatible with number of observations." );
        }
        return residuals_.at( index );
    }


    void setConstantWeight( const double weight )
    {
        for ( int k = 0 ; k < numberOfObservations_ ; k++ )
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
        if( weightsVector.rows( ) != singleObservationSize_ * static_cast< int >( observations_.size( ) ) )
        {
            throw std::runtime_error( "Error when setting weights in single observation set, sizes are incompatible." );
        }

        for ( int k = 0 ; k < numberOfObservations_ ; k++ )
        {
            for ( int i = 0 ; i < singleObservationSize_ ; i++ )
            {
                weights_.at( k )[ i ] = weightsVector[ k * singleObservationSize_ + i ];
            }
        }
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
        for( int i = 0; i < numberOfObservations_ ; i++ )
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

    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > getFilteredObservationSet( ) const
    {
        return filteredObservationSet_;
    }


    void removeSingleObservation( unsigned int indexToRemove )
    {
        if ( static_cast< int >( indexToRemove ) >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when removing single observation from SingleObservationSet, index incompatible with number of observations." );
        }

        // Update observations
        observations_.erase( observations_.begin( ) + indexToRemove );
        observationTimes_.erase( observationTimes_.begin( ) + indexToRemove );
        residuals_.erase( residuals_.begin( ) + indexToRemove );
        weights_.erase( weights_.begin( ) + indexToRemove );

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

    void filterObservations( const std::shared_ptr< ObservationFilterBase > observationFilter, const bool saveFilteredObservations = true )
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

        int nbObservationsToTest = ( observationFilter->filterOut( ) ? numberOfObservations_ : filteredObservationSet_->getNumberOfObservables( ) );
        bool useOppositeCondition = observationFilter->useOppositeCondition( );

        std::vector< unsigned int > indicesToRemove;
        switch ( observationFilter->getFilterType( ) )
        {
            case residual_filtering:
            {
                Eigen::VectorXd residualCutOff = Eigen::VectorXd::Zero( singleObservationSize_ );
                if ( std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter ) != nullptr )
                {
                    residualCutOff = std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter )->getFilterValue( ) * Eigen::VectorXd::Ones( singleObservationSize_ );
                }
                else if ( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter ) != nullptr )
                {
                    if ( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( ).size( ) != singleObservationSize_ )
                    {
                        throw std::runtime_error( "Error when performing residual filtering, size of the residual cut off vector inconsistent with observable size." );
                    }
                    residualCutOff = std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( );
                }

                for ( int j = 0 ; j < nbObservationsToTest ; j++ )
                {
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservationResidual =
                            ( observationFilter->filterOut( ) ? residuals_.at( j ) : filteredObservationSet_->getResidual( j ) );
                    bool removeObservation = false;
                    for( int k = 0 ; k < singleObservationSize_ ; k++ )
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
            case absolute_value_filtering:
            {
                Eigen::VectorXd absoluteValueCutOff = Eigen::VectorXd::Zero( singleObservationSize_ );
                if ( std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter ) != nullptr )
                {
                    absoluteValueCutOff = std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter )->getFilterValue( ) * Eigen::VectorXd::Ones( singleObservationSize_ );
                }
                else if ( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter ) != nullptr )
                {
                    if ( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( ).size( ) != singleObservationSize_ )
                    {
                        throw std::runtime_error( "Error when performing observation value filtering, size of the filter value inconsistent with observable size." );
                    }
                    absoluteValueCutOff = std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( );
                }

                for ( int j = 0 ; j < nbObservationsToTest ; j++ )
                {
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservation =
                            ( observationFilter->filterOut( ) ? observations_.at( j ) : filteredObservationSet_->getObservation( j ) );
                    bool removeObservation = false;
                    for( int k = 0 ; k < singleObservationSize_ ; k++ )
                    {
                        if ( ( !useOppositeCondition && ( singleObservation[ k ] > absoluteValueCutOff[ k ] ) ) ||
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
            case epochs_filtering:
            {
                std::vector< double > filterEpochs = std::dynamic_pointer_cast< ObservationFilter< std::vector< double > > >( observationFilter )->getFilterValue( );
                for ( int j = 0 ; j < nbObservationsToTest ; j++ )
                {
                    TimeType singleObservationTime = ( observationFilter->filterOut( ) ? observationTimes_.at( j ) : filteredObservationSet_->getObservationTime( j ) );
                    if ( ( !useOppositeCondition && ( std::count( filterEpochs.begin( ), filterEpochs.end( ), singleObservationTime ) > 0 ) ) ||
                    ( useOppositeCondition && ( std::count( filterEpochs.begin( ), filterEpochs.end( ), singleObservationTime ) == 0 ) ) )
                    {
                        indicesToRemove.push_back( j );
                    }
                }
                break;
            }
            case time_bounds_filtering:
            {
                std::pair< double, double > timeBounds = std::dynamic_pointer_cast< ObservationFilter< std::pair< double, double > > >( observationFilter )->getFilterValue( );
                for ( int j = 0 ; j < nbObservationsToTest ; j++ )
                {
                    TimeType singleObservationTime = ( observationFilter->filterOut( ) ? observationTimes_.at( j ) : filteredObservationSet_->getObservationTime( j ) );
                    if ( ( !useOppositeCondition && ( ( singleObservationTime >= timeBounds.first ) && ( singleObservationTime <= timeBounds.second ) ) ) ||
                            ( useOppositeCondition && ( ( singleObservationTime < timeBounds.first ) || ( singleObservationTime > timeBounds.second ) ) ) )
                    {
                        indicesToRemove.push_back( j );
                    }
                }
                break;
            }
            case dependent_variable_filtering:
            {
                throw std::runtime_error( "Error when filtering observations based on dependent variable, filter type not yet implemented." );
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
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType >& times,
            const std::vector< Eigen::VectorXd >& dependentVariables = { },
            const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights = { },
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals = { },
            const bool sortObservations = true )
    {
        if ( ( observations.size( ) != times.size( ) ) ||
                ( weights.size( ) > 0 && ( observations.size( ) != weights.size( ) ) ) ||
                ( residuals.size( ) > 0 && ( observations.size( ) != residuals.size( ) ) ) ||
                ( dependentVariables.size( ) > 0 && ( observations.size( ) != dependentVariables.size( ) ) ) )
        {
            throw std::runtime_error( "Error when adding observations to SingleObservationSet, input sizes are inconsistent." );
        }

        for ( unsigned int k = 0 ; k < observations.size( ) ; k++ )
        {
            if ( observations.at( k ).size( ) != singleObservationSize_ )
            {
                throw std::runtime_error( "Error when adding observations to SingleObservationSet, new observation size is inconsistent." );
            }

            observations_.push_back( observations.at( k ) );
            observationTimes_.push_back( times.at( k ) );

            // If residuals are provided as inputs
            if ( residuals.size( ) > 0 )
            {
                if ( residuals.at( k ).size( ) != singleObservationSize_ )
                {
                    throw std::runtime_error( "Error when adding observations to SingleObservationSet, new residual size is inconsistent." );
                }
                residuals_.push_back( residuals.at( k ) );
            }
            else // Otherwise, set to zero by default
            {
                residuals_.push_back( Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_, 1 )  );
            }

            // If weights are provided as inputs
            if ( weights.size( ) > 0 )
            {
                if ( weights.at( k ).size( ) != singleObservationSize_ )
                {
                    throw std::runtime_error( "Error when adding observations to SingleObservationSet, new weight size is inconsistent." );
                }
                weights_.push_back( weights.at( k ) );
            }
            else // Otherwise, set to one by default
            {
                weights_.push_back( Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( singleObservationSize_, 1 ) );
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

        // Update time bounds
        updateTimeBounds( );
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
                std::multimap< TimeType, Eigen::VectorXd > observationsDependentVariablesMap;
                for( unsigned int i = 0; i < observationsDependentVariables_.size( ); i++ )
                {
                    observationsDependentVariablesMap.insert( { observationTimes_.at( i ), observationsDependentVariables_.at( i ) } );
                }
                observationsDependentVariables_ = utilities::createVectorFromMultiMapValues( observationsDependentVariablesMap );
            }


            if( static_cast< int >( weights_.size( ) ) != numberOfObservations_ )
            {
                throw std::runtime_error( "Error when making SingleObservationSet, weights size is incompatible after time ordering" );
            }
            std::multimap< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > weightsMap;
            for( unsigned int i = 0; i < weights_.size( ); i++ )
            {
                weightsMap.insert( { observationTimes_.at( i ), weights_.at( i ) } );
            }
            weights_ = utilities::createVectorFromMultiMapValues( weightsMap );

            if( static_cast< int >( residuals_.size( ) ) != numberOfObservations_ )
            {
                throw std::runtime_error( "Error when making SingleObservationSet, residuals size is incompatible after time ordering" );
            }
            std::multimap< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residualsMap;
            for( unsigned int i = 0; i < residuals_.size( ); i++ )
            {
                residualsMap.insert( {  observationTimes_.at( i ), residuals_.at( i ) } );
            }
            residuals_ = utilities::createVectorFromMultiMapValues( residualsMap );
        }
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
                weights.push_back( weights_.at( index ) );
                residuals.push_back( residuals_.at( index ) );

                // If dependent variables not empty
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
                if ( static_cast< int >( index ) >= getNumberOfFilteredObservations( ) )
                {
                    throw std::runtime_error( "Error when moving observation back from filtered observation set, index incompatible with number of observations." );
                }

                observations.push_back( filteredObservationSet_->getObservations( ).at( index ) );
                times.push_back( filteredObservationSet_->getObservationTimes( ).at( index ) );
                weights.push_back( filteredObservationSet_->getWeights( ).at( index ) );
                residuals.push_back( filteredObservationSet_->getResiduals( ).at( index ) );

                // If dependent variables are set
                if ( filteredObservationSet_->getObservationsDependentVariables( ).size( ) > 0 )
                {
                    dependentVariables.push_back( filteredObservationSet_->getObservationsDependentVariables( ).at( index ) );
                }

            }
            addObservations( observations, times, dependentVariables, weights, residuals, true );
            filteredObservationSet_->removeObservations( indices );
        }
    }

    const ObservableType observableType_;

    LinkDefinition linkEnds_;

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

template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > filterObservations(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > singleObservationSet,
        const std::shared_ptr< ObservationFilterBase > observationFilter, const bool saveFilteredObservations = false )
{
    if ( !observationFilter->filterOut( ) )
    {
        throw std::runtime_error( "Error when creating new single observation set post-filtering, the filterOut option should be set to true" );
    }
    // Create new observation set
    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > newObservationSet =
            std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                    singleObservationSet->getObservableType( ), singleObservationSet->getLinkEnds( ), singleObservationSet->getObservationsReference( ),
                    singleObservationSet->getObservationTimesReference( ), singleObservationSet->getReferenceLinkEnd( ), singleObservationSet->getObservationsDependentVariablesReference( ),
                    singleObservationSet->getDependentVariableCalculator( ), singleObservationSet->getAncilliarySettings( ) );
    newObservationSet->setTabulatedWeights( singleObservationSet->getWeightsVector( ) );
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals = singleObservationSet->getResidualsReference( );
    newObservationSet->setResiduals( singleObservationSet->getResidualsReference( ) /*residuals*/ );

    // Filter observations from new observation set
    newObservationSet->filterObservations( observationFilter, saveFilteredObservations );

    return newObservationSet;
}

template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > splitObservationSet(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > observationSet,
        const std::shared_ptr< ObservationSetSplitterBase > observationSetSplitter )
{
    if ( observationSet->getFilteredObservationSet( ) != nullptr )
    {
        std::cerr << "Warning when splitting single observation set, the filtered observation set pointer is not empty and "
                     " any filtered observation will be lost after splitting." ;
    }

    std::vector< unsigned int > rawStartIndicesNewSets = { 0 };
    std::vector< TimeType > observationTimes = observationSet->getObservationTimes( );
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations = observationSet->getObservations( );
    std::vector< Eigen::VectorXd > dependentVariables = observationSet->getObservationsDependentVariables( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > weightsVector = observationSet->getWeightsVector( );
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals = observationSet->getResiduals( );

    switch ( observationSetSplitter->getSplitterType( ) )
    {
        case time_tags_splitter:
        {
            std::vector< double > timeTags = std::dynamic_pointer_cast< ObservationSetSplitter< std::vector< double > > >( observationSetSplitter )->getSplitterValue( );
            for ( auto currentTimeTag : timeTags )
            {
                if ( currentTimeTag > observationSet->getTimeBounds( ).first )
                {
                    bool detectedStartSet = false;
                    unsigned int indexObs = rawStartIndicesNewSets.at( rawStartIndicesNewSets.size( ) - 1 );
                    while ( !detectedStartSet && indexObs < observationTimes.size( ) )
                    {
                        if ( observationTimes.at( indexObs ) > currentTimeTag )
                        {
                            rawStartIndicesNewSets.push_back( indexObs );
                            detectedStartSet = true;
                        }
                        indexObs++;
                    }
                }
            }
            rawStartIndicesNewSets.push_back( observationTimes.size( ) );
            break;
        }
        case time_interval_splitter:
        {
            double maxTimeInterval = std::dynamic_pointer_cast< ObservationSetSplitter< double > >( observationSetSplitter )->getSplitterValue( );
            for( unsigned int i = 1; i < observationTimes.size( ); i++ )
            {
                if( ( observationTimes.at( i ) - observationTimes.at( i - 1 ) ) > maxTimeInterval )
                {
                    rawStartIndicesNewSets.push_back( i );
                }
            }
            rawStartIndicesNewSets.push_back( observationTimes.size( ) );
            break;
        }
        case time_span_splitter:
        {
            double maxTimeSpan = std::dynamic_pointer_cast< ObservationSetSplitter< double > >( observationSetSplitter )->getSplitterValue( );
            if ( observationSet->getNumberOfObservables( ) > 0 )
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
                rawStartIndicesNewSets.push_back( observationTimes.size( ) );
            }
            break;
        }
        case nb_observations_splitter:
        {
            int maxNbObs = std::dynamic_pointer_cast< ObservationSetSplitter< int > >( observationSetSplitter )->getSplitterValue( );
            if ( maxNbObs < observationSetSplitter->getMinNumberObservations( ) )
            {
                throw std::runtime_error( "Error when splitting observation sets, the maximum number of observations cannot be smaller than the minimum number of observations." );
            }
            for ( int ind = maxNbObs ; ind < observationSet->getNumberOfObservables( ) ; ind+=maxNbObs )
            {
                rawStartIndicesNewSets.push_back( ind );
            }
            rawStartIndicesNewSets.push_back( observationTimes.size( ) );
            break;
        }
        case dependent_variables_splitter:
        {
            throw std::runtime_error( "Error when splitting observation set based on dependent variable, splitter type not yet implemented." );
            break;
        }
        default:
            throw std::runtime_error( "Observation set splitter type not recognised." );
    }

    // Check that the minimum number of observations is met
    std::vector< std::pair< unsigned int, unsigned int > > indicesNewSets;
    for( unsigned int j = 1; j < rawStartIndicesNewSets.size( ); j++ )
    {
        if( static_cast< int >( rawStartIndicesNewSets.at( j ) - rawStartIndicesNewSets.at( j - 1 ) ) >= observationSetSplitter->getMinNumberObservations( ) )
        {
            indicesNewSets.push_back( std::make_pair( rawStartIndicesNewSets.at( j-1 ), rawStartIndicesNewSets.at( j ) - rawStartIndicesNewSets.at( j-1 ) ) );
        }
    }

    // Split current observation set based on indices
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > newObsSets;
    for ( unsigned int k = 0 ; k < indicesNewSets.size( ) ; k++ )
    {
        unsigned int startIndex = indicesNewSets.at( k ).first;
        unsigned int sizeCurrentSet = indicesNewSets.at( k ).second;

        std::vector< Eigen::VectorXd > newDependentVariables;
        if ( !dependentVariables.empty( ) )
        {
            newDependentVariables = utilities::getStlVectorSegment(
                    observationSet->getObservationsDependentVariablesReference( ), startIndex, sizeCurrentSet );
        }

        std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > newSet = std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                observationSet->getObservableType( ), observationSet->getLinkEnds( ),
                utilities::getStlVectorSegment( observationSet->getObservationsReference( ), startIndex, sizeCurrentSet ),
                utilities::getStlVectorSegment( observationSet->getObservationTimesReference( ), startIndex, sizeCurrentSet ),
                observationSet->getReferenceLinkEnd( ), newDependentVariables, observationSet->getDependentVariableCalculator( ),
                observationSet->getAncilliarySettings( ) );

        Eigen::Matrix< double, Eigen::Dynamic, 1 > newWeightsVector = weightsVector.segment( startIndex, sizeCurrentSet * observationSet->getSingleObservableSize( ) );
        newSet->setTabulatedWeights( newWeightsVector );

        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > newResiduals =
                utilities::getStlVectorSegment( observationSet->getResidualsReference( ), startIndex, sizeCurrentSet );
        newSet->setResiduals( newResiduals );

        newObsSets.push_back( newSet );
    }

    return newObsSets;
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

    void setObservations( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& newObservations )
    {
        if ( newObservations.size( ) != totalObservableSize_ )
        {
            throw std::runtime_error( "Error when resetting observations in ObservationCollection, size of input observation vector is inconsistent." );
        }

        unsigned int startIndexObsSet = 0;
        for ( auto observableIt : observationSetList_ )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto set : linkEndsIt.second )
                {
                    unsigned int sizeCurrentSet = set->getTotalObservationSetSize( );
                    set->setObservations( newObservations.segment( startIndexObsSet, sizeCurrentSet ) );
                    startIndexObsSet += sizeCurrentSet;
                }
            }
        }
    }

    void setResiduals( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& newResiduals )
    {
        if ( newResiduals.size( ) != totalObservableSize_ )
        {
            throw std::runtime_error( "Error when resetting observations in ObservationCollection, size of input observation vector is inconsistent." );
        }

        unsigned int startIndexObsSet = 0;
        for ( auto observableIt : observationSetList_ )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto set : linkEndsIt.second )
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
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObsSets = getSingleObservationSets( observationParser );
        if ( singleObsSets.empty( ) )
        {
            std::cerr << "Warning when setting observations, no single observation set found for specified observation parser. Weights not set";
        }
        int totalSizeAllObsSets = 0;
        for ( unsigned int k = 0 ; k < singleObsSets.size( ) ; k++ )
        {
            totalSizeAllObsSets += singleObsSets.at( k )->getTotalObservationSetSize( );
        }

        if ( observations.size( ) == totalSizeAllObsSets )
        {
            unsigned int startObsSet = 0;
            for ( auto obsSet : singleObsSets )
            {
                obsSet->setObservations( observations.segment( startObsSet, obsSet->getTotalObservationSetSize( ) ) );
                startObsSet += obsSet->getTotalObservationSetSize( );
            }
        }
        else
        {
            throw std::runtime_error( "Error when setting observations, the size of the input observation vector should be consistent with "
                                      "the combined size of all required observation sets." );
        }
    }

    void setObservations(
            const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observationsPerParser )
    {
        for ( auto parserIt : observationsPerParser )
        {
            setObservations( parserIt.second, parserIt.first );
        }
    }

    void setResiduals(
            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals,
            const std::shared_ptr< ObservationCollectionParser > observationParser )
    {
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObsSets = getSingleObservationSets( observationParser );
        if ( singleObsSets.empty( ) )
        {
            std::cerr << "Warning when setting residuals, no single observation set found for specified observation parser. Weights not set";
        }
        int totalSizeAllObsSets = 0;
        for ( unsigned int k = 0 ; k < singleObsSets.size( ) ; k++ )
        {
            totalSizeAllObsSets += singleObsSets.at( k )->getTotalObservationSetSize( );
        }

        if ( residuals.size( ) == totalSizeAllObsSets )
        {
            unsigned int startObsSet = 0;
            for ( auto obsSet : singleObsSets )
            {
                obsSet->setResiduals( residuals.segment( startObsSet, obsSet->getTotalObservationSetSize( ) ) );
                startObsSet += obsSet->getTotalObservationSetSize( );
            }
        }
        else
        {
            throw std::runtime_error( "Error when setting residuals, the size of the input residual vector should be consistent with "
                                      "the combined size of all required observation sets." );
        }
    }

    void setResiduals(
            const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residualsPerParser )
    {
        for ( auto parserIt : residualsPerParser )
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

    SortedObservationSets getObservationsSets( )
    {
        return observationSetList_;
    }

    const SortedObservationSets& getObservationsReference ( ) const
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

    Eigen::VectorXd getUnparsedConcatenatedWeights( ) const
    {
        return getConcatenatedWeights( );
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getResiduals(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;

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

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getConcatenatedResiduals(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals = getResiduals( observationParser );

        unsigned int totalObsSize = 0;
        for ( auto residual : residuals )
        {
            totalObsSize += residual.size( );
        }

        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > concatenatedResiduals = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObsSize, 1 );
        unsigned int obsIndex = 0;
        for ( unsigned int k = 0 ; k < residuals.size( ) ; k++ )
        {
            concatenatedResiduals.block( obsIndex, 0, residuals.at( k ).size( ), 1 ) = residuals.at( k );
            obsIndex += residuals.at( k ).size( );
        }

        return concatenatedResiduals;
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getComputedObservations(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > computedObservations;

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

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getConcatenatedComputedObservations(
            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > computedObservations = getComputedObservations( observationParser );

        unsigned int totalObsSize = 0;
        for ( auto obs : computedObservations )
        {
            totalObsSize += obs.size( );
        }

        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > concatenatedComputedObservations = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObsSize );
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

//                    bool addDependentVariables = false;
                    if ( dependentVariableCalculator != nullptr )
                    {
                        std::pair<int, int> variableIndices = dependentVariableCalculator->getDependentVariableIndices(
                            dependentVariable );

                        if( variableIndices.second > 0 )
                        {
//                            addDependentVariables = true;
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

    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > getSingleLinkAndTypeObservationSets(
            const ObservableType observableType,
            const LinkDefinition linkEndsDefinition )
    {
        std::vector< std::shared_ptr< ObservationCollectionParser > > multiTypeParserList;
        multiTypeParserList.push_back( observationParser( observableType ) );
        multiTypeParserList.push_back( observationParser( linkEndsDefinition.linkEnds_ ) );
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > observationSets =
                getSingleObservationSets( observationParser( multiTypeParserList, true ) );
        if ( observationSets.empty( ) )
        {
            throw std::runtime_error( " Error when getting single link observation sets, not observations of type "
                                      + std::to_string( observableType ) + " for given link ends." );
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
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observations = getConcatenatedObservations( observationParser( multiTypeParserList, true ) );
        if ( observations.size( ) == 0 )
        {
            throw std::runtime_error( " Error when getting single link observations, not observations of type "
                                      + std::to_string( observableType ) + " for given link ends." );
        }
        return observations;
    }

    std::vector< TimeType > getSingleLinkTimes(
            const ObservableType observableType,
            const LinkDefinition& linkEndsDefinition )
    {
        std::vector< std::shared_ptr< ObservationCollectionParser > > multiTypeParserList;
        multiTypeParserList.push_back( observationParser( observableType ) );
        multiTypeParserList.push_back( observationParser( linkEndsDefinition.linkEnds_ ) );
        std::vector< TimeType > observationTimes = getConcatenatedObservationTimes( observationParser( multiTypeParserList, true ) );
        if ( observationTimes.empty( ) )
        {
            throw std::runtime_error( " Error when getting single link observation times, not observations of type "
                                      + std::to_string( observableType ) + " for given link ends." );
        }
        return observationTimes;
    }

    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > getSingleLinkObservationsAndTimes(
            const ObservableType observableType,
            const LinkDefinition& linkEndsDefinition )
    {
        std::vector< std::shared_ptr< ObservationCollectionParser > > multiTypeParserList;
        multiTypeParserList.push_back( observationParser( observableType ) );
        multiTypeParserList.push_back( observationParser( linkEndsDefinition.linkEnds_ ) );
        std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > observationsAndTimes =
                getConcatenatedObservationsAndTimes( observationParser( multiTypeParserList, true ) );
        if ( observationsAndTimes.second.empty( ) )
        {
            throw std::runtime_error( " Error when getting single link observations and times, not observation set of type "
                                      + std::to_string( observableType ) + " for given link ends." );
        }
        return observationsAndTimes;
    }


    void setConstantWeight(
            const double weight = 1.0,
            const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
    {
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObsSets = getSingleObservationSets( observationParser );
        if ( singleObsSets.empty( ) )
        {
            std::cerr << "Warning when setting constant weights, no single observation set found for specified observation parser. Weights not set";
        }
        for ( auto singleObsSet : singleObsSets )
        {
            singleObsSet->setConstantWeight( weight );
        }
    }

    void setConstantWeight(
            const Eigen::VectorXd weight,
            const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
    {
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObsSets = getSingleObservationSets( observationParser );
        if ( singleObsSets.empty( ) )
        {
            std::cerr << "Warning when setting constant weights, no single observation set found for specified observation parser. Weights not set";
        }
        for ( auto singleObsSet : singleObsSets )
        {
            singleObsSet->setConstantWeight( weight );
        }
    }

    void setConstantWeightPerObservable(
            const std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser )
    {
        for ( auto parserIt : weightsPerObservationParser )
        {
            setConstantWeight( parserIt.second, parserIt.first );
        }
    }

    void setConstantWeightPerObservable(
            const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd > weightsPerObservationParser )
    {
        for ( auto parserIt : weightsPerObservationParser )
        {
            setConstantWeight( parserIt.second, parserIt.first );
        }
    }

    void setTabulatedWeights(
            const Eigen::VectorXd tabulatedWeights,
            const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
    {
        std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObsSets = getSingleObservationSets( observationParser );
        if ( singleObsSets.empty( ) )
        {
            std::cerr << "Warning when setting tabulated weights, no single observation set found for specified observation parser. Weights not set";
        }

        bool areObsSetsSameSize = true;
        int totalSizeAllObsSets = singleObsSets.at( 0 )->getTotalObservationSetSize( );
        for ( unsigned int k = 1 ; k < singleObsSets.size( ) ; k++ )
        {
            int currentObsSetSize = singleObsSets.at( k )->getTotalObservationSetSize( );
            totalSizeAllObsSets += currentObsSetSize;
            if ( currentObsSetSize != static_cast< int >( singleObsSets.at( 0 )->getTotalObservationSetSize( ) ) )
            {
                areObsSetsSameSize = false;
            }
        }

        unsigned int startObsSet = 0;
        for ( auto obsSet : singleObsSets )
        {
            if ( tabulatedWeights.size( ) == totalSizeAllObsSets )
            {
                Eigen::VectorXd singleSetWeights = tabulatedWeights.segment( startObsSet, obsSet->getTotalObservationSetSize( ) );
                startObsSet += obsSet->getTotalObservationSetSize( );
                obsSet->setTabulatedWeights( singleSetWeights );
            }
            else if ( areObsSetsSameSize && ( tabulatedWeights.size( ) == singleObsSets.at( 0 )->getTotalObservationSetSize( ) ) )
            {
                obsSet->setTabulatedWeights( tabulatedWeights );
            }
            else
            {
                throw std::runtime_error( "Error when setting tabulated weights, the size of the input weight vector should be consistent with "
                                          "either the size of each individual observation set, or the combined size of all required observation sets." );
            }
        }
    }

    void setTabulatedWeights(
            const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd > weightsPerObservationParser )
    {
        for ( auto parserIt : weightsPerObservationParser )
        {
            setTabulatedWeights( parserIt.second, parserIt.first );
        }
    }

    void filterObservations( const std::map< std::shared_ptr< ObservationCollectionParser >, std::shared_ptr< ObservationFilterBase > >& observationFilters,
                             const bool saveFilteredObservations = true )
    {
        // Parse all observation filters
        for ( auto filterIt : observationFilters )
        {
            std::shared_ptr< ObservationCollectionParser > observationParser = filterIt.first;
            std::shared_ptr< ObservationFilterBase > filter = filterIt.second;

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

   void filterObservations( std::shared_ptr< ObservationFilterBase > observationFilter,
                            std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ),
                            const bool saveFilteredObservations = true )
    {
        std::map< std::shared_ptr< ObservationCollectionParser >, std::shared_ptr< ObservationFilterBase > > observationFilters
                { { observationParser, observationFilter } };

        filterObservations( observationFilters, saveFilteredObservations );
    }

    void splitObservationSets( std::shared_ptr< ObservationSetSplitterBase > observationSetSplitter,
                               std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
    {
        // Retrieve single observation sets based on observation parser
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > singleObsSetsIndices = getSingleObservationSetsIndices( observationParser );

        for ( auto observableIt : singleObsSetsIndices )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > singleObsSets = observationSetList_.at( observableIt.first ).at( linkEndsIt.first );

                unsigned int obsSetCounter = 0;
                for ( auto indexSetToSplit : linkEndsIt.second )
                {
                    // Get new observation sets after splitting
                    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > newObsSets =
                            splitObservationSet( singleObsSets.at( indexSetToSplit + obsSetCounter ), observationSetSplitter );

                    // Remove original observation set
                    singleObsSets.erase( singleObsSets.begin( ) + ( indexSetToSplit + obsSetCounter ) );

                    // Replace by new sets
                    for ( auto newSet : newObsSets )
                    {
                        singleObsSets.insert( singleObsSets.begin( ) + ( indexSetToSplit + obsSetCounter ), newSet );
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

    void replaceSingleObservationSet( const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >& newSet,
                                      const unsigned int setIndex )
    {
        if ( observationSetList_.count( newSet->getObservableType( ) ) == 0 )
        {
            throw std::runtime_error( "Error when replacing single observation set in observation collection, observable type not found." );
        }
        if ( observationSetList_.at( newSet->getObservableType( ) ).count( newSet->getLinkEnds( ).linkEnds_ ) == 0 )
        {
            throw std::runtime_error( "Error when replacing single observation set in observation collection, link ends not found for given observable type." );
        }
        if ( setIndex > observationSetList_.at( newSet->getObservableType( ) ).at( newSet->getLinkEnds( ).linkEnds_ ).size( ) - 1 )
        {
            throw std::runtime_error( "Error when replacing single observation set in observation collection, "
                                      "set index exceeds number of sets found for given observable and link ends." );
        }
        observationSetList_.at( newSet->getObservableType( ) ).at( newSet->getLinkEnds( ).linkEnds_ ).at( setIndex ) = newSet;

        // Reset observation set indices and concatenated observations and times
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    void removeSingleObservationSets( const std::shared_ptr< ObservationCollectionParser > observationParser )
    {
        // Retrieve single observation sets based on observation parser
        removeSingleObservationSets( getSingleObservationSetsIndices( observationParser ) );
    }

    void removeEmptySingleObservationSets( )
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > indicesSetsToRemove;

        for ( auto observableIt : observationSetList_ )
        {
            std::map< LinkEnds, std::vector< unsigned int > > indicesToRemovePerLinkEnds;
            for ( auto linkEndsIt : observableIt.second )
            {
                // Identify empty observation sets
                for ( unsigned int k = 0 ; k < linkEndsIt.second.size( ) ; k++ )
                {
                    if ( linkEndsIt.second.at( k )->getNumberOfObservables( ) == 0 )
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

    void removeSingleObservationSets( const std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > >& indicesSetsToRemove )
    {
        // Parse observation set list and remove selected sets
        for ( auto observableIt : indicesSetsToRemove )
        {
            if ( observationSetList_.count( observableIt.first ) == 0 )
            {
                throw std::runtime_error( "Error when removing single observation sets, given observable type not found." );
            }
            for ( auto linkEndsIt : observableIt.second )
            {
                if ( observationSetList_.at( observableIt.first ).count( linkEndsIt.first ) == 0 )
                {
                    throw std::runtime_error( "Error when removing single observation sets, link ends not found for given observable type." );
                }
                unsigned int counterRemovedSets = 0;
                for ( auto indexToRemove : linkEndsIt.second )
                {
                    if ( indexToRemove > observationSetList_.at( observableIt.first ).at( linkEndsIt.first ).size( ) - 1 )
                    {
                        throw std::runtime_error( "Error when removing single observation sets, set index exceeds number of observation sets"
                                                  " for given observation type and link ends." );
                    }
                    observationSetList_.at( observableIt.first ).at( linkEndsIt.first ).erase(
                            observationSetList_.at( observableIt.first ).at( linkEndsIt.first ).begin( ) + indexToRemove - counterRemovedSets );
                    counterRemovedSets += 1;
                }
            }
        }

        // Reset observation set indices and concatenated observations and times
        setObservationSetIndices( );
        setConcatenatedObservationsAndTimes( );
    }

    std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > getSingleObservationSetsIndices(
            const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) ) const
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
                            if ( indicesPerObservable.size( ) > 0 )
                            {
                                observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                            }
                        }
                    }
                }
                break;
            }
            case link_end_string_parser:
            {
                std::shared_ptr< ObservationCollectionLinkEndStringParser > linkEndIdObservationParser =
                        std::dynamic_pointer_cast< ObservationCollectionLinkEndStringParser >( observationParser );
                std::vector< std::string > linkEndsNames = linkEndIdObservationParser->getLinkEndNames( );
                bool isReferencePoint = linkEndIdObservationParser->isReferencePoint( );

                // Retrieve all body names
                std::vector< std::string > allLinkEndsNames;
                if ( !isReferencePoint )
                {
                    allLinkEndsNames = getBodiesInLinkEnds( );
                }
                else
                {
                    allLinkEndsNames = getReferencePointsInLinkEnds( );
                }

                for ( auto name : allLinkEndsNames )
                {
                    if ( ( ( std::count( linkEndsNames.begin( ), linkEndsNames.end( ), name ) ) && ( !observationParser->useOppositeCondition( ) ) ) ||
                         ( ( std::count( linkEndsNames.begin( ), linkEndsNames.end( ), name ) == 0 ) && ( observationParser->useOppositeCondition( ) ) ) )
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

                                if ( ( !isReferencePoint && isBodyInLinkEnds/*( linkEndsIt.first, name )*/ ) ||
                                     ( isReferencePoint && isGroundStationInLinkEnds/*( linkEndsIt.first, name )*/ ) )
                                {
                                    for ( unsigned int k = 0 ; k < linkEndsIt.second.size( ) ; k++ )
                                    {
                                        indicesPerObservable[ linkEndsIt.first ].push_back( k );
                                    }
                                }
                            }
                            if ( indicesPerObservable.size( ) > 0 )
                            {
                                observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                            }
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
                        if ( indicesPerObservable.size( ) > 0 )
                        {
                            observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                        }
                    }
                }
                break;
            }
            case ancillary_settings_parser:
            {
                std::vector< std::shared_ptr< ObservationAncilliarySimulationSettings > > ancillarySettings =
                        std::dynamic_pointer_cast< ObservationCollectionAncillarySettingsParser >( observationParser )->getAncillarySettings( );

                for ( auto setting : ancillarySettings )
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
                                bool identicalAncillarySettings = false;
                                if ( ( linkEndsIt.second.at( k )->getAncilliarySettings( )->getDoubleData( ) == setting->getDoubleData( ) ) &&
                                    ( linkEndsIt.second.at( k )->getAncilliarySettings( )->getDoubleVectorData( ) == setting->getDoubleVectorData( ) ) )
                                {
                                    identicalAncillarySettings = true;
                                }
                                if ( ( identicalAncillarySettings && !observationParser->useOppositeCondition( ) ) ||
                                ( !identicalAncillarySettings && observationParser->useOppositeCondition( ) ) )
                                {
                                    indicesPerObservable[ linkEndsIt.first ].push_back( k );
                                }
                            }
                        }
                        if ( indicesPerObservable.size( ) > 0 )
                        {
                            observationSetsIndices[ observableIt.first ] = indicesPerObservable;
                        }
                    }
                }
                break;
            }
            case multi_type_parser:
            {
                std::shared_ptr< ObservationCollectionMultiTypeParser > multiTypeParser = std::dynamic_pointer_cast< ObservationCollectionMultiTypeParser >( observationParser );
                std::vector< std::shared_ptr< ObservationCollectionParser > > observationParsers = multiTypeParser->getObservationParsers_( );

                bool areConditionsCombined = multiTypeParser->areConditionsCombined( );

                if ( !areConditionsCombined )
                {
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
                }

                else
                {
                    // First retrieve all observation sets
                    std::vector< std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > > allObservationSetsIndices;
                    for ( auto parser : observationParsers )
                    {
                        allObservationSetsIndices.push_back( getSingleObservationSetsIndices( parser ) );
                    }

                    if ( allObservationSetsIndices.size( ) > 0 )
                    {
                        std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > originalSetsIndices = allObservationSetsIndices.at( 0 );

                        // Retrieve common observable types
                        std::vector< ObservableType > originalObsTypes = utilities::createVectorFromMapKeys( originalSetsIndices );
                        std::vector< ObservableType > commonObsTypes;
                        for ( unsigned int k = 1 ; k < allObservationSetsIndices.size( ) ; k++ )
                        {
                            std::vector< ObservableType > currentObsTypes = utilities::createVectorFromMapKeys( allObservationSetsIndices.at( k ) );
                            std::set_intersection( originalObsTypes.begin( ), originalObsTypes.end( ), currentObsTypes.begin( ), currentObsTypes.end( ),
                                                   std::back_inserter( commonObsTypes ) );
                        }

                        for ( auto obsType : commonObsTypes )
                        {
                            std::map< LinkEnds, std::vector< unsigned int > > indicesPerObservable;

                            // For given observable type, retrieve common link ends
                            std::vector< LinkEnds > originalLinkEndsList = utilities::createVectorFromMapKeys( originalSetsIndices.at( obsType ) );
                            std::vector< LinkEnds > commonLinkEndsList;
                            for ( unsigned int k = 1 ; k < allObservationSetsIndices.size( ) ; k++ )
                            {
                                std::vector< LinkEnds > currentLinkEndsList = utilities::createVectorFromMapKeys( allObservationSetsIndices.at( k ).at( obsType ) );
                                std::set_intersection( originalLinkEndsList.begin( ), originalLinkEndsList.end( ), currentLinkEndsList.begin( ), currentLinkEndsList.end( ),
                                                       std::back_inserter( commonLinkEndsList ) );
                            }

                            for ( auto linkEnds : commonLinkEndsList )
                            {
                                // For given observable type and link ends, retrieve common observation set indices
                                std::vector< unsigned int > originalIndices = originalSetsIndices.at( obsType ).at( linkEnds );
                                std::vector< unsigned int > commonIndices;
                                for ( unsigned int k = 1 ; k < allObservationSetsIndices.size( ) ; k++ )
                                {
                                    std::vector< unsigned int > currentIndices = allObservationSetsIndices.at( k ).at( obsType ).at( linkEnds );
                                    std::set_intersection( originalIndices.begin( ), originalIndices.end( ), currentIndices.begin( ), currentIndices.end( ),
                                                           std::back_inserter( commonIndices ) );
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

    const std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > getSingleObservationSets(
            const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
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

    void printObservationSetsStartAndSize( ) const
    {
        std::cout << "Observation collection structure: start and size of individual observation sets" << std::endl;
        for ( auto observableIt : observationSetStartAndSize_ )
        {
            for ( auto linkEndsIt : observableIt.second )
            {
                for ( auto indices : linkEndsIt.second )
                {
                    std::cout << "start index " << indices.first << " - size " << indices.second << " : " << "type " << observableIt.first << ", link ends " << getLinkEndsString( linkEndsIt.first ) << std::endl;
                }
            }
        }
        std::cout << "\n\n";
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


    SortedObservationSets observationSetList_;

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
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > createNewObservationCollection(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
{
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > oldObservationSets =
            observationCollection->getSingleObservationSets( observationParser );

    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > newSingleObservationSets;
    for ( auto oldObsSet : oldObservationSets )
    {
        std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > newObsSet = std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                oldObsSet->getObservableType( ), oldObsSet->getLinkEnds( ), oldObsSet->getObservations( ), oldObsSet->getObservationTimes( ), oldObsSet->getReferenceLinkEnd( ),
                oldObsSet->getObservationsDependentVariablesReference( ), oldObsSet->getDependentVariableCalculator( ), oldObsSet->getAncilliarySettings( ) );
        newObsSet->setTabulatedWeights( oldObsSet->getWeightsVector( ) );
        newObsSet->setResiduals( oldObsSet->getResiduals( ) );

        newSingleObservationSets.push_back( newObsSet );
    }
    return std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >( newSingleObservationSets );
}


template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > filterObservations(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::map< std::shared_ptr< ObservationCollectionParser >, std::shared_ptr< ObservationFilterBase > >& observationFilters )
{
    // Create new observation collection
    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > newObservationCollection = createNewObservationCollection( observationCollection );

    // Apply filtering to observation collection
    newObservationCollection->filterObservations( observationFilters, false );

    return newObservationCollection;
}

template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > filterObservations(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< ObservationFilterBase > observationFilter,
        const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
{
    // Create new observation collection
    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > newObservationCollection = createNewObservationCollection( observationCollection );

    // Apply filtering to observation collection
    newObservationCollection->filterObservations( observationFilter, observationParser, false );

    return newObservationCollection;
}

template< typename ObservationScalarType = double, typename TimeType = double,
        typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > splitObservationSets(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< ObservationSetSplitterBase > observationSetSplitter,
        const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
{
    // Create new observation collection
    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > newObservationCollection = createNewObservationCollection( observationCollection );

    // Split observation sets within new observation collection
    newObservationCollection->splitObservationSets( observationSetSplitter, observationParser );

    return newObservationCollection;
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
    if( observationCollection->getObservationsSets( ).count( observableType ) != 0 )
    {
        std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > >
                observationsOfGivenType = observationCollection->getObservationsSets( ).at( observableType );

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
    if( observationCollection->getObservationsSets( ).count( observableType ) != 0 )
    {
        if( observationCollection->getObservationsSets( ).at( observableType ).count( linkEnds ) != 0 )
        {
            std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > fullObservationList =
                    observationCollection->getObservationsSets( ).at( observableType ).at( linkEnds );
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
                std::vector< Eigen::VectorXd >( ), nullptr, ancilliarySettings );
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


} // namespace observation_models

} // namespace tudat

#endif // TUDAT_OBSERVATIONS_H

















