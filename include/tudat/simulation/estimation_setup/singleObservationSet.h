/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SINGLE_OBSERVATION_SET_H
#define TUDAT_SINGLE_OBSERVATION_SET_H

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
class SingleObservationSet
{
   public:
    SingleObservationSet(
            const ObservableType observableType,
            const LinkDefinition& linkEnds,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
                    observations,
            const std::vector< TimeType > observationTimes,
            const LinkEndType referenceLinkEnd,
            const std::vector< Eigen::VectorXd >& observationsDependentVariables =
                    std::vector< Eigen::VectorXd >( ),
            const std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator >
                    dependentVariableCalculator = nullptr,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >
                    ancilliarySettings = nullptr ) :
        observableType_( observableType ), linkEnds_( linkEnds ), observations_( observations ),
        observationTimes_( observationTimes ), referenceLinkEnd_( referenceLinkEnd ),
        observationsDependentVariables_( observationsDependentVariables ),
        dependentVariableCalculator_( dependentVariableCalculator ),
        ancilliarySettings_( ancilliarySettings ), numberOfObservations_( observations_.size( ) )
    {
        if( dependentVariableCalculator_ != nullptr )
        {
            if( dependentVariableCalculator_->getObservableType( ) != observableType_ )
            {
                std::cout << dependentVariableCalculator_->getObservableType( ) << " "
                          << observableType_ << std::endl;
                throw std::runtime_error(
                        "Error when creating SingleObservationSet, "
                        "ObservationDependentVariableCalculator has incompatible type " );
            }

            if( !( dependentVariableCalculator_->getLinkEnds( ) == linkEnds ) )
            {
                throw std::runtime_error(
                        "Error when creating SingleObservationSet, "
                        "ObservationDependentVariableCalculator has incompatible link ends " );
            }
        }

        if( observations_.size( ) != observationTimes_.size( ) )
        {
            throw std::runtime_error(
                    "Error when making SingleObservationSet, input sizes are inconsistent." +
                    std::to_string( observations_.size( ) ) + ", " +
                    std::to_string( observationTimes_.size( ) ) );
        }

        for( unsigned int i = 1; i < observations.size( ); i++ )
        {
            if( observations.at( i ).rows( ) != observations.at( i - 1 ).rows( ) )
            {
                throw std::runtime_error(
                        "Error when making SingleObservationSet, input observables not of "
                        "consistent size." );
            }
        }

        singleObservationSize_ = getObservableSize( observableType );

        // Initialise weights
        for( unsigned int k = 0; k < numberOfObservations_; k++ )
        {
            weights_.push_back(
                    Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( singleObservationSize_, 1 ) );
        }

        // Initialise residuals
        for( unsigned int k = 0; k < numberOfObservations_; k++ )
        {
            residuals_.push_back( Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero(
                    singleObservationSize_, 1 ) );
        }

        // Check observation dependent variables size
        if( observationsDependentVariables_.size( ) > 0 )
        {
            if( observationsDependentVariables_.size( ) != numberOfObservations_ )
            {
                throw std::runtime_error(
                        "Error when creating SingleObservationSet, the size of the observation "
                        "dependent variables input should be consistent "
                        "with the number of observations." );
            }
            if( observationsDependentVariables_[ 0 ].size( ) !=
                dependentVariableCalculator_->getTotalDependentVariableSize( ) )
            {
                throw std::runtime_error(
                        "Error when creating SingleObservationSet, the size of the observation "
                        "dependent variables input "
                        "should be consistent with the total dependent variable size." );
            }
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

    void setObservations(
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
                    observations )
    {
        if( observations.size( ) != observations_.size( ) )
        {
            throw std::runtime_error(
                    "Error when resetting observations, number of observations is incompatible." );
        }
        observations_ = observations;
        if( !observations_.empty( ) )
        {
            singleObservationSize_ = observations_.at( 0 ).size( );
        }
    }

    void setObservations(
            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationsVector )
    {
        if( observationsVector.size( ) != numberOfObservations_ * singleObservationSize_ )
        {
            throw std::runtime_error(
                    "Error when resetting observations, number of observations is incompatible." );
        }

        observations_.clear( );
        for( unsigned int k = 0; k < numberOfObservations_; k++ )
        {
            observations_.push_back( observationsVector.segment( k * singleObservationSize_,
                                                                 singleObservationSize_ ) );
        }
        if( !observations_.empty( ) )
        {
            singleObservationSize_ = observations_.at( 0 ).size( );
        }
    }

    void setResiduals(
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
                    residuals )
    {
        if( residuals.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when setting residuals, number of observations is inconsistent." );
        }
        residuals_ = residuals;
    }

    void setResiduals(
            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualsVector )
    {
        if( residualsVector.size( ) != numberOfObservations_ * singleObservationSize_ )
        {
            throw std::runtime_error(
                    "Error when setting residuals, number of observations is inconsistent." );
        }

        residuals_.clear( );
        for( unsigned int k = 0; k < numberOfObservations_; k++ )
        {
            residuals_.push_back(
                    residualsVector.segment( k * singleObservationSize_, singleObservationSize_ ) );
        }
    }

    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
    getObservationsReference( )
    {
        return observations_;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservation(
            const unsigned int index )
    {
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when retrieving single observation, index is out of bounds" );
        }
        return observations_.at( index );
    }

    void setObservation( const unsigned int index,
                         Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observation )
    {
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when setting single observation value, index is out of bounds" );
        }
        if( observation.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error(
                    "Error when setting single observation value, the observation size is "
                    "inconsistent." );
        }
        observations_.at( index ) = observation;
    }

    std::vector< TimeType > getObservationTimes( )
    {
        return observationTimes_;
    }

    TimeType getObservationTime( unsigned int index ) const
    {
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when retrieving single observation time, required index incompatible "
                    "with number of observations." );
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

    unsigned int getNumberOfObservables( )
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
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero(
                        singleObservationSize_ * numberOfObservations_, 1 );
        for( unsigned int i = 0; i < observations_.size( ); i++ )
        {
            observationsVector.segment( i * singleObservationSize_, singleObservationSize_ ) =
                    observations_.at( i );
        }
        return observationsVector;
    }

    std::pair< TimeType, TimeType > getTimeBounds( )
    {
        if( observationTimes_.size( ) == 0 )
        {
            throw std::runtime_error(
                    "Error when getting time bounds of observation set, no observations found" );
        }
        return timeBounds_;
    }

    std::map< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
    getObservationsHistory( )
    {
        return utilities::createMapFromVectors<
                TimeType,
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( observationTimes_,
                                                                             observations_ );
    }

    std::vector< Eigen::VectorXd > getObservationsDependentVariables( )
    {
        return observationsDependentVariables_;
    }

    Eigen::MatrixXd getObservationsDependentVariablesMatrix( )
    {
        Eigen::MatrixXd dependentVariablesMatrix = Eigen::MatrixXd::Zero(
                numberOfObservations_,
                dependentVariableCalculator_->getTotalDependentVariableSize( ) );
        for( unsigned int i = 0; i < observationsDependentVariables_.size( ); i++ )
        {
            dependentVariablesMatrix.block(
                    i, 0, 1, dependentVariableCalculator_->getTotalDependentVariableSize( ) ) =
                    observationsDependentVariables_[ i ].transpose( );
        }
        return dependentVariablesMatrix;
    }

    //! Function returning the dependent variable values for a single observation (indicated by index)
    Eigen::VectorXd getDependentVariablesForSingleObservation( unsigned int index ) const
    {
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when retrieving observation dependent variables for single observation, "
                    "required index incompatible with number of observations." );
        }
        return observationsDependentVariables_.at( index );
    }

    //! Function returning the values of a single dependent variable (specified by dependent variable settings)
    Eigen::MatrixXd getSingleDependentVariable(
            std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings,
            const bool returnFirstCompatibleSettings = false )
    {
        // Retrieve full map of dependent variables start indices and sizes based on settings
        std::map< std::pair< int, int >,
                  std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > >
                settingsIndicesAndSizes =
                        dependentVariableCalculator_->getSettingsIndicesAndSizes( );

        // Get the start indices and sizes of all dependent variables that would be compatible with
        // the settings provided as inputs.
        std::vector< std::pair< int, int > > indicesAndSizes;
        for( auto it: settingsIndicesAndSizes )
        {
            if( dependentVariableSettings->areSettingsCompatible( it.second ) )
            {
                indicesAndSizes.push_back( it.first );
            }
        }

        // Check that a single settings is identified
        if( indicesAndSizes.size( ) == 0 )
        {
            throw std::runtime_error(
                    "Error when getting dependent variable, no dependent variable values found for "
                    "given settings." );
        }
        else if( indicesAndSizes.size( ) > 1 && !returnFirstCompatibleSettings )
        {
            throw std::runtime_error(
                    "Error when getting dependent variable, multiple dependent variables found for "
                    "given settings." );
        }

        // Return the dependent variable values for the first (and/or single) compatible settings
        // identified
        return getSingleDependentVariable( indicesAndSizes.at( 0 ) );
    }

    //! Function returning the list of all dependent variable settings compatible with the settings provided as inputs
    //! (which might not be fully defined, i.e. with missing link ends information, etc.)
    std::vector< std::shared_ptr< ObservationDependentVariableSettings > >
    getCompatibleDependentVariablesSettingsList(
            std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings )
            const
    {
        // Retrieve all dependent variables settings for this observation set
        std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > >
                allDependentVariablesSettings =
                        dependentVariableCalculator_->getDependentVariableSettings( );

        // Check which settings are compatible with the input settings object
        std::vector< std::shared_ptr< ObservationDependentVariableSettings > > compatibleSettings;
        for( auto it: allDependentVariablesSettings )
        {
            if( dependentVariableSettings->areSettingsCompatible( it ) )
            {
                compatibleSettings.push_back( it );
            }
        }
        return compatibleSettings;
    }

    //! Function returning a vector containing the values of all dependent variables compatible with the settings provided as input
    //! The order in which they are provided matches the list of compatible settings given by the getCompatibleDependentVariablesSettingsList function
    std::vector< Eigen::MatrixXd > getAllCompatibleDependentVariables(
            std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings )
            const
    {
        // Retrieve start indices and sizes for each dependent variable settings
        std::map< std::pair< int, int >,
                  std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > >
                settingsIndicesAndSizes =
                        dependentVariableCalculator_->getSettingsIndicesAndSizes( );

        // Retrieve all relevant dependent variables
        std::vector< Eigen::MatrixXd > dependentVariablesList;
        for( auto it: settingsIndicesAndSizes )
        {
            if( dependentVariableSettings->areSettingsCompatible( it.second ) )
            {
                dependentVariablesList.push_back( getSingleDependentVariable( it.first ) );
            }
        }

        return dependentVariablesList;
    }

    std::vector< Eigen::VectorXd >& getObservationsDependentVariablesReference( )
    {
        return observationsDependentVariables_;
    }

    //! Function to reset the observation dependent variable values
    void setObservationsDependentVariables( std::vector< Eigen::VectorXd >& dependentVariables )
    {
        if( observationsDependentVariables_.size( ) > 0 )
        {
            if( observationsDependentVariables_.size( ) != numberOfObservations_ )
            {
                throw std::runtime_error(
                        "Error when resetting observation dependent variables in "
                        "SingleObservationSet, the input size should be consistent "
                        "with the number of observations." );
            }
            if( ( dependentVariableCalculator_ != nullptr ) &&
                ( dependentVariables[ 0 ].size( ) !=
                  dependentVariableCalculator_->getTotalDependentVariableSize( ) ) )
            {
                throw std::runtime_error(
                        "Error when resetting observation dependent variables in "
                        "SingleObservationSet, the size of the observation dependent variables "
                        "input "
                        "should be consistent with the total dependent variable size." );
            }
        }
        observationsDependentVariables_ = dependentVariables;
    }

    std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator >
    getDependentVariableCalculator( )
    {
        return dependentVariableCalculator_;
    }

    //! Function that returns the time history of all observation dependent variables. It must be noted that the reported epochs are the times at which the
    //! observations are computed/acquired, which might differ from the times at which the dependent variables are evaluated.
    std::map< TimeType, Eigen::VectorXd > getDependentVariableHistory( )
    {
        return utilities::createMapFromVectors< TimeType, Eigen::VectorXd >(
                observationTimes_, observationsDependentVariables_ );
    }

    //! Function that returns the time history of a single dependent variables (specified by settings). It must be noted that the reported epochs are the times
    //! at which the observations are computed/acquired, which might differ from the times at which the dependent variables are evaluated.
    std::map< TimeType, Eigen::VectorXd > getSingleDependentVariableHistory(
            std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings,
            const bool returnFirstCompatibleSettings = false )
    {
        Eigen::MatrixXd singleDependentVariableValues = getSingleDependentVariable(
                dependentVariableSettings, returnFirstCompatibleSettings );
        std::map< TimeType, Eigen::VectorXd > singleDependentVariableMap;
        for( unsigned int i = 0; i < numberOfObservations_; i++ )
        {
            Eigen::VectorXd dependentVariableCurrentTime =
                    singleDependentVariableValues
                            .block( i, 0, 1, singleDependentVariableValues.cols( ) )
                            .transpose( );
            singleDependentVariableMap[ observationTimes_[ i ] ] = dependentVariableCurrentTime;
        }
        return singleDependentVariableMap;
    }

    std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >
    getAncilliarySettings( )
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
                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero(
                        singleObservationSize_ * numberOfObservations_, 1 );
        for( unsigned int i = 0; i < numberOfObservations_; i++ )
        {
            weightsVector.block( i * singleObservationSize_, 0, singleObservationSize_, 1 ) =
                    weights_.at( i );
        }
        return weightsVector;
    }

    Eigen::Matrix< double, Eigen::Dynamic, 1 > getWeight( unsigned int index ) const
    {
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when retrieving single observation weight, required index incompatible "
                    "with number of observations." );
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
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero(
                        singleObservationSize_ * numberOfObservations_, 1 );
        for( unsigned int i = 0; i < numberOfObservations_; i++ )
        {
            residualsVector.block( i * singleObservationSize_, 0, singleObservationSize_, 1 ) =
                    residuals_.at( i );
        }
        return residualsVector;
    }

    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
    getResidualsReference( )
    {
        return residuals_;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidual(
            unsigned int index ) const
    {
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when retrieving single observation residual, required index "
                    "incompatible with number of observations." );
        }
        return residuals_.at( index );
    }

    Eigen::VectorXd getRmsResiduals( )
    {
        Eigen::VectorXd rmsResiduals = Eigen::VectorXd::Zero( singleObservationSize_ );
        for( unsigned int i = 0; i < singleObservationSize_; i++ )
        {
            // Calculate RMS of the residuals for each observation component
            for( int j = 0; j < numberOfObservations_; j++ )
            {
                rmsResiduals[ i ] += residuals_[ j ]( i, 0 ) * residuals_[ j ]( i, 0 );
            }
            rmsResiduals[ i ] = std::sqrt( rmsResiduals[ i ] / numberOfObservations_ );
        }

        return rmsResiduals;
    }

    Eigen::VectorXd getMeanResiduals( )
    {
        Eigen::VectorXd meanResiduals = Eigen::VectorXd::Zero( singleObservationSize_ );
        for( unsigned int i = 0; i < singleObservationSize_; i++ )
        {
            // Calculate mean residual for each observation component
            for( unsigned int j = 0; j < numberOfObservations_; j++ )
            {
                meanResiduals[ i ] += residuals_[ j ]( i, 0 );
            }
            meanResiduals[ i ] /= numberOfObservations_;
        }
        return meanResiduals;
    }

    void setConstantWeight( const double weight )
    {
        for( unsigned int k = 0; k < numberOfObservations_; k++ )
        {
            weights_.at( k ) = weight *
                    Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( singleObservationSize_, 1 );
        }
    }

    void setConstantWeight( const Eigen::Matrix< double, Eigen::Dynamic, 1 >& weight )
    {
        if( weight.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error(
                    "Error when setting constant weight in single observation set, weight size is "
                    "inconsistent with single observation size." );
        }
        for( unsigned int k = 0; k < weights_.size( ); k++ )
        {
            weights_.at( k ) = weight;
        }
    }

    void setTabulatedWeights( const Eigen::VectorXd& weightsVector )
    {
        if( weightsVector.rows( ) !=
            static_cast< int >( singleObservationSize_ * observations_.size( ) ) )
        {
            throw std::runtime_error(
                    "Error when setting weights in single observation set, sizes are "
                    "incompatible." );
        }
        for( unsigned int k = 0; k < numberOfObservations_; k++ )
        {
            for( unsigned int i = 0; i < singleObservationSize_; i++ )
            {
                weights_.at( k )[ i ] = weightsVector[ k * singleObservationSize_ + i ];
            }
        }
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
    getComputedObservations( ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                computedObservations;
        for( unsigned int k = 0; k < observations_.size( ); k++ )
        {
            computedObservations.push_back( observations_.at( k ) - residuals_.at( k ) );
        }
        return computedObservations;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getComputedObservationsVector( ) const
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > computedObservationsVector =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero(
                        singleObservationSize_ * numberOfObservations_, 1 );

        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                computedObservations = getComputedObservations( );
        for( unsigned int i = 0; i < numberOfObservations_; i++ )
        {
            computedObservationsVector.segment( i * singleObservationSize_,
                                                singleObservationSize_ ) =
                    computedObservations.at( i );
        }
        return computedObservationsVector;
    }

    unsigned int getNumberOfFilteredObservations( ) const
    {
        int numberFilteredObservations = 0;
        if( filteredObservationSet_ != nullptr )
        {
            numberFilteredObservations = filteredObservationSet_->getNumberOfObservables( );
        }
        return numberFilteredObservations;
    }

    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >
    getFilteredObservationSet( ) const
    {
        return filteredObservationSet_;
    }

    void removeSingleObservation( unsigned int indexToRemove )
    {
        if( indexToRemove >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when removing single observation from SingleObservationSet, index "
                    "incompatible with number of observations." );
        }

        // Update observations
        observations_.erase( observations_.begin( ) + indexToRemove );
        observationTimes_.erase( observationTimes_.begin( ) + indexToRemove );
        residuals_.erase( residuals_.begin( ) + indexToRemove );
        weights_.erase( weights_.begin( ) + indexToRemove );

        if( observationsDependentVariables_.size( ) > 0 )
        {
            observationsDependentVariables_.erase( observationsDependentVariables_.begin( ) +
                                                   indexToRemove );
        }

        // Update number of observations and time bounds
        numberOfObservations_ = observations_.size( );
        updateTimeBounds( );
    }

    void removeObservations( const std::vector< unsigned int >& indicesToRemove )
    {
        unsigned int counter = 0;
        for( auto ind: indicesToRemove )
        {
            removeSingleObservation( ind -
                                     counter );  // observations are already filtered and sorted
            counter += 1;
        }
    }

    void filterObservations( const std::shared_ptr< ObservationFilterBase > observationFilter,
                             const bool saveFilteredObservations = true )
    {
        if( observationFilter->filterOut( ) && filteredObservationSet_ == nullptr )
        {
            // Initialise empty filtered observation set
            filteredObservationSet_ = std::make_shared<
                    SingleObservationSet< ObservationScalarType, TimeType > >(
                    observableType_,
                    linkEnds_,
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                    std::vector< TimeType >( ),
                    referenceLinkEnd_,
                    std::vector< Eigen::VectorXd >( ),
                    dependentVariableCalculator_,
                    ancilliarySettings_ );
        }
        if( !observationFilter->filterOut( ) && filteredObservationSet_ == nullptr )
        {
            throw std::runtime_error(
                    "Error when attempting to un-filter observations, filtered observation set is "
                    "empty." );
        }

        unsigned int nbObservationsToTest =
                ( observationFilter->filterOut( )
                          ? numberOfObservations_
                          : filteredObservationSet_->getNumberOfObservables( ) );
        bool useOppositeCondition = observationFilter->useOppositeCondition( );

        std::vector< unsigned int > indicesToRemove;
        switch( observationFilter->getFilterType( ) )
        {
            case residual_filtering: {
                Eigen::VectorXd residualCutOff = Eigen::VectorXd::Zero( singleObservationSize_ );
                if( std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter ) !=
                    nullptr )
                {
                    residualCutOff = std::dynamic_pointer_cast< ObservationFilter< double > >(
                                             observationFilter )
                                             ->getFilterValue( ) *
                            Eigen::VectorXd::Ones( singleObservationSize_ );
                }
                else if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >(
                                 observationFilter ) != nullptr )
                {
                    if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >(
                                observationFilter )
                                ->getFilterValue( )
                                .size( ) != singleObservationSize_ )
                    {
                        throw std::runtime_error(
                                "Error when performing residual filtering, size of the residual "
                                "cut off vector inconsistent with observable size." );
                    }
                    residualCutOff =
                            std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >(
                                    observationFilter )
                                    ->getFilterValue( );
                }

                for( unsigned int j = 0; j < nbObservationsToTest; j++ )
                {
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >
                            singleObservationResidual =
                                    ( observationFilter->filterOut( )
                                              ? residuals_.at( j )
                                              : filteredObservationSet_->getResidual( j ) );
                    bool removeObservation = false;
                    for( unsigned int k = 0; k < singleObservationSize_; k++ )
                    {
                        if( ( !useOppositeCondition &&
                              ( std::fabs( singleObservationResidual[ k ] ) >
                                residualCutOff[ k ] ) ) ||
                            ( useOppositeCondition &&
                              ( std::fabs( singleObservationResidual[ k ] ) <=
                                residualCutOff[ k ] ) ) )
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
                Eigen::VectorXd absoluteValueCutOff =
                        Eigen::VectorXd::Zero( singleObservationSize_ );
                if( std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter ) !=
                    nullptr )
                {
                    absoluteValueCutOff = std::dynamic_pointer_cast< ObservationFilter< double > >(
                                                  observationFilter )
                                                  ->getFilterValue( ) *
                            Eigen::VectorXd::Ones( singleObservationSize_ );
                }
                else if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >(
                                 observationFilter ) != nullptr )
                {
                    if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >(
                                observationFilter )
                                ->getFilterValue( )
                                .size( ) != singleObservationSize_ )
                    {
                        throw std::runtime_error(
                                "Error when performing observation value filtering, size of the "
                                "filter value inconsistent with observable size." );
                    }
                    absoluteValueCutOff =
                            std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >(
                                    observationFilter )
                                    ->getFilterValue( );
                }

                for( unsigned int j = 0; j < nbObservationsToTest; j++ )
                {
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservation =
                            ( observationFilter->filterOut( )
                                      ? observations_.at( j )
                                      : filteredObservationSet_->getObservation( j ) );
                    bool removeObservation = false;
                    for( unsigned int k = 0; k < singleObservationSize_; k++ )
                    {
                        if( ( !useOppositeCondition &&
                              ( singleObservation[ k ] > absoluteValueCutOff[ k ] ) ) ||
                            ( useOppositeCondition &&
                              ( singleObservation[ k ] <= absoluteValueCutOff[ k ] ) ) )
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
                std::vector< double > filterEpochs =
                        std::dynamic_pointer_cast< ObservationFilter< std::vector< double > > >(
                                observationFilter )
                                ->getFilterValue( );
                for( unsigned int j = 0; j < nbObservationsToTest; j++ )
                {
                    TimeType singleObservationTime =
                            ( observationFilter->filterOut( )
                                      ? observationTimes_.at( j )
                                      : filteredObservationSet_->getObservationTime( j ) );
                    if( ( !useOppositeCondition &&
                          ( std::count( filterEpochs.begin( ),
                                        filterEpochs.end( ),
                                        singleObservationTime ) > 0 ) ) ||
                        ( useOppositeCondition &&
                          ( std::count( filterEpochs.begin( ),
                                        filterEpochs.end( ),
                                        singleObservationTime ) == 0 ) ) )
                    {
                        indicesToRemove.push_back( j );
                    }
                }
                break;
            }
            case time_bounds_filtering: {
                std::pair< double, double > timeBounds =
                        std::dynamic_pointer_cast<
                                ObservationFilter< std::pair< double, double > > >(
                                observationFilter )
                                ->getFilterValue( );
                for( unsigned int j = 0; j < nbObservationsToTest; j++ )
                {
                    TimeType singleObservationTime =
                            ( observationFilter->filterOut( )
                                      ? observationTimes_.at( j )
                                      : filteredObservationSet_->getObservationTime( j ) );
                    if( ( !useOppositeCondition &&
                          ( ( singleObservationTime >= timeBounds.first ) &&
                            ( singleObservationTime <= timeBounds.second ) ) ) ||
                        ( useOppositeCondition &&
                          ( ( singleObservationTime < timeBounds.first ) ||
                            ( singleObservationTime > timeBounds.second ) ) ) )
                    {
                        indicesToRemove.push_back( j );
                    }
                }
                break;
            }
            case dependent_variable_filtering: {
                if( std::dynamic_pointer_cast< ObservationDependentVariableFilter >(
                            observationFilter ) == nullptr )
                {
                    throw std::runtime_error(
                            "Error when performing dependent variable filtering, inconsistent "
                            "filter input (should be ObservationDependentVariableFilter object)." );
                }

                // Retrieve dependent variable settings and size
                std::shared_ptr< ObservationDependentVariableSettings > settings =
                        std::dynamic_pointer_cast< ObservationDependentVariableFilter >(
                                observationFilter )
                                ->getDependentVariableSettings( );
                unsigned int dependentVariableSize =
                        getObservationDependentVariableSize( settings, linkEnds_.linkEnds_ );

                // Retrieve dependent variable cut-off value
                Eigen::VectorXd dependentVariableCutOff =
                        std::dynamic_pointer_cast< ObservationDependentVariableFilter >(
                                observationFilter )
                                ->getFilterValue( );
                if( dependentVariableCutOff.size( ) != dependentVariableSize )
                {
                    throw std::runtime_error(
                            "Error when performing dependent variable filtering, size of the "
                            "dependent variable cut off vector inconsistent with dependent "
                            "variable size." );
                }

                // Retrieve dependent variable values
                Eigen::MatrixXd singleDependentVariableValues =
                        ( observationFilter->filterOut( )
                                  ? getSingleDependentVariable( settings )
                                  : filteredObservationSet_->getSingleDependentVariable(
                                            settings ) );
                if( ( singleDependentVariableValues.rows( ) != nbObservationsToTest ) ||
                    ( singleDependentVariableValues.cols( ) != dependentVariableSize ) )
                {
                    throw std::runtime_error(
                            "Error when performing dependent variable filtering, size of "
                            "observation dependent variables is inconsistent with the "
                            "number of observations and presupposed dependent variable size." );
                }

                // Check dependent variable values against cut-off value
                for( unsigned int j = 0; j < nbObservationsToTest; j++ )
                {
                    bool removeObservation = false;
                    for( unsigned int k = 0; k < dependentVariableSize; k++ )
                    {
                        if( ( !useOppositeCondition &&
                              ( singleDependentVariableValues( j, k ) ) >
                                      dependentVariableCutOff[ k ] ) ||
                            ( useOppositeCondition &&
                              ( singleDependentVariableValues( j, k ) ) <=
                                      dependentVariableCutOff[ k ] ) )
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

        if( observationFilter->filterOut( ) )
        {
            moveObservationsInOutFilteredSet( indicesToRemove, true, saveFilteredObservations );
        }
        else
        {
            moveObservationsInOutFilteredSet( indicesToRemove, false, true );
        }
    }

    void addObservations(
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
                    observations,
            const std::vector< TimeType >& times,
            const std::vector< Eigen::VectorXd >& dependentVariables = { },
            const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights = { },
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
                    residuals = { },
            const bool sortObservations = true )
    {
        if( ( observations.size( ) != times.size( ) ) ||
            ( weights.size( ) > 0 && ( observations.size( ) != weights.size( ) ) ) ||
            ( residuals.size( ) > 0 && ( observations.size( ) != residuals.size( ) ) ) ||
            ( dependentVariables.size( ) > 0 &&
              ( observations.size( ) != dependentVariables.size( ) ) ) )
        {
            throw std::runtime_error(
                    "Error when adding observations to SingleObservationSet, input sizes are "
                    "inconsistent." );
        }

        for( unsigned int k = 0; k < observations.size( ); k++ )
        {
            if( observations.at( k ).size( ) != singleObservationSize_ )
            {
                throw std::runtime_error(
                        "Error when adding observations to SingleObservationSet, new observation "
                        "size is inconsistent." );
            }

            observations_.push_back( observations.at( k ) );
            observationTimes_.push_back( times.at( k ) );

            // If residuals are provided as inputs
            if( residuals.size( ) > 0 )
            {
                if( residuals.at( k ).size( ) != singleObservationSize_ )
                {
                    throw std::runtime_error(
                            "Error when adding observations to SingleObservationSet, new residual "
                            "size is inconsistent." );
                }
                residuals_.push_back( residuals.at( k ) );
            }
            else  // Otherwise, set to zero by default
            {
                residuals_.push_back(
                        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero(
                                singleObservationSize_, 1 ) );
            }

            // If weights are provided as inputs
            if( weights.size( ) > 0 )
            {
                if( weights.at( k ).size( ) != singleObservationSize_ )
                {
                    throw std::runtime_error(
                            "Error when adding observations to SingleObservationSet, new weight "
                            "size is inconsistent." );
                }
                weights_.push_back( weights.at( k ) );
            }
            else  // Otherwise, set to one by default
            {
                weights_.push_back( Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones(
                        singleObservationSize_, 1 ) );
            }

            // if dependent variables are set
            if( ( observationsDependentVariables_.size( ) > 0 || numberOfObservations_ == 0 ) &&
                dependentVariables.size( ) > 0 )
            {
                observationsDependentVariables_.push_back( dependentVariables.at( k ) );
            }

            numberOfObservations_ += 1;
        }

        // Sort observations
        if( sortObservations )
        {
            orderObservationsAndMetadata( );
        }

        // Update time bounds
        updateTimeBounds( );
    }

    void addDependentVariables(
            const std::vector<
                    std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > >
                    dependentVariableSettings,
            const simulation_setup::SystemOfBodies& bodies )
    {
        if( dependentVariableCalculator_ == nullptr )
        {
            dependentVariableCalculator_ =
                    std::make_shared< ObservationDependentVariableCalculator >(
                            observableType_, linkEnds_.linkEnds_ );
        }
        dependentVariableCalculator_->addDependentVariables( dependentVariableSettings, bodies );
    }

   private:
    void orderObservationsAndMetadata( )
    {
        if( !std::is_sorted( observationTimes_.begin( ), observationTimes_.end( ) ) )
        {
            if( observations_.size( ) != numberOfObservations_ )
            {
                throw std::runtime_error(
                        "Error when making SingleObservationSet, number of observations is "
                        "incompatible after time ordering" );
            }

            std::multimap< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                    observationsMap;
            for( unsigned int i = 0; i < observations_.size( ); i++ )
            {
                observationsMap.insert( { observationTimes_.at( i ), observations_.at( i ) } );
            }
            observationTimes_ = utilities::createVectorFromMultiMapKeys( observationsMap );
            observations_ = utilities::createVectorFromMultiMapValues( observationsMap );

            if( observationsDependentVariables_.size( ) > 0 )
            {
                if( observationsDependentVariables_.size( ) != numberOfObservations_ )
                {
                    throw std::runtime_error(
                            "Error when making SingleObservationSet, dependent variables vector "
                            "size is incompatible after time ordering" );
                }
                std::multimap< TimeType, Eigen::VectorXd > observationsDependentVariablesMap;
                for( unsigned int i = 0; i < observationsDependentVariables_.size( ); i++ )
                {
                    observationsDependentVariablesMap.insert(
                            { observationTimes_.at( i ),
                              observationsDependentVariables_.at( i ) } );
                }
                observationsDependentVariables_ = utilities::createVectorFromMultiMapValues(
                        observationsDependentVariablesMap );
            }

            if( weights_.size( ) != numberOfObservations_ )
            {
                throw std::runtime_error(
                        "Error when making SingleObservationSet, weights size is incompatible "
                        "after time ordering" );
            }
            std::multimap< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > weightsMap;
            for( unsigned int i = 0; i < weights_.size( ); i++ )
            {
                weightsMap.insert( { observationTimes_.at( i ), weights_.at( i ) } );
            }
            weights_ = utilities::createVectorFromMultiMapValues( weightsMap );

            if( residuals_.size( ) != numberOfObservations_ )
            {
                throw std::runtime_error(
                        "Error when making SingleObservationSet, residuals size is incompatible "
                        "after time ordering" );
            }
            std::multimap< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                    residualsMap;
            for( unsigned int i = 0; i < residuals_.size( ); i++ )
            {
                residualsMap.insert( { observationTimes_.at( i ), residuals_.at( i ) } );
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
            timeBounds_ = std::make_pair(
                    *std::min_element( observationTimes_.begin( ), observationTimes_.end( ) ),
                    *std::max_element( observationTimes_.begin( ), observationTimes_.end( ) ) );
        }
    }

    void moveObservationsInOutFilteredSet( const std::vector< unsigned int >& indices,
                                           const bool moveInFilteredSet = true,
                                           const bool saveFilteredObservations = true )
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
        std::vector< TimeType > times;
        std::vector< Eigen::VectorXd > dependentVariables;
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights;
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;

        if( moveInFilteredSet )
        {
            for( auto index: indices )
            {
                if( index >= numberOfObservations_ )
                {
                    throw std::runtime_error(
                            "Error when moving observation to filtered observation set, index "
                            "incompatible with number of observations." );
                }

                observations.push_back( observations_.at( index ) );
                times.push_back( observationTimes_.at( index ) );
                weights.push_back( weights_.at( index ) );
                residuals.push_back( residuals_.at( index ) );

                // If dependent variables not empty
                if( observationsDependentVariables_.size( ) > 0 )
                {
                    dependentVariables.push_back( observationsDependentVariables_.at( index ) );
                }
            }

            if( saveFilteredObservations )
            {
                filteredObservationSet_->addObservations(
                        observations, times, dependentVariables, weights, residuals, true );
            }
            removeObservations( indices );
        }
        else
        {
            for( auto index: indices )
            {
                if( getNumberOfFilteredObservations( ) == 0 )
                {
                    throw std::runtime_error(
                            "Error when moving observation back from filtered observation set, "
                            "filtered observation set is empty." );
                }
                if( index >= getNumberOfFilteredObservations( ) )
                {
                    throw std::runtime_error(
                            "Error when moving observation back from filtered observation set, "
                            "index incompatible with number of observations." );
                }

                observations.push_back( filteredObservationSet_->getObservations( ).at( index ) );
                times.push_back( filteredObservationSet_->getObservationTimes( ).at( index ) );
                weights.push_back( filteredObservationSet_->getWeights( ).at( index ) );
                residuals.push_back( filteredObservationSet_->getResiduals( ).at( index ) );

                // If dependent variables are set
                if( filteredObservationSet_->getObservationsDependentVariables( ).size( ) > 0 )
                {
                    dependentVariables.push_back(
                            filteredObservationSet_->getObservationsDependentVariables( ).at(
                                    index ) );
                }
            }
            addObservations( observations, times, dependentVariables, weights, residuals, true );
            filteredObservationSet_->removeObservations( indices );
        }
    }

    //! Function extracting the values of a single dependent variable
    Eigen::MatrixXd getSingleDependentVariable(
            std::pair< int, int > dependentVariableIndexAndSize ) const
    {
        Eigen::MatrixXd singleDependentVariable = Eigen::MatrixXd::Zero(
                numberOfObservations_, dependentVariableIndexAndSize.second );
        for( unsigned int i = 0; i < observationsDependentVariables_.size( ); i++ )
        {
            if( dependentVariableIndexAndSize.first + dependentVariableIndexAndSize.second >
                observationsDependentVariables_.at( i ).size( ) )
            {
                throw std::runtime_error(
                        "Error when retrieving single observation dependent variable, required "
                        "index and size incompatible with "
                        "dependent variables size." );
            }
            else
            {
                Eigen::VectorXd singleDependentVariableVector =
                        observationsDependentVariables_.at( i ).segment(
                                dependentVariableIndexAndSize.first,
                                dependentVariableIndexAndSize.second );
                singleDependentVariable.block( i, 0, 1, dependentVariableIndexAndSize.second ) =
                        singleDependentVariableVector.transpose( );
            }
        }
        return singleDependentVariable;
    }

    const ObservableType observableType_;

    LinkDefinition linkEnds_;

    std::pair< TimeType, TimeType > timeBounds_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations_;

    std::vector< TimeType > observationTimes_;

    const LinkEndType referenceLinkEnd_;

    std::vector< Eigen::VectorXd > observationsDependentVariables_;

    std::shared_ptr< simulation_setup::ObservationDependentVariableCalculator >
            dependentVariableCalculator_;

    const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >
            ancilliarySettings_;

    unsigned int numberOfObservations_;

    unsigned int singleObservationSize_;

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights_;

    //    Eigen::VectorXd weightsVector_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals_;

    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >
            filteredObservationSet_;
};

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if<
                  is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value,
                  int >::type = 0 >
std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > filterObservations(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >
                singleObservationSet,
        const std::shared_ptr< ObservationFilterBase > observationFilter,
        const bool saveFilteredObservations = false )
{
    if( !observationFilter->filterOut( ) )
    {
        throw std::runtime_error(
                "Error when creating new single observation set post-filtering, the filterOut "
                "option should be set to true" );
    }
    // Create new observation set
    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > newObservationSet =
            std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                    singleObservationSet->getObservableType( ),
                    singleObservationSet->getLinkEnds( ),
                    singleObservationSet->getObservationsReference( ),
                    singleObservationSet->getObservationTimesReference( ),
                    singleObservationSet->getReferenceLinkEnd( ),
                    singleObservationSet->getObservationsDependentVariablesReference( ),
                    singleObservationSet->getDependentVariableCalculator( ),
                    singleObservationSet->getAncilliarySettings( ) );
    newObservationSet->setTabulatedWeights( singleObservationSet->getWeightsVector( ) );
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals =
            singleObservationSet->getResidualsReference( );
    newObservationSet->setResiduals( singleObservationSet->getResidualsReference( ) /*residuals*/ );

    // Filter observations from new observation set
    newObservationSet->filterObservations( observationFilter, saveFilteredObservations );

    return newObservationSet;
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if<
                  is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value,
                  int >::type = 0 >
std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
splitObservationSet(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >
                observationSet,
        const std::shared_ptr< ObservationSetSplitterBase > observationSetSplitter,
        const bool printWarning = true )
{
    if( printWarning && observationSet->getFilteredObservationSet( ) != nullptr )
    {
        std::cerr << "Warning when splitting single observation set, the filtered observation set "
                     "pointer is not empty and "
                     " any filtered observation will be lost after splitting."
                  << std::endl;
    }

    std::vector< int > rawStartIndicesNewSets = { 0 };
    std::vector< TimeType > observationTimes = observationSet->getObservationTimes( );
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations =
            observationSet->getObservations( );
    std::vector< Eigen::VectorXd > dependentVariables =
            observationSet->getObservationsDependentVariables( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > weightsVector = observationSet->getWeightsVector( );
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals =
            observationSet->getResiduals( );

    switch( observationSetSplitter->getSplitterType( ) )
    {
        case time_tags_splitter: {
            std::vector< double > timeTags =
                    std::dynamic_pointer_cast< ObservationSetSplitter< std::vector< double > > >(
                            observationSetSplitter )
                            ->getSplitterValue( );
            for( auto currentTimeTag: timeTags )
            {
                if( currentTimeTag > observationSet->getTimeBounds( ).first )
                {
                    bool detectedStartSet = false;
                    int indexObs = rawStartIndicesNewSets.at( rawStartIndicesNewSets.size( ) - 1 );
                    while( !detectedStartSet &&
                           indexObs < static_cast< int >( observationTimes.size( ) ) )
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
            double maxTimeInterval = std::dynamic_pointer_cast< ObservationSetSplitter< double > >(
                                             observationSetSplitter )
                                             ->getSplitterValue( );
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
            double maxTimeSpan = std::dynamic_pointer_cast< ObservationSetSplitter< double > >(
                                         observationSetSplitter )
                                         ->getSplitterValue( );
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
            int maxNbObs = std::dynamic_pointer_cast< ObservationSetSplitter< int > >(
                                   observationSetSplitter )
                                   ->getSplitterValue( );
            if( maxNbObs < observationSetSplitter->getMinNumberObservations( ) )
            {
                throw std::runtime_error(
                        "Error when splitting observation sets, the maximum number of observations "
                        "cannot be smaller than the minimum number of observations." );
            }
            for( int ind = maxNbObs;
                 ind < static_cast< int >( observationSet->getNumberOfObservables( ) );
                 ind += maxNbObs )
            {
                rawStartIndicesNewSets.push_back( ind );
            }
            rawStartIndicesNewSets.push_back( static_cast< int >( observationTimes.size( ) ) );
            break;
        }
        default:
            throw std::runtime_error( "Observation set splitter type not recognised." );
    }

    // Check that the minimum number of observations is met
    std::vector< std::pair< int, int > > indicesNewSets;
    for( unsigned int j = 1; j < rawStartIndicesNewSets.size( ); j++ )
    {
        if( ( rawStartIndicesNewSets.at( j ) - rawStartIndicesNewSets.at( j - 1 ) ) >=
            observationSetSplitter->getMinNumberObservations( ) )
        {
            indicesNewSets.push_back( std::make_pair(
                    rawStartIndicesNewSets.at( j - 1 ),
                    rawStartIndicesNewSets.at( j ) - rawStartIndicesNewSets.at( j - 1 ) ) );
        }
    }

    // Split current observation set based on indices
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > >
            newObsSets;
    for( unsigned int k = 0; k < indicesNewSets.size( ); k++ )
    {
        int startIndex = indicesNewSets.at( k ).first;
        int sizeCurrentSet = indicesNewSets.at( k ).second;

        std::vector< Eigen::VectorXd > newDependentVariables;
        if( !dependentVariables.empty( ) )
        {
            newDependentVariables = utilities::getStlVectorSegment(
                    observationSet->getObservationsDependentVariablesReference( ),
                    startIndex,
                    sizeCurrentSet );
        }

        std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > newSet =
                std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                        observationSet->getObservableType( ),
                        observationSet->getLinkEnds( ),
                        utilities::getStlVectorSegment( observationSet->getObservationsReference( ),
                                                        startIndex,
                                                        sizeCurrentSet ),
                        utilities::getStlVectorSegment(
                                observationSet->getObservationTimesReference( ),
                                startIndex,
                                sizeCurrentSet ),
                        observationSet->getReferenceLinkEnd( ),
                        newDependentVariables,
                        observationSet->getDependentVariableCalculator( ),
                        observationSet->getAncilliarySettings( ) );

        Eigen::Matrix< double, Eigen::Dynamic, 1 > newWeightsVector = weightsVector.segment(
                startIndex, sizeCurrentSet * observationSet->getSingleObservableSize( ) );
        newSet->setTabulatedWeights( newWeightsVector );

        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > newResiduals =
                utilities::getStlVectorSegment(
                        observationSet->getResidualsReference( ), startIndex, sizeCurrentSet );
        newSet->setResiduals( newResiduals );

        newObsSets.push_back( newSet );
    }

    return newObsSets;
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_SINGLE_OBSERVATION_SET_H
