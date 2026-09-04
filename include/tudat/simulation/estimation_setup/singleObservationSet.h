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
#include <set>
#include <vector>

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/simulation/estimation_setup/observationsProcessing.h"
#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{

namespace observation_models
{

using namespace simulation_setup;

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class SingleObservationSet
{
public:
    SingleObservationSet( const ObservableType observableType,
                          const LinkDefinition& linkEnds,
                          const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                          const std::vector< TimeType > observationTimes,
                          const LinkEndType referenceLinkEnd,
                          const std::vector< Eigen::VectorXd >& observationsDependentVariables = std::vector< Eigen::VectorXd >( ),
                          const std::shared_ptr< ObservationDependentVariableBookkeeping > dependentVariableBookkeeping = nullptr,
                          const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr,
                          const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights = {},
                          const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals = {},
                          const bool eraseDuplicates = false ):

        observableType_( observableType ), referenceLinkEnd_( referenceLinkEnd ),
        dependentVariableBookkeeping_( dependentVariableBookkeeping ), ancillarySettings_( ancillarySettings ),
        filteredObservationSet_( nullptr ), dataset_( std::make_shared< ObservationDataset< ObservationScalarType, TimeType > >( ) ),
        setId_( 0 )
    {
        if( dependentVariableBookkeeping_ != nullptr )
        {
            if( dependentVariableBookkeeping_->getObservableType( ) != observableType_ )
            {
                std::cout << dependentVariableBookkeeping_->getObservableType( ) << " " << observableType_ << std::endl;
                throw std::runtime_error(
                        "Error when creating SingleObservationSet, "
                        "ObservationDependentVariableBookkeeping has incompatible type " );
            }

            if( !( dependentVariableBookkeeping_->getLinkEnds( ) == linkEnds ) )
            {
                throw std::runtime_error(
                        "Error when creating SingleObservationSet, "
                        "ObservationDependentVariableBookkeeping has incompatible link ends " );
            }
        }

        if( observations.size( ) != observationTimes.size( ) )
        {
            throw std::runtime_error( "Error when making SingleObservationSet, input sizes are inconsistent." +
                                      std::to_string( observations.size( ) ) + ", " + std::to_string( observationTimes.size( ) ) );
        }

        // Check observation dependent variables size
        if( observationsDependentVariables.size( ) > 0 )
        {
            if( observationsDependentVariables.size( ) != observations.size( ) )
            {
                throw std::runtime_error(
                        "Error when creating SingleObservationSet, the size of the observation "
                        "dependent variables input should be consistent "
                        "with the number of observations." );
            }
            if( dependentVariableBookkeeping_ != nullptr &&
                observationsDependentVariables[ 0 ].size( ) != dependentVariableBookkeeping_->getTotalDependentVariableSize( ) )
            {
                throw std::runtime_error(
                        "Error when creating SingleObservationSet, the size of the observation "
                        "dependent variables input "
                        "should be consistent with the total dependent variable size." );
            }
        }

        setId_ = dataset_->addObservationSet( observableType_,
                                              linkEnds,
                                              observations,
                                              observationTimes,
                                              referenceLinkEnd_,
                                              observationsDependentVariables,
                                              dependentVariableBookkeeping_,
                                              ancillarySettings_,
                                              weights,
                                              residuals,
                                              true,
                                              eraseDuplicates );
    }

    SingleObservationSet( const std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > >& dataset, const int setId ):
        observableType_( dataset->getObservationSetMetadata( setId ).observableType_ ),
        referenceLinkEnd_( dataset->getObservationSetMetadata( setId ).referenceLinkEnd_ ),
        dependentVariableBookkeeping_(
                dataset->getDependentVariableBookkeeping( dataset->getObservationSetMetadata( setId ).dependentVariableLayoutId_ ) ),
        ancillarySettings_( dataset->getAncillarySettings( dataset->getObservationSetMetadata( setId ).ancillarySettingsId_ ) ),
        filteredObservationSet_( nullptr ), dataset_( dataset ), setId_( setId )
    {
        if( dataset_ == nullptr )
        {
            throw std::runtime_error( "Error when creating SingleObservationSet wrapper, input dataset is null." );
        }
    }

    ObservableType getObservableType( )
    {
        return dataset_->getObservationSetMetadata( setId_ ).observableType_;
    }

    LinkDefinition getLinkEnds( )
    {
        return dataset_->getLinkDefinition( dataset_->getObservationSetMetadata( setId_ ).linkDefinitionId_ );
    }

    void setLinkEnds( LinkDefinition& linkEnds )
    {
        dataset_->resetLinkDefinitionForSet( setId_, linkEnds );
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservations( )
    {
        return dataset_->getObservationsForSet( setId_ );
    }

    void setObservations( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations )
    {
        dataset_->setObservationsForSet( setId_, observations );
    }

    void setObservations( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationsVector )
    {
        dataset_->setObservationVectorForSet( setId_, observationsVector );
    }

    void setResiduals( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals )
    {
        dataset_->setResidualsForSet( setId_, residuals );
    }

    void setResiduals( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualsVector )
    {
        dataset_->setResidualVectorForSet( setId_, residualsVector );
    }

    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& getObservationsReference( )
    {
        observations_ = dataset_->getObservationsForSet( setId_ );
        return observations_;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservation( const unsigned int index )
    {
        if( index >= dataset_->getNumberOfObservationsForSet( setId_ ) )
        {
            throw std::runtime_error( "Error when retrieving single observation, index is out of bounds" );
        }
        return dataset_->getObservationValue( dataset_->getObservationIdsForSet( setId_ ).at( index ) );
    }

    void setObservation( const unsigned int index, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observation )
    {
        if( index >= dataset_->getNumberOfObservationsForSet( setId_ ) )
        {
            throw std::runtime_error( "Error when setting single observation value, index is out of bounds" );
        }
        if( observation.size( ) != dataset_->getObservationSetMetadata( setId_ ).observableSize_ )
        {
            throw std::runtime_error(
                    "Error when setting single observation value, the observation size is "
                    "inconsistent." );
        }
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations = getObservations( );
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > updatedObservations = observations;
        updatedObservations.at( index ) = observation;
        setObservations( updatedObservations );
    }

    std::vector< TimeType > getObservationTimes( )
    {
        return dataset_->getObservationTimesForSet( setId_ );
    }

    TimeType getObservationTime( unsigned int index ) const
    {
        if( index >= dataset_->getNumberOfObservationsForSet( setId_ ) )
        {
            throw std::runtime_error(
                    "Error when retrieving single observation time, required index incompatible "
                    "with number of observations." );
        }
        return dataset_->getObservationTime( dataset_->getObservationIdsForSet( setId_ ).at( index ) );
    }

    const std::vector< TimeType >& getObservationTimesReference( )
    {
        observationTimes_ = dataset_->getObservationTimesForSet( setId_ );
        return observationTimes_;
    }

    LinkEndType getReferenceLinkEnd( )
    {
        return dataset_->getObservationSetMetadata( setId_ ).referenceLinkEnd_;
    }

    unsigned int getNumberOfObservables( )
    {
        return static_cast< unsigned int >( dataset_->getNumberOfObservationsForSet( setId_ ) );
    }

    unsigned int getSingleObservableSize( ) const
    {
        return dataset_->getObservationSetMetadata( setId_ ).observableSize_;
    }

    unsigned int getTotalObservationSetSize( ) const
    {
        return dataset_->getTotalScalarSizeForSet( setId_ );
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationsVector( )
    {
        return dataset_->getObservationVectorForSet( setId_ );
    }

    std::pair< TimeType, TimeType > getTimeBounds( )
    {
        if( dataset_->getNumberOfObservationsForSet( setId_ ) == 0 )
        {
            throw std::runtime_error( "Error when getting time bounds of observation set, no observations found" );
        }
        return dataset_->getTimeBoundsForSet( setId_ );
    }

    std::map< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservationsHistory( )
    {
        return utilities::createMapFromVectors< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >(
                getObservationTimes( ), getObservations( ) );
    }

    std::vector< Eigen::VectorXd > getObservationsDependentVariables( )
    {
        return dataset_->getDependentVariablesForSet( setId_ );
    }

    Eigen::MatrixXd getObservationsDependentVariablesMatrix( )
    {
        const std::vector< Eigen::VectorXd > observationsDependentVariables = getObservationsDependentVariables( );
        Eigen::MatrixXd dependentVariablesMatrix = Eigen::MatrixXd::Zero( dataset_->getNumberOfObservationsForSet( setId_ ),
                                                                          dependentVariableBookkeeping_->getTotalDependentVariableSize( ) );
        for( unsigned int i = 0; i < observationsDependentVariables.size( ); i++ )
        {
            dependentVariablesMatrix.block( i, 0, 1, dependentVariableBookkeeping_->getTotalDependentVariableSize( ) ) =
                    observationsDependentVariables.at( i ).transpose( );
        }
        return dependentVariablesMatrix;
    }

    //! Function returning the dependent variable values for a single observation (indicated by index)
    Eigen::VectorXd getDependentVariablesForSingleObservation( unsigned int index ) const
    {
        if( index >= dataset_->getNumberOfObservationsForSet( setId_ ) )
        {
            throw std::runtime_error(
                    "Error when retrieving observation dependent variables for single observation, "
                    "required index incompatible with number of observations." );
        }
        return dataset_->getDependentVariables( dataset_->getObservationIdsForSet( setId_ ).at( index ) );
    }

    //! Function returning the values of a single dependent variable (specified by dependent variable settings)
    Eigen::MatrixXd getSingleDependentVariable( std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings,
                                                const bool returnFirstCompatibleSettings = false )
    {
        return dataset_->getSingleDependentVariableForSet( setId_, dependentVariableSettings, returnFirstCompatibleSettings );
    }

    //! Function returning the list of all dependent variable settings compatible with the settings provided as inputs
    //! (which might not be fully defined, i.e. with missing link ends information, etc.)
    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > getCompatibleDependentVariablesSettingsList(
            std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings ) const
    {
        return dataset_->getCompatibleDependentVariableSettingsForSet( setId_, dependentVariableSettings );
    }

    //! Function returning a vector containing the values of all dependent variables compatible with the settings provided as input
    //! The order in which they are provided matches the list of compatible settings given by the getCompatibleDependentVariablesSettingsList function
    std::vector< Eigen::MatrixXd > getAllCompatibleDependentVariables(
            std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings ) const
    {
        return dataset_->getAllCompatibleDependentVariablesForSet( setId_, dependentVariableSettings );
    }

    std::vector< Eigen::VectorXd >& getObservationsDependentVariablesReference( )
    {
        observationsDependentVariables_ = dataset_->getDependentVariablesForSet( setId_ );
        return observationsDependentVariables_;
    }

    //! Function to reset the observation dependent variable values
    void setObservationsDependentVariables( std::vector< Eigen::VectorXd >& dependentVariables )
    {
        if( dependentVariables.size( ) > 0 )
        {
            if( dependentVariables.size( ) != dataset_->getNumberOfObservationsForSet( setId_ ) )
            {
                throw std::runtime_error(
                        "Error when resetting observation dependent variables in "
                        "SingleObservationSet, the input size should be consistent "
                        "with the number of observations." );
            }
            if( ( dependentVariableBookkeeping_ != nullptr ) &&
                ( dependentVariables[ 0 ].size( ) != dependentVariableBookkeeping_->getTotalDependentVariableSize( ) ) )
            {
                throw std::runtime_error(
                        "Error when resetting observation dependent variables in "
                        "SingleObservationSet, the size of the observation dependent variables "
                        "input "
                        "should be consistent with the total dependent variable size." );
            }
            dataset_->setDependentVariablesForSet( setId_, dependentVariables );
        }
        else
        {
            dataset_->clearDependentVariablesForSet( setId_ );
        }
    }

    std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > getDependentVariableBookkeeping( )
    {
        return dependentVariableBookkeeping_;
    }

    //! Function that returns the time history of all observation dependent variables. It must be noted that the reported epochs are the times at which the
    //! observations are computed/acquired, which might differ from the times at which the dependent variables are evaluated.
    std::map< TimeType, Eigen::VectorXd > getDependentVariableHistory( )
    {
        return utilities::createMapFromVectors< TimeType, Eigen::VectorXd >( getObservationTimes( ), getObservationsDependentVariables( ) );
    }

    //! Function that returns the time history of a single dependent variables (specified by settings). It must be noted that the reported epochs are the times
    //! at which the observations are computed/acquired, which might differ from the times at which the dependent variables are evaluated.
    std::map< TimeType, Eigen::VectorXd > getSingleDependentVariableHistory(
            std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings,
            const bool returnFirstCompatibleSettings = false )
    {
        Eigen::MatrixXd singleDependentVariableValues =
                getSingleDependentVariable( dependentVariableSettings, returnFirstCompatibleSettings );
        std::map< TimeType, Eigen::VectorXd > singleDependentVariableMap;
        for( unsigned int i = 0; i < dataset_->getNumberOfObservationsForSet( setId_ ); i++ )
        {
            Eigen::VectorXd dependentVariableCurrentTime =
                    singleDependentVariableValues.block( i, 0, 1, singleDependentVariableValues.cols( ) ).transpose( );
            singleDependentVariableMap[ getObservationTime( i ) ] = dependentVariableCurrentTime;
        }
        return singleDependentVariableMap;
    }

    std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > getAncillarySettings( )
    {
        return ancillarySettings_;
    }

    std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > getObservationDataset( ) const
    {
        return dataset_;
    }

    int getObservationSetId( ) const
    {
        return setId_;
    }

    void resetObservationDatasetReference( const std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > >& dataset,
                                           const int setId )
    {
        if( dataset == nullptr )
        {
            throw std::runtime_error( "Error when resetting SingleObservationSet dataset reference, dataset is null." );
        }
        dataset_ = dataset;
        setId_ = setId;
    }

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getWeights( ) const
    {
        return dataset_->getWeightsForSet( setId_ );
    }

    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& getWeightsReference( )
    {
        weights_ = dataset_->getWeightsForSet( setId_ );
        return weights_;
    }

    Eigen::VectorXd getWeightsVector( ) const
    {
        return dataset_->getWeightVectorForSet( setId_ );
    }

    Eigen::Matrix< double, Eigen::Dynamic, 1 > getWeight( unsigned int index ) const
    {
        if( index >= dataset_->getNumberOfObservationsForSet( setId_ ) )
        {
            throw std::runtime_error(
                    "Error when retrieving single observation weight, required index incompatible "
                    "with number of observations." );
        }
        return dataset_->getWeightValue( dataset_->getObservationIdsForSet( setId_ ).at( index ) );
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getResiduals( ) const
    {
        return dataset_->getResidualsForSet( setId_ );
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidualsVector( ) const
    {
        return dataset_->getResidualVectorForSet( setId_ );
    }

    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& getResidualsReference( )
    {
        residuals_ = dataset_->getResidualsForSet( setId_ );
        return residuals_;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidual( unsigned int index ) const
    {
        if( index >= dataset_->getNumberOfObservationsForSet( setId_ ) )
        {
            throw std::runtime_error(
                    "Error when retrieving single observation residual, required index "
                    "incompatible with number of observations." );
        }
        return dataset_->getResidualValue( dataset_->getObservationIdsForSet( setId_ ).at( index ) );
    }

    Eigen::VectorXd getRmsResiduals( )
    {
        return dataset_->getRmsResidualsForSet( setId_ );
    }

    Eigen::VectorXd getMeanResiduals( )
    {
        return dataset_->getMeanResidualsForSet( setId_ );
    }

    void setConstantWeight( const double weight )
    {
        dataset_->setConstantSingleObservationScalarWeightForSet( setId_, weight );
    }

    void setConstantWeight( const Eigen::Matrix< double, Eigen::Dynamic, 1 >& weight )
    {
        dataset_->setConstantSingleObservationDiagonalWeightForSet( setId_, weight );
    }

    void setWeights( const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights )
    {
        const unsigned int numberOfObservations = dataset_->getNumberOfObservationsForSet( setId_ );
        const unsigned int singleObservationSize = dataset_->getObservationSetMetadata( setId_ ).observableSize_;
        if( weights.size( ) != numberOfObservations )
        {
            throw std::runtime_error(
                    "Error when settings weights in single observation set, numbers of weights and observations are inconsistent." );
        }

        Eigen::VectorXd weightsVector( numberOfObservations * singleObservationSize );
        for( unsigned int k = 0; k < weights.size( ); k++ )
        {
            if( weights[ k ].size( ) != singleObservationSize )
            {
                throw std::runtime_error(
                        "Error when settings weights in single observation set, size of single weight entry is inconsistent with single "
                        "observation size." );
            }
            weightsVector.segment( k * singleObservationSize, singleObservationSize ) = weights[ k ];
        }
        dataset_->setWeightVectorForSet( setId_, weightsVector );
    }

    void setTabulatedWeights( const Eigen::VectorXd& weightsVector )
    {
        dataset_->setWeightVectorForSet( setId_, weightsVector );
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getComputedObservations( ) const
    {
        return dataset_->getComputedObservationsForSet( setId_ );
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getComputedObservationsVector( ) const
    {
        return dataset_->getComputedObservationVectorForSet( setId_ );
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

    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > getFilteredObservationSet( ) const
    {
        return filteredObservationSet_;
    }

    void removeSingleObservation( unsigned int indexToRemove )
    {
        if( indexToRemove >= dataset_->getNumberOfObservationsForSet( setId_ ) )
        {
            throw std::runtime_error(
                    "Error when removing single observation from SingleObservationSet, index "
                    "incompatible with number of observations." );
        }

        dataset_->removeObservationsFromSet( setId_, std::vector< unsigned int >( { indexToRemove } ) );
    }

    void removeObservations( const std::vector< unsigned int >& indicesToRemove )
    {
        dataset_->removeObservationsFromSet( setId_, indicesToRemove );
    }

    void eraseDuplicateObservations( )
    {
        dataset_->eraseDuplicateObservationsFromSet( setId_ );
    }

    void filterObservations( const std::shared_ptr< ObservationFilterBase > observationFilter, const bool saveFilteredObservations = true )
    {
        if( observationFilter->filterOut( ) && filteredObservationSet_ == nullptr )
        {
            filteredObservationSet_ = std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                    observableType_,
                    getLinkEnds( ),
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                    std::vector< TimeType >( ),
                    referenceLinkEnd_,
                    std::vector< Eigen::VectorXd >( ),
                    dependentVariableBookkeeping_,
                    ancillarySettings_ );
        }
        if( !observationFilter->filterOut( ) && filteredObservationSet_ == nullptr )
        {
            throw std::runtime_error(
                    "Error when attempting to un-filter observations, filtered observation set is "
                    "empty." );
        }

        if( observationFilter->filterOut( ) )
        {
            const std::vector< unsigned int > indicesToRemove = dataset_->getFilteredObservationIndices( setId_, observationFilter );
            if( saveFilteredObservations )
            {
                dataset_->moveObservationsToSet( setId_, *filteredObservationSet_->getObservationDataset( ), 0, indicesToRemove, true );
            }
            else
            {
                dataset_->removeObservationsFromSet( setId_, indicesToRemove );
            }
        }
        else
        {
            const std::vector< unsigned int > indicesToRemove =
                    filteredObservationSet_->getObservationDataset( )->getFilteredObservationIndices( 0, observationFilter );
            filteredObservationSet_->getObservationDataset( )->moveObservationsToSet( 0, *dataset_, setId_, indicesToRemove, true );
        }
    }

    void addObservations( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                          const std::vector< TimeType >& times,
                          const std::vector< Eigen::VectorXd >& dependentVariables = {},
                          const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights = {},
                          const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals = {},
                          const bool sortObservations = true )
    {
        dataset_->addObservationsToSet( setId_, observations, times, dependentVariables, weights, residuals, sortObservations );
    }

    void addDependentVariables(
            const std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > dependentVariableSettings )
    {
        if( dependentVariableBookkeeping_ == nullptr )
        {
            dependentVariableBookkeeping_ =
                    std::make_shared< ObservationDependentVariableBookkeeping >( observableType_, getLinkEnds( ).linkEnds_ );
            dataset_->resetDependentVariableBookkeepingForSet( setId_, dependentVariableBookkeeping_ );
        }
        else if( dataset_->getDependentVariablesForSet( setId_ ).size( ) != 0 )
        {
            throw std::runtime_error(
                    "Error, cannot add dependent variable settings to SingleObservationSet that has dependent variables calculated "
                    "already" );
        }
        else
        {
            dependentVariableBookkeeping_ = std::make_shared< ObservationDependentVariableBookkeeping >( *dependentVariableBookkeeping_ );
            dataset_->resetDependentVariableBookkeepingForSet( setId_, dependentVariableBookkeeping_ );
        }
        dependentVariableBookkeeping_->addDependentVariables( dependentVariableSettings );
    }

    void clearDependentVariableValues( )
    {
        dataset_->clearDependentVariablesForSet( setId_ );
    }

private:
    //! Function extracting the values of a single dependent variable
    Eigen::MatrixXd getSingleDependentVariable( std::pair< int, int > dependentVariableIndexAndSize ) const
    {
        return dataset_->getSingleDependentVariableForSet( setId_, dependentVariableIndexAndSize );
    }

    const ObservableType observableType_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations_;

    std::vector< TimeType > observationTimes_;

    const LinkEndType referenceLinkEnd_;

    std::vector< Eigen::VectorXd > observationsDependentVariables_;

    std::shared_ptr< ObservationDependentVariableBookkeeping > dependentVariableBookkeeping_;

    const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings_;

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights_;

    //    Eigen::VectorXd weightsVector_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals_;

    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > filteredObservationSet_;

    std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > dataset_;

    int setId_;
};

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createObservationDataset(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > >& observationSet )
{
    if( observationSet == nullptr )
    {
        throw std::runtime_error( "Error when creating observation dataset, input set is null." );
    }

    std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > dataset =
            std::make_shared< ObservationDataset< ObservationScalarType, TimeType > >( );
    dataset->addObservationSetFromDataset( *observationSet->getObservationDataset( ), observationSet->getObservationSetId( ) );
    return dataset;
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > createSingleObservationSet(
        const std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > >& dataset,
        const int setId = 0 )
{
    if( dataset == nullptr )
    {
        throw std::runtime_error( "Error when creating single observation set wrapper, input dataset is null." );
    }
    return std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >( dataset, setId );
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > filterObservations(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > singleObservationSet,
        const std::shared_ptr< ObservationFilterBase > observationFilter,
        const bool saveFilteredObservations = false )
{
    if( !observationFilter->filterOut( ) )
    {
        throw std::runtime_error(
                "Error when creating new single observation set post-filtering, the filterOut "
                "option should be set to true" );
    }
    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > newObservationSet =
            createSingleObservationSet( createObservationDataset( singleObservationSet ) );

    // Filter observations from new observation set
    newObservationSet->filterObservations( observationFilter, saveFilteredObservations );

    return newObservationSet;
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > splitObservationSet(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > observationSet,
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

    switch( observationSetSplitter->getSplitterType( ) )
    {
        case time_tags_splitter: {
            std::vector< double > timeTags =
                    std::dynamic_pointer_cast< ObservationSetSplitter< std::vector< double > > >( observationSetSplitter )
                            ->getSplitterValue( );
            for( auto currentTimeTag : timeTags )
            {
                if( currentTimeTag > observationSet->getTimeBounds( ).first )
                {
                    bool detectedStartSet = false;
                    int indexObs = rawStartIndicesNewSets.at( rawStartIndicesNewSets.size( ) - 1 );
                    while( !detectedStartSet && indexObs < static_cast< int >( observationTimes.size( ) ) )
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
            double maxTimeInterval =
                    std::dynamic_pointer_cast< ObservationSetSplitter< double > >( observationSetSplitter )->getSplitterValue( );
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
            double maxTimeSpan =
                    std::dynamic_pointer_cast< ObservationSetSplitter< double > >( observationSetSplitter )->getSplitterValue( );
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
            int maxNbObs = std::dynamic_pointer_cast< ObservationSetSplitter< int > >( observationSetSplitter )->getSplitterValue( );
            if( maxNbObs < observationSetSplitter->getMinNumberObservations( ) )
            {
                throw std::runtime_error(
                        "Error when splitting observation sets, the maximum number of observations "
                        "cannot be smaller than the minimum number of observations." );
            }
            for( int ind = maxNbObs; ind < static_cast< int >( observationSet->getNumberOfObservables( ) ); ind += maxNbObs )
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
        if( ( rawStartIndicesNewSets.at( j ) - rawStartIndicesNewSets.at( j - 1 ) ) >= observationSetSplitter->getMinNumberObservations( ) )
        {
            indicesNewSets.push_back( std::make_pair( rawStartIndicesNewSets.at( j - 1 ),
                                                      rawStartIndicesNewSets.at( j ) - rawStartIndicesNewSets.at( j - 1 ) ) );
        }
    }

    // Split current observation set based on indices
    std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > newObsSets;
    for( unsigned int k = 0; k < indicesNewSets.size( ); k++ )
    {
        int startIndex = indicesNewSets.at( k ).first;
        int sizeCurrentSet = indicesNewSets.at( k ).second;

        const std::vector< unsigned int >& sourceObservationIds =
                observationSet->getObservationDataset( )->getObservationIdsForSet( observationSet->getObservationSetId( ) );
        std::set< unsigned int > selectedObservationIds;
        for( int i = startIndex; i < startIndex + sizeCurrentSet; ++i )
        {
            selectedObservationIds.insert( sourceObservationIds.at( i ) );
        }
        const ObservationSelectionCondition< ObservationScalarType, TimeType > selectedRowsCondition(
                [ selectedObservationIds ]( const ObservationDataset< ObservationScalarType, TimeType >&, const int observationId ) {
                    return selectedObservationIds.count( observationId ) > 0;
                } );
        const std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > splitDataset =
                observationSet->getObservationDataset( )->createNewAndKeep( selectedRowsCondition );
        std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > newSet = createSingleObservationSet( splitDataset );

        newObsSets.push_back( newSet );
    }

    return newObsSets;
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_SINGLE_OBSERVATION_SET_H
