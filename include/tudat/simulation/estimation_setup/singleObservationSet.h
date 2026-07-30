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

#include <cmath>
#include <algorithm>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <functional>
#include <memory>
#include <string>
#include <utility>
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

//! Enumeration of supported base weight representations in a single observation set.
enum ObservationWeightsMatrixType {
    //! Independent per-component weights per observation.
    diagonal_weights_matrix = 0,
    //! Dense per-observation blocks, with no coupling between different observations.
    block_diagonal_weights_matrix = 1
};

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

        observableType_( observableType ), linkEnds_( linkEnds ), observations_( observations ), observationTimes_( observationTimes ),
        referenceLinkEnd_( referenceLinkEnd ), observationsDependentVariables_( observationsDependentVariables ),
        dependentVariableBookkeeping_( dependentVariableBookkeeping ), ancillarySettings_( ancillarySettings ),
        numberOfObservations_( observations_.size( ) ), residuals_( residuals )
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

        validateObservationsAndTimesSize( observations_.size( ), observationTimes_.size( ), "making SingleObservationSet" );
        validateConsistentObservationDimensions( observations_, "making SingleObservationSet" );

        singleObservationSize_ = getObservableSize( observableType );
        if( getObservableSize( observableType ) > 0 && !observations_.empty( ) )
        {
            validateObservationDimensionsAgainstSingleSize( observations_, "making SingleObservationSet" );
        }
        // Initialise weights
        if( weights.size( ) == 0 )
        {
            weightState_.diagonalWeights.reserve( numberOfObservations_ );
            for( unsigned int k = 0; k < numberOfObservations_; k++ )
            {
                weightState_.diagonalWeights.push_back( Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( singleObservationSize_, 1 ) );
            }
        }
        else
        {
            validateOptionalBatchSize( weights, observationTimes.size( ), "weights", "creating observation set with weights" );
            validatePerObservationVectorSize( weights, "weight", "creating observation set with weights" );
            weightState_.diagonalWeights = weights;
        }

        if( residuals.size( ) == 0 )
        {
            // Initialise residuals
            for( unsigned int k = 0; k < numberOfObservations_; k++ )
            {
                residuals_.push_back( Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_, 1 ) );
            }
        }
        else
        {
            validateOptionalBatchSize( residuals, observationTimes.size( ), "residuals", "creating observation set with residuals" );
            validatePerObservationVectorSize( residuals, "residual", "creating observation set with residuals" );
        }

        validateDependentVariableBatch( observationsDependentVariables_, numberOfObservations_, true, "creating SingleObservationSet" );

        // Sort observations and metadata per observation time
        orderObservationsAndMetadata( );

        // Erase duplicate observations if requested
        if( eraseDuplicates )
        {
            eraseDuplicateObservations( );
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
        if( observations.size( ) != observations_.size( ) )
        {
            throw std::runtime_error( "Error when resetting observations, number of observations is incompatible." );
        }
        observations_ = observations;
        if( !observations_.empty( ) )
        {
            singleObservationSize_ = observations_.at( 0 ).size( );
        }
    }

    void setObservations( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationsVector )
    {
        if( observationsVector.size( ) != numberOfObservations_ * singleObservationSize_ )
        {
            throw std::runtime_error( "Error when resetting observations, number of observations is incompatible." );
        }

        observations_.clear( );
        for( unsigned int k = 0; k < numberOfObservations_; k++ )
        {
            observations_.push_back( observationsVector.segment( k * singleObservationSize_, singleObservationSize_ ) );
        }
        if( !observations_.empty( ) )
        {
            singleObservationSize_ = observations_.at( 0 ).size( );
        }
    }

    void setResiduals( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals )
    {
        if( residuals.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Error when setting residuals, number of observations is inconsistent." );
        }
        residuals_ = residuals;
    }

    void setResiduals( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualsVector )
    {
        if( residualsVector.size( ) != numberOfObservations_ * singleObservationSize_ )
        {
            throw std::runtime_error( "Error when setting residuals, number of observations is inconsistent." );
        }

        residuals_.clear( );
        for( unsigned int k = 0; k < numberOfObservations_; k++ )
        {
            residuals_.push_back( residualsVector.segment( k * singleObservationSize_, singleObservationSize_ ) );
        }
    }

    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& getObservationsReference( )
    {
        return observations_;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservation( const unsigned int index )
    {
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when retrieving single observation, index is out of bounds" );
        }
        return observations_.at( index );
    }

    void setObservation( const unsigned int index, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observation )
    {
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when setting single observation value, index is out of bounds" );
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
        return utilities::createMapFromVectors< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( observationTimes_,
                                                                                                                       observations_ );
    }

    std::vector< Eigen::VectorXd > getObservationsDependentVariables( )
    {
        return observationsDependentVariables_;
    }

    Eigen::MatrixXd getObservationsDependentVariablesMatrix( )
    {
        Eigen::MatrixXd dependentVariablesMatrix =
                Eigen::MatrixXd::Zero( numberOfObservations_, dependentVariableBookkeeping_->getTotalDependentVariableSize( ) );
        for( unsigned int i = 0; i < observationsDependentVariables_.size( ); i++ )
        {
            dependentVariablesMatrix.block( i, 0, 1, dependentVariableBookkeeping_->getTotalDependentVariableSize( ) ) =
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
    Eigen::MatrixXd getSingleDependentVariable( std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings,
                                                const bool returnFirstCompatibleSettings = false )
    {
        // Retrieve full map of dependent variables start indices and sizes based on settings
        std::map< std::pair< int, int >, std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > >
                settingsIndicesAndSizes = dependentVariableBookkeeping_->getSettingsIndicesAndSizes( );

        // Get the start indices and sizes of all dependent variables that would be compatible with
        // the settings provided as inputs.
        std::vector< std::pair< int, int > > indicesAndSizes;
        for( auto it : settingsIndicesAndSizes )
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
    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > getCompatibleDependentVariablesSettingsList(
            std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings ) const
    {
        // Check which settings are compatible with the input settings object
        std::vector< std::shared_ptr< ObservationDependentVariableSettings > > compatibleSettings;
        for( auto it : dependentVariableBookkeeping_->getDependentVariableSettings( ) )
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
            std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings ) const
    {
        // Retrieve start indices and sizes for each dependent variable settings
        std::map< std::pair< int, int >, std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > >
                settingsIndicesAndSizes = dependentVariableBookkeeping_->getSettingsIndicesAndSizes( );

        // Retrieve all relevant dependent variables
        std::vector< Eigen::MatrixXd > dependentVariablesList;
        for( auto it : settingsIndicesAndSizes )
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
            if( ( dependentVariableBookkeeping_ != nullptr ) &&
                ( dependentVariables[ 0 ].size( ) != dependentVariableBookkeeping_->getTotalDependentVariableSize( ) ) )
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

    std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > getDependentVariableBookkeeping( )
    {
        return dependentVariableBookkeeping_;
    }

    //! Function that returns the time history of all observation dependent variables. It must be noted that the reported epochs are the times at which the
    //! observations are computed/acquired, which might differ from the times at which the dependent variables are evaluated.
    std::map< TimeType, Eigen::VectorXd > getDependentVariableHistory( )
    {
        return utilities::createMapFromVectors< TimeType, Eigen::VectorXd >( observationTimes_, observationsDependentVariables_ );
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
        for( unsigned int i = 0; i < numberOfObservations_; i++ )
        {
            Eigen::VectorXd dependentVariableCurrentTime =
                    singleDependentVariableValues.block( i, 0, 1, singleDependentVariableValues.cols( ) ).transpose( );
            singleDependentVariableMap[ observationTimes_[ i ] ] = dependentVariableCurrentTime;
        }
        return singleDependentVariableMap;
    }

    std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > getAncillarySettings( )
    {
        return ancillarySettings_;
    }

    //! Returns per-observation diagonal weights including any full-matrix diagonal contribution.
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getWeights( ) const
    {
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > combinedWeights;
        Eigen::VectorXd combinedDiagonal = getWeightsDiagonalVector( );
        combinedWeights.reserve( numberOfObservations_ );
        for( unsigned int i = 0; i < numberOfObservations_; i++ )
        {
            combinedWeights.push_back( combinedDiagonal.segment( i * singleObservationSize_, singleObservationSize_ ) );
        }
        return combinedWeights;
    }

    //! Returns a reference to stored base diagonal weights (diagonal base mode only).
    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& getWeightsReference( )
    {
        if( weightState_.matrixType != diagonal_weights_matrix )
        {
            throw std::runtime_error(
                    "Error when retrieving diagonal weight vectors by reference, base weight matrix type is not diagonal." );
        }
        return weightState_.diagonalWeights;
    }

    //! Returns the diagonal vector of the base weights (without full-matrix contribution).
    Eigen::VectorXd getBaseWeightsDiagonalVector( ) const
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > weightsVector =
                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( singleObservationSize_ * numberOfObservations_, 1 );

        if( weightState_.matrixType == diagonal_weights_matrix )
        {
            if( weightState_.diagonalWeights.size( ) != numberOfObservations_ )
            {
                throw std::runtime_error(
                        "Error when retrieving base diagonal weights vector, number of stored diagonal weights is inconsistent." );
            }
            for( unsigned int i = 0; i < numberOfObservations_; i++ )
            {
                weightsVector.block( i * singleObservationSize_, 0, singleObservationSize_, 1 ) = weightState_.diagonalWeights.at( i );
            }
        }
        else if( weightState_.matrixType == block_diagonal_weights_matrix )
        {
            for( unsigned int i = 0; i < numberOfObservations_; i++ )
            {
                weightsVector.block( i * singleObservationSize_, 0, singleObservationSize_, 1 ) =
                        weightState_.blockWeights.at( i ).diagonal( );
            }
        }

        return weightsVector;
    }

    //! Returns the diagonal vector of the combined weights matrix.
    Eigen::VectorXd getWeightsDiagonalVector( ) const
    {
        Eigen::VectorXd weightsVector = getBaseWeightsDiagonalVector( );
        if( weightState_.fullWeights.rows( ) > 0 )
        {
            weightsVector += weightState_.fullWeights.diagonal( );
        }
        return weightsVector;
    }

    //! Returns the diagonal-weight segment for a single observation index.
    Eigen::Matrix< double, Eigen::Dynamic, 1 > getWeight( unsigned int index ) const
    {
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when retrieving single observation weight, required index incompatible "
                    "with number of observations." );
        }
        return getWeightsDiagonalVector( ).segment( index * singleObservationSize_, singleObservationSize_ );
    }

    //! Returns which base representation is currently used for the weights.
    ObservationWeightsMatrixType getWeightsMatrixType( ) const
    {
        return weightState_.matrixType;
    }

    //! Checks whether the combined weights matrix currently contains off-diagonal terms.
    bool hasOffDiagonalWeights( ) const
    {
        if( weightState_.matrixType == block_diagonal_weights_matrix )
        {
            return true;
        }

        if( weightState_.fullWeights.rows( ) > 0 )
        {
            Eigen::MatrixXd offDiagonalContribution = weightState_.fullWeights;
            offDiagonalContribution.diagonal( ).setZero( );
            return !( offDiagonalContribution.isZero( 0.0 ) );
        }

        return false;
    }

    //! Returns per-observation dense blocks for block-diagonal base mode.
    std::vector< Eigen::MatrixXd > getBlockDiagonalWeightMatrices( ) const
    {
        if( weightState_.matrixType != block_diagonal_weights_matrix )
        {
            throw std::runtime_error( "Error when retrieving block-diagonal weights, current weight matrix type is not block-diagonal." );
        }
        return weightState_.blockWeights;
    }

    //! Returns the optional full-matrix contribution added to the base weights.
    Eigen::MatrixXd getFullWeightMatrix( ) const
    {
        if( weightState_.fullWeights.rows( ) == 0 )
        {
            throw std::runtime_error( "Error when retrieving full weights matrix contribution, contribution is not set." );
        }
        return weightState_.fullWeights;
    }

    //! Indicates whether an additional full-matrix weight contribution is present.
    bool hasFullWeightMatrixContribution( ) const
    {
        return ( weightState_.fullWeights.rows( ) > 0 );
    }

    //! Builds and returns the full combined weight matrix for active observations.
    Eigen::MatrixXd getWeightMatrix( ) const
    {
        const int totalObservationSetSize = static_cast< int >( getTotalObservationSetSize( ) );
        Eigen::MatrixXd weightMatrix = Eigen::MatrixXd::Zero( totalObservationSetSize, totalObservationSetSize );

        if( weightState_.matrixType == diagonal_weights_matrix )
        {
            if( weightState_.diagonalWeights.size( ) != numberOfObservations_ )
            {
                throw std::runtime_error(
                        "Error when retrieving weights matrix in single observation set, number of diagonal weights is inconsistent." );
            }
            weightMatrix.diagonal( ) = Eigen::VectorXd::Zero( totalObservationSetSize );
            for( unsigned int i = 0; i < numberOfObservations_; i++ )
            {
                weightMatrix.block( i * singleObservationSize_, i * singleObservationSize_, singleObservationSize_, singleObservationSize_ )
                        .diagonal( ) = weightState_.diagonalWeights.at( i );
            }
        }
        else if( weightState_.matrixType == block_diagonal_weights_matrix )
        {
            for( unsigned int i = 0; i < numberOfObservations_; i++ )
            {
                const int startIndex = static_cast< int >( i * singleObservationSize_ );
                weightMatrix.block( startIndex, startIndex, singleObservationSize_, singleObservationSize_ ) =
                        weightState_.blockWeights.at( i );
            }
        }

        if( weightState_.fullWeights.rows( ) > 0 )
        {
            weightMatrix += weightState_.fullWeights;
        }
        return weightMatrix;
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getResiduals( ) const
    {
        return residuals_;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidualsVector( ) const
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualsVector =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_ * numberOfObservations_, 1 );
        for( unsigned int i = 0; i < numberOfObservations_; i++ )
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
            for( unsigned int j = 0; j < numberOfObservations_; j++ )
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

    //! Sets one constant diagonal-weight vector for all observations.
    void setConstantWeight( const double weight )
    {
        setDiagonalWeightState(
                createUniformDiagonalWeights( weight * Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( singleObservationSize_, 1 ) ) );
        validateCombinedWeightsMatrix( );
    }

    //! Sets one component-wise diagonal-weight vector for all observations.
    void setConstantWeight( const Eigen::Matrix< double, Eigen::Dynamic, 1 >& weight )
    {
        if( weight.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error(
                    "Error when setting constant weight in single observation set, weight size is "
                    "inconsistent with single observation size." );
        }
        setDiagonalWeightState( createUniformDiagonalWeights( weight ) );
        validateCombinedWeightsMatrix( );
    }

    //! Sets per-observation diagonal weights from one concatenated vector.
    void setTabulatedWeights( const Eigen::VectorXd& weightsVector )
    {
        if( weightsVector.rows( ) != static_cast< int >( singleObservationSize_ * observations_.size( ) ) )
        {
            throw std::runtime_error(
                    "Error when setting weights in single observation set, sizes are "
                    "incompatible." );
        }
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > diagonalWeights;
        diagonalWeights.reserve( numberOfObservations_ );
        for( unsigned int k = 0; k < numberOfObservations_; k++ )
        {
            Eigen::Matrix< double, Eigen::Dynamic, 1 > currentWeights =
                    Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( singleObservationSize_, 1 );
            for( unsigned int i = 0; i < singleObservationSize_; i++ )
            {
                currentWeights[ i ] = weightsVector[ k * singleObservationSize_ + i ];
            }
            diagonalWeights.push_back( currentWeights );
        }
        setDiagonalWeightState( std::move( diagonalWeights ) );
        validateCombinedWeightsMatrix( );
    }

    //! Sets per-observation dense blocks as the base block-diagonal weight model.
    void setBlockDiagonalWeights( const std::vector< Eigen::MatrixXd >& blockDiagonalWeights )
    {
        if( blockDiagonalWeights.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when setting block-diagonal weights in single observation set, number of blocks is inconsistent." );
        }

        for( unsigned int i = 0; i < numberOfObservations_; i++ )
        {
            if( blockDiagonalWeights.at( i ).rows( ) != static_cast< int >( singleObservationSize_ ) ||
                blockDiagonalWeights.at( i ).cols( ) != static_cast< int >( singleObservationSize_ ) )
            {
                throw std::runtime_error(
                        "Error when setting block-diagonal weights in single observation set, block size is inconsistent." );
            }
            validateWeightMatrixPositiveDefinite(
                    blockDiagonalWeights.at( i ),
                    "Error when setting block-diagonal weights in single observation set, block matrix is invalid." );
        }
        weightState_.diagonalWeights.clear( );
        weightState_.blockWeights = blockDiagonalWeights;
        weightState_.matrixType = block_diagonal_weights_matrix;
        validateCombinedWeightsMatrix( );
    }

    //! Sets or clears the additional full-matrix weight contribution.
    void setFullWeightMatrix( const Eigen::MatrixXd& fullWeightMatrix )
    {
        if( fullWeightMatrix.rows( ) == 0 && fullWeightMatrix.cols( ) == 0 )
        {
            weightState_.fullWeights.resize( 0, 0 );
            return;
        }

        const int totalObservationSetSize = static_cast< int >( getTotalObservationSetSize( ) );
        if( fullWeightMatrix.rows( ) != totalObservationSetSize || fullWeightMatrix.cols( ) != totalObservationSetSize )
        {
            throw std::runtime_error( "Error when setting full weights matrix in single observation set, matrix size is inconsistent." );
        }

        validateSymmetricWeightMatrix( fullWeightMatrix,
                                       "Error when setting full weights matrix in single observation set, matrix is invalid." );
        weightState_.fullWeights = fullWeightMatrix;
        validateCombinedWeightsMatrix( );
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getComputedObservations( ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > computedObservations;
        for( unsigned int k = 0; k < observations_.size( ); k++ )
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
        for( unsigned int i = 0; i < numberOfObservations_; i++ )
        {
            computedObservationsVector.segment( i * singleObservationSize_, singleObservationSize_ ) = computedObservations.at( i );
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

    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > getFilteredObservationSet( ) const
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
        if( weightState_.matrixType == diagonal_weights_matrix )
        {
            weightState_.diagonalWeights.erase( weightState_.diagonalWeights.begin( ) + indexToRemove );
        }
        if( weightState_.matrixType == block_diagonal_weights_matrix )
        {
            weightState_.blockWeights.erase( weightState_.blockWeights.begin( ) + indexToRemove );
        }
        if( weightState_.fullWeights.rows( ) > 0 )
        {
            const int startIndex = static_cast< int >( indexToRemove * singleObservationSize_ );
            const int reducedSize = weightState_.fullWeights.rows( ) - static_cast< int >( singleObservationSize_ );
            Eigen::MatrixXd reducedWeights = Eigen::MatrixXd::Zero( reducedSize, reducedSize );

            if( startIndex > 0 )
            {
                reducedWeights.block( 0, 0, startIndex, startIndex ) = weightState_.fullWeights.block( 0, 0, startIndex, startIndex );
            }

            const int trailingSize = reducedSize - startIndex;
            if( trailingSize > 0 )
            {
                reducedWeights.block( startIndex, startIndex, trailingSize, trailingSize ) = weightState_.fullWeights.block(
                        startIndex + singleObservationSize_, startIndex + singleObservationSize_, trailingSize, trailingSize );
            }
            if( startIndex > 0 && trailingSize > 0 )
            {
                reducedWeights.block( 0, startIndex, startIndex, trailingSize ) =
                        weightState_.fullWeights.block( 0, startIndex + singleObservationSize_, startIndex, trailingSize );
                reducedWeights.block( startIndex, 0, trailingSize, startIndex ) =
                        weightState_.fullWeights.block( startIndex + singleObservationSize_, 0, trailingSize, startIndex );
            }

            weightState_.fullWeights = reducedWeights;
        }
        if( fullWeightCrossCorrelationWithFilteredSet_.rows( ) > 0 )
        {
            const int startIndex = static_cast< int >( indexToRemove * singleObservationSize_ );
            const int reducedRows = fullWeightCrossCorrelationWithFilteredSet_.rows( ) - static_cast< int >( singleObservationSize_ );
            Eigen::MatrixXd reducedCrossWeights = Eigen::MatrixXd::Zero( reducedRows, fullWeightCrossCorrelationWithFilteredSet_.cols( ) );

            if( startIndex > 0 )
            {
                reducedCrossWeights.block( 0, 0, startIndex, fullWeightCrossCorrelationWithFilteredSet_.cols( ) ) =
                        fullWeightCrossCorrelationWithFilteredSet_.block(
                                0, 0, startIndex, fullWeightCrossCorrelationWithFilteredSet_.cols( ) );
            }

            const int trailingRows = reducedRows - startIndex;
            if( trailingRows > 0 )
            {
                reducedCrossWeights.block( startIndex, 0, trailingRows, fullWeightCrossCorrelationWithFilteredSet_.cols( ) ) =
                        fullWeightCrossCorrelationWithFilteredSet_.block(
                                startIndex + singleObservationSize_, 0, trailingRows, fullWeightCrossCorrelationWithFilteredSet_.cols( ) );
            }

            fullWeightCrossCorrelationWithFilteredSet_ = reducedCrossWeights;
        }

        if( observationsDependentVariables_.size( ) > 0 )
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
        for( auto ind : indicesToRemove )
        {
            removeSingleObservation( ind - counter );  // observations are already filtered and sorted
            counter += 1;
        }
    }

    void eraseDuplicateObservations( )
    {
        std::vector< unsigned int > indicesToRemove;

        // Single pass through sorted observations
        for( unsigned int i = 1; i < numberOfObservations_; i++ )
        {
            // Check if current observation time equals previous observation time
            if( observationTimes_[ i ] == observationTimes_[ i - 1 ] )
            {
                const double currentObsValue = observationTimes_[ i ];
                const double previousObsValue = observationTimes_[ i - 1 ];

                // Check if observation values are also identical (with relative tolerance)
                if( std::abs( currentObsValue - previousObsValue ) <=
                    1e-12 * std::max( std::abs( currentObsValue ), std::abs( previousObsValue ) ) )
                {
                    // Mark current observation for removal
                    indicesToRemove.push_back( i );
                }
            }
        }

        // Remove duplicates if any were found
        if( indicesToRemove.size( ) > 0 )
        {
            int beforeCount = numberOfObservations_;
            removeObservations( indicesToRemove );
            std::cerr << "[WARNING] Detected and removed " << beforeCount - numberOfObservations_
                      << "duplicate observations when creating instance of SingleObservationSet" << std::endl;
        }
    }

    void filterObservations( const std::shared_ptr< ObservationFilterBase > observationFilter, const bool saveFilteredObservations = true );

    void addObservations( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                          const std::vector< TimeType >& times,
                          const std::vector< Eigen::VectorXd >& dependentVariables = {},
                          const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights = {},
                          const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals = {},
                          const bool sortObservations = true );

    void addDependentVariables(
            const std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > dependentVariableSettings )
    {
        if( dependentVariableBookkeeping_ == nullptr )
        {
            dependentVariableBookkeeping_ =
                    std::make_shared< ObservationDependentVariableBookkeeping >( observableType_, linkEnds_.linkEnds_ );
        }
        else if( observationsDependentVariables_.size( ) != 0 )
        {
            throw std::runtime_error(
                    "Error, cannot add dependent variable settings to SingleObservationSet that has dependent variables calculated "
                    "already" );
        }
        dependentVariableBookkeeping_->addDependentVariables( dependentVariableSettings );
    }

    void clearDependentVariableValues( )
    {
        observationsDependentVariables_.clear( );
    }

private:
    //! Canonical weight state stored for the currently active observations.
    struct WeightState {
        //! Base matrix representation used for per-observation weights.
        ObservationWeightsMatrixType matrixType = diagonal_weights_matrix;
        //! Per-observation diagonal weight vectors (used when matrixType is diagonal).
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > diagonalWeights;
        //! Per-observation dense blocks (used when matrixType is block-diagonal).
        std::vector< Eigen::MatrixXd > blockWeights;
        //! Optional additional full contribution over the complete active observation set.
        Eigen::MatrixXd fullWeights;
    };

    //! Validates that the number of observations matches the number of epochs.
    void validateObservationsAndTimesSize( const std::size_t numberOfObservations,
                                           const std::size_t numberOfTimes,
                                           const std::string& context ) const
    {
        if( numberOfObservations != numberOfTimes )
        {
            throw std::runtime_error( "Error when " + context + ", input sizes are inconsistent." );
        }
    }

    //! Validates that an optional per-observation batch is either empty or has the expected size.
    template< typename DataType >
    void validateOptionalBatchSize( const std::vector< DataType >& data,
                                    const std::size_t expectedSize,
                                    const std::string& dataLabel,
                                    const std::string& context ) const
    {
        if( !data.empty( ) && data.size( ) != expectedSize )
        {
            throw std::runtime_error( "Error when " + context + ", number of " + dataLabel + " is inconsistent." );
        }
    }

    //! Validates that all observation vectors in the batch have the same dimension.
    void validateConsistentObservationDimensions(
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::string& context ) const
    {
        for( unsigned int i = 1; i < observations.size( ); i++ )
        {
            if( observations.at( i ).rows( ) != observations.at( i - 1 ).rows( ) )
            {
                throw std::runtime_error( "Error when " + context + ", input observables are not of consistent size." );
            }
        }
    }

    //! Validates that each observation vector matches the configured single-observation size.
    void validateObservationDimensionsAgainstSingleSize(
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::string& context ) const
    {
        for( unsigned int i = 0; i < observations.size( ); i++ )
        {
            if( observations.at( i ).size( ) != static_cast< int >( singleObservationSize_ ) )
            {
                throw std::runtime_error( "Error when " + context + ", observation size is inconsistent." );
            }
        }
    }

    //! Validates that each per-observation vector in a batch has the expected component size.
    template< typename ScalarType >
    void validatePerObservationVectorSize( const std::vector< Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > >& data,
                                           const std::string& dataLabel,
                                           const std::string& context ) const
    {
        for( unsigned int i = 0; i < data.size( ); i++ )
        {
            if( data.at( i ).size( ) != static_cast< int >( singleObservationSize_ ) )
            {
                throw std::runtime_error( "Error when " + context + ", " + dataLabel + " size is inconsistent." );
            }
        }
    }

    //! Validates cross-batch consistency for observations, times, residuals and dependent variables.
    void validateObservationBatchSizes( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                                        const std::vector< TimeType >& times,
                                        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals,
                                        const std::vector< Eigen::VectorXd >& dependentVariables,
                                        const std::string& context ) const
    {
        validateObservationsAndTimesSize( observations.size( ), times.size( ), context );
        validateOptionalBatchSize( residuals, observations.size( ), "residuals", context );
        validateOptionalBatchSize( dependentVariables, observations.size( ), "dependent variable entries", context );
    }

    //! Validates dependent-variable batch size and bookkeeping compatibility.
    void validateDependentVariableBatch( const std::vector< Eigen::VectorXd >& dependentVariables,
                                         const std::size_t expectedObservations,
                                         const bool requireBookkeeping,
                                         const std::string& context ) const
    {
        validateOptionalBatchSize( dependentVariables, expectedObservations, "dependent variable entries", context );
        if( dependentVariables.empty( ) )
        {
            return;
        }

        if( dependentVariableBookkeeping_ == nullptr )
        {
            if( requireBookkeeping )
            {
                throw std::runtime_error( "Error when " + context +
                                          ", dependent variable bookkeeping is required when dependent variables are provided." );
            }
            return;
        }

        const int expectedDependentVariableSize = dependentVariableBookkeeping_->getTotalDependentVariableSize( );
        for( unsigned int i = 0; i < dependentVariables.size( ); i++ )
        {
            if( dependentVariables.at( i ).size( ) != expectedDependentVariableSize )
            {
                throw std::runtime_error( "Error when " + context + ", dependent variable vector size is inconsistent." );
            }
        }
    }

    //! Resets base weights to unit diagonal and clears the optional full-matrix contribution.
    void resetWeightsToUnitDiagonal( )
    {
        setDiagonalWeightState(
                createUniformDiagonalWeights( Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( singleObservationSize_, 1 ) ) );
        weightState_.fullWeights.resize( 0, 0 );
    }

    //! Creates one identical diagonal-weight vector per currently active observation.
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > createUniformDiagonalWeights(
            const Eigen::Matrix< double, Eigen::Dynamic, 1 >& singleWeight ) const
    {
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > diagonalWeights;
        diagonalWeights.assign( numberOfObservations_, singleWeight );
        return diagonalWeights;
    }

    //! Stores diagonal base weights and clears block-diagonal storage.
    void setDiagonalWeightState( std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > diagonalWeights )
    {
        weightState_.matrixType = diagonal_weights_matrix;
        weightState_.diagonalWeights = std::move( diagonalWeights );
        weightState_.blockWeights.clear( );
    }

    //! Builds the permutation that sorts observation indices by ascending epoch.
    static std::vector< std::size_t > getTimeSortingPermutation( const std::vector< TimeType >& observationTimes )
    {
        const std::size_t numberOfObservations = observationTimes.size( );

        std::vector< std::size_t > permutation( numberOfObservations );
        for( std::size_t i = 0; i < numberOfObservations; ++i )
        {
            permutation[ i ] = i;
        }

        std::sort( permutation.begin( ), permutation.end( ), [ &observationTimes ]( const std::size_t i, const std::size_t j ) {
            return observationTimes.at( i ) < observationTimes.at( j );
        } );

        return permutation;
    }

    //! Applies a precomputed permutation to a vector in-place.
    template< typename T >
    void reorderVectorInPlace( std::vector< T >& data, const std::vector< std::size_t >& permutation )
    {
        const std::size_t numberOfElements = data.size( );

        std::vector< T > reorderedData( numberOfElements );
        for( std::size_t i = 0; i < numberOfElements; ++i )
        {
            reorderedData[ i ] = data.at( permutation.at( i ) );
        }

        data.swap( reorderedData );
    }

    //! Reorders a square matrix interpreted as observation-sized blocks.
    void reorderSquareMatrixInObservationBlocksInPlace( Eigen::MatrixXd& matrix, const std::vector< std::size_t >& permutation )
    {
        if( matrix.rows( ) == 0 )
        {
            return;
        }

        const int blockSize = static_cast< int >( singleObservationSize_ );
        const int numberOfObservations = static_cast< int >( permutation.size( ) );
        const int matrixSize = numberOfObservations * blockSize;
        if( matrix.rows( ) != matrixSize || matrix.cols( ) != matrixSize )
        {
            throw std::runtime_error(
                    "Error when reordering full weights matrix in SingleObservationSet, matrix size is incompatible with "
                    "observation metadata." );
        }

        Eigen::MatrixXd reorderedMatrix = Eigen::MatrixXd::Zero( matrixSize, matrixSize );
        for( int i = 0; i < numberOfObservations; i++ )
        {
            const int sourceI = static_cast< int >( permutation.at( i ) );
            for( int j = 0; j < numberOfObservations; j++ )
            {
                const int sourceJ = static_cast< int >( permutation.at( j ) );
                reorderedMatrix.block( i * blockSize, j * blockSize, blockSize, blockSize ) =
                        matrix.block( sourceI * blockSize, sourceJ * blockSize, blockSize, blockSize );
            }
        }
        matrix.swap( reorderedMatrix );
    }

    //! Sorts all active observation-aligned containers by ascending observation epoch.
    void orderObservationsAndMetadata( )
    {
        const std::size_t numberOfObservations = observationTimes_.size( );

        if( numberOfObservations < 2 )
        {
            return;
        }

        const std::vector< std::size_t > permutation = getTimeSortingPermutation( observationTimes_ );

        reorderVectorInPlace( observationTimes_, permutation );
        reorderVectorInPlace( observations_, permutation );
        if( weightState_.matrixType == diagonal_weights_matrix )
        {
            reorderVectorInPlace( weightState_.diagonalWeights, permutation );
        }
        else if( weightState_.matrixType == block_diagonal_weights_matrix )
        {
            reorderVectorInPlace( weightState_.blockWeights, permutation );
        }
        if( weightState_.fullWeights.rows( ) > 0 )
        {
            reorderSquareMatrixInObservationBlocksInPlace( weightState_.fullWeights, permutation );
        }
        reorderVectorInPlace( residuals_, permutation );

        if( !observationsDependentVariables_.empty( ) )
        {
            reorderVectorInPlace( observationsDependentVariables_, permutation );
        }
    }

    //! Updates cached minimum/maximum epochs for the active observation set.
    void updateTimeBounds( )
    {
        if( observationTimes_.size( ) == 0 )
        {
            timeBounds_ = std::make_pair( TUDAT_NAN, TUDAT_NAN );
        }
        else
        {
            timeBounds_ = std::make_pair( *std::min_element( observationTimes_.begin( ), observationTimes_.end( ) ),
                                          *std::max_element( observationTimes_.begin( ), observationTimes_.end( ) ) );
        }
    }

    //! Creates consecutive observation indices [0, ..., N-1].
    std::vector< unsigned int > createSequentialObservationIndices( const unsigned int numberOfObservations ) const
    {
        std::vector< unsigned int > indices( numberOfObservations );
        for( unsigned int i = 0; i < numberOfObservations; i++ )
        {
            indices.at( i ) = i;
        }
        return indices;
    }

    //! Creates consecutive observation indices with a constant offset.
    std::vector< unsigned int > createSequentialObservationIndicesWithOffset( const unsigned int numberOfObservations,
                                                                              const unsigned int offset ) const
    {
        std::vector< unsigned int > indices( numberOfObservations );
        for( unsigned int i = 0; i < numberOfObservations; i++ )
        {
            indices.at( i ) = i + offset;
        }
        return indices;
    }

    //! Applies a constant offset to each observation index in the input list.
    std::vector< unsigned int > applyOffsetToObservationIndices( const std::vector< unsigned int >& indices,
                                                                 const unsigned int offset ) const
    {
        std::vector< unsigned int > offsetIndices;
        offsetIndices.reserve( indices.size( ) );
        for( const unsigned int index : indices )
        {
            offsetIndices.push_back( index + offset );
        }
        return offsetIndices;
    }

    //! Concatenates two observation-index lists in-order.
    std::vector< unsigned int > concatenateObservationIndices( const std::vector< unsigned int >& firstIndices,
                                                               const std::vector< unsigned int >& secondIndices ) const
    {
        std::vector< unsigned int > concatenatedIndices;
        concatenatedIndices.reserve( firstIndices.size( ) + secondIndices.size( ) );
        concatenatedIndices.insert( concatenatedIndices.end( ), firstIndices.begin( ), firstIndices.end( ) );
        concatenatedIndices.insert( concatenatedIndices.end( ), secondIndices.begin( ), secondIndices.end( ) );
        return concatenatedIndices;
    }

    //! Applies a permutation to an index list and returns the permuted indices.
    std::vector< unsigned int > applyPermutationToObservationIndices( const std::vector< unsigned int >& indices,
                                                                      const std::vector< std::size_t >& permutation ) const
    {
        if( indices.size( ) != permutation.size( ) )
        {
            throw std::runtime_error( "Error when applying permutation to observation indices, sizes are inconsistent." );
        }

        std::vector< unsigned int > permutedIndices( indices.size( ) );
        for( unsigned int i = 0; i < indices.size( ); i++ )
        {
            permutedIndices.at( i ) = indices.at( permutation.at( i ) );
        }
        return permutedIndices;
    }

    //! Returns all observation indices not present in the selected set.
    std::vector< unsigned int > getComplementObservationIndices( const unsigned int totalNumberOfObservations,
                                                                 const std::vector< unsigned int >& selectedIndices ) const
    {
        std::vector< bool > isSelected( totalNumberOfObservations, false );
        for( const unsigned int index : selectedIndices )
        {
            if( index >= totalNumberOfObservations )
            {
                throw std::runtime_error( "Error when creating complement observation index set, input index is out of bounds." );
            }
            isSelected.at( index ) = true;
        }

        std::vector< unsigned int > complementIndices;
        complementIndices.reserve( totalNumberOfObservations - selectedIndices.size( ) );
        for( unsigned int i = 0; i < totalNumberOfObservations; i++ )
        {
            if( !isSelected.at( i ) )
            {
                complementIndices.push_back( i );
            }
        }
        return complementIndices;
    }

    //! Extracts a block matrix using observation indices as block selectors.
    Eigen::MatrixXd extractObservationBlockMatrix( const Eigen::MatrixXd& matrix,
                                                   const std::vector< unsigned int >& rowObservationIndices,
                                                   const std::vector< unsigned int >& columnObservationIndices ) const
    {
        const int blockSize = static_cast< int >( singleObservationSize_ );
        const int extractedRows = static_cast< int >( rowObservationIndices.size( ) ) * blockSize;
        const int extractedColumns = static_cast< int >( columnObservationIndices.size( ) ) * blockSize;
        Eigen::MatrixXd extractedMatrix = Eigen::MatrixXd::Zero( extractedRows, extractedColumns );
        if( extractedRows == 0 || extractedColumns == 0 || matrix.rows( ) == 0 || matrix.cols( ) == 0 )
        {
            return extractedMatrix;
        }

        for( unsigned int i = 0; i < rowObservationIndices.size( ); i++ )
        {
            const int sourceRow = static_cast< int >( rowObservationIndices.at( i ) ) * blockSize;
            for( unsigned int j = 0; j < columnObservationIndices.size( ); j++ )
            {
                const int sourceColumn = static_cast< int >( columnObservationIndices.at( j ) ) * blockSize;
                extractedMatrix.block( static_cast< int >( i ) * blockSize, static_cast< int >( j ) * blockSize, blockSize, blockSize ) =
                        matrix.block( sourceRow, sourceColumn, blockSize, blockSize );
            }
        }
        return extractedMatrix;
    }

    //! Returns the full-weight contribution for one state (or a size-compatible zero matrix if absent).
    Eigen::MatrixXd getEffectiveFullWeightContributionMatrix( const WeightState& weightState,
                                                              const unsigned int numberOfObservations ) const
    {
        const int totalObservationSize = static_cast< int >( numberOfObservations * singleObservationSize_ );
        if( totalObservationSize == 0 )
        {
            return Eigen::MatrixXd( 0, 0 );
        }
        if( weightState.fullWeights.rows( ) == 0 )
        {
            return Eigen::MatrixXd::Zero( totalObservationSize, totalObservationSize );
        }
        if( weightState.fullWeights.rows( ) != totalObservationSize || weightState.fullWeights.cols( ) != totalObservationSize )
        {
            throw std::runtime_error(
                    "Error when retrieving full weight matrix contribution, size is inconsistent with observation metadata." );
        }
        return weightState.fullWeights;
    }

    //! Returns the active-filtered cross-correlation block (or zeros if absent).
    Eigen::MatrixXd getEffectiveCrossWeightContributionMatrix( const unsigned int currentNumberOfObservations,
                                                               const unsigned int filteredNumberOfObservations ) const
    {
        const int currentObservationSize = static_cast< int >( currentNumberOfObservations * singleObservationSize_ );
        const int filteredObservationSize = static_cast< int >( filteredNumberOfObservations * singleObservationSize_ );
        if( currentObservationSize == 0 || filteredObservationSize == 0 )
        {
            return Eigen::MatrixXd::Zero( currentObservationSize, filteredObservationSize );
        }
        if( fullWeightCrossCorrelationWithFilteredSet_.rows( ) == 0 )
        {
            return Eigen::MatrixXd::Zero( currentObservationSize, filteredObservationSize );
        }
        if( fullWeightCrossCorrelationWithFilteredSet_.rows( ) != currentObservationSize ||
            fullWeightCrossCorrelationWithFilteredSet_.cols( ) != filteredObservationSize )
        {
            throw std::runtime_error(
                    "Error when retrieving cross weight matrix contribution with filtered set, size is inconsistent "
                    "with observation metadata." );
        }
        return fullWeightCrossCorrelationWithFilteredSet_;
    }

    //! Builds the combined full-weight matrix over [active; filtered] observations.
    Eigen::MatrixXd getCombinedFullWeightContributionMatrix( const unsigned int currentNumberOfObservations,
                                                             const unsigned int filteredNumberOfObservations ) const
    {
        const int currentObservationSize = static_cast< int >( currentNumberOfObservations * singleObservationSize_ );
        const int filteredObservationSize = static_cast< int >( filteredNumberOfObservations * singleObservationSize_ );
        Eigen::MatrixXd combinedFullWeights =
                Eigen::MatrixXd::Zero( currentObservationSize + filteredObservationSize, currentObservationSize + filteredObservationSize );

        const Eigen::MatrixXd currentFullWeights = getEffectiveFullWeightContributionMatrix( weightState_, currentNumberOfObservations );
        const Eigen::MatrixXd filteredFullWeights =
                getEffectiveFullWeightContributionMatrix( filteredObservationSet_->weightState_, filteredNumberOfObservations );
        const Eigen::MatrixXd crossWeights =
                getEffectiveCrossWeightContributionMatrix( currentNumberOfObservations, filteredNumberOfObservations );

        if( currentObservationSize > 0 )
        {
            combinedFullWeights.block( 0, 0, currentObservationSize, currentObservationSize ) = currentFullWeights;
        }
        if( filteredObservationSize > 0 )
        {
            combinedFullWeights.block( currentObservationSize, currentObservationSize, filteredObservationSize, filteredObservationSize ) =
                    filteredFullWeights;
        }
        if( currentObservationSize > 0 && filteredObservationSize > 0 )
        {
            combinedFullWeights.block( 0, currentObservationSize, currentObservationSize, filteredObservationSize ) = crossWeights;
            combinedFullWeights.block( currentObservationSize, 0, filteredObservationSize, currentObservationSize ) =
                    crossWeights.transpose( );
        }
        return combinedFullWeights;
    }

    //! Splits a combined [active; filtered] full-weight matrix back into internal storage blocks.
    void setCombinedFullWeightContributionMatrix( const Eigen::MatrixXd& combinedFullWeights,
                                                  const unsigned int currentNumberOfObservations,
                                                  const unsigned int filteredNumberOfObservations )
    {
        const int currentObservationSize = static_cast< int >( currentNumberOfObservations * singleObservationSize_ );
        const int filteredObservationSize = static_cast< int >( filteredNumberOfObservations * singleObservationSize_ );
        const int combinedObservationSize = currentObservationSize + filteredObservationSize;
        if( combinedFullWeights.rows( ) != combinedObservationSize || combinedFullWeights.cols( ) != combinedObservationSize )
        {
            throw std::runtime_error(
                    "Error when setting combined full weight contribution matrix, size is inconsistent with "
                    "current/filtered observation metadata." );
        }

        if( currentObservationSize > 0 )
        {
            weightState_.fullWeights =
                    compressWeightContributionMatrix( combinedFullWeights.block( 0, 0, currentObservationSize, currentObservationSize ) );
        }
        else
        {
            weightState_.fullWeights.resize( 0, 0 );
        }

        if( filteredObservationSize > 0 )
        {
            filteredObservationSet_->weightState_.fullWeights = compressWeightContributionMatrix( combinedFullWeights.block(
                    currentObservationSize, currentObservationSize, filteredObservationSize, filteredObservationSize ) );
        }
        else
        {
            filteredObservationSet_->weightState_.fullWeights.resize( 0, 0 );
        }

        if( currentObservationSize > 0 && filteredObservationSize > 0 )
        {
            fullWeightCrossCorrelationWithFilteredSet_ = compressWeightContributionMatrix(
                    combinedFullWeights.block( 0, currentObservationSize, currentObservationSize, filteredObservationSize ) );
        }
        else
        {
            fullWeightCrossCorrelationWithFilteredSet_.resize( 0, 0 );
        }

        validateCombinedWeightsMatrix( );
        filteredObservationSet_->validateCombinedWeightsMatrix( );
    }

    //! Stores an empty matrix when the contribution is numerically zero to keep sparse semantics.
    Eigen::MatrixXd compressWeightContributionMatrix( const Eigen::MatrixXd& matrix ) const
    {
        if( matrix.rows( ) == 0 || matrix.cols( ) == 0 || matrix.isZero( 0.0 ) )
        {
            return Eigen::MatrixXd( 0, 0 );
        }
        return matrix;
    }

    //! Validates that a weight matrix is square and symmetric.
    void validateSymmetricWeightMatrix( const Eigen::MatrixXd& weightMatrix, const std::string& errorPrefix ) const
    {
        if( weightMatrix.rows( ) != weightMatrix.cols( ) )
        {
            throw std::runtime_error( errorPrefix + " Weights matrix is not square." );
        }
        if( !weightMatrix.isApprox( weightMatrix.transpose( ), 1.0E-14 ) )
        {
            throw std::runtime_error( errorPrefix + " Weights matrix is not symmetric." );
        }
    }

    //! Validates that a weight matrix is symmetric positive definite.
    void validateWeightMatrixPositiveDefinite( const Eigen::MatrixXd& weightMatrix, const std::string& errorPrefix ) const
    {
        validateSymmetricWeightMatrix( weightMatrix, errorPrefix );
        Eigen::LLT< Eigen::MatrixXd > llt( weightMatrix );
        if( llt.info( ) != Eigen::Success )
        {
            throw std::runtime_error( errorPrefix + " Weights matrix is not positive definite." );
        }
    }

    //! Validates internal consistency and definiteness of the currently active combined weight model.
    void validateCombinedWeightsMatrix( ) const
    {
        if( numberOfObservations_ == 0 )
        {
            return;
        }

        if( weightState_.fullWeights.rows( ) == 0 )
        {
            if( weightState_.matrixType == diagonal_weights_matrix )
            {
                if( weightState_.diagonalWeights.size( ) != numberOfObservations_ )
                {
                    throw std::runtime_error(
                            "Error when setting weights in single observation set, number of diagonal weights is inconsistent." );
                }
                for( unsigned int i = 0; i < numberOfObservations_; i++ )
                {
                    if( weightState_.diagonalWeights.at( i ).rows( ) != static_cast< int >( singleObservationSize_ ) )
                    {
                        throw std::runtime_error(
                                "Error when setting weights in single observation set, diagonal weights size is inconsistent." );
                    }

                    for( unsigned int j = 0; j < singleObservationSize_; j++ )
                    {
                        const double currentWeight = weightState_.diagonalWeights.at( i )( j );
                        if( !std::isfinite( currentWeight ) || currentWeight <= 0.0 )
                        {
                            throw std::runtime_error(
                                    "Error when setting weights in single observation set, diagonal weights must be finite and strictly "
                                    "positive." );
                        }
                    }
                }
                return;
            }
            else if( weightState_.matrixType == block_diagonal_weights_matrix )
            {
                if( weightState_.blockWeights.size( ) != numberOfObservations_ )
                {
                    throw std::runtime_error(
                            "Error when setting weights in single observation set, number of block-diagonal weights is inconsistent." );
                }

                for( unsigned int i = 0; i < numberOfObservations_; i++ )
                {
                    validateWeightMatrixPositiveDefinite(
                            weightState_.blockWeights.at( i ),
                            "Error when setting weights in single observation set, block-diagonal weights are invalid." );
                }
                return;
            }
        }

        validateWeightMatrixPositiveDefinite( getWeightMatrix( ),
                                              "Error when setting weights in single observation set, combined weights matrix is invalid." );
    }

    //! Extracts the full-weight submatrix corresponding to the selected observation indices.
    Eigen::MatrixXd extractFullWeightMatrixSubset( const std::vector< unsigned int >& indices ) const;

    //! Extracts base/full weight data for a selected observation subset.
    WeightState extractWeightSubsetData( const std::vector< unsigned int >& indices ) const;

    //! Appends a batch of observations together with explicit pre-validated weight payload data.
    void addObservationsWithWeightData( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                                        const std::vector< TimeType >& times,
                                        const std::vector< Eigen::VectorXd >& dependentVariables,
                                        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals,
                                        const WeightState& weightSubsetData,
                                        const bool sortObservations );

    //! Moves observations between active and filtered sets while preserving weight/residual/dependent-variable consistency.
    void moveObservationsInOutFilteredSet( const std::vector< unsigned int >& indices,
                                           const bool moveInFilteredSet = true,
                                           const bool saveFilteredObservations = true );

    //! Function extracting the values of a single dependent variable
    Eigen::MatrixXd getSingleDependentVariable( std::pair< int, int > dependentVariableIndexAndSize ) const
    {
        Eigen::MatrixXd singleDependentVariable = Eigen::MatrixXd::Zero( numberOfObservations_, dependentVariableIndexAndSize.second );
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
                Eigen::VectorXd singleDependentVariableVector = observationsDependentVariables_.at( i ).segment(
                        dependentVariableIndexAndSize.first, dependentVariableIndexAndSize.second );
                singleDependentVariable.block( i, 0, 1, dependentVariableIndexAndSize.second ) = singleDependentVariableVector.transpose( );
            }
        }
        return singleDependentVariable;
    }

    //! Observable type shared by all observations in this set.
    const ObservableType observableType_;

    //! Link-end definition shared by all observations in this set.
    LinkDefinition linkEnds_;

    //! Cached minimum/maximum epoch over the active observations.
    std::pair< TimeType, TimeType > timeBounds_;

    //! Active observed values, one vector per observation epoch.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations_;

    //! Active observation epochs, aligned with observations_.
    std::vector< TimeType > observationTimes_;

    //! Reference link end used when interpreting the observable.
    const LinkEndType referenceLinkEnd_;

    //! Optional dependent-variable vectors, aligned with observations_.
    std::vector< Eigen::VectorXd > observationsDependentVariables_;

    //! Bookkeeping object describing dependent-variable layout per observation.
    std::shared_ptr< ObservationDependentVariableBookkeeping > dependentVariableBookkeeping_;

    //! Ancillary settings shared by all observations in this set.
    const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings_;

    //! Number of active observations currently stored.
    unsigned int numberOfObservations_;

    //! Observable size (number of components) per observation.
    unsigned int singleObservationSize_;

    //! Weight state for active observations only.
    WeightState weightState_;

    //    Eigen::VectorXd weightsVector_;

    //! Active residual vectors, aligned with observations_.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals_;

    //! Storage for observations currently filtered out of this set.
    std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > filteredObservationSet_;

    //! Full-weight cross block coupling active observations to filtered observations.
    Eigen::MatrixXd fullWeightCrossCorrelationWithFilteredSet_;
};

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
//! Returns a filtered copy of the input single observation set.
std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > filterObservations(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > singleObservationSet,
        const std::shared_ptr< ObservationFilterBase > observationFilter,
        const bool saveFilteredObservations = false );

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
//! Splits one single observation set into multiple sets according to a splitter configuration.
std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > splitObservationSet(
        const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > observationSet,
        const std::shared_ptr< ObservationSetSplitterBase > observationSetSplitter,
        const bool printWarning = true );

}  // namespace observation_models

}  // namespace tudat

#include "tudat/simulation/estimation_setup/singleObservationSetFilteringAndAppending.h"
#include "tudat/simulation/estimation_setup/singleObservationSetUtilities.h"

#endif  // TUDAT_SINGLE_OBSERVATION_SET_H
