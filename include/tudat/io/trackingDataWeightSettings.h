/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRACKING_DATA_WEIGHT_SETTINGS_H
#define TUDAT_TRACKING_DATA_WEIGHT_SETTINGS_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/io/observationWeightSettings.h"

namespace tudat
{

namespace data
{

//! Named methods to compute observation weights when creating an observation dataset
enum NamedWeightingMethod {
    no_named_weighting,  //!< No named weights (default)
    vfcc17               //!< Weights of Veres et al. (2017), for all observations that have no explicit weights
};

//! Definitions of blocks of observations for which the weights are reduced
enum ObservationBlockDefinition {
    local_night  //!< Observations taken in the same local night, from the same station
};

//! Rules to reduce the weights of the observations in a block
enum BlockDeweightingRule {
    divide_by_number_of_observations  //!< Weights divided by the number of observations N in the block
};

//! Settings to reduce the weights of all observations that belong to the same block (e.g. the same local night)
class BlockDeweightingSettings
{
public:
    BlockDeweightingSettings( const ObservationBlockDefinition blockDefinition, const BlockDeweightingRule rule ):
        blockDefinition_( blockDefinition ), rule_( rule )
    {}

    ObservationBlockDefinition blockDefinition_;

    BlockDeweightingRule rule_;
};

//! Base class for uncertainties of the observation epochs, that are converted to observation uncertainties when creating an
//! observation dataset (using the time derivatives of the observations)
class TimingBasedUncertaintySettings
{
public:
    TimingBasedUncertaintySettings( const std::vector< double >& timeUncertainties ): timeUncertainties_( timeUncertainties ) {}

    void removeObservation( const unsigned int observationIndex )
    {
        timeUncertainties_.erase( timeUncertainties_.begin( ) + observationIndex );
    }

    //! (1-sigma) uncertainties of the observation epochs (in seconds), one per observation
    std::vector< double > timeUncertainties_;
};

//! Function to check that settings that are only available for angular position observations are not used for other observables
inline void checkAngularPositionObservable( const std::string& observableType, const std::string& settingsName )
{
    if( observableType != "AngularPosition" )
    {
        throw std::runtime_error( "Error in TrackingData settings, " + settingsName +
                                  " are only available for angular position observations." );
    }
}

//! Settings that define the weights of the observations in a TrackingData object
class TrackingDataWeightSettings
{
public:
    //! Constructor, for the observations of a TrackingData object
    /*!
     * \param numberOfObservations Number of observations
     * \param singleObservationSize Size of a single observation
     * \param observableType Observable type of the observations (as in TrackingData)
     */
    TrackingDataWeightSettings( const unsigned int numberOfObservations,
                                const unsigned int singleObservationSize,
                                const std::string& observableType ):
        numberOfObservations_( numberOfObservations ), singleObservationSize_( singleObservationSize ), observableType_( observableType )
    {}

    unsigned int getNumberOfObservations( ) const
    {
        return numberOfObservations_;
    }

    unsigned int getSingleObservationSize( ) const
    {
        return singleObservationSize_;
    }

    const std::string& getObservableType( ) const
    {
        return observableType_;
    }

    //! Backwards-compatible interface for setting one diagonal weight vector per observation (stored as explicit weights).
    //! New code should use setObservationWeightSettings, which supports every weight representation.
    void setObservationWeights( const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& observationWeights )
    {
        // Check size consistency (for the total number of observations)
        if( observationWeights.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Error when adding observation weights to tracking data object, size of weights (" +
                                      std::to_string( observationWeights.size( ) ) + ") does not match number of observations (" +
                                      std::to_string( numberOfObservations_ ) + ")." );
        }

        // Check size consistency (for a single observable)
        for( auto weight : observationWeights )
        {
            if( weight.size( ) != singleObservationSize_ )
            {
                throw std::runtime_error( "Error when adding observation weights to tracking data object, size of single weight (" +
                                          std::to_string( weight.size( ) ) + ") does not match single observable size (should be " +
                                          std::to_string( singleObservationSize_ ) + ")." );
            }
        }

        const std::shared_ptr< observation_models::ObservationWeightSettings > diagonalWeights =
                observation_models::diagonalPerObservationWeights( observationWeights );
        checkWeightsConsistency( diagonalWeights );

        // Warn only after the complete replacement has passed validation. Rejected
        // replacements must leave the existing weights untouched.
        if( explicitWeights_ != nullptr )
        {
            std::cerr << "Warning when adding observation weights to tracking data object, weights already existed and are overwritten ."
                      << std::endl;
        }

        explicitWeights_ = diagonalWeights;
    }

    //! Set the explicit weights (known before creating an observation dataset), for all observations or for the observations given by
    //! the observation indices of the weight settings (in which case the other observations require named weights, see
    //! setNamedWeights).
    void setObservationWeightSettings( const std::shared_ptr< observation_models::ObservationWeightSettings > observationWeightSettings )
    {
        if( observationWeightSettings == nullptr )
        {
            throw std::runtime_error( "Error when setting TrackingData weights, no weight settings provided." );
        }
        observationWeightSettings->validateDimensions( numberOfObservations_, singleObservationSize_ );
        checkWeightsConsistency( observationWeightSettings );
        explicitWeights_ = observationWeightSettings;
    }

    //! Return the explicit weights of the observations (nullptr if none)
    std::shared_ptr< observation_models::ObservationWeightSettings > getObservationWeightSettings( ) const
    {
        return explicitWeights_;
    }

    //! Return whether explicit weights have been defined.
    bool hasObservationWeightSettings( ) const
    {
        return explicitWeights_ != nullptr;
    }

    //! Set weights computed with a named method when creating an observation dataset, for all observations without explicit weights
    void setNamedWeights( const NamedWeightingMethod method )
    {
        if( method != no_named_weighting )
        {
            checkAngularPositionObservable( observableType_, "named weights" );
        }
        namedWeightingMethod_ = method;
    }

    //! Return the named weighting method (no_named_weighting if none)
    NamedWeightingMethod getNamedWeightingMethod( ) const
    {
        return namedWeightingMethod_;
    }

    //! Set a reduction of the weights of all observations in the same block (e.g. local night), applied when creating an
    //! observation dataset
    void setBlockDeweighting( const ObservationBlockDefinition blockDefinition, const BlockDeweightingRule rule )
    {
        blockDeweighting_ = std::make_shared< BlockDeweightingSettings >( blockDefinition, rule );
    }

    //! Return the reduction of the weights per block (nullptr if none)
    std::shared_ptr< BlockDeweightingSettings > getBlockDeweighting( ) const
    {
        return blockDeweighting_;
    }

    //! Set the (1-sigma) uncertainties of the observation epochs (in seconds), one per observation, which are added to the observation
    //! covariance when creating an observation dataset (currently for angular position observations only)
    void setTimeUncertainties( const std::vector< double >& timeUncertainties )
    {
        if( timeUncertainties.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Error when setting observation time uncertainties, one entry per observation is required." );
        }
        for( const double timeUncertainty : timeUncertainties )
        {
            if( !std::isfinite( timeUncertainty ) || timeUncertainty < 0.0 )
            {
                throw std::runtime_error(
                        "Error when setting observation time uncertainties, values must be finite and non-negative, found " +
                        std::to_string( timeUncertainty ) + "." );
            }
        }
        timingUncertainties_ = std::make_shared< TimingBasedUncertaintySettings >( timeUncertainties );
    }

    //! Return the uncertainties of the observation epochs (empty if not provided)
    std::vector< double > getTimeUncertainties( ) const
    {
        return timingUncertainties_ != nullptr ? timingUncertainties_->timeUncertainties_ : std::vector< double >( );
    }

    //! Return the settings for the uncertainties of the observation epochs (nullptr if none)
    std::shared_ptr< TimingBasedUncertaintySettings > getTimingBasedUncertaintySettings( ) const
    {
        return timingUncertainties_;
    }

    //! Backwards-compatible interface for resetting one diagonal weight vector set through setObservationWeights.
    //! It cannot modify scalar, matrix-block, set-level or partial weight representations.
    void setSingleObservationWeight( const unsigned int index, const Eigen::Matrix< double, Eigen::Dynamic, 1 >& observationWeight )
    {
        // Check if weights are already available
        if( !hasDiagonalWeightsPerObservation( ) )
        {
            throw std::runtime_error(
                    "Error when resetting single observation weight in TrackingData object, per-observation diagonal weights not yet "
                    "defined." );
        }

        // Check that the observation index for which the weight needs resetting does not exceed the size of the observation vector
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when resetting single observation weight in TrackingData object, index exceeds number of observations contained "
                    "in TrackingData object." );
        }

        // Check size consistency
        if( observationWeight.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error( "Error when resetting single observation weight in TrackingData object, weight size (" +
                                      std::to_string( observationWeight.size( ) ) + ") inconsistent with observation size (" +
                                      std::to_string( singleObservationSize_ ) + ")." );
        }

        checkWeightsConsistency( observation_models::diagonalPerObservationWeights( std::vector< Eigen::VectorXd >( 1, observationWeight ),
                                                                                    std::vector< unsigned int >( 1, index ) ) );

        // Overwrites specific weight entry
        std::dynamic_pointer_cast< observation_models::DiagonalPerObservationWeightSettings >( explicitWeights_ )
                ->setWeight( index, observationWeight );
    }

    //! Backwards-compatible accessor for diagonal weight vectors set through setObservationWeights.
    //! Use getObservationWeightSettings for all representations; this function returns an empty vector for every other
    //! representation.
    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& getObservationWeights( ) const
    {
        if( hasDiagonalWeightsPerObservation( ) )
        {
            return std::dynamic_pointer_cast< observation_models::DiagonalPerObservationWeightSettings >( explicitWeights_ )->getWeights( );
        }
        static const std::vector< Eigen::VectorXd > emptyWeights;
        return emptyWeights;
    }

    //! Backwards-compatible accessor that concatenates diagonal weights set through setObservationWeights.
    //! Use getObservationWeightSettings for all representations; this function returns a size-zero vector for every other
    //! representation.
    Eigen::VectorXd getObservationWeightsVector( ) const
    {
        const std::vector< Eigen::VectorXd >& observationWeights = getObservationWeights( );
        Eigen::Matrix< double, Eigen::Dynamic, 1 > weightsVector =
                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( observationWeights.size( ) * singleObservationSize_, 1 );
        for( unsigned int i = 0; i < observationWeights.size( ); i++ )
        {
            weightsVector.block( i * singleObservationSize_, 0, singleObservationSize_, 1 ) = observationWeights.at( i );
        }
        return weightsVector;
    }

    //! Remove all weight settings
    void clear( )
    {
        explicitWeights_ = nullptr;
        namedWeightingMethod_ = no_named_weighting;
        blockDeweighting_ = nullptr;
        timingUncertainties_ = nullptr;
    }

    //! Remove a single observation from all settings
    void removeObservation( const unsigned int observationIndex )
    {
        if( explicitWeights_ != nullptr && !explicitWeights_->removeObservation( observationIndex, singleObservationSize_ ) )
        {
            explicitWeights_ = nullptr;
        }
        if( timingUncertainties_ != nullptr )
        {
            timingUncertainties_->removeObservation( observationIndex );
        }
        numberOfObservations_--;
    }

private:
    //! Check the consistency of weights before they are stored. For angular position observations, each observation must have a weight
    //! for both right ascension and declination, or for neither of the two
    void checkWeightsConsistency( const std::shared_ptr< observation_models::ObservationWeightSettings >& weightSettings ) const
    {
        if( observableType_ != "AngularPosition" || singleObservationSize_ != 2 )
        {
            return;
        }
        const std::vector< unsigned int >& observationIndices = weightSettings->getObservationIndices( );
        const std::shared_ptr< observation_models::SetBlockWeightSettings > setLevelWeights =
                std::dynamic_pointer_cast< observation_models::SetBlockWeightSettings >( weightSettings );
        const unsigned int numberOfWeights = observationIndices.empty( ) ? numberOfObservations_ : observationIndices.size( );
        for( unsigned int i = 0; i < numberOfWeights; ++i )
        {
            const Eigen::Vector2d weights = setLevelWeights != nullptr
                    ? Eigen::Vector2d( setLevelWeights->getWeightMatrix( ).diagonal( ).segment< 2 >( 2 * i ) )
                    : Eigen::Vector2d( weightSettings->getWeightBlock( i, 2 ).diagonal( ) );
            if( ( weights( 0 ) == 0.0 ) != ( weights( 1 ) == 0.0 ) )
            {
                throw std::runtime_error( "Error when setting TrackingData weights, observation " +
                                          std::to_string( observationIndices.empty( ) ? i : observationIndices.at( i ) ) +
                                          " has a weight for only one of right ascension and declination." );
            }
        }
    }

    //! Check whether the explicit weights are diagonal weights given for every observation
    bool hasDiagonalWeightsPerObservation( ) const
    {
        return explicitWeights_ != nullptr && explicitWeights_->getObservationIndices( ).empty( ) &&
                std::dynamic_pointer_cast< observation_models::DiagonalPerObservationWeightSettings >( explicitWeights_ ) != nullptr;
    }

    //! Number of observations
    unsigned int numberOfObservations_;

    //! Size of a single observation
    unsigned int singleObservationSize_;

    //! Observable type of the observations (as in TrackingData)
    std::string observableType_;

    //! Weights that are known before creating the observation dataset (nullptr if none)
    std::shared_ptr< observation_models::ObservationWeightSettings > explicitWeights_;

    //! Weights computed with a named method, for the observations without explicit weights (nullptr if none)
    NamedWeightingMethod namedWeightingMethod_ = no_named_weighting;

    //! Reduction of the weights of observations in the same block (nullptr if none)
    std::shared_ptr< BlockDeweightingSettings > blockDeweighting_;

    //! Uncertainties of the observation epochs, added to the observation covariance (nullptr if none)
    std::shared_ptr< TimingBasedUncertaintySettings > timingUncertainties_;
};

}  // namespace data

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_WEIGHT_SETTINGS_H
