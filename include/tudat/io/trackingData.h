/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRACKING_DATA_H
#define TUDAT_TRACKING_DATA_H

#include <Eigen/Core>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"
#include "tudat/io/observationWeightSettings.h"

namespace tudat
{

namespace data
{

using PlainLinkDefinition = std::vector< std::pair< std::pair< std::string, std::string >, std::string > >;

inline bool isOpticalObservationMetadata( const std::string& key )
{
    return key == "note2" || key == "catalog" || key == "band" || key == "phottype" || key == "custom_name" || key == "mag" ||
            key == "discovery" || key == "number";
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class TrackingData
{
public:
    TrackingData( const std::string observableType,
                  const PlainLinkDefinition& linkEnds,
                  const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                  const std::vector< TimeType > epochs,
                  const std::string referenceLinkEnd,
                  const std::string timeScale = "TDB",
                  const std::string weighingScheme = "" ):
        observableType_( observableType ), linkEnds_( linkEnds ), observations_( observations ), epochs_( epochs ),
        referenceLinkEnd_( referenceLinkEnd ), timeScale_( timeScale ), weighingScheme_( weighingScheme ),
        numberOfObservations_( observations.size( ) ),
        singleObservationSize_( !observations.empty( ) ? static_cast< unsigned int >( observations[ 0 ].size( ) ) : 0 )
    {
        // Check inputs size consistency
        if( observations_.size( ) != epochs_.size( ) )
        {
            throw std::runtime_error( "Error when creating TrackingData object, numbers of epochs (" + std::to_string( epochs_.size( ) ) +
                                      ") and observations (" + std::to_string( observations_.size( ) ) + ") are inconsistent." );
        }

        for( unsigned int i = 1; i < observations.size( ); i++ )
        {
            if( observations.at( i ).rows( ) != observations.at( i - 1 ).rows( ) )
            {
                throw std::runtime_error( "Error when making TrackingData, input observables not of consistent size." );
            }
        }
    }

    //! Function that returns the number of observations
    unsigned int getNumberOfObservations( ) const
    {
        return numberOfObservations_;
    }

    //! Function that returns the size of a single observation
    unsigned int getSingleObservationSize( ) const
    {
        return singleObservationSize_;
    }

    //! Function that returns the total size of the observations contained in the TrackingData object
    unsigned int getTotalObservationSetSize( ) const
    {
        return numberOfObservations_ * singleObservationSize_;
    }

    //! Function that returns the observable type
    const std::string& getObservableType( ) const
    {
        return observableType_;
    }

    //! Function that returns the link ends description
    const PlainLinkDefinition& getLinkEnds( ) const
    {
        return linkEnds_;
    }

    //! Function that returns the link end used as reference
    const std::string& getReferenceLinkEnd( ) const
    {
        return referenceLinkEnd_;
    }

    std::string getReferencePointName( ) const
    {
        for( auto const& [ linkEndPair, referencePoint ] : linkEnds_ )
        {
            if( referencePoint == referenceLinkEnd_ )
            {
                return linkEndPair.second;
            }
        }
        throw std::runtime_error( "Reference point name not found." );
    }

    //! Function that returns the time scale of the observation epochs
    const std::string& getTimeScale( ) const
    {
        return timeScale_;
    }

    //! Function that returns the requested weighing scheme
    const std::string& getWeighingScheme( ) const
    {
        return weighingScheme_;
    }

    //! Function that returns observation values
    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& getObservations( ) const
    {
        return observations_;
    }

    //! Function that returns a concatenated vector of observations
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationsVector( ) const
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_ * numberOfObservations_, 1 );
        for( unsigned int i = 0; i < observations_.size( ); i++ )
        {
            observationsVector.segment( i * singleObservationSize_, singleObservationSize_ ) = observations_.at( i );
        }
        return observationsVector;
    }

    //! Function that overwrites the full observation vector
    void resetObservations( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& newObservations )
    {
        // Check that the number of observations is consistent
        if( newObservations.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when resetting the observation values in TrackingData object, the new number of observations (" +
                    std::to_string( newObservations.size( ) ) + ") is inconsistent with original number of observations (" +
                    std::to_string( numberOfObservations_ ) + ")." );
        }

        // Check that the size of each single observation is consistent
        for( auto obs : newObservations )
        {
            if( obs.size( ) != singleObservationSize_ )
            {
                throw std::runtime_error(
                        "Error when resetting the observation values in TrackingData object, the size of a single new observation (" +
                        std::to_string( newObservations.size( ) ) + ") is inconsistent with the original single observation size (" +
                        std::to_string( singleObservationSize_ ) + ")." );
            }
        }

        // Overwrite observation values
        observations_ = newObservations;
    }

    //! Function that overwrites a single observation value
    void resetObservation( const unsigned int index, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observation )
    {
        // Check that the index of the observation that needs overwriting does not exceed the size of the observation vector
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when resetting single observation value in TrackingData object, index exceeds number of observations contained "
                    "in TrackingData object." );
        }

        // Check size consistency
        if( observation.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error(
                    "Error when resetting single observation value in TrackingData object, the observation size is inconsistent." );
        }

        // Overwrites specific observation entry
        observations_.at( index ) = observation;
    }

    //! Function that returns the observation epochs
    const std::vector< TimeType >& getObservationEpochs( ) const
    {
        return epochs_;
    }

    //! Function that returns a concatenated vector of observation epochs
    Eigen::Matrix< TimeType, Eigen::Dynamic, 1 > getObservationEpochsVector( ) const
    {
        Eigen::Matrix< TimeType, Eigen::Dynamic, 1 > epochsVector =
                Eigen::Matrix< TimeType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_ * numberOfObservations_, 1 );
        for( unsigned int i = 0; i < epochs_.size( ); i++ )
        {
            epochsVector.segment( i * singleObservationSize_, singleObservationSize_ ).setConstant( epochs_.at( i ) );
        }
        return epochsVector;
    }

    // Add ancillary settings (string type)
    void addAncillarySettings( const std::string ancillarySettingsType, const std::string ancillarySettingsValue )
    {
        ancillarySettingsString_[ ancillarySettingsType ] = ancillarySettingsValue;
    }

    //! Function that returns map of ancillary settings (string type)
    const std::map< std::string, std::string >& getAncillarySettingsString( ) const
    {
        return ancillarySettingsString_;
    }

    // Add ancillary settings (string type)
    void addAncillarySettings( const std::string ancillarySettingsType, const std::vector< std::string > ancillarySettingsValue )
    {
        // Retain the original optical metadata entry point, with explicit row semantics.
        if( isOpticalObservationMetadata( ancillarySettingsType ) || isObservationMetadata( ancillarySettingsType ) )
        {
            addObservationMetadata( ancillarySettingsType, ancillarySettingsValue );
            return;
        }
        ancillarySettingsStringVector_[ ancillarySettingsType ] = ancillarySettingsValue;
    }

    //! Register one string value per observation; link-level ancillary vectors use addAncillarySettings.
    void addObservationMetadata( const std::string& key, const std::vector< std::string >& values )
    {
        if( values.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Observation metadata '" + key + "' must have one entry per observation." );
        }
        ancillarySettingsStringVector_[ key ] = values;
        observationMetadataKeys_.insert( key );
    }

    bool isObservationMetadata( const std::string& key ) const
    {
        return observationMetadataKeys_.count( key ) != 0;
    }

    //! Function that returns map of ancillary settings (string type)
    const std::map< std::string, std::vector< std::string > >& getAncillarySettingsStringVector( ) const
    {
        return ancillarySettingsStringVector_;
    }

    // Add ancillary settings (double type)
    void addAncillarySettings( const std::string ancillarySettingsType, const double ancillarySettingsValue )
    {
        ancillarySettingsDouble_[ ancillarySettingsType ] = ancillarySettingsValue;
    }

    //! Function that returns map of ancillary settings (double type)
    const std::map< std::string, double >& getAncillarySettingsDouble( ) const
    {
        return ancillarySettingsDouble_;
    }

    // Add ancillary settings (vector type)
    void addAncillarySettings( const std::string ancillarySettingsType, const std::vector< double > ancillarySettingsValue )
    {
        ancillarySettingsDoubleVector_[ ancillarySettingsType ] = ancillarySettingsValue;
    }

    //! Function that returns map of ancillary settings (vector type)
    const std::map< std::string, std::vector< double > >& getAncillarySettingsDoubleVector( ) const
    {
        return ancillarySettingsDoubleVector_;
    }

    //! Set one diagonal weight vector for each observation (optional).
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

        // Warn only after the complete replacement has passed validation. Rejected
        // replacements must leave the existing weights untouched.
        if( observationWeightSettings_.has_value( ) )
        {
            std::cerr << "Warning when adding observation weights to tracking data object, weights already existed and are overwritten ."
                      << std::endl;
        }

        // Preserve the existing diagonal-weight interface while storing all
        // weight definitions through ObservationWeightSettings.
        observationWeightSettings_ = observation_models::ObservationWeightSettings::diagonalPerObservation( observationWeights );
    }

    //! Set the representation and values used to weight these observations.
    void setObservationWeightSettings( const observation_models::ObservationWeightSettings& observationWeightSettings )
    {
        observationWeightSettings_ = observationWeightSettings;
    }

    //! Return the settings that define the observation weights.
    const observation_models::ObservationWeightSettings& getObservationWeightSettings( ) const
    {
        if( !observationWeightSettings_.has_value( ) )
        {
            throw std::runtime_error( "Error when retrieving TrackingData weight settings, no weights have been defined." );
        }
        return observationWeightSettings_.value( );
    }

    //! Return whether observation weights have been defined.
    bool hasObservationWeightSettings( ) const
    {
        return observationWeightSettings_.has_value( );
    }

    //! Remove any observation weights currently stored.
    void clearObservationWeightSettings( )
    {
        observationWeightSettings_.reset( );
    }

    //! Reset one stored per-observation diagonal weight vector.
    void setSingleObservationWeight( const unsigned int index, const Eigen::Matrix< double, Eigen::Dynamic, 1 >& observationWeight )
    {
        // Check if weights are already available
        if( !observationWeightSettings_.has_value( ) ||
            observationWeightSettings_->type_ != observation_models::ObservationWeightSettings::WeightsBlockType::diagonal_per_observation )
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

        // Overwrites specific weight entry
        observationWeightSettings_->diagonalWeights_.at( index ) = observationWeight;
    }

    //! Return stored per-observation diagonal weight vectors, or an empty vector for any other weight representation.
    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& getObservationWeights( ) const
    {
        if( observationWeightSettings_.has_value( ) &&
            observationWeightSettings_->type_ == observation_models::ObservationWeightSettings::WeightsBlockType::diagonal_per_observation )
        {
            return observationWeightSettings_->diagonalWeights_;
        }
        static const std::vector< Eigen::VectorXd > emptyWeights;
        return emptyWeights;
    }

    //! Return the concatenated per-observation diagonal weights, or a size-zero vector for any other representation.
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

    //! Set corrections to the observations (optional)
    void setObservationCorrections( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observationCorrections )
    {
        // Check size consistency (for the total number of observations)
        if( observationCorrections.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Error when adding observation corrections to tracking data object, size of corrections (" +
                                      std::to_string( observationCorrections.size( ) ) + ") does not match number of observations (" +
                                      std::to_string( numberOfObservations_ ) + ")." );
        }

        // Check size consistency (for a single observable)
        for( auto correction : observationCorrections )
        {
            if( correction.size( ) != singleObservationSize_ )
            {
                throw std::runtime_error( "Error when adding observation corrections to tracking data object, size of single correction (" +
                                          std::to_string( correction.size( ) ) + ") does not match single observable size (should be " +
                                          std::to_string( singleObservationSize_ ) + ")." );
            }
        }

        // Warn only after the complete replacement has passed validation. Rejected
        // replacements must leave the existing corrections untouched.
        if( !observationCorrections_.empty( ) )
        {
            std::cerr << "Warning when adding observation corrections to tracking data object, corrections already existed and are "
                         "overwritten ."
                      << std::endl;
        }

        // If all sizes are consistent, store observation corrections
        observationCorrections_ = observationCorrections;
    }

    //! Function that (re-)sets a single observation correction
    void setSingleObservationCorrection( const unsigned int index,
                                         Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationCorrection )
    {
        // Check that observation corrections are already defined
        if( observationCorrections_.empty( ) )
        {
            throw std::runtime_error(
                    "Error when resetting single observation correction in TrackingData object, corrections not yet defined." );
        }

        // Check that the index of the observation correction that needs overwriting does not exceed the size of the observation vector
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when resetting single observation correction in TrackingData object, index exceeds number of observations "
                    "contained in TrackingData object." );
        }

        // Check size consistency
        if( observationCorrection.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error( "Error when resetting single observation correction in TrackingData object, new correction size (" +
                                      std::to_string( observationCorrection.size( ) ) + ") inconsistent with observation size (" +
                                      std::to_string( singleObservationSize_ ) + ")." );
        }

        // Overwrites specific correction entry
        observationCorrections_.at( index ) = observationCorrection;
    }

    //! Function that return a vector of observation corrections (empty of observation corrections are not provided)
    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& getObservationCorrections( ) const
    {
        return observationCorrections_;
    }

    //! Function that remove a single observation entry
    void removeSingleObservationEntry( const unsigned int index )
    {
        // Check that the index of the observation that needs deleting does not exceed the size of the observation vector
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when removing single observation in TrackingData object, index exceeds number of observations." );
        }

        // Remove observation value and associated epoch
        observations_.erase( observations_.begin( ) + index );
        epochs_.erase( epochs_.begin( ) + index );

        // Keep every per-observation metadata vector aligned with the data.
        for( const auto& key : observationMetadataKeys_ )
        {
            auto& values = ancillarySettingsStringVector_.at( key );
            values.erase( values.begin( ) + index );
        }

        // Keep any per-observation weight values aligned with the remaining observations.
        if( observationWeightSettings_.has_value( ) )
        {
            using WeightsBlockType = observation_models::ObservationWeightSettings::WeightsBlockType;
            switch( observationWeightSettings_->type_ )
            {
                case WeightsBlockType::scalar_per_observation:
                    observationWeightSettings_->scalarWeights_.erase( observationWeightSettings_->scalarWeights_.begin( ) + index );
                    break;
                case WeightsBlockType::diagonal_per_observation:
                    observationWeightSettings_->diagonalWeights_.erase( observationWeightSettings_->diagonalWeights_.begin( ) + index );
                    break;
                case WeightsBlockType::block_per_observation:
                    observationWeightSettings_->weightBlocks_.erase( observationWeightSettings_->weightBlocks_.begin( ) + index );
                    break;
                case WeightsBlockType::set_block: {
                    // A set-level matrix is ordered by observation and then by
                    // scalar component. Remove the rows and columns belonging
                    // to the deleted observation.
                    const Eigen::MatrixXd originalMatrix = observationWeightSettings_->weightBlock_;
                    const Eigen::Index removedScalarStart = static_cast< Eigen::Index >( index * singleObservationSize_ );
                    Eigen::MatrixXd reducedMatrix( originalMatrix.rows( ) - singleObservationSize_,
                                                   originalMatrix.cols( ) - singleObservationSize_ );
                    Eigen::Index newRow = 0;
                    for( Eigen::Index oldRow = 0; oldRow < originalMatrix.rows( ); ++oldRow )
                    {
                        if( oldRow >= removedScalarStart && oldRow < removedScalarStart + singleObservationSize_ )
                        {
                            continue;
                        }
                        Eigen::Index newColumn = 0;
                        for( Eigen::Index oldColumn = 0; oldColumn < originalMatrix.cols( ); ++oldColumn )
                        {
                            if( oldColumn >= removedScalarStart && oldColumn < removedScalarStart + singleObservationSize_ )
                            {
                                continue;
                            }
                            reducedMatrix( newRow, newColumn++ ) = originalMatrix( oldRow, oldColumn );
                        }
                        ++newRow;
                    }
                    observationWeightSettings_->weightBlock_ = reducedMatrix;
                    break;
                }
                case WeightsBlockType::default_weights:
                case WeightsBlockType::constant_scalar:
                case WeightsBlockType::constant_block:
                    // These representations do not store values separately for each observation.
                    break;
                default: {
                    throw std::runtime_error( "Error when removing TrackingData observation, weight type is not recognized." );
                }
            }
        }

        // Remove associated correction (if it exists)
        if( !observationCorrections_.empty( ) )
        {
            observationCorrections_.erase( observationCorrections_.begin( ) + index );
        }

        // Update total number of observations contained in TrackingData object
        numberOfObservations_ = observations_.size( );
    }

private:
    const std::string observableType_;

    const PlainLinkDefinition linkEnds_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations_;

    std::vector< TimeType > epochs_;

    const std::string referenceLinkEnd_;

    const std::string timeScale_;

    const std::string weighingScheme_;

    unsigned int numberOfObservations_;

    const unsigned int singleObservationSize_;

    std::optional< observation_models::ObservationWeightSettings > observationWeightSettings_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observationCorrections_;

    std::map< std::string, std::string > ancillarySettingsString_;

    std::map< std::string, std::vector< std::string > > ancillarySettingsStringVector_;

    std::set< std::string > observationMetadataKeys_;

    std::map< std::string, double > ancillarySettingsDouble_;

    std::map< std::string, std::vector< double > > ancillarySettingsDoubleVector_;
};

}  // namespace data

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_H
