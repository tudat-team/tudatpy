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
#include <cmath>
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
#include "tudat/io/trackingDataSettings.h"

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
                  const std::string timeScale = "TDB" ):
        observableType_( observableType ), linkEnds_( linkEnds ), observations_( observations ), epochs_( epochs ),
        referenceLinkEnd_( referenceLinkEnd ), timeScale_( timeScale ), numberOfObservations_( observations.size( ) ),
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

    //! Attach one numerical value per observation under a user-defined key.
    void addNumericalObservationMetadata( const std::string& key, const std::vector< double >& values )
    {
        if( values.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Numerical observation metadata '" + key + "' must have one entry per observation." );
        }
        numericalObservationMetadata_[ key ] = values;
    }

    const std::map< std::string, std::vector< double > >& getNumericalObservationMetadata( ) const
    {
        return numericalObservationMetadata_;
    }

    const std::vector< double >& getNumericalObservationMetadata( const std::string& key ) const
    {
        return numericalObservationMetadata_.at( key );
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

    //! Function that returns the settings that define the observation weights
    const std::shared_ptr< TrackingDataWeightSettings >& getWeightSettings( ) const
    {
        return weightSettings_;
    }

    //! Function that returns the settings that define the observation corrections
    const std::shared_ptr< TrackingDataCorrectionSettings< ObservationScalarType > >& getCorrectionSettings( ) const
    {
        return correctionSettings_;
    }

    //! Backwards-compatible interface for setting one diagonal weight vector per observation (stored as positional weights).
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

        const observation_models::ObservationWeightSettings diagonalWeights =
                observation_models::ObservationWeightSettings::diagonalPerObservation( observationWeights );
        checkAngularPositionWeights( diagonalWeights );

        // Warn only after the complete replacement has passed validation. Rejected
        // replacements must leave the existing weights untouched.
        if( weightSettings_->positionalWeights_ != nullptr )
        {
            std::cerr << "Warning when adding observation weights to tracking data object, weights already existed and are overwritten ."
                      << std::endl;
        }

        weightSettings_->positionalWeights_ = std::make_shared< PositionalWeightSettings >( diagonalWeights );
    }

    //! Set the positional weights (known before creating an observation dataset), for all observations or for the observations
    //! given by observationIndices (in which case the other observations require named weights, see setNamedWeights).
    void setObservationWeightSettings( const observation_models::ObservationWeightSettings& observationWeightSettings,
                                       const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) )
    {
        if( !observationIndices.empty( ) )
        {
            if( observationWeightSettings.type_ == observation_models::ObservationWeightSettings::WeightsBlockType::set_block )
            {
                throw std::runtime_error(
                        "Error when setting TrackingData weights, set-level weights must be given for all observations." );
            }
            const std::set< unsigned int > uniqueIndices( observationIndices.begin( ), observationIndices.end( ) );
            if( uniqueIndices.size( ) != observationIndices.size( ) || *uniqueIndices.rbegin( ) >= numberOfObservations_ )
            {
                throw std::runtime_error(
                        "Error when setting TrackingData weights, observation indices must be unique and smaller than the number of "
                        "observations (" +
                        std::to_string( numberOfObservations_ ) + ")." );
            }
        }
        observationWeightSettings.validateDimensions( observationIndices.empty( ) ? numberOfObservations_ : observationIndices.size( ),
                                                      singleObservationSize_ );
        checkAngularPositionWeights( observationWeightSettings, observationIndices );
        weightSettings_->positionalWeights_ = std::make_shared< PositionalWeightSettings >( observationWeightSettings, observationIndices );
    }

    //! Return the positional weights of the observations (for all observations, or for those given by the indices of the positional
    //! weight settings, see getWeightSettings).
    const observation_models::ObservationWeightSettings& getObservationWeightSettings( ) const
    {
        if( weightSettings_->positionalWeights_ == nullptr )
        {
            throw std::runtime_error( "Error when retrieving TrackingData weight settings, no positional weights have been defined." );
        }
        return weightSettings_->positionalWeights_->weights_;
    }

    //! Return whether positional weights have been defined.
    bool hasObservationWeightSettings( ) const
    {
        return weightSettings_->positionalWeights_ != nullptr;
    }

    //! Remove all weight settings currently stored.
    void clearObservationWeightSettings( )
    {
        weightSettings_ = std::make_shared< TrackingDataWeightSettings >( );
    }

    //! Set weights computed with a named method when creating an observation dataset, for all observations without positional weights
    void setNamedWeights( const NamedWeightingMethod method )
    {
        checkAngularPositionObservable( "named weights" );
        weightSettings_->namedWeights_ = std::make_shared< NamedWeightSettings >( method );
    }

    //! Set a reduction of the weights of all observations in the same block (e.g. local night), applied when creating an
    //! observation dataset
    void setBlockDeweighting( const ObservationBlockDefinition blockDefinition, const BlockDeweightingRule rule )
    {
        weightSettings_->blockDeweighting_ = std::make_shared< BlockDeweightingSettings >( blockDefinition, rule );
    }

    //! Set the uncertainties of the observation epochs, which are added to the observation covariance when creating an observation
    //! dataset
    void setTimingUncertainties( const std::shared_ptr< TimingBasedUncertaintySettings > timingUncertainties )
    {
        if( timingUncertainties == nullptr )
        {
            throw std::runtime_error( "Error when setting observation time uncertainties, no settings provided." );
        }
        if( timingUncertainties->timeUncertainties_.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Error when setting observation time uncertainties, one entry per observation is required." );
        }
        for( const double timeUncertainty : timingUncertainties->timeUncertainties_ )
        {
            if( !std::isfinite( timeUncertainty ) || timeUncertainty < 0.0 )
            {
                throw std::runtime_error(
                        "Error when setting observation time uncertainties, values must be finite and non-negative, found " +
                        std::to_string( timeUncertainty ) + "." );
            }
        }
        if( std::dynamic_pointer_cast< AstrometricTimeUncertaintySettings >( timingUncertainties ) != nullptr )
        {
            checkAngularPositionObservable( "astrometric time uncertainties" );
        }
        weightSettings_->timingUncertainties_ = timingUncertainties;
    }

    //! Set the (1-sigma) uncertainties of the epochs of angular position observations (in seconds), one per observation
    void setAstrometricTimeUncertainties( const std::vector< double >& timeUncertainties )
    {
        setTimingUncertainties( std::make_shared< AstrometricTimeUncertaintySettings >( timeUncertainties ) );
    }

    //! Return the uncertainties of the epochs of angular position observations (empty if not provided)
    std::vector< double > getAstrometricTimeUncertainties( ) const
    {
        std::shared_ptr< AstrometricTimeUncertaintySettings > timingUncertainties =
                std::dynamic_pointer_cast< AstrometricTimeUncertaintySettings >( weightSettings_->timingUncertainties_ );
        return timingUncertainties != nullptr ? timingUncertainties->timeUncertainties_ : std::vector< double >( );
    }

    //! Set the biases of the epochs of angular position observations (in seconds), one per observation. An observation with time
    //! bias c is taken to be h(t+c), with t its epoch. When creating an observation dataset with corrections applied, each
    //! observation is corrected to h(t) using its time derivative.
    void setAstrometricTimeBiases( const std::vector< double >& timeBiases )
    {
        if( timeBiases.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Error when setting observation time biases, one entry per observation is required." );
        }
        for( const double timeBias : timeBiases )
        {
            if( !std::isfinite( timeBias ) )
            {
                throw std::runtime_error( "Error when setting observation time biases, values must be finite, found " +
                                          std::to_string( timeBias ) + "." );
            }
        }
        checkAngularPositionObservable( "astrometric time biases" );
        correctionSettings_->timeBiases_ = std::make_shared< AstrometricTimeBiasSettings >( timeBiases );
    }

    //! Return the biases of the epochs of angular position observations (empty if not provided)
    std::vector< double > getAstrometricTimeBiases( ) const
    {
        std::shared_ptr< AstrometricTimeBiasSettings > timeBiases =
                std::dynamic_pointer_cast< AstrometricTimeBiasSettings >( correctionSettings_->timeBiases_ );
        return timeBiases != nullptr ? timeBiases->timeBiases_ : std::vector< double >( );
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

        checkAngularPositionWeights( observation_models::ObservationWeightSettings::diagonalPerObservation(
                                             std::vector< Eigen::VectorXd >( 1, observationWeight ) ),
                                     std::vector< unsigned int >( 1, index ) );

        // Overwrites specific weight entry
        weightSettings_->positionalWeights_->weights_.diagonalWeights_.at( index ) = observationWeight;
    }

    //! Backwards-compatible accessor for diagonal weight vectors set through setObservationWeights.
    //! Use getWeightSettings for all representations; this function returns an empty vector for every other representation.
    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& getObservationWeights( ) const
    {
        if( hasDiagonalWeightsPerObservation( ) )
        {
            return weightSettings_->positionalWeights_->weights_.diagonalWeights_;
        }
        static const std::vector< Eigen::VectorXd > emptyWeights;
        return emptyWeights;
    }

    //! Backwards-compatible accessor that concatenates diagonal weights set through setObservationWeights.
    //! Use getWeightSettings for all representations; this function returns a size-zero vector for every other representation.
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

    //! Store corrections to be evaluated when constructing an observation dataset.
    void setObservationCorrectionSettings( const AngularObservationCorrectionSettings& settings )
    {
        if( observableType_ != "AngularPosition" || singleObservationSize_ != 2 || referenceLinkEnd_ != "receiver" )
        {
            throw std::runtime_error( "Angular correction settings require two-component reception-time angular observations." );
        }
        correctionSettings_->modelCorrections_ = std::make_shared< AngularObservationCorrectionSettings >( settings );
    }

    //! Return the angular correction settings (nullptr if not provided)
    std::shared_ptr< AngularObservationCorrectionSettings > getObservationCorrectionSettings( ) const
    {
        return std::dynamic_pointer_cast< AngularObservationCorrectionSettings >( correctionSettings_->modelCorrections_ );
    }

    bool hasObservationCorrectionSettings( ) const
    {
        const std::shared_ptr< AngularObservationCorrectionSettings > settings = getObservationCorrectionSettings( );
        return settings != nullptr && settings->hasCorrections( );
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
        if( correctionSettings_->explicitCorrections_ != nullptr )
        {
            std::cerr << "Warning when adding observation corrections to tracking data object, corrections already existed and are "
                         "overwritten ."
                      << std::endl;
        }

        // If all sizes are consistent, store observation corrections
        correctionSettings_->explicitCorrections_ =
                std::make_shared< ExplicitCorrectionSettings< ObservationScalarType > >( observationCorrections );
    }

    //! Function that (re-)sets a single observation correction
    void setSingleObservationCorrection( const unsigned int index,
                                         Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationCorrection )
    {
        // Check that observation corrections are already defined
        if( correctionSettings_->explicitCorrections_ == nullptr )
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
        correctionSettings_->explicitCorrections_->corrections_.at( index ) = observationCorrection;
    }

    //! Function that return a vector of observation corrections (empty of observation corrections are not provided)
    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& getObservationCorrections( ) const
    {
        if( correctionSettings_->explicitCorrections_ != nullptr )
        {
            return correctionSettings_->explicitCorrections_->corrections_;
        }
        static const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > emptyCorrections;
        return emptyCorrections;
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
        for( auto& metadata : numericalObservationMetadata_ )
        {
            metadata.second.erase( metadata.second.begin( ) + index );
        }

        // Keep the per-observation weight and correction settings aligned with the remaining observations.
        weightSettings_->removeObservation( index, singleObservationSize_ );
        correctionSettings_->removeObservation( index );

        // Update total number of observations contained in TrackingData object
        numberOfObservations_ = observations_.size( );
    }

private:
    //! Check that the observations are angular positions, for settings that are only available for these
    void checkAngularPositionObservable( const std::string& settingsName ) const
    {
        if( observableType_ != "AngularPosition" )
        {
            throw std::runtime_error( "Error in TrackingData, " + settingsName + " are only available for angular position observations." );
        }
    }

    //! Check that, for angular position observations, each observation has a weight for both right ascension and declination, or for
    //! neither of the two
    void checkAngularPositionWeights( const observation_models::ObservationWeightSettings& weightSettings,
                                      const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) ) const
    {
        if( observableType_ != "AngularPosition" || singleObservationSize_ != 2 )
        {
            return;
        }
        const unsigned int numberOfWeights = observationIndices.empty( ) ? numberOfObservations_ : observationIndices.size( );
        for( unsigned int i = 0; i < numberOfWeights; ++i )
        {
            const Eigen::Vector2d weights =
                    weightSettings.type_ == observation_models::ObservationWeightSettings::WeightsBlockType::set_block
                    ? Eigen::Vector2d( weightSettings.weightBlock_.diagonal( ).segment< 2 >( 2 * i ) )
                    : Eigen::Vector2d( weightSettings.getObservationWeightBlock( i, 2 ).diagonal( ) );
            if( ( weights( 0 ) == 0.0 ) != ( weights( 1 ) == 0.0 ) )
            {
                throw std::runtime_error( "Error when setting TrackingData weights, observation " +
                                          std::to_string( observationIndices.empty( ) ? i : observationIndices.at( i ) ) +
                                          " has a weight for only one of right ascension and declination." );
            }
        }
    }

    //! Check whether the positional weights are diagonal weights given for every observation
    bool hasDiagonalWeightsPerObservation( ) const
    {
        return weightSettings_->positionalWeights_ != nullptr && weightSettings_->positionalWeights_->observationIndices_.empty( ) &&
                weightSettings_->positionalWeights_->weights_.type_ ==
                observation_models::ObservationWeightSettings::WeightsBlockType::diagonal_per_observation;
    }

    const std::string observableType_;

    const PlainLinkDefinition linkEnds_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations_;

    std::vector< TimeType > epochs_;

    const std::string referenceLinkEnd_;

    const std::string timeScale_;

    unsigned int numberOfObservations_;

    const unsigned int singleObservationSize_;

    //! Settings that define the weights of the observations
    std::shared_ptr< TrackingDataWeightSettings > weightSettings_ = std::make_shared< TrackingDataWeightSettings >( );

    //! Settings that define the corrections to the observations
    std::shared_ptr< TrackingDataCorrectionSettings< ObservationScalarType > > correctionSettings_ =
            std::make_shared< TrackingDataCorrectionSettings< ObservationScalarType > >( );

    std::map< std::string, std::string > ancillarySettingsString_;

    std::map< std::string, std::vector< std::string > > ancillarySettingsStringVector_;

    std::set< std::string > observationMetadataKeys_;

    std::map< std::string, std::vector< double > > numericalObservationMetadata_;

    std::map< std::string, double > ancillarySettingsDouble_;

    std::map< std::string, std::vector< double > > ancillarySettingsDoubleVector_;
};

}  // namespace data

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_H
