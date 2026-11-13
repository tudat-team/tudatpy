/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRACKING_DATA_SETTINGS_H
#define TUDAT_TRACKING_DATA_SETTINGS_H

#include <algorithm>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "tudat/io/angularObservationCorrectionSettings.h"
#include "tudat/io/observationWeightSettings.h"

namespace tudat
{

namespace data
{

//! Weights of (a subset of) the observations in a TrackingData object, that are known before creating an observation dataset
class PositionalWeightSettings
{
public:
    //! Constructor
    /*!
     * \param weights Weights of the observations given by observationIndices (of all observations if observationIndices is empty)
     * \param observationIndices Indices of the observations for which weights are given (empty for all observations)
     */
    PositionalWeightSettings( const observation_models::ObservationWeightSettings& weights,
                              const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) ):
        weights_( weights ), observationIndices_( observationIndices )
    {}

    //! Remove a single observation, return whether weights for any observation remain
    bool removeObservation( const unsigned int observationIndex, const unsigned int singleObservationSize )
    {
        if( observationIndices_.empty( ) )
        {
            weights_.removeObservation( observationIndex, singleObservationSize );
            return true;
        }

        auto indexIterator = std::find( observationIndices_.begin( ), observationIndices_.end( ), observationIndex );
        if( indexIterator != observationIndices_.end( ) )
        {
            weights_.removeObservation( std::distance( observationIndices_.begin( ), indexIterator ), singleObservationSize );
            observationIndices_.erase( indexIterator );
        }
        for( unsigned int& index : observationIndices_ )
        {
            if( index > observationIndex )
            {
                index--;
            }
        }
        return !observationIndices_.empty( );
    }

    //! Weights of the observations given by observationIndices_ (of all observations if observationIndices_ is empty)
    observation_models::ObservationWeightSettings weights_;

    //! Indices of the observations for which weights are given (empty for all observations)
    std::vector< unsigned int > observationIndices_;
};

//! Named methods to compute observation weights when creating an observation dataset
enum class NamedWeightingMethod {
    vfcc17  //!< Weights of Veres et al. (2017), for all observations that have no positional weights
};

//! Settings for observation weights computed with a named method, for all observations that have no positional weights
class NamedWeightSettings
{
public:
    NamedWeightSettings( const NamedWeightingMethod method ): method_( method ) {}

    NamedWeightingMethod method_;
};

//! Definitions of blocks of observations for which the weights are reduced
enum class ObservationBlockDefinition {
    local_night  //!< Observations taken in the same local night, from the same station
};

//! Rules to reduce the weights of the observations in a block
enum class BlockDeweightingRule {
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

    virtual ~TimingBasedUncertaintySettings( ) {}

    void removeObservation( const unsigned int observationIndex )
    {
        timeUncertainties_.erase( timeUncertainties_.begin( ) + observationIndex );
    }

    //! (1-sigma) uncertainties of the observation epochs (in seconds), one per observation
    std::vector< double > timeUncertainties_;
};

//! Uncertainties of the epochs of angular position observations. For each observation, the covariance is augmented by
//! sigma_T^2 r r^T (Farnocchia et al., 2022), with r the rates of right ascension and declination
class AstrometricTimeUncertaintySettings : public TimingBasedUncertaintySettings
{
public:
    using TimingBasedUncertaintySettings::TimingBasedUncertaintySettings;
};

//! Settings that define the weights of the observations in a TrackingData object
class TrackingDataWeightSettings
{
public:
    //! Remove a single observation from all settings
    void removeObservation( const unsigned int observationIndex, const unsigned int singleObservationSize )
    {
        if( positionalWeights_ != nullptr && !positionalWeights_->removeObservation( observationIndex, singleObservationSize ) )
        {
            positionalWeights_ = nullptr;
        }
        if( timingUncertainties_ != nullptr )
        {
            timingUncertainties_->removeObservation( observationIndex );
        }
    }

    //! Weights that are known before creating the observation dataset (nullptr if none)
    std::shared_ptr< PositionalWeightSettings > positionalWeights_;

    //! Weights computed with a named method, for the observations without positional weights (nullptr if none)
    std::shared_ptr< NamedWeightSettings > namedWeights_;

    //! Reduction of the weights of observations in the same block (nullptr if none)
    std::shared_ptr< BlockDeweightingSettings > blockDeweighting_;

    //! Uncertainties of the observation epochs, added to the observation covariance (nullptr if none)
    std::shared_ptr< TimingBasedUncertaintySettings > timingUncertainties_;
};

//! Corrections that are known before creating an observation dataset, one per observation
template< typename ObservationScalarType = double >
class ExplicitCorrectionSettings
{
public:
    ExplicitCorrectionSettings( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& corrections ):
        corrections_( corrections )
    {}

    void removeObservation( const unsigned int observationIndex )
    {
        corrections_.erase( corrections_.begin( ) + observationIndex );
    }

    //! Corrections that are added to the observations, one per observation
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > corrections_;
};

//! Base class for biases of the observation epochs, that are converted to observation corrections when creating an observation
//! dataset (using the time derivatives of the observations)
class TimingBasedCorrectionSettings
{
public:
    TimingBasedCorrectionSettings( const std::vector< double >& timeBiases ): timeBiases_( timeBiases ) {}

    virtual ~TimingBasedCorrectionSettings( ) {}

    void removeObservation( const unsigned int observationIndex )
    {
        timeBiases_.erase( timeBiases_.begin( ) + observationIndex );
    }

    //! Biases of the observation epochs (in seconds), one per observation. An observation with epoch t and time bias c is taken to
    //! be h(t+c), and is corrected to h(t) using its time derivative
    std::vector< double > timeBiases_;
};

//! Biases of the epochs of angular position observations, converted to corrections using the rates of right ascension and declination
class AstrometricTimeBiasSettings : public TimingBasedCorrectionSettings
{
public:
    using TimingBasedCorrectionSettings::TimingBasedCorrectionSettings;
};

//! Settings that define the corrections to the observations in a TrackingData object (applied when creating an observation dataset
//! with corrections applied)
template< typename ObservationScalarType = double >
class TrackingDataCorrectionSettings
{
public:
    //! Remove a single observation from all settings
    void removeObservation( const unsigned int observationIndex )
    {
        if( explicitCorrections_ != nullptr )
        {
            explicitCorrections_->removeObservation( observationIndex );
        }
        if( timeBiases_ != nullptr )
        {
            timeBiases_->removeObservation( observationIndex );
        }
    }

    //! Corrections that are known before creating the observation dataset (nullptr if none)
    std::shared_ptr< ExplicitCorrectionSettings< ObservationScalarType > > explicitCorrections_;

    //! Biases of the observation epochs (nullptr if none)
    std::shared_ptr< TimingBasedCorrectionSettings > timeBiases_;

    //! Corrections computed from models of the environment (nullptr if none)
    std::shared_ptr< ModelBasedCorrectionSettings > modelCorrections_;
};

}  // namespace data

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_SETTINGS_H
