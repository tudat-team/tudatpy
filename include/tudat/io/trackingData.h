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
#include "tudat/io/trackingDataCorrectionSettings.h"
#include "tudat/io/trackingDataWeightSettings.h"

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

    //! Function that returns the settings that define the observation weights (set and retrieve the weights through these settings)
    const std::shared_ptr< TrackingDataWeightSettings >& getWeightSettings( ) const
    {
        return weightSettings_;
    }

    //! Function that replaces the settings that define the observation weights
    void setWeightSettings( const std::shared_ptr< TrackingDataWeightSettings > weightSettings )
    {
        if( weightSettings == nullptr || weightSettings->getNumberOfObservations( ) != numberOfObservations_ ||
            weightSettings->getSingleObservationSize( ) != singleObservationSize_ ||
            weightSettings->getObservableType( ) != observableType_ )
        {
            throw std::runtime_error( "Error when setting TrackingData weight settings, settings are not defined for these observations." );
        }
        weightSettings_ = weightSettings;
    }

    //! Function that returns the settings that define the observation corrections (set and retrieve the corrections through these
    //! settings)
    const std::shared_ptr< TrackingDataCorrectionSettings< ObservationScalarType > >& getCorrectionSettings( ) const
    {
        return correctionSettings_;
    }

    //! Function that replaces the settings that define the observation corrections
    void setCorrectionSettings( const std::shared_ptr< TrackingDataCorrectionSettings< ObservationScalarType > > correctionSettings )
    {
        if( correctionSettings == nullptr || correctionSettings->getNumberOfObservations( ) != numberOfObservations_ ||
            correctionSettings->getSingleObservationSize( ) != singleObservationSize_ ||
            correctionSettings->getObservableType( ) != observableType_ || correctionSettings->getReferenceLinkEnd( ) != referenceLinkEnd_ )
        {
            throw std::runtime_error(
                    "Error when setting TrackingData correction settings, settings are not defined for these observations." );
        }
        correctionSettings_ = correctionSettings;
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
        weightSettings_->removeObservation( index );
        correctionSettings_->removeObservation( index );

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

    unsigned int numberOfObservations_;

    const unsigned int singleObservationSize_;

    //! Settings that define the weights of the observations
    std::shared_ptr< TrackingDataWeightSettings > weightSettings_ =
            std::make_shared< TrackingDataWeightSettings >( numberOfObservations_, singleObservationSize_, observableType_ );

    //! Settings that define the corrections to the observations
    std::shared_ptr< TrackingDataCorrectionSettings< ObservationScalarType > > correctionSettings_ =
            std::make_shared< TrackingDataCorrectionSettings< ObservationScalarType > >( numberOfObservations_,
                                                                                         singleObservationSize_,
                                                                                         observableType_,
                                                                                         referenceLinkEnd_ );

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
