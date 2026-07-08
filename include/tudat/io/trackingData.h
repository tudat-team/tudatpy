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
#include <memory>
#include <vector>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace data
{

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class TrackingData
{
public:
    TrackingData( const std::string observableType,
                  const std::vector< std::pair< std::pair< std::string, std::string >, std::string > >& linkEnds,
                  const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                  const std::vector< TimeType > epochs,
                  const std::string referenceLinkEnd ):
        observableType_( observableType ), linkEnds_( linkEnds ), observations_( observations ), epochs_( epochs ),
        referenceLinkEnd_( referenceLinkEnd ), numberOfObservations_( observations.size( ) )
    {
        // Check inputs size consistency
        if( observations_.size( ) != epochs_.size( ) )
        {
            throw std::runtime_error( "Error when creating TrackingData object, 
                numbers of epochs ("+std::to_string( epochs_.size( ) )+") 
                and observations ("+std::to_string( observations_.size( ) )+") are inconsistent." );
        }

        for( unsigned int i = 1; i < observations.size( ); i++ )
        {
            if( observations.at( i ).rows( ) != observations.at( i - 1 ).rows( ) )
            {
                throw std::runtime_error( "Error when making TrackingData, input observables not of consistent size." );
            }
        }

        // Set size of a single observation (zero is observations vector is empty)
        singleObservationSize_ = ( !observations.empty( ) ? observations[ 0 ].size( ) : 0 );
    }

    //! Function that returns the number of observations
    unsigned int getNumberOfObservations( )
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
    std::string getObservableType( )
    {
        return observableType_;
    }

    //! Function that returns the link ends description
    std::vector< std::pair< std::pair< std::string, std::string >, std::string > > getLinkEnds( )
    {
        return linkEnds_;
    }

    //! Function that returns the link end used as reference
    std::string getReferenceLinkEnd( )
    {
        return referenceLinkEnd_;
    }

    //! Function that returns observation values
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservations( )
    {
        return observations_;
    }

    //! Function that returns a concatenated vector of observations
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

    //! Function that overwrites the full observation vector
    void resetObservations( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& newObservations )
    {
        // Check that the number of observations is consistent
        if( newObservations.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error("Error when resetting the observation values in TrackingData object, the new number of 
                observations (" + std::to_string(newObservations.size( )) + ") is inconsistent with original number of observations
                (" + std::to_string(numberOfObservations_) +")." );
        }

        // Check that the size of each single observation is consistent
        for( auto obs : newObservations )
        {
            if( obs.size( ) != singleObservationSize_ )
            {
                throw std::runtime_error("Error when resetting the observation values in TrackingData object, the size of a 
                    single new observation (" + std::to_string(newObservations.size( )) + ") is inconsistent with the 
                    original single observation size (" + std::to_string(singleObservationSize_) +")." );
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
            throw std::runtime_error( "Error when resetting single observation value in TrackingData object, index exceeds 
                number of observations contained in TrackingData object." );
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
    std::vector< TimeType > getObservationEpochs( )
    {
        return epochs_;
    }

    //! Function that returns a concatenated vector of observation epochs
    Eigen::Matrix< TimeType, Eigen::Dynamic, 1 > getObservationEpochsVector( )
    {
        Eigen::Matrix< TimeType, Eigen::Dynamic, 1 > epochsVector =
                Eigen::Matrix< TimeType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_ * numberOfObservations_, 1 );
        for( unsigned int i = 0; i < epochs_.size( ); i++ )
        {
            epochsVector.segment( i * singleObservationSize_, singleObservationSize_ ) = epochs_.at( i );
        }
        return epochsVector;
    }

    // Add ancillary settings (string type)
    void addAncillarySettings( const std::string ancillarySettingsType, const std::string ancillarySettingsValue )
    {
        ancillarySettingsString_[ ancillarySettingsType ] = ancillarySettingsValue;
    }

    //! Function that returns map of ancillary settings (string type)
    std::map< std::string, std::string > getAncillarySettingsString( ) const
    {
        return ancillarySettingsString_;
    }

    // Add ancillary settings (double type)
    void addAncillarySettings( const std::string ancillarySettingsType, const double ancillarySettingsValue )
    {
        ancillarySettingsDouble_[ ancillarySettingsType ] = ancillarySettingsValue;
    }

    //! Function that returns map of ancillary settings (double type)
    std::map< std::string, double > getAncillarySettingsDouble( ) const
    {
        return ancillarySettingsDouble_;
    }

    // Add ancillary settings (vector type)
    void addAncillarySettings( const std::string ancillarySettingsType, const std::vector< double > ancillarySettingsValue )
    {
        ancillarySettingsVector_[ ancillarySettingsType ] = ancillarySettingsValue;
    }

    //! Function that returns map of ancillary settings (vector type)
    std::map< std::string, vector< double > > getAncillarySettingsVector( ) const
    {
        return ancillarySettingsVector_;
    }

    //! Set observation weights to the tracking data object (optional)
    void setObservationWeights( const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& observationWeights )
    {
        // Check if observation weights already existed and overwrite them if they did (+throw a warning)
        if( !weights_.empty( ) )
        {
            std::cerr << "Warning when adding observation weights to tracking data object, weights already existed and 
                         are overwritten
                                 ." << std::endl;
                         weights_.clear( );
        }

        // Check size consistency (for the total number of observations)
        if( observationWeights.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Error when adding observation weights to tracking data object, size of weights (" +
                                      std::to_string( observationWeights.size( ) ) + ") does not match number of observations (" +
                                      std::to_string( numberOfObservations_ ) ")." );
        }

        // Check size consistency (for a single observable)
        for( auto weight : observationWeights )
        {
            if( weight.size( ) != singleObservationSize_ )
            {
                throw std::runtime_error( "Error when adding observation weights to tracking data object, size of single weight (" +
                                          std::to_string( weight.size( ) ) + ") does not match single observable size (should be " +
                                          std::to_string( singleObservationSize_ ) ")." );
            }
        }

        // If all sizes are consistent, store observation weights
        weights_ = observationWeights;
    }

    //! Function that reset a specific observation weight (only possible if weights are already available)
    void setSingleObservationWeight( const unsigned int index, const Eigen::Matrix< double, Eigen::Dynamic, 1 >& observationWeight )
    {
        // Check if weights are already available
        if( weights.empty( ) )
        {
            throw std::runtime_error( "Error when resetting single observation weight in TrackingData object, weights not 
                yet defined." );
        }

        // Check that the observation index for which the weight needs resetting does not exceed the size of the observation vector
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when resetting single observation weight in TrackingData object, index exceeds 
                number of observations contained in TrackingData object." );
        }

        // Check size consistency
        if( observationWeight.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error( "Error when resetting single observation weight in TrackingData object, weight size (" +
                                      std::to_string( observationWeight.size( ) ) + ") inconsistent with observation size (" +
                                      std::to_string( singleObservationSize_ ) + ")." );
        }

        // Overwrites specific weight entry
        weights_.at( index ) = observationWeight;
    }

    //! Function that returns a vector of observation weights (empty if observation weights are not provided)
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservationWeights( ) const
    {
        return weights_;
    }

    //! Function that returns a concatenated vector of observation weights (size zero if no observation weights provided)
    Eigen::VectorXd getWeightsVector( ) const
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > weightsVector =
                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( weights_.size( ) * singleObservationSize_, 1 );
        for( unsigned int i = 0; i < weights_; i++ )
        {
            weightsVector.block( i * singleObservationSize_, 0, singleObservationSize_, 1 ) = weights_.at( i );
        }
        return weightsVector;
    }

    //! Set corrections to the observations (optional)
    void setObservationCorrections( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observationCorrections )
    {
        // Check if observation corrections already existed and clear them if they did + throw a warning (overwritten)
        if( !observationCorrections_.empty( ) )
        {
            std::cerr << "Warning when adding observation corrections to tracking data object, corrections already existed and 
                         are overwritten
                                 ." << std::endl;
                         observationCorrections_.clear( );
        }

        // Check size consistency (for the total number of observations)
        if( observationCorrections.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Error when adding observation corrections to tracking data object, size of corrections (" +
                                      std::to_string( observationCorrections.size( ) ) + ") does not match number of observations (" +
                                      std::to_string( numberOfObservations_ ) ")." );
        }

        // Check size consistency (for a single observable)
        for( auto correction : observationCorrections )
        {
            if( correction.size( ) != singleObservationSize_ )
            {
                throw std::runtime_error( "Error when adding observation corrections to tracking data object, size of single correction (" +
                                          std::to_string( correction.size( ) ) + ") does not match single observable size (should be " +
                                          std::to_string( singleObservationSize_ ) ")." );
            }
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
            throw std::runtime_error( "Error when resetting single observation correction in TrackingData object, corrections not 
                yet defined." );
        }

        // Check that the index of the observation correction that needs overwriting does not exceed the size of the observation vector
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error( "Error when resetting single observation correction in TrackingData object, index exceeds 
                number of observations contained in TrackingData object." );
        }

        // Check size consistency
        if( observationCorrection.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error( "Error when resetting single observation correction in TrackingData object, new correction size (" +
                                      std::to_string( observationCorrection.size( ) ) + " inconsistent with observation size (" +
                                      std::to_string( singleObservationSize_ ) + ")." );
        }

        // Overwrites specific correction entry
        observationCorrections_.at( index ) = observationCorrection;
    }

    //! Function that return a vector of observation corrections (empty of observation corrections are not provided)
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservationCorrections( ) const
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

        // Remove associated weight (if it exists)
        if( !weights_.empty( ) )
        {
            weights_.erase( weights_.begin( ) + index );
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

    const std::vector< std::pair< std::pair< std::string, std::string >, std::string > > linkEnds_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations_;

    std::vector< TimeType > epochs_;

    const std::string referenceLinkEnd_;

    unsigned int numberOfObservations_;

    const unsigned int singleObservationSize_;

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observationCorrections_;

    std::map< std::string, std::string > ancillarySettingsString_;

    std::map< std::string, double > ancillarySettingsDouble_;

    std::map< std::string, vector< double > > ancillarySettingsVector_;
};

}  // namespace data

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_H
