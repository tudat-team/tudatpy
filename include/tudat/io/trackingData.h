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

// #include "tudat/astro/observation_models/linkTypeDefs.h"
// #include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"
// #include "tudat/simulation/estimation_setup/observationOutput.h"
// #include "tudat/simulation/estimation_setup/observationsProcessing.h"

namespace tudat
{

namespace data
{

std::vector< std::string > observableTypeStrings = {};

// using namespace simulation_setup;

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
                  const std::string referenceLinkEnd,
                  const unsigned int singleObservationSize,
                  const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights = {} ):
        observableType_( observableType ), linkEnds_( linkEnds ), observations_( observations ), epochs_( epochs ),
        referenceLinkEnd_( referenceLinkEnd ), numberOfObservations_( observations.size( ) ),
        singleObservationSize_( singleObservationSize ), weights_( weights )
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

        // Check weight size consistency (if it exists)
        if( weights.size( ) != observationTimes.size( ) )
        {
            throw std::runtime_error( "Error when creating TrackingData set with weights; size is incompatible" );
        }

        for( std::size_t k = 0; k < weights.size( ); ++k )
        {
            if( weights.at( k ).size( ) != static_cast< int >( singleObservationSize_ ) )
            {
                throw std::runtime_error( "Error when creating observation set with weights; individual weight size is incompatible" );
            }
        }
    }

    bool checkTrackingDataValidity( ) {}

    std::string getObservableType( )
    {
        return observableType_;
    }

    std::vector< std::pair< std::pair< std::string, std::string >, std::string > > getLinkEnds( )
    {
        return linkEnds_;
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservations( )
    {
        return observations_;
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
            throw std::runtime_error( "Error when setting single observation value in tracking data, index is out of bounds" );
        }
        if( observation.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error(
                    "Error when setting single observation value in tracking data, the observation size is inconsistent." );
        }
        observations_.at( index ) = observation;
    }

    std::vector< TimeType > getObservationEpochs( )
    {
        return epochs_;
    }

    const std::vector< TimeType >& getObservationEpochs( )
    {
        return epochs_;
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

    std::map< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservationsHistory( )
    {
        return utilities::createMapFromVectors< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( observationTimes_,
                                                                                                                       observations_ );
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
        for( unsigned int i = 0; i < numberOfObservations_; i++ )
        {
            weightsVector.block( i * singleObservationSize_, 0, singleObservationSize_, 1 ) = weights_.at( i );
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

    // Add ancillary settings (string type)
    void addAncillarySettings( const std::string ancillarySettingsType, const std::string ancillarySettingsValue )
    {
        ancillarySettingsString_[ ancillarySettingsType ] = ancillarySettingsValue;
    }

    // Add ancillary settings (double type)
    void addAncillarySettings( const std::string ancillarySettingsType, const double ancillarySettingsValue )
    {
        ancillarySettingsDouble_[ ancillarySettingsType ] = ancillarySettingsValue;
    }

    // Add ancillary settings (vector type)
    void addAncillarySettings( const std::string ancillarySettingsType, const std::vector< double > ancillarySettingsValue )
    {
        ancillarySettingsVector_[ ancillarySettingsType ] = ancillarySettingsValue;
    }

    // Retrieve ancillary settings (string type)
    std::map< std::string, std::string > getAncillarySettingsString( ) const
    {
        return ancillarySettingsString_;
    }

    // Retrieve ancillary settings (double type)
    std::map< std::string, double > getAncillarySettingsDouble( ) const
    {
        return ancillarySettingsDouble_;
    }

    // Retrieve ancillary settings (vector type)
    std::map< std::string, vector< double > > getAncillarySettingsVector( ) const
    {
        return ancillarySettingsVector_;
    }

private:
    const std::string observableType_;

    const std::vector< std::pair< std::pair< std::string, std::string >, std::string > > linkEnds_;

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations_;

    const std::vector< TimeType > epochs_;

    const std::string referenceLinkEndIndex_;

    const unsigned int numberOfObservations_;

    const unsigned int singleObservationSize_;

    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights_;

    const std::map< std::string, std::string > ancillarySettingsString_;

    const std::map< std::string, double > ancillarySettingsDouble_;

    const std::map< std::string, vector< double > > ancillarySettingsVector_;
};

}  // namespace data

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_H
