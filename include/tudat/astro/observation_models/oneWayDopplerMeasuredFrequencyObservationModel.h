/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 */

#ifndef TUDAT_ONEWAYDOPPLERMEASUREDFREQUENCYOBSERVATIONMODEL_H
#define TUDAT_ONEWAYDOPPLERMEASUREDFREQUENCYOBSERVATIONMODEL_H

#include <stdexcept>
#include <string>

#include "tudat/astro/observation_models/nWayRangeObservationModel.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationFrequencies.h"
#include "tudat/astro/observation_models/oneWayDopplerObservationModel.h"  // Include this header
#include "tudat/astro/observation_models/twoWayDopplerObservationModel.h"  // Include this header
#include "tudat/simulation/simulation.h"

namespace tudat
{

namespace observation_models
{

// inline double getMeasuredFrequencyDopplerScalingFactor(
//         const std::function< double( std::vector< FrequencyBands > frequencyBands, double time ) > receivedFrequencyFunction,
//         const observation_models::LinkEndType referenceLinkEnd,
//         const std::vector< Eigen::Vector6d >& linkEndStates,
//         const std::vector< double >& linkEndTimes,
//         const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings )
// {
//     double integrationTime;  // Not used?
//     std::vector< FrequencyBands > frequencyBands;
//     try
//     {
//         frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncilliaryDoubleVectorData( frequency_bands ) );
//     }
//     catch( std::runtime_error& caughtException )
//     {
//         throw std::runtime_error( "Error when retrieving ancillary settings for measured frequency observable: " +
//                                   std::string( caughtException.what( ) ) );
//     }

//     double transmissionTime = linkEndTimes.at( 0 );
//     double frequency = receivedFrequencyFunction( frequencyBands, transmissionTime );

//     // Moyer (2000), eq. 13-59
//     return frequency / physical_constants::getSpeedOfLight< double >( );
// }

template< typename ObservationScalarType = double, typename TimeType = Time >
class OneWayDopplerMeasuredFrequencyObservationModel : public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;

    using ObservationModel< 1, ObservationScalarType, TimeType >::timeScaleConverter_;
    using ObservationModel< 1, ObservationScalarType, TimeType >::frequencyInterpolator_;

    /**
     * @brief Constructor for the OneWayDopplerMeasuredFrequencyObservationModel.
     *
     * Initializes a one-way Doppler observation model that computes measured frequency observations
     * based on an underlying Doppler model and ground station frequency interpolation.
     *
     * @tparam ObservationScalarType The scalar type for observation values.
     * @tparam TimeType The scalar type for time values.
     *
     * @param linkEnds The link ends defining the transmitter and receiver stations. Must contain exactly 2 link ends.
     * @param oneWayDopplerModel A shared pointer to the underlying OneWayDopplerObservationModel used for
     *        Doppler shift calculations.
     * @param transmittingFrequencyCalculator A shared pointer to the StationFrequencyInterpolator for
     *        computing transmitting station frequencies.
     * @param observationBiasCalculator Optional shared pointer to an ObservationBias calculator for
     *        applying systematic biases to observations. Defaults to nullptr if not provided.
     * @param groundStationStates Optional map of ground station states indexed by LinkEndType.
     *        Defaults to an empty map if not provided.
     *
     * @throws std::runtime_error If the number of link ends is not exactly 2.
     *
     * @note The light-time calculator is obtained directly from the underlying one-way Doppler model
     *       and stored as a single LightTimeCalculator instance for use in this observation model.
     */
    OneWayDopplerMeasuredFrequencyObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< OneWayDopplerObservationModel< ObservationScalarType, TimeType > > oneWayDopplerModel,
            const std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyCalculator,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
            const std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > groundStationStates =
                    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > >( ) ):
        ObservationModel< 1, ObservationScalarType, TimeType >( one_way_doppler_measured_frequency, linkEnds, observationBiasCalculator ),
        oneWayDopplerModel_( oneWayDopplerModel ), numberOfLinkEnds_( linkEnds.size( ) ), stationStates_( groundStationStates )
    {
        // Check if OneWayDopplerModel is not nullptr
        if( oneWayDopplerModel_ == nullptr )
        {
            throw std::runtime_error(
                    "Error when creating OneWayDopplerMeasuredFrequencyObservationModel: input one-way Doppler model is nullptr." );
        }
        // Check if transmitting frequency calculator is not nullptr
        if( transmittingFrequencyCalculator == nullptr )
        {
            throw std::runtime_error(
                    "Error when creating OneWayDopplerMeasuredFrequencyObservationModel: input transmitting frequency calculator is "
                    "nullptr." );
        }

        this->setFrequencyInterpolator( transmittingFrequencyCalculator );

        if( numberOfLinkEnds_ != 2 )
        {
            throw std::runtime_error(
                    "Error when defining One Way Doppler Measured Frequency Model: model allows exactly 2 "
                    "link ends, " +
                    std::to_string( numberOfLinkEnds_ ) + " were selected." );
        }

        lightTimeCalculator_ = oneWayDopplerModel_->getLightTimeCalculator( );
    }

    //! Destructor
    ~OneWayDopplerMeasuredFrequencyObservationModel( ) {}

    /*!
     * Function to compute Measured Frequency for a doppler observation model
     * \param time Time at which observable is to be evaluated.
     * \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for
     * which associated time is kept constant (to input value) \param linkEndTimes List of times at
     * each link end during observation. Set in this method. \param linkEndStates List of states at
     * each link end during observation. Set in this method. \param ancillarySettings Ancillary
     * settings for the observation model \return Measured Frequency for a doppler observation model
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings = nullptr )
    {
        // Check if selected reference link end is valid
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error(
                    "Error when computing Doppler measured frequency observables: the selected "
                    "reference link end (" +
                    getLinkEndTypeString( linkEndAssociatedWithTime ) + ") is not valid. Must be the receiver." );
        }

        if( lightTimeCalculator_->doCorrectionsNeedFrequency( ) )
        {
            setTransmissionFrequency(
                    lightTimeCalculator_, timeScaleConverter_, frequencyInterpolator_, time, linkEndAssociatedWithTime, ancillarySettings );
        }

        // Calculate the light time using correct signature
        StateType receiverState, transmitterState;
        bool isTimeAtReception = ( linkEndAssociatedWithTime == receiver );
        TimeType lightTime = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                receiverState, transmitterState, time, isTimeAtReception, ancillarySettings );

        // Update link end states and times
        linkEndStates.resize( 2 );
        linkEndTimes.resize( 2 );
        linkEndStates[ 0 ] = transmitterState.template cast< double >( );
        linkEndStates[ 1 ] = receiverState.template cast< double >( );
        if( isTimeAtReception )
        {
            linkEndTimes[ 1 ] = static_cast< double >( time );
            linkEndTimes[ 0 ] = static_cast< double >( time - lightTime );
        }
        else
        {
            linkEndTimes[ 0 ] = static_cast< double >( time );
            linkEndTimes[ 1 ] = static_cast< double >( time + lightTime );
        }

        // Get the frequency of the transmitter
        TimeType transmitterTime = time - lightTime;

        Eigen::Vector3d transmitterPosition = transmitterState.template segment< 3 >( 0 ).template cast< double >( );

        TimeType transmitterUtcTime = timeScaleConverter_->template getCurrentTime< TimeType >(
                basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, transmitterTime, transmitterPosition );

        ObservationScalarType transmittedFrequency =
                frequencyInterpolator_->template getTemplatedCurrentFrequency< ObservationScalarType, TimeType >( transmitterUtcTime );

        // Calculate the Doppler observable
        ObservationScalarType dopplerMultiplicationTerm = oneWayDopplerModel_->getMultiplicationTerm( );
        ObservationScalarType oneWayDoppler =
                oneWayDopplerModel_->computeIdealObservationsWithLinkEndData(
                        time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancillarySettings )( 0, 0 ) /
                dopplerMultiplicationTerm;

        ObservationScalarType receivedFrequency =
                ( transmittedFrequency * ( mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 ) + oneWayDoppler ) );

        Eigen::Matrix< ObservationScalarType, 1, 1 > observation =
                ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << receivedFrequency ).finished( );
        return observation;
    }

    std::shared_ptr< OneWayDopplerObservationModel< ObservationScalarType, TimeType > > getOneWayDopplerModel( )
    {
        return oneWayDopplerModel_;
    }

    //! Function to get the frequency interpolator for the transmitting station
    /*!
     * Function to get the frequency interpolator for the transmitting station
     * \return Frequency interpolator for the transmitting station
     */
    std::shared_ptr< ground_stations::StationFrequencyInterpolator > getFrequencyInterpolator( )
    {
        return frequencyInterpolator_;
    }

    std::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculator( )
    {
        return lightTimeCalculator_;
    }

private:
    // Doppler observation model associated with the measurement
    std::shared_ptr< OneWayDopplerObservationModel< ObservationScalarType, TimeType > > oneWayDopplerModel_;

    // Number of link ends
    unsigned int numberOfLinkEnds_;

    // Light time calculator
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator_;

    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > stationStates_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_ONEWAYDOPPLERMEASUREDFREQUENCYOBSERVATIONMODEL_H
