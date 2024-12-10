/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 */

#ifndef TUDAT_DOPPLERMEASUREDFREQUENCYOBSERVATIONMODEL_H
#define TUDAT_DOPPLERMEASUREDFREQUENCYOBSERVATIONMODEL_H

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


inline double getMeasuredFrequencyDopplerScalingFactor(
    const std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > receivedFrequencyFunction,
    const observation_models::LinkEndType referenceLinkEnd,
    const std::vector< Eigen::Vector6d >& linkEndStates,
    const std::vector< double >& linkEndTimes,
    const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings )
{
    double integrationTime;
    std::vector< FrequencyBands > frequencyBands;
    try
    {
        frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncilliaryDoubleVectorData( frequency_bands ) );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error(
            "Error when retrieving ancillary settings for measured frequency observable: " +
            std::string( caughtException.what( ) ) );
    }

    double transmissionTime = linkEndTimes.at( 0 );
    double frequency = receivedFrequencyFunction( frequencyBands, transmissionTime );

    // Moyer (2000), eq. 13-59
    return frequency / physical_constants::getSpeedOfLight< double >( );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
class DopplerMeasuredFrequencyObservationModel
    : public ObservationModel< 1, ObservationScalarType, TimeType >
{
   public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;

    /*! Constructor
     * Constructor for Doppler observation model.
     * \param linkEnds Link ends for observation model.
     * \param rangeObservationModel N-way range observation model associated with the Doppler
     * measurement. \param transmittingFrequencyCalculator Object returning the transmitted
     * frequency as the transmitting link end \param turnaroundRatio Function returning the
     * turnaround ratio for given uplink and downlink bands \param observationBiasCalculator Object
     * for calculating system-dependent errors in the observation
     */
    DopplerMeasuredFrequencyObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< TwoWayDopplerObservationModel< ObservationScalarType,
                                                                  TimeType > > twoWayDopplerModel,
            const std::shared_ptr< ground_stations::StationFrequencyInterpolator >
                    transmittingFrequencyCalculator,
            const std::function< double( observation_models::FrequencyBands uplinkBand,
                                         observation_models::FrequencyBands downlinkBand ) >&
                    turnaroundRatio,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
            const std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > >
                    groundStationStates = std::map<
                            LinkEndType,
                            std::shared_ptr< ground_stations::GroundStationState > >( ) ) :
        ObservationModel< 1, ObservationScalarType, TimeType >( doppler_measured_frequency,
                                                                linkEnds,
                                                                observationBiasCalculator ),
        twoWayDopplerModel_( twoWayDopplerModel ), numberOfLinkEnds_( linkEnds.size( ) ),
        transmittingFrequencyCalculator_( transmittingFrequencyCalculator ),
        turnaroundRatio_( turnaroundRatio ), stationStates_( groundStationStates )
    {
        if( numberOfLinkEnds_ != 3 )
        {
            throw std::runtime_error(
                    "Error when defining Doppler Measured Frequency Model: model allows exactly 3 "
                    "link ends, " +
                    std::to_string( numberOfLinkEnds_ ) + "were selected." );
        }

        uplinkDopplerModel_ = twoWayDopplerModel_->getUplinkDopplerCalculator( );
        downlinkDopplerModel_ = twoWayDopplerModel_->getDownlinkDopplerCalculator( );

        std::shared_ptr<
                observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                uplinkLightTimeCalculator = uplinkDopplerModel_->getLightTimeCalculator( );
        std::shared_ptr<
                observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                downlinkLightTimeCalculator = downlinkDopplerModel_->getLightTimeCalculator( );

        std::vector< std::shared_ptr<
                observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >
                lightTimeCalculators = { uplinkLightTimeCalculator, downlinkLightTimeCalculator };
        lighTimeCalculator_ = std::make_shared<
                observation_models::MultiLegLightTimeCalculator< ObservationScalarType,
                                                                 TimeType > >(
                lightTimeCalculators );

        terrestrialTimeScaleConverter_ = earth_orientation::createDefaultTimeConverter( );
    }

    //! Destructor
    ~DopplerMeasuredFrequencyObservationModel( ) { }

    /*!
     * Function to compute Measured Frequency for a doppler observation model
     * \param time Time at which observable is to be evaluated.
     * \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for
     * which associated time is kept constant (to input value) \param linkEndTimes List of times at
     * each link end during observation. Set in this mothod. \param linkEndStates List of states at
     * each link end during observation. Set in this method. \param ancillarySettings Ancillary
     * settings for the observation model \return Measured Frequency for a doppler observation model
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings =
                    nullptr )
    {
        // Check if selected reference link end is valid
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error(
                    "Error when computing Doppler measured frequency observables: the selected "
                    "reference link end (" +
                    getLinkEndTypeString( linkEndAssociatedWithTime ) +
                    ") is not valid. Must be the receiver." );
        }

        // Check if ancillary settings were provided
        if( ancillarySettings == nullptr )
        {
            throw std::runtime_error(
                    "Error when simulating Doppler Measured Frequency observable; no ancillary "
                    "settings found. " );
        }

        std::vector< FrequencyBands > frequencyBands;

        try
        {
            frequencyBands = convertDoubleVectorToFrequencyBands(
                    ancillarySettings->getAncilliaryDoubleVectorData( frequency_bands ) );
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error(
                    "Error when retrieving ancillary settings for Doppler measured frequency "
                    "observable: " +
                    std::string( caughtException.what( ) ) );
        }

        if( frequencyBands.size( ) != numberOfLinkEnds_ - 1 )
        {
            throw std::runtime_error(
                    "Error when retrieving frequency bands ancillary settings: "
                    "size (" +
                    std::to_string( frequencyBands.size( ) ) +
                    ") is inconsistent with number of links (" +
                    std::to_string( numberOfLinkEnds_ - 1 ) + ")." );
        }

        FrequencyBands uplinkBand = frequencyBands.at( 0 );
        FrequencyBands downlinkBand = frequencyBands.at( 1 );

        // Calculate the light time
        TimeType lightTime = lighTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancillarySettings );

        // Get the time when the signal left the transmitter
        Eigen::Vector3d nominalTransmittingStationState =
                ( stationStates_.count( transmitter ) == 0 )
                ? Eigen::Vector3d::Zero( )
                : stationStates_.at( transmitter )->getNominalCartesianPosition( );
        TimeType transmitterTime = time - lightTime;

        TimeType transmitterUtcTime = terrestrialTimeScaleConverter_->getCurrentTime< TimeType >(
                basic_astrodynamics::tdb_scale,
                basic_astrodynamics::utc_scale,
                transmitterTime,
                nominalTransmittingStationState );

        // Get the frequency of the transmitter
        ObservationScalarType transmittedFrequency =
                transmittingFrequencyCalculator_
                        ->getTemplatedCurrentFrequency< ObservationScalarType, TimeType >(
                                transmitterUtcTime );

        // Calculate the Doppler observable
        ObservationScalarType dopplerMultiplicationTerm =
                twoWayDopplerModel_->getMultiplicationTerm( );
        ObservationScalarType twoWayDoppler =
                twoWayDopplerModel_->computeIdealObservationsWithLinkEndData(
                        time,
                        linkEndAssociatedWithTime,
                        linkEndTimes,
                        linkEndStates,
                        ancillarySettings )( 0, 0 ) /
                dopplerMultiplicationTerm;

        ObservationScalarType receivedFrequency =
                ( transmittedFrequency *
                  ( mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 ) +
                    twoWayDoppler ) ) *
                turnaroundRatio_( uplinkBand, downlinkBand );

        Eigen::Matrix< ObservationScalarType, 1, 1 > observation =
                ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << receivedFrequency )
                        .finished( );
        return observation;
    }


    std::shared_ptr< TwoWayDopplerObservationModel< ObservationScalarType, TimeType > > getTwoWayDopplerModel( )
    {
        return twoWayDopplerModel_;
    }


    // Doppler observation model associated with the measurement
    std::shared_ptr< TwoWayDopplerObservationModel< ObservationScalarType, TimeType > >
            twoWayDopplerModel_;

    // Number of link ends
    unsigned int numberOfLinkEnds_;

    // Object returning the transmitted frequency as the transmitting link end
    std::shared_ptr< ground_stations::StationFrequencyInterpolator >
            transmittingFrequencyCalculator_;

    // Function returning the turnaround ratio for given uplink and downlink bands
    std::function< double( FrequencyBands uplinkBand, FrequencyBands downlinkBand ) >
            turnaroundRatio_;

    // Individual doppler observation models
    std::shared_ptr<
            observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > >
            uplinkDopplerModel_;

    std::shared_ptr<
            observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > >
            downlinkDopplerModel_;

    // Light time calculator
    std::shared_ptr<
            observation_models::MultiLegLightTimeCalculator< ObservationScalarType, TimeType > >
            lighTimeCalculator_;

    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >
            terrestrialTimeScaleConverter_;

    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > stationStates_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_DOPPLERMEASUREDFREQUENCYOBSERVATIONMODEL_H
