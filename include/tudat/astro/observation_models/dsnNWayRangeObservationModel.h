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
 *          T. Moyer (2000), Formulation for Observed and Computed Values of Deep Space Network Data
 * Types for Navigation, DEEP SPACE COMMUNICATIONS AND NAVIGATION SERIES, JPL/NASA
 */

#ifndef TUDAT_DSNNWAYRANGEOBSERVATIONMODEL_H
#define TUDAT_DSNNWAYRANGEOBSERVATIONMODEL_H

#include <stdexcept>
#include <string>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/nWayRangeObservationModel.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationFrequencies.h"
#include "tudat/simulation/simulation.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = Time >
class DsnNWayRangeObservationModel : public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;

    /*! Constructor.
     *
     * @param linkEnds Map of the linkEnds defining the observation model
     * @param arcStartObservationModel N-way range observation model associated with the start of the Doppler integration time.
     * @param arcEndObservationModel N-way range observation model associated with the end of the Doppler integration time.
     * @param bodyWithGroundStations Body object where the ground stations are located.
     * @param observationBiasCalculator Object for calculating (system-dependent) errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points
     * (default none).
     */
    DsnNWayRangeObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< observation_models::MultiLegLightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator,
            const std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyCalculator,
            const std::function< double( observation_models::FrequencyBands uplinkBand, observation_models::FrequencyBands downlinkBand ) >&
                    turnaroundRatio,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
            const std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > groundStationStates =
                    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > >( ) ):
        ObservationModel< 1, ObservationScalarType, TimeType >( dsn_n_way_range, linkEnds, observationBiasCalculator ),
        lightTimeCalculator_( lightTimeCalculator ), numberOfLinkEnds_( linkEnds.size( ) ),
        transmittingFrequencyCalculator_( transmittingFrequencyCalculator ), turnaroundRatio_( turnaroundRatio ),
        stationStates_( groundStationStates )
    {
        if( !std::is_same< Time, TimeType >::value )
        {
            //            std::cerr<<
            //                    "Warning when defining DSN N-way range observation model: the
            //                    selected time type " "is not valid, using it would lead to large
            //                    numerical errors."<<std::endl;
        }

        if( numberOfLinkEnds_ != 3 )
        {
            throw std::runtime_error(
                    "Error when defining DSN N-way range observation model: model allows exactly 3 "
                    "link ends, " +
                    std::to_string( numberOfLinkEnds_ ) + "were selected." );
        }
        terrestrialTimeScaleConverter_ = earth_orientation::createDefaultTimeConverter( );
    }

    //! Destructor
    ~DsnNWayRangeObservationModel( ) { }

    /*! Function to compute DSN n-way Doppler observation at given time.
     *
     * Function to compute DSN n-way Doppler observation at given time. Only implemented for
     * receiver as the linkEndAssociatedWithTime. Computes the observable according to
     * section 13.3.2.2 of Moyer (2000).
     *
     * @param time Time at which observable is to be evaluated.
     * @param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     * @param linkEndTimes List of times at each link end during observation.
     * @param linkEndStates List of states at each link end during observation.
     * @param ancillarySettings Observation ancillary simulation settings.
     * @return Observation value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings = nullptr )
    {
        // Check if selected reference link end is valid
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error(
                    "Error when computing DSN N-way range observables: the selected reference link "
                    "end (" +
                    getLinkEndTypeString( linkEndAssociatedWithTime ) + ") is not valid." );
        }
        // Check if ancillary settings were provided
        if( ancillarySettings == nullptr )
        {
            throw std::runtime_error(
                    "Error when simulating n-way DSN range observable; no ancillary settings "
                    "found. " );
        }

        ObservationScalarType lowestRangingComponent;
        std::vector< FrequencyBands > frequencyBands;
        try
        {
            lowestRangingComponent = ancillarySettings->getAncilliaryDoubleData( sequential_range_lowest_ranging_component );
            //            referenceFrequency = ancillarySettings->getAncilliaryDoubleData(
            //            sequential_range_reference_frequency );
            frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncilliaryDoubleVectorData( frequency_bands ) );
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error( "Error when retrieving ancillary settings for DSN N-way range observable: " +
                                      std::string( caughtException.what( ) ) );
        }

        if( frequencyBands.size( ) != numberOfLinkEnds_ - 1 )
        {
            throw std::runtime_error(
                    "Error when retrieving frequency bands ancillary settings for DSN N-way range "
                    "observable: "
                    "size (" +
                    std::to_string( frequencyBands.size( ) ) + ") is inconsistent with number of links (" +
                    std::to_string( numberOfLinkEnds_ - 1 ) + ")." );
        }
        FrequencyBands uplinkBand = frequencyBands.at( 0 );
        FrequencyBands downlinkBand = frequencyBands.at( 1 );

        // Define the conversion factor based on the value of the uplink frequency band
        double conversionFactor;
        if( uplinkBand == 0 )  // S-band
        {
            conversionFactor = 0.5;
        }
        else if( uplinkBand == 1 )  // X-band
        {
            conversionFactor = 221.0 / ( 749.0 * 2.0 );
        }
        else if( uplinkBand == 2 )  // Ka-band
        {
            conversionFactor = 221 / ( 3599 * 2.0 );
        }
        else
        {
            throw std::runtime_error( "Unsupported uplink frequency band" );
        }

        // Set approximate up- and down-link frequencies.
        double currentTurnAroundRatio = static_cast< ObservationScalarType >( turnaroundRatio_( uplinkBand, downlinkBand ) );

        if( true )
        {
            setTransmissionReceptionFrequencies( lightTimeCalculator_,
                                                 terrestrialTimeScaleConverter_,
                                                 transmittingFrequencyCalculator_,
                                                 time,
                                                 ancillarySettings,
                                                 currentTurnAroundRatio );
        }

        TimeType lightTime = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancillarySettings );

        Eigen::Vector3d nominalReceivingStationState = ( stationStates_.count( receiver ) == 0 )
                ? Eigen::Vector3d::Zero( )
                : stationStates_.at( receiver )->getNominalCartesianPosition( );
        TimeType utcReceptionTime = terrestrialTimeScaleConverter_->getCurrentTime< TimeType >(
                basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, time, nominalReceivingStationState );
        TimeType utcTransmissionTime = terrestrialTimeScaleConverter_->getCurrentTime< TimeType >(
                basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, time - lightTime, nominalReceivingStationState );

        ObservationScalarType uplinkFrequency =
                transmittingFrequencyCalculator_->template getTemplatedCurrentFrequency< ObservationScalarType, TimeType >(
                        utcTransmissionTime );
        ancillarySettings->setAncilliaryDoubleData( observation_models::range_conversion_factor,
                                                    physical_constants::SPEED_OF_LIGHT / ( uplinkFrequency * conversionFactor ) );

        ObservationScalarType transmitterFrequencyIntegral =
                transmittingFrequencyCalculator_->template getTemplatedFrequencyIntegral< ObservationScalarType, TimeType >(
                        utcTransmissionTime, utcReceptionTime );
        ObservationScalarType rangeUnitIntegral = conversionFactor * transmitterFrequencyIntegral;

        // Moyer (2000), eq. 13-54
        Eigen::Matrix< ObservationScalarType, 1, 1 > observation =
                ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << static_cast< ObservationScalarType >(
                          basic_mathematics::computeModulo( rangeUnitIntegral, std::pow( 2.0, lowestRangingComponent + 6.0 ) ) ) )
                        .finished( );

        return observation;
    }

    std::shared_ptr< observation_models::MultiLegLightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculator( )
    {
        return lightTimeCalculator_;
    }

private:
    std::shared_ptr< observation_models::MultiLegLightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator_;

    // Number of link ends
    unsigned int numberOfLinkEnds_;

    // Object returning the transmitted frequency as the transmitting link end
    std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyCalculator_;

    // Function returning the turnaround ratio for given uplink and downlink bands
    std::function< double( FrequencyBands uplinkBand, FrequencyBands downlinkBand ) > turnaroundRatio_;

    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter_;

    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > stationStates_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_DSNNWAYRANGEOBSERVATIONMODEL_H
