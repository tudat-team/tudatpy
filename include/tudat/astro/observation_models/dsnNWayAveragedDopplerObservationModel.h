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
 *          T. Moyer (2000), Formulation for Observed and Computed Values of Deep Space Network Data Types for Navigation,
 *              DEEP SPACE COMMUNICATIONS AND NAVIGATION SERIES, JPL/NASA
 */

#ifndef TUDAT_DSNNWAYAVERAGEDDOPPLEROBSERVATIONMODEL_H
#define TUDAT_DSNNWAYAVERAGEDDOPPLEROBSERVATIONMODEL_H

#include <stdexcept>
#include <string>

#include "tudat/simulation/simulation.h"

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationFrequencies.h"
#include "tudat/astro/observation_models/nWayRangeObservationModel.h"

namespace tudat
{

namespace observation_models
{

/*! Calculate the scaling factor for computing partials via DifferencedObservablePartial.
 *
 * Calculate the scaling factor for computing partials via DifferencedObservablePartial, for DSN n-way Doppler
 * observations. The scaling factor are selected according to eq. 13-59 of Moyer(2000).
 * The scaling factor depends on whether it is the first or the second partial being calculated. In
 * DifferencedObservablePartial, the first partial is multiplied by -1, hence here corresponds to the start of the
 * integration interval. The second partial is multiplied by +1, hence here corresponds to the end of the
 * integration interval.
 *
 * @param bodies System of bodies
 * @param linkEnds Map of the linkEnds defining the observation model
 * @param referenceLinkEnd Link end at which given time is valid, i.e. link end for which associated time
 *      is kept constant (to input value)
 * @param linkEndStates List of states at each link end during observation.
 * @param linkEndTimes List of times at each link end during observation.
 * @param ancillarySettings Observation ancillary simulation settings.
 * @param isFirstPartial Boolean indicating whether the scaling factor should be computed for the first (true) or
 *      second (false) partial.
 * @return Scaling factor
 */
inline double getDsnNWayAveragedDopplerScalingFactor(
        const std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > receivedFrequencyFunction,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings,
        const bool isFirstPartial )
{
    double integrationTime;
    std::vector< FrequencyBands > frequencyBands;
    try
    {
        integrationTime = ancillarySettings->getAncilliaryDoubleData( doppler_integration_time );
        frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncilliaryDoubleVectorData( frequency_bands ) );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error(
                "Error when retrieving integration ancillary settings for DSN N-way averaged Doppler observable: " +
                std::string( caughtException.what( ) ) );
    }

    double transmissionTime;
    if ( referenceLinkEnd == receiver )
    {
        if ( isFirstPartial )
        {
            transmissionTime = linkEndTimes.at( 0 );
        }
        else
        {
            transmissionTime = linkEndTimes.at( 4 );
        }
    }
//    else if ( referenceLinkEnd == transmitter )
//    {
//        if ( isFirstPartial )
//        {
//            transmissionTime = linkEndTimes.at( 3 );
//        }
//        else
//        {
//            transmissionTime = linkEndTimes.at( 7 );
//        }
//    }
    else
    {
        throw std::runtime_error(
                "Error when getting DSN N-way Doppler partials scaling factor: the selected reference link end (" +
                getLinkEndTypeString( referenceLinkEnd ) + ") is not valid." );
    }

    double frequency = receivedFrequencyFunction( frequencyBands, transmissionTime );

    // Moyer (2000), eq. 13-59
    return frequency / integrationTime / physical_constants::getSpeedOfLight< double >( );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
class DsnNWayAveragedDopplerObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
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
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    DsnNWayAveragedDopplerObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel,
            const std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel,
            const std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyCalculator,
            const std::function< double ( observation_models::FrequencyBands uplinkBand,
                    observation_models::FrequencyBands downlinkBand ) >& turnaroundRatio,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( dsn_n_way_averaged_doppler , linkEnds, observationBiasCalculator),
        arcStartObservationModel_( arcStartObservationModel ),
        arcEndObservationModel_( arcEndObservationModel ),
        numberOfLinkEnds_( linkEnds.size( ) ),
        transmittingFrequencyCalculator_( transmittingFrequencyCalculator ),
        turnaroundRatio_( turnaroundRatio )
    {
        if( !std::is_same< Time, TimeType >::value )
        {
            std::cerr<<
                    "Warning when defining DSN N-way averaged Doppler observation model: the selected time type "
                    "is not valid, using it would lead to large numerical errors."<<std::endl;
        }

        if ( numberOfLinkEnds_ != 3 )
        {
            throw std::runtime_error(
                    "Error when defining DSN N-way averaged Doppler observation model: model allows exactly 3 link ends, " +
                    std::to_string( numberOfLinkEnds_ ) + "were selected.");
        }
    }

    //! Destructor
    ~DsnNWayAveragedDopplerObservationModel( ) { }

    /*! Function to compute DSN n-way Doppler observation at given time.
     *
     * Function to compute DSN n-way Doppler observation at given time. Only implemented for receiver as the
     * linkEndAssociatedWithTime. Computes the observable according to section 13.3.2.2 of Moyer (2000).
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
        if ( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error(
                "Error when computing DSN N-way Doppler observables: the selected reference link end (" +
                getLinkEndTypeString( linkEndAssociatedWithTime ) + ") is not valid." );
        }
        // Check if ancillary settings were provided
        if( ancillarySettings == nullptr )
        {
            throw std::runtime_error(
                    "Error when simulating n-way DSN averaged Doppler observable; no ancillary settings found. " );
        }

        std::vector< double > arcStartLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > arcStartLinkEndStates;
        std::vector< double > arcEndLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > arcEndLinkEndStates;

        TimeType integrationTime;
        ObservationScalarType referenceFrequency;
        std::vector< FrequencyBands > frequencyBands;
        FrequencyBands referenceUplinkBand;
        try
        {
            integrationTime = ancillarySettings->getAncilliaryDoubleData( doppler_integration_time );
            referenceFrequency = ancillarySettings->getAncilliaryDoubleData( doppler_reference_frequency );
            frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncilliaryDoubleVectorData( frequency_bands ) );
            referenceUplinkBand = convertDoubleToFrequencyBand( ancillarySettings->getAncilliaryDoubleData( reception_reference_frequency_band ) );
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error(
                    "Error when retrieving ancillary settings for DSN N-way averaged Doppler observable: " +
                    std::string( caughtException.what( ) ) );
        }

        if ( frequencyBands.size( ) != numberOfLinkEnds_ - 1 )
        {
            throw std::runtime_error(
                    "Error when retrieving frequency bands ancillary settings for DSN N-way averaged Doppler observable: "
                    "size (" + std::to_string( frequencyBands.size( ) ) + ") is inconsistent with number of links (" +
                    std::to_string( numberOfLinkEnds_ - 1 ) + ")." );
        }
        FrequencyBands uplinkBand = frequencyBands.at( 0 );
        FrequencyBands downlinkBand = frequencyBands.at( 1 );
        TimeType receptionStartTime = time - integrationTime / 2.0;
        TimeType receptionEndTime = time + integrationTime / 2.0;

        TimeType startLightTime = arcStartObservationModel_->computeIdealObservationsWithLinkEndData(
                receptionStartTime, linkEndAssociatedWithTime, arcStartLinkEndTimes, arcStartLinkEndStates,
                ancillarySettings )( 0, 0 ) / physical_constants::getSpeedOfLight< ObservationScalarType >( );
        TimeType endLightTime = arcEndObservationModel_->computeIdealObservationsWithLinkEndData(
                receptionEndTime, linkEndAssociatedWithTime, arcEndLinkEndTimes, arcEndLinkEndStates,
                ancillarySettings )( 0, 0 ) / physical_constants::getSpeedOfLight< ObservationScalarType >( );

        // Moyer (2000), eqs. 13-52 and 13-53
        TimeType transmissionStartTime = receptionStartTime - startLightTime;
        TimeType transmissionEndTime = receptionEndTime - endLightTime;

        ObservationScalarType transmitterFrequencyIntegral =
                transmittingFrequencyCalculator_->template getTemplatedFrequencyIntegral< ObservationScalarType, TimeType >(
                        transmissionStartTime, transmissionEndTime );

        // Moyer (2000), eq. 13-54
        Eigen::Matrix< ObservationScalarType, 1, 1 > observation = ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) <<
                turnaroundRatio_( referenceUplinkBand, downlinkBand ) * referenceFrequency -
                turnaroundRatio_( uplinkBand, downlinkBand ) / static_cast< ObservationScalarType >( integrationTime ) *
                transmitterFrequencyIntegral ).finished( );

        linkEndTimes.clear( );
        linkEndStates.clear( );
        linkEndTimes.resize( 4 * ( numberOfLinkEnds_ - 1 ) );
        linkEndStates.resize( 4 * ( numberOfLinkEnds_ - 1 ) );

        for( unsigned int i = 0; i < 2 * ( numberOfLinkEnds_ - 1 ) ; i++ )
        {
            linkEndTimes[ i ] = arcStartLinkEndTimes[ i ];
            linkEndTimes[ i + 2 * ( numberOfLinkEnds_ - 1 ) ] = arcEndLinkEndTimes[ i ];

            linkEndStates[ i ] = arcStartLinkEndStates[ i ];
            linkEndStates[ i + 2 * ( numberOfLinkEnds_ - 1 ) ] = arcEndLinkEndStates[ i ];
        }

        return observation;
    }

    // Function to retrieve the arc end observation model
    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > getArcEndObservationModel( )
    {
        return arcEndObservationModel_;
    }

    // Function to retrieve the arc start observation model
    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > getArcStartObservationModel( )
    {
        return arcStartObservationModel_;
    }

private:

    // N-way range observation model associated with the start of the Doppler integration time.
    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel_;

    // N-way range observation model associated with the end of the Doppler integration time.
    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel_;

    // Number of link ends
    unsigned int numberOfLinkEnds_;

    // Object returning the transmitted frequency as the transmitting link end
    std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyCalculator_;

    // Function returning the turnaround ratio for given uplink and downlink bands
    std::function< double ( FrequencyBands uplinkBand, FrequencyBands downlinkBand ) > turnaroundRatio_;
};



} // namespace observation_models

} // namespace tudat


#endif //TUDAT_DSNNWAYAVERAGEDDOPPLEROBSERVATIONMODEL_H
