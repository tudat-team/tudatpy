/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ONEWAYDOPPLERMEASUREDFREQUENCYPARTIAL_H
#define TUDAT_ONEWAYDOPPLERMEASUREDFREQUENCYPARTIAL_H

#include <functional>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"

namespace tudat
{

namespace observation_partials
{

/*! Calculate the factor scaling one-way Doppler partials to one-way Doppler measured frequency partials.
 *
 * Calculate the factor by which partials of the (unnormalized) one-way Doppler observable are scaled to obtain partials of the
 * one-way Doppler measured frequency observable f_R = f_T( t_T ) * ( 1 + D ). The scaling factor is the transmitted
 * frequency f_T at the transmission time t_T (converted from TDB to UTC for the frequency lookup), in the frequency band of
 * the transmitter as retrieved from the ancillary settings (see getTransmitterFrequencyBand).
 *
 * @param frequencyInterpolator Transmitted frequency interpolator of the transmitter.
 * @param timeScaleConverter Time scale converter, used to convert the transmission time from TDB to UTC.
 * @param linkEndStates List of states at each link end during observation (first entry is the transmitter).
 * @param linkEndTimes List of times (TDB) at each link end during observation (first entry is the transmission time).
 * @param ancillarySettings Observation ancillary simulation settings, from which the transmitter frequency band is retrieved.
 * @return Scaling factor from one-way Doppler partials to one-way Doppler measured frequency partials.
 */
double getOneWayDopplerMeasuredFrequencyScalingFactor(
        const std::shared_ptr< ground_stations::StationFrequencyInterpolator > frequencyInterpolator,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter,
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings );

//! Class to compute the partial derivatives of a one-way Doppler measured frequency observation.
/*!
 *  Class to compute the partial derivatives of a one-way Doppler measured frequency observation. The partials are
 *  computed from the partials of the one-way Doppler observable defined on the same link ends, scaled by the transmitted
 *  frequency of the current observation (see getOneWayDopplerMeasuredFrequencyScalingFactor).
 */
class OneWayDopplerMeasuredFrequencyPartial : public ObservationPartial< 1 >
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param oneWayDopplerPartial Partial object for the one-way Doppler observable defined on the same link ends.
     * \param scalingFactorFunction Function computing the factor scaling one-way Doppler partials to one-way Doppler measured
     *      frequency partials, as a function of the link end states, link end times and ancillary settings.
     */
    OneWayDopplerMeasuredFrequencyPartial(
            const std::shared_ptr< ObservationPartial< 1 > > oneWayDopplerPartial,
            const std::function< double( const std::vector< Eigen::Vector6d >&,
                                         const std::vector< double >&,
                                         const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ) >
                    scalingFactorFunction );

    //! Destructor
    ~OneWayDopplerMeasuredFrequencyPartial( ) = default;

    //! Function to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param ancillarySettings Observation ancillary simulation settings.
     *  \param currentObservation Value of the observation for which the partial is to be computed.
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Vector1d& currentObservation = Eigen::Vector1d::Constant( TUDAT_NAN ) ) override;

    //! Function to retrieve the partial object for the one-way Doppler observable defined on the same link ends.
    std::shared_ptr< ObservationPartial< 1 > > getOneWayDopplerPartial( )
    {
        return oneWayDopplerPartial_;
    }

protected:
    //! Partial object for the one-way Doppler observable defined on the same link ends.
    std::shared_ptr< ObservationPartial< 1 > > oneWayDopplerPartial_;

    //! Function computing the factor scaling one-way Doppler partials to one-way Doppler measured frequency partials.
    const std::function< double( const std::vector< Eigen::Vector6d >&,
                                 const std::vector< double >&,
                                 const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ) >
            scalingFactorFunction_;
};

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_ONEWAYDOPPLERMEASUREDFREQUENCYPARTIAL_H
