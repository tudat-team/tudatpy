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
 *          T. Moyer (2003), Formulation for Observed and Computed Values of Deep Space Network Data
 * Types for Navigation, DEEP SPACE COMMUNICATIONS AND NAVIGATION SERIES, JPL/NASA
 */

#ifndef TUDAT_DSNNWAYRANGEPARTIAL_H
#define TUDAT_DSNNWAYRANGEPARTIAL_H

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/astro/observation_models/observationFrequencies.h"
#include "tudat/astro/orbit_determination/observation_partials/nWayRangePartial.h"

namespace tudat
{

namespace observation_partials
{

/*! Calculate the factor scaling n-way range partials to DSN n-way range partials.
 *
 * Calculate the factor by which partials of the n-way range observable (in meters) are scaled to obtain
 * partials of the DSN n-way (sequential) range observable (in range units), for a reception-time-referenced
 * observation. Following from Moyer (2003), eq. 13-128, the scaling factor is
 *
 *      C_b * f_T( t_1 ) / c,
 *
 * with C_b the uplink-band-dependent range unit conversion factor, f_T the frequency transmitted by the
 * transmitting ground station, t_1 the transmission time (converted to UTC for the frequency lookup) and c the
 * speed of light.
 *
 * @param transmittedFrequencyFunction Function returning the frequency transmitted by the transmitting ground
 *      station as a function of frequency bands (unused for the transmitted frequency) and time (UTC).
 * @param timeScaleConverter Time scale converter, used to convert the transmission time from TDB to UTC.
 * @param linkEndTimes List of times (TDB) at each link end during observation, as returned by the DSN n-way
 *      range observation model (first entry is the transmission time).
 * @param ancillarySettings Observation ancillary simulation settings, from which the frequency bands are retrieved.
 * @return Scaling factor from n-way range partials (meters) to DSN n-way range partials (range units).
 */
double getDsnNWayRangeScalingFactor(
        const std::function< double( std::vector< observation_models::FrequencyBands >, double ) > transmittedFrequencyFunction,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings );

//! Class to compute the partial derivatives of a DSN n-way (sequential) range observation.
/*!
 *  Class to compute the partial derivatives of a DSN n-way (sequential) range observation. The partials are
 *  computed as the partials of the corresponding n-way range observable (defined on the same link ends and
 *  light time calculator), scaled by the factor converting a range change in meters to a change of the
 *  observable in range units (see getDsnNWayRangeScalingFactor).
 */
class DsnNWayRangePartial : public ObservationPartial< 1 >
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param nWayRangePartial Partial object for the n-way range observable defined on the same link ends.
     * \param scalingFactorFunction Function computing the factor scaling n-way range partials (in meters) to DSN
     *      n-way range partials (in range units), as a function of the link end times and ancillary settings.
     */
    DsnNWayRangePartial(
            const std::shared_ptr< NWayRangePartial > nWayRangePartial,
            const std::function< double( const std::vector< double >&,
                                         const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ) >
                    scalingFactorFunction );

    //! Destructor
    ~DsnNWayRangePartial( ) = default;

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

protected:
    //! Partial object for the n-way range observable defined on the same link ends.
    std::shared_ptr< NWayRangePartial > nWayRangePartial_;

    //! Function computing the factor scaling n-way range partials to DSN n-way range partials.
    const std::function< double( const std::vector< double >&,
                                 const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ) >
            scalingFactorFunction_;
};

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_DSNNWAYRANGEPARTIAL_H
