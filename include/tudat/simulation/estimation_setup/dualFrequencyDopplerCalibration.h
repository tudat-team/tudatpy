/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DUALFREQUENCYDOPPLERCALIBRATION_H
#define TUDAT_DUALFREQUENCYDOPPLERCALIBRATION_H

#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/basics/timeType.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrection.h"

namespace tudat
{
namespace observation_models
{

//! Function to compute the quantities needed for a dual-frequency (dispersive) calibration of two-way DSN Doppler data
/*!
 * Function to compute the quantities needed for a dual-frequency (dispersive) calibration of two-way DSN Doppler data
 * (receiver-tagged UTC epochs, ramped transmitter, constant receiver reference frequency). The emission interval of each count is
 * computed from the reference ephemerides in the bodies, using the non-dispersive light-time corrections provided. The ionospheric
 * correction is evaluated independently on uplink and downlink.
 * \param linkEnds Link ends of observations (transmitter, reflector1, receiver)
 * \param bodies System of bodies
 * \param receptionUtcTimes Count midpoint reception times (UTC)
 * \param observations Doppler observations [Hz], including any corrections already applied to the data
 * \param ancillarySettings Ancillary settings of observations
 * \param transmittingFrequency Ramp tables of transmitting station
 * \param nonDispersiveCorrections Non-dispersive light-time corrections (used to compute light time)
 * \param ionosphericCorrection Ionospheric light-time correction
 * \return Matrix with, for each observation (row): received frequency normalized by the carrier normalization (next column),
 * carrier normalization (mean transmitted frequency times turnaround ratio) [Hz], contribution of the ionospheric correction to
 * the normalized received frequency, squared inverse mean uplink and downlink frequencies (relative to 1 GHz), and duration
 * of emission interval [s]
 */
Eigen::Matrix< long double, Eigen::Dynamic, 6 > getDsnDopplerFrequencyCalibration(
        const LinkEnds& linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< Time >& receptionUtcTimes,
        const std::vector< long double >& observations,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings,
        const std::shared_ptr< ground_stations::StationFrequencyInterpolator >& transmittingFrequency,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& nonDispersiveCorrections,
        const std::shared_ptr< LightTimeCorrectionSettings >& ionosphericCorrection );

}  // namespace observation_models
}  // namespace tudat
#endif  // TUDAT_DUALFREQUENCYDOPPLERCALIBRATION_H
