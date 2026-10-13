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
 *    T. Moyer (2000), Formulation for Observed and Computed Values of Deep Space Network Data Types for Navigation,
 *              DEEP SPACE COMMUNICATIONS AND NAVIGATION SERIES, JPL/NASA
 */

#ifndef TUDAT_TRANSPONDERTURNAROUNDRATIOS_H
#define TUDAT_TRANSPONDERTURNAROUNDRATIOS_H

#include "tudat/astro/basic_astro/frequencyBands.h"

namespace tudat
{

namespace system_models
{

/*!
 * Returns the default turnaround ratios used for DSN spacecraft, according to Moyer (2000), table 13-1.
 *
 * @param uplinkBand Uplink band
 * @param downlinkBand Downlink band
 * @return Turnaround ratio
 */
double getDsnDefaultTurnaroundRatios( basic_astrodynamics::FrequencyBands uplinkBand, basic_astrodynamics::FrequencyBands downlinkBand );

/*!
 * Returns the Cassini turnaround ratio for Ka uplink and downlink bands, according to Moyer (2000), section 13.2.2.
 */
double getCassiniKaBandTurnaroundRatio( );

/*!
 * Returns the Cassini turnaround ratio as a function of the uplink and downlink band. Returns the Cassini-specific
 * turnaround ratio for Ka band uplink and downlink, and the default DSN turnaround ratios for other bands, according
 * to Moyer (2000), section 13.2.2.
 *
 * @param uplinkBand Uplink band
 * @param downlinkBand Downlink band
 * @return Turnaround ratio
 */
double getCassiniTurnaroundRatio( basic_astrodynamics::FrequencyBands uplinkBand, basic_astrodynamics::FrequencyBands downlinkBand );

}  // namespace system_models

}  // namespace tudat

#endif  // TUDAT_TRANSPONDERTURNAROUNDRATIOS_H
