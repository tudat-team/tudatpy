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


#ifndef TUDAT_OBSERVATIONFREQUENCIES_H
#define TUDAT_OBSERVATIONFREQUENCIES_H

#include <string>
#include <vector>

namespace tudat
{

namespace observation_models
{

enum FrequencyBands
{
    s_band,
    x_band,
    ka_band,
    ku_band
};

/*!
 * Function returning a string with the name of the frequency band.
 *
 * @param frequencyBand Frequency band
 * @return Name of the frequency band.
 */
std::string getFrequencyBandString( FrequencyBands frequencyBand );

/*!
 * Returns the default turnaround ratios used for DSN spacecraft, according to Moyer (2000), table 13-1.
 *
 * @param uplinkBand Uplink band
 * @param downlinkBand Downlink band
 * @return Turnaround ratio
 */
double getDsnDefaultTurnaroundRatios( FrequencyBands uplinkBand, FrequencyBands downlinkBand );

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
double getCassiniTurnaroundRatio( FrequencyBands uplinkBand, FrequencyBands downlinkBand );

/*!
 * Converts a vector of frequency bands to the corresponding vector of doubles, using the correspondence between each
 * enum and an integer.
 *
 * @param frequencyBands Vector of frequency bands
 * @return Vector of doubles
 */
std::vector< double > convertFrequencyBandsToDoubleVector( const std::vector< FrequencyBands >& frequencyBands );

inline double convertFrequencyBandToDouble( const FrequencyBands& frequencyBand )
{
    return convertFrequencyBandsToDoubleVector( { frequencyBand } ).front( );
}

/*!
 * Converts a vector of doubles to the corresponding vector of frequency bands, using the correspondence between each
 * enum and an integer.
 *
 * @param frequencyBands Vector of doubles
 * @return Vector of frequency bands
 */
std::vector< FrequencyBands > convertDoubleVectorToFrequencyBands( const std::vector< double >& frequencyBands );

inline FrequencyBands convertDoubleToFrequencyBand( const double frequencyBand )
{
    return convertDoubleVectorToFrequencyBands( { frequencyBand } ).front( );
}

} // namespace observation_models

} // namespace tudat

#endif //TUDAT_OBSERVATIONFREQUENCIES_H
