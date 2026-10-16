/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FREQUENCYBANDS_H
#define TUDAT_FREQUENCYBANDS_H

#include <string>
#include <vector>

namespace tudat
{

namespace basic_astrodynamics
{

enum class FrequencyBands { s_band, x_band, ka_band, ku_band };

/*!
 * Function returning a string with the name of the frequency band.
 *
 * @param frequencyBand Frequency band
 * @return Name of the frequency band.
 */
std::string getFrequencyBandString( FrequencyBands frequencyBand );

/*!
 * Function returning a frequency band from a string identifier
 *
 * @param frequencyBand String identifier for frequency band
 * @return Frequency band
 */
FrequencyBands getFrequencyBandFromString( const std::string& frequencyBand );

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

/*!
 * Returns the name of the DSN uplink frequency band that contains a given (sky-level) uplink frequency, according to the
 * DSN Earth-to-space frequency allocations for deep space (Category B) and near-Earth (Category A) missions (DSN
 * Telecommunications Link Design Handbook 810-005, module 201): S-band (2025-2120 MHz), X-band (7145-7235 MHz) and
 * Ka-band (34200-34700 MHz).
 *
 * @param uplinkFrequency Uplink frequency in Hz
 * @return Name of the frequency band (as used by getFrequencyBandFromString), or an empty string if the frequency is not
 * in any of the DSN uplink allocations.
 */
inline std::string getDsnUplinkFrequencyBandName( const double uplinkFrequency )
{
    if( uplinkFrequency >= 2025.0E6 && uplinkFrequency <= 2120.0E6 )
    {
        return "S-band";
    }
    else if( uplinkFrequency >= 7145.0E6 && uplinkFrequency <= 7235.0E6 )
    {
        return "X-band";
    }
    else if( uplinkFrequency >= 34200.0E6 && uplinkFrequency <= 34700.0E6 )
    {
        return "Ka-band";
    }
    return "";
}

}  // namespace basic_astrodynamics

}  // namespace tudat

#endif  // TUDAT_FREQUENCYBANDS_H
