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

}  // namespace basic_astrodynamics

}  // namespace tudat

#endif  // TUDAT_FREQUENCYBANDS_H
