/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

std::string getFrequencyBandString( FrequencyBands frequencyBand );

double getDsnDefaultTurnaroundRatios( FrequencyBands uplinkBand, FrequencyBands downlinkBand );

double getCassiniTurnaroundRatio( );

double getCassiniTurnaroundRatio( FrequencyBands uplinkBand, FrequencyBands downlinkBand );

std::vector< double > convertFrequencyBandsToDoubleVector( const std::vector< FrequencyBands >& frequencyBands );

std::vector< FrequencyBands > convertDoubleVectorToFrequencyBands( const std::vector< double >& frequencyBands );

} // namespace observation_models

} // namespace tudat

#endif //TUDAT_OBSERVATIONFREQUENCIES_H
