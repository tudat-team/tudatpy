/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/system_models/transponderTurnaroundRatios.h"

#include <stdexcept>

namespace tudat
{

namespace system_models
{

using basic_astrodynamics::FrequencyBands;

double getDsnDefaultTurnaroundRatios( FrequencyBands uplinkBand, FrequencyBands downlinkBand )
{
    double numerator;
    double denominator;

    if( uplinkBand == FrequencyBands::s_band )
    {
        denominator = 221.0;
    }
    else if( uplinkBand == FrequencyBands::x_band )
    {
        denominator = 749.0;
    }
    else if( uplinkBand == FrequencyBands::ka_band )
    {
        denominator = 3599.0;
    }
    else
    {
        throw std::runtime_error( "Error when retrieving default turnaround ratios: uplink frequency band " +
                                  basic_astrodynamics::getFrequencyBandString( uplinkBand ) + " is not recognized." );
    }

    if( downlinkBand == FrequencyBands::s_band )
    {
        numerator = 240.0;
    }
    else if( downlinkBand == FrequencyBands::x_band )
    {
        numerator = 880.0;
    }
    else if( downlinkBand == FrequencyBands::ka_band )
    {
        numerator = 3344.0;
    }
    else
    {
        throw std::runtime_error( "Error when retrieving default turnaround ratios: downlink frequency band" +
                                  basic_astrodynamics::getFrequencyBandString( downlinkBand ) + "is not recognized." );
    }

    return numerator / denominator;
}

double getCassiniKaBandTurnaroundRatio( )
{
    return 14.0 / 15.0;
}

double getCassiniTurnaroundRatio( FrequencyBands uplinkBand, FrequencyBands downlinkBand )
{
    if( uplinkBand == FrequencyBands::ka_band && downlinkBand == FrequencyBands::ka_band )
    {
        return getCassiniKaBandTurnaroundRatio( );
    }
    else
    {
        return getDsnDefaultTurnaroundRatios( uplinkBand, downlinkBand );
    }
}

}  // namespace system_models

}  // namespace tudat
