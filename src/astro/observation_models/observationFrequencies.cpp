/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include "tudat/astro/observation_models/observationFrequencies.h"

#include <iostream>

namespace tudat
{

namespace observation_models
{

std::string getFrequencyBandString( FrequencyBands frequencyBand )
{
    std::string frequencyBandString = "";
    switch( frequencyBand )
    {
    case s_band:
        frequencyBandString = "S-band";
        break;
    case x_band:
        frequencyBandString = "X-band";
        break;
    case ka_band:
        frequencyBandString = "Ka-band";
        break;
    case ku_band:
        frequencyBandString = "Ku-band";
        break;
    default:
        std::string errorMessage = "Error when getting link end string for type " +
                std::to_string( frequencyBand ) + ", type not found.";
        throw std::runtime_error( errorMessage );
    }
    return frequencyBandString;
}

double getDsnDefaultTurnaroundRatios( FrequencyBands uplinkBand, FrequencyBands downlinkBand )
{
    double numerator, denominator;

    if ( uplinkBand == s_band )
    {
        denominator = 221.0;
    }
    else if ( uplinkBand == x_band )
    {
        denominator = 749.0;
    }
    else if ( uplinkBand == ka_band )
    {
        denominator = 3599.0;
    }
    else
    {
        throw std::runtime_error("Error when retrieving default turnaround ratios: uplink frequency band " +
            getFrequencyBandString( uplinkBand ) + " is not recognized." );
    }

    if ( downlinkBand == s_band )
    {
        numerator = 240.0;
    }
    else if ( downlinkBand == x_band )
    {
        numerator = 880.0;
    }
    else if ( downlinkBand == ka_band )
    {
        numerator = 3344.0;
    }
    else
    {
        throw std::runtime_error("Error when retrieving default turnaround ratios: downlink frequency band" +
            getFrequencyBandString( downlinkBand ) + "is not recognized." );
    }

    return numerator / denominator;
}

double getCassiniTurnaroundRatio( )
{
    return 14.0 / 15.0;
}

double getCassiniTurnaroundRatio( FrequencyBands uplinkBand, FrequencyBands downlinkBand )
{
    return getCassiniTurnaroundRatio( );
}

std::vector< double > convertFrequencyBandsToDoubleVector( const std::vector< FrequencyBands >& frequencyBands )
{
    std::vector< double > doubleFrequencyBands;

    for ( FrequencyBands frequencyBand : frequencyBands )
    {
        doubleFrequencyBands.push_back( frequencyBand );
    }

    for ( unsigned int i = 0; i < doubleFrequencyBands.size(); ++i )
    {
        std::cerr << doubleFrequencyBands.at(i) << " ";
    }

    return doubleFrequencyBands;
}

std::vector< FrequencyBands > convertDoubleVectorToFrequencyBands( const std::vector< double >& frequencyBands )
{
    std::vector< FrequencyBands > frequencyBandFrequencyBands;

    for ( double frequencyBand : frequencyBands )
    {
        frequencyBandFrequencyBands.push_back( static_cast< FrequencyBands >( frequencyBand ) );
    }

    return frequencyBandFrequencyBands;
}

} // namespace observation_models

} // namespace tudat