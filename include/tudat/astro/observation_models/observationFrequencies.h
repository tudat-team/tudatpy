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

inline double getDsnDefaultTurnaroundRatios( FrequencyBands uplinkBand, FrequencyBands downlinkBand )
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
        throw std::runtime_error("Error when retrieving default turnaround ratios: uplink frequency band" +
            std::to_string( uplinkBand ) + "is not recognized." );
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
            std::to_string( downlinkBand ) + "is not recognized." );
    }

    return numerator / denominator;
}

inline double getCassiniTurnaroundRatio( )
{
    return 14.0 / 15.0;
}

} // namespace observation_models

} // namespace tudat

#endif //TUDAT_OBSERVATIONFREQUENCIES_H
