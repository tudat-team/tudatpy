/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/basic_astro/frequencyBands.h"

#include <map>
#include <stdexcept>

namespace tudat
{

namespace basic_astrodynamics
{

std::string getFrequencyBandString( FrequencyBands frequencyBand )
{
    std::string frequencyBandString = "";
    switch( frequencyBand )
    {
        case FrequencyBands::s_band:
            frequencyBandString = "S-band";
            break;
        case FrequencyBands::x_band:
            frequencyBandString = "X-band";
            break;
        case FrequencyBands::ka_band:
            frequencyBandString = "Ka-band";
            break;
        case FrequencyBands::ku_band:
            frequencyBandString = "Ku-band";
            break;
        default:
            std::string errorMessage = "Error when getting frequency band string for " +
                    std::to_string( static_cast< int >( frequencyBand ) ) + ", frequency band not found.";
            throw std::runtime_error( errorMessage );
    }
    return frequencyBandString;
}

FrequencyBands getFrequencyBandFromString( const std::string& frequencyBand )
{
    static const std::map< std::string, FrequencyBands > frequencyBandFromStringMap = { { "S-band", FrequencyBands::s_band },
                                                                                        { "X-band", FrequencyBands::x_band },
                                                                                        { "Ka-band", FrequencyBands::ka_band },
                                                                                        { "Ku-band", FrequencyBands::ku_band } };
    const auto it = frequencyBandFromStringMap.find( frequencyBand );
    if( it == frequencyBandFromStringMap.end( ) )
        throw std::runtime_error( "Error when getting frequency band from string identifier, " + frequencyBand + " not recognised." );
    return it->second;
}

std::vector< double > convertFrequencyBandsToDoubleVector( const std::vector< FrequencyBands >& frequencyBands )
{
    std::vector< double > doubleFrequencyBands;

    for( FrequencyBands frequencyBand : frequencyBands )
    {
        doubleFrequencyBands.push_back( static_cast< double >( frequencyBand ) );
    }

    return doubleFrequencyBands;
}

std::vector< FrequencyBands > convertDoubleVectorToFrequencyBands( const std::vector< double >& frequencyBands )
{
    std::vector< FrequencyBands > frequencyBandFrequencyBands;

    for( double frequencyBand : frequencyBands )
    {
        frequencyBandFrequencyBands.push_back( static_cast< FrequencyBands >( frequencyBand ) );
    }

    return frequencyBandFrequencyBands;
}

}  // namespace basic_astrodynamics

}  // namespace tudat
