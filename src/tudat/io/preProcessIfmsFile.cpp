/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/preProcessIfmsFile.h"

namespace tudat
{

namespace input_output
{

data::PlainLinkDefinition getIfmsLinkEnds( const std::string& spacecraftName,
                                           const std::string& earthName,
                                           const std::string& groundStationName )
{
    return { { { earthName, groundStationName }, "transmitter" },
             { { spacecraftName, "" }, "reflector_1" },
             { { earthName, groundStationName }, "receiver" } };
}

double getIfmsFrequencyRampRate( const double frequencyRampRate )
{
    return frequencyRampRate == -99999.999999 ? 0.0 : frequencyRampRate;
}

void addIfmsFrequencyRamps( const std::shared_ptr< TrackingTxtFileContents >& rawIfmsData,
                            const std::shared_ptr< data::RampedFrequencySupplementaryData >& frequencySupplementaryData )
{
    const std::vector< double > rampUtcTimes = rawIfmsData->getDoubleDataColumn( TrackingDataType::utc_ramp_referencee_j2000 );
    const std::vector< double > frequencyValues =
            rawIfmsData->getDoubleDataColumn( TrackingDataType::transmission_frequency_constant_term );
    const std::vector< double > frequencyRampRates =
            rawIfmsData->getDoubleDataColumn( TrackingDataType::transmission_frequency_linear_term );

    if( rampUtcTimes.empty( ) )
    {
        return;
    }

    unsigned int currentRampStartIndex = 0;
    double currentFrequency = frequencyValues.front( );
    double currentFrequencyRampRate = getIfmsFrequencyRampRate( frequencyRampRates.front( ) );

    for( unsigned int i = 1; i < rampUtcTimes.size( ); ++i )
    {
        const double testFrequency = frequencyValues.at( i );
        const double testFrequencyRampRate = getIfmsFrequencyRampRate( frequencyRampRates.at( i ) );

        if( testFrequency != currentFrequency || testFrequencyRampRate != currentFrequencyRampRate )
        {
            frequencySupplementaryData->addFrequencyRamp(
                    rampUtcTimes.at( currentRampStartIndex ), rampUtcTimes.at( i - 1 ), currentFrequency, currentFrequencyRampRate );
            currentRampStartIndex = i;
            currentFrequency = testFrequency;
            currentFrequencyRampRate = testFrequencyRampRate;
        }
    }

    frequencySupplementaryData->addFrequencyRamp(
            rampUtcTimes.at( currentRampStartIndex ), rampUtcTimes.back( ), currentFrequency, currentFrequencyRampRate );
}

}  // namespace input_output

}  // namespace tudat
