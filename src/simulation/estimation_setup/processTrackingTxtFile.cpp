/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>

#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"

namespace tudat
{
namespace observation_models
{

std::vector<ObservableType> findAvailableObservableTypes(const std::vector<input_output::TrackingDataType> availableDataTypes)
{
    // Initialise container for available types
    std::vector<ObservableType> availableObservableTypes;

    // Loop over map with observables and their required data types. Add observabletype to vector if those data types are present
    for (const auto& pair : observableRequiredDataTypesMap)
    {
        std::vector<input_output::TrackingDataType> requiredDataTypeSet = pair.second;
        if (utilities::containsAll(availableDataTypes, requiredDataTypeSet))
        {
            availableObservableTypes.push_back(pair.first);
        }
    }
    return availableObservableTypes;
}

void setStationFrequenciesFromTrackingData(
    const std::map< std::string, std::vector< std::tuple< std::vector< double >, std::vector< double >, std::vector< double > > > >& rampInformation,
    simulation_setup::SystemOfBodies& bodies )
{
    std::map< std::string, std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > > rampInterpolators;

    for( auto it : rampInformation )
    {
        std::vector< Time > rampStartTimes;
        std::vector< Time > rampEndTimes;
        std::vector< double > rampRates;
        std::vector< double > rampStartFrequencies;

        for( unsigned int i = 0; i < it.second.size( ); i++ )
        {
            std::vector< double > currentRampUtcTimes = std::get< 0 >( it.second.at( i ) );
            std::vector< double > currentFrequencyValues = std::get< 1 >( it.second.at( i ) );
            std::vector< double > currentFrequencyRampRates = std::get< 2 >( it.second.at( i ) );

            rampStartTimes.push_back( currentRampUtcTimes.at( 0 ) );
            rampEndTimes.push_back( currentRampUtcTimes.at( currentRampUtcTimes.size( ) - 1 ) );

            if ( std::adjacent_find( currentFrequencyValues.begin(), currentFrequencyValues.end(), std::not_equal_to<>() ) == currentFrequencyValues.end() )
            {
                double constantTransmitterFrequency = currentFrequencyValues.at( 0 );
                if( !( std::adjacent_find( currentFrequencyRampRates.begin(), currentFrequencyRampRates.end(), std::not_equal_to<>() ) == currentFrequencyRampRates.end() ) )
                {
                    throw std::runtime_error( "Error when reading IFMS transmitter frequencies, frequency is constant, but ramp is not constant" );
                }
                else if( currentFrequencyRampRates.at( 0 ) != 0.0 && currentFrequencyRampRates.at( 0 ) != -99999.999999 )
                {
                    throw std::runtime_error( "Error when reading IFMS transmitter frequencies, frequency is constant, but ramp is not zero" + std::to_string( currentFrequencyRampRates.at( 0 ) ) );
                }
                rampRates.push_back( 0.0 );
                rampStartFrequencies.push_back( constantTransmitterFrequency );
            }
            else
            {
                std::cout<<utilities::convertStlVectorToEigenVector( currentFrequencyValues ).transpose( )<<std::endl;
                throw std::runtime_error( "Error when reading IFMS transmitter frequencies, only unramped data currently supported." );
            }
        }
        rampStartTimes[ 0 ] -= 1.0;
        for( unsigned int i = 0; i < rampStartTimes.size( ) - 1 ; i++ )
        {
            double timeDifference = rampStartTimes.at( i + 1 ) - rampEndTimes.at( i );
            rampStartTimes[ i + 1 ] -= timeDifference/2.0;
            rampEndTimes[ i ] += timeDifference/2.0;
        }
        rampEndTimes[ rampEndTimes.size( ) - 1 ] += 1.0;
        rampInterpolators[ it.first ] = std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >(
            rampStartTimes, rampEndTimes, rampRates, rampStartFrequencies );
    }

    for( auto it = rampInterpolators.begin( ); it != rampInterpolators.end( ); it++ )
    {
        if( bodies.at( "Earth" )->getGroundStationMap( ).count( it->first ) == 0 )
        {
            throw std::runtime_error( "Error when setting frequencies for station " + it->first + ", station not found." );
        }
        bodies.at( "Earth" )->getGroundStation( it->first )->setTransmittingFrequencyCalculator( it->second );
    }

}




} // namespace observation_models
} // namespace tudat
