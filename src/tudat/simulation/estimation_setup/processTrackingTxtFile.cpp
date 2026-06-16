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
#include "tudat/astro/ground_stations/transmittingFrequencies.h"

namespace tudat
{
namespace observation_models
{

std::vector< ObservableType > findAvailableObservableTypes( const std::vector< input_output::TrackingDataType > availableDataTypes )
{
    // Initialise container for available types
    std::vector< ObservableType > availableObservableTypes;

    // Loop over map with observables and their required data types. Add observabletype to vector if those data types are present
    for( const auto& pair : observableRequiredDataTypesMap )
    {
        std::vector< input_output::TrackingDataType > requiredDataTypeSet = pair.second;
        if( utilities::containsAll( availableDataTypes, requiredDataTypeSet ) )
        {
            availableObservableTypes.push_back( pair.first );
        }
    }
    return availableObservableTypes;
}

void setStationFrequenciesFromTrackingData( const StationRampInformation& rampInformation,
                                            simulation_setup::SystemOfBodies& bodies,
                                            const std::string& stationReferenceBodyName )
{
    std::map< std::string, std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > > rampInterpolators;

    for( auto const& [ stationName, frequencyRampDataVector ] : rampInformation )
    {
        std::vector< Time > rampStartTimes;
        std::vector< Time > rampEndTimes;
        std::vector< double > rampRates;
        std::vector< double > rampStartFrequencies;

        for( auto const& currentRampData : frequencyRampDataVector )
        {
            rampStartTimes.emplace_back( currentRampData.rampUtcTimes.at( 0 ) );
            rampEndTimes.emplace_back( currentRampData.rampUtcTimes.at( currentRampData.rampUtcTimes.size( ) - 1 ) );

            if( std::adjacent_find( currentRampData.frequencyValues.begin( ),
                                    currentRampData.frequencyValues.end( ),
                                    std::not_equal_to<>( ) ) == currentRampData.frequencyValues.end( ) )
            {
                double constantTransmitterFrequency = currentRampData.frequencyValues.at( 0 );
                if( !( std::adjacent_find( currentRampData.frequencyRampRates.begin( ),
                                           currentRampData.frequencyRampRates.end( ),
                                           std::not_equal_to<>( ) ) == currentRampData.frequencyRampRates.end( ) ) )
                {
                    throw std::runtime_error(
                            "Error when reading IFMS transmitter frequencies, frequency is constant, but ramp is not constant" );
                }
                else if( currentRampData.frequencyRampRates.at( 0 ) != 0.0 && currentRampData.frequencyRampRates.at( 0 ) != -99999.999999 )
                {
                    throw std::runtime_error(
                            "Error when reading IFMS transmitter frequencies, frequency is constant, but ramp is not zero" +
                            std::to_string( currentRampData.frequencyRampRates.at( 0 ) ) );
                }
                rampRates.push_back( 0.0 );
                rampStartFrequencies.push_back( constantTransmitterFrequency );
            }
            else
            {
                std::cout << utilities::convertStlVectorToEigenVector( currentRampData.frequencyValues ).transpose( ) << std::endl;
                throw std::runtime_error( "Error when reading IFMS transmitter frequencies, only unramped data currently supported." );
            }
        }

        rampInterpolators[ stationName ] = std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >(
                rampStartTimes, rampEndTimes, rampRates, rampStartFrequencies );
    }

    if( bodies.count( stationReferenceBodyName ) == 0 )
    {
        throw std::runtime_error( "Error when setting station frequencies, reference body \"" + stationReferenceBodyName +
                                  "\" not found in system of bodies." );
    }
    auto stationReferenceBody = bodies.at( stationReferenceBodyName );

    for( auto const& [ stationName, frequencyInterpolator ] : rampInterpolators )
    {
        if( stationReferenceBody->getGroundStationMap( ).count( stationName ) == 0 )
        {
            throw std::runtime_error( "Error when setting frequencies for station " + stationName + ", station not found." );
        }

        auto groundStation = stationReferenceBody->getGroundStation( stationName );
        if( !groundStation->hasFrequencyCalculator( ) )
        {
            groundStation->setTransmittingFrequencyCalculator( frequencyInterpolator );
        }
        else if( auto existingFrequencyInterpolator = std::dynamic_pointer_cast< ground_stations::PiecewiseLinearFrequencyInterpolator >(
                         groundStation->getTransmittingFrequencyCalculator( ) ) )
        {
            existingFrequencyInterpolator->addFrequencyInterpolator( frequencyInterpolator );
        }
        else
        {
            throw std::runtime_error( "Error when adding ramp tables for station " + stationName +
                                      ", existing frequency calculator implemented, but not of correct type" );
        }
    }
}

}  // namespace observation_models
}  // namespace tudat
