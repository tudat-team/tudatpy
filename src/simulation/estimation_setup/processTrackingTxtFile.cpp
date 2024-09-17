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

void ProcessedTrackingTxtFileContents::updateObservations()
{
    // Update the observableTypes that one can expect to process
    observationMap_.clear();

    for (const ObservableType observableType : observableTypes_)
    {
        std::vector<double> observableValues;

        std::cout<<"Getting observable type "<<observableType<<std::endl;
        // Convert the raw data to required observables
        // TODO: This function should set up the ancilliary data
        switch (observableType)
        {
        case n_way_range:
        {
            // Conversion function for n-way range
            auto lightTimeRangeConversion = [](double lightTime, double lightTimeDelay) {
            return (lightTime - lightTimeDelay) * physical_constants::SPEED_OF_LIGHT; };

            // Extract columns from raw data and convert to observable values
            std::vector<double> lightTimes = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::n_way_light_time);
            std::vector<double> lightTimeDelays = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::light_time_measurement_delay, 0.0);
            observableValues = utilities::convertVectors(lightTimeRangeConversion, lightTimes, lightTimeDelays);
            break;
        }
        case doppler_measured_frequency:
        {
            // Conversion function for doppler measured frequency
            auto dopplerFrequencyConversion = [](double dopplerFrequency, double dopplerBaseFrequency) {
            return dopplerFrequency + dopplerBaseFrequency; };

            // Extract columns from raw data and convert to observable values
            std::vector<double> dopplerFrequencies = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::doppler_measured_frequency);
            std::vector<double> dopplerBaseFrequencies = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::doppler_base_frequency);

            // Check if any of the base frequencies are zero
            if (std::any_of(dopplerBaseFrequencies.begin(), dopplerBaseFrequencies.end(), [](double baseFrequency) { return baseFrequency == 0.0; }))
            {
                std::cerr << "Warning when processing doppler_measured_frequency. Doppler base frequency is zero." << std::endl;
            }

            observableValues = utilities::convertVectors(dopplerFrequencyConversion, dopplerFrequencies, dopplerBaseFrequencies);
            break;
        }
        case dsn_n_way_averaged_doppler:
        {
            observableValues = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::doppler_averaged_frequency );
            break;
        }
        default: {
        throw std::runtime_error("Error while processing tracking txt file. ObservableType conversion not implemented");
        }
    }

    // Store observables
    observationMap_[observableType] = observableValues;
    }

}

void ProcessedTrackingTxtFileContents::updateObservationTimes()
{
    // Clear any previous values
    observationTimes_.clear();

    // Get data map and time representation
    const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();
    TimeRepresentation timeRepresentation = getTimeRepresentation();

    // Depending on the time representation, convert further to tdb seconds since j2000
    switch (timeRepresentation)
    {
        case tdb_seconds_j2000:
        {
            observationTimes_ = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::tdb_reception_time_j2000);
            break;
        }
        case utc_seconds_j2000:
        {
            std::vector<double> observationTimesUtc = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::utc_reception_time_j2000);
            observationTimes_ = convertTimesTdbFromJ2000(observationTimesUtc, observation_models::receiver );
            break;
        }
        case calendar_day_time:
        {
            // Convert dates to Julian days since J2000
            auto years = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::year);
            auto months = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::month);
            auto days = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::day);
            auto hours = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::hour);
            auto minutes = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::minute);
            auto seconds = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::second);
            // Convert to seconds and add to utc times
            std::vector<double> observationTimesUtc;
            for( unsigned int i = 0; i < years.size( ); i++ )
            {
                observationTimesUtc.push_back(basic_astrodynamics::DateTime(
                years.at( i ), months.at( i ), days.at( i ), hours.at( i ), minutes.at( i ), seconds.at( i ) ).epoch< double >() );
            }
            // Convert to TDB
            observationTimes_ = convertTimesTdbFromJ2000(observationTimesUtc, observation_models::receiver );

            break;
        }
        // Throw error if representation not implemented
        default:
        {
            throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
        }
    }

    // Get the delays in the time tag (or set to 0.0 if not specified)
    std::vector<double> timeTagDelays = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::time_tag_delay, 0.0);
    for (size_t idx = 0; idx < observationTimes_.size(); ++idx) {
    observationTimes_[idx] -= timeTagDelays[idx];
    }
}

std::vector<double> ProcessedTrackingTxtFileContents::convertTimesTdbFromJ2000(
    std::vector<double> observationTimesUtc, const observation_models::LinkEndType referenceLinkEnd )
{
    // Get the timescale converter

    // Check if there is one LinkEnds per observation time
    if (linkEndsVector_.size() != observationTimesUtc.size())
    {
        throw std::runtime_error("Error while processing tracking data: vector of linkEnds and observationTimes not of equal size");
    }

    // Ge a vector of ground station positions
    std::vector<Eigen::Vector3d> groundStationPositions;
    for (const auto& linkEnds : linkEndsVector_)
    {
        try
        {
            std::string currentGroundStation = linkEnds.at( referenceLinkEnd ).getStationName( );
            groundStationPositions.push_back( earthFixedGroundStationPositions_.at( currentGroundStation ));
        }
        catch( const std::runtime_error& error )
        {
            throw std::runtime_error( "Error when creating getting link ends for time conversion in tracking file processing: "
                + std::string( error.what( ) ) + ". For reference link end " + std::to_string( static_cast< int >( referenceLinkEnd ) ) );
        }
    }

    // Convert to TDB using the GS positions
    std::vector<double> observationTimesTdb = timeScaleConverter_->getCurrentTimes(basic_astrodynamics::utc_scale,
                                                                                 basic_astrodynamics::tdb_scale,
                                                                                 observationTimesUtc,
                                                                                 groundStationPositions);
    return observationTimesTdb;
}

void ProcessedTrackingTxtFileContents::updateLinkEnds()
{
    // Clear any previous values
    linkEndsVector_.clear();

    // Get information from raw data file
    const auto& metaDataStrMap = rawTrackingTxtFileContents_->getMetaDataStrMap();
    const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();

    // Deduce linkends representation
    LinkEndsRepresentation linkEndsRepresentation = getLinkEndsRepresentation();

    // Create a vector of LinkEnds based on how they are represented
    // This currently only implements the DSN transmitter and receiver
    switch (linkEndsRepresentation)
    {
        case dsn_transmitting_receiving_station_nr:
        {
            const auto& dsnTransmitterIds = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::dsn_transmitting_station_nr);
            const auto& dsnReceiverIds = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::dsn_receiving_station_nr);

            for (size_t i = 0; i < numDataRows; ++i)
            {
                std::string transmitterName = getStationNameFromStationId(dsnTransmitterIds.at(i));
                std::string receiverName = getStationNameFromStationId(dsnReceiverIds.at(i));
                LinkEnds currentLinkEnds{
                    {transmitter, LinkEndId("Earth", transmitterName)},
                    {reflector, LinkEndId(spacecraftName_, "")},
                    {receiver, LinkEndId("Earth", receiverName)} };
                linkEndsVector_.push_back(currentLinkEnds);
            }
            break;
        }

        case transmitting_receiving_station_name:
        {
            std::string receivingStationName = rawTrackingTxtFileContents_->getMetaDataStrMap().at(input_output::TrackingDataType::receiving_station_name);
            std::string transmittingStationName = rawTrackingTxtFileContents_->getMetaDataStrMap().at(input_output::TrackingDataType::transmitting_station_name);

            for (size_t i = 0; i < numDataRows; ++i)
            {
                LinkEnds currentLinkEnds{
                    {transmitter, LinkEndId("Earth", transmittingStationName)},
                    {reflector, LinkEndId(spacecraftName_, "")},
                    {receiver, LinkEndId("Earth", receivingStationName)},
                    };
                linkEndsVector_.push_back(currentLinkEnds);
            }
            break;
        }

        // Throw error if representation not implemented
        default: {
        throw std::runtime_error("Error while processing tracking txt file: LinkEnds representation not recognised or implemented.");
        }
    }

    // Creating a set with all the distinct LinkEnds
    linkEndsSet_ = utilities::vectorToSet(linkEndsVector_);
}

ProcessedTrackingTxtFileContents::TimeRepresentation ProcessedTrackingTxtFileContents::getTimeRepresentation()
{

    // Get all the data types from the raw file contents
    auto const& availableDataTypes = rawTrackingTxtFileContents_->getDataColumnTypes();

    // Return representation based on available data types

    if (utilities::containsAll(availableDataTypes, std::vector<input_output::TrackingDataType>{input_output::TrackingDataType::tdb_reception_time_j2000}))
    {
        return tdb_seconds_j2000;
    }
    if (utilities::containsAll(availableDataTypes, std::vector<input_output::TrackingDataType>{input_output::TrackingDataType::utc_reception_time_j2000}))
    {
        return utc_seconds_j2000;
    }

    if (utilities::containsAll(availableDataTypes,
                 std::vector<input_output::TrackingDataType>{
                     input_output::TrackingDataType::year,
                     input_output::TrackingDataType::month,
                     input_output::TrackingDataType::day,
                     input_output::TrackingDataType::hour,
                     input_output::TrackingDataType::minute,
                     input_output::TrackingDataType::second }))
    {
        return calendar_day_time;
    }

    // Throw an error if no match is found
    throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
}

ProcessedTrackingTxtFileContents::LinkEndsRepresentation ProcessedTrackingTxtFileContents::getLinkEndsRepresentation( )
{
    // Get all the available data columns
    auto const& availableDataTypes = rawTrackingTxtFileContents_->getAllAvailableDataTypes();

    // Porvide link Ends representation based on available columns

    if (utilities::containsAll(availableDataTypes,
                         std::vector<input_output::TrackingDataType>{
                             input_output::TrackingDataType::dsn_transmitting_station_nr,
                             input_output::TrackingDataType::dsn_receiving_station_nr } ) )
    {
        return dsn_transmitting_receiving_station_nr;
    }

    if (utilities::containsAll(availableDataTypes,
                           std::vector<input_output::TrackingDataType>{
                               input_output::TrackingDataType::transmitting_station_name,
                               input_output::TrackingDataType::receiving_station_name } ) )
    {
        return transmitting_receiving_station_name;
    }

    // Throw error if no match is found
    throw std::runtime_error("Error while processing tracking txt file: Link Ends representation not recognised or implemented.");
}


double ProcessedTrackingTxtFileContents::getObservationTimeStep( )
{
    if( observationTimes_.size( ) < 2 )
    {
        throw std::runtime_error( "Error when getting integration time for processed file contents, size is < 2" );
    }
    double observationTimeStep = observationTimes_.at( 1 ) - observationTimes_.at( 0 );
    for( unsigned int i = 1; i < observationTimes_.size( ); i++ )
    {
        double testObservationTimeStep = observationTimes_.at( i ) - observationTimes_.at( i - 1 );
        if( std::fabs( observationTimeStep - testObservationTimeStep ) > 50.0 * std::numeric_limits< double >::epsilon( ) * observationTimes_.at( i - 1 )  )
        {
            std::cout<<i<<" "<<testObservationTimeStep<<" "<<observationTimeStep<<" "<<testObservationTimeStep - observationTimeStep<<" "<<
                                                                                                                                         50.0 * std::numeric_limits< double >::epsilon( ) * observationTimes_.at( i - 1 )<<std::endl;
            throw std::runtime_error( "Error when getting integration time for processed file contents, step is not equal" );
        }
    }

    return observationTimeStep;
}

void setStationFrequenciesFromTrackingData(
    const std::map< std::string, std::vector< std::tuple< std::vector< double >, std::vector< double >, std::vector< double > > > >& rampInformation,
    simulation_setup::SystemOfBodies& bodies )
{
    std::map< std::string, std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > > rampInterpolators;

    for( auto it : rampInformation )
    {
        std::vector< double > rampStartTimes;
        std::vector< double > rampEndTimes;
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
                std::cout<<utilities::convertStlVectorToEigenVector( currentFrequencyValues )<<std::endl;
                throw std::runtime_error( "Error when reading IFMS transmitter frequencies, only unramped data currently supported." );
            }
        }
        for( unsigned int i = 0; i < rampStartTimes.size( ) - 1 ; i++ )
        {
            double timeDifference = rampStartTimes.at( i + 1 ) - rampEndTimes.at( i );
            rampStartTimes[ i + 1 ] -= timeDifference/2.0;
            rampEndTimes[ i ] += timeDifference/2.0;

        }
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

void setStationFrequenciesFromTrackingData(
    const std::vector< std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents > > processedTrackingTxtFileContents,
    simulation_setup::SystemOfBodies& bodies )
{
    std::map< std::string, std::vector< std::tuple< std::vector< double >, std::vector< double >, std::vector< double > > > > rampInformation;
    for( unsigned int i = 0; i < processedTrackingTxtFileContents.size( ); i++ )
    {
        std::vector<double> observationTimes = processedTrackingTxtFileContents.at( i )->getObservationTimes( );
        std::shared_ptr<input_output::TrackingTxtFileContents> fileContents = processedTrackingTxtFileContents.at( i )->getRawTrackingTxtFileContents( );
        std::vector< double > rampUtcTimes = fileContents->getDoubleDataColumn( input_output::TrackingDataType::utc_ramp_referencee_j2000 );
        std::vector< double > frequencyRampRates = fileContents->getDoubleDataColumn( input_output::TrackingDataType::transmission_frequency_linear_term );
        std::vector< double > frequencyValues = fileContents->getDoubleDataColumn( input_output::TrackingDataType::transmission_frequency_constant_term );

        std::string transmitterName;
        if( processedTrackingTxtFileContents.at( i )->getLinkEndsSet( ).size( ) != 1 )
        {
            throw std::runtime_error( "Error when getting link ends from IFMS file, found multiple link ends sets." +
                                          std::to_string( processedTrackingTxtFileContents.at( i )->getLinkEndsSet( ).size( ) ) );
        }
        else
        {
            LinkEnds currentLinkEnds = *(processedTrackingTxtFileContents.at( i )->getLinkEndsSet( ).begin( ) );
            transmitterName = currentLinkEnds.at( transmitter ).stationName_;
        }
        rampInformation[ transmitterName ].push_back( std::make_tuple( rampUtcTimes, frequencyValues, frequencyRampRates ) );
    }
    setStationFrequenciesFromTrackingData( rampInformation, bodies );
}


void setTrackingDataInformationInBodies(
    const std::vector< std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents> > processedTrackingTxtFileContents,
    simulation_setup::SystemOfBodies& bodies,
    ObservableType observableType )
{
    for( unsigned int i = 0; i < processedTrackingTxtFileContents.size( ); i++ )
    {
        std::vector<ObservableType> availableObservableTypes = processedTrackingTxtFileContents.at( i )->getObservableTypes( );
        if ( !utilities::containsAll( availableObservableTypes, std::vector< ObservableType >( { observableType } ) ) )
        {
            throw std::runtime_error(
                "Error while processing Tracking txt file for body properties. Observable not found: " + getObservableName( observableType ) );
        }
    }


    switch( observableType )
    {
    case dsn_n_way_averaged_doppler:
    {
        setStationFrequenciesFromTrackingData( processedTrackingTxtFileContents, bodies );
    }
    default:
        break;
    }

}


} // namespace observation_models
} // namespace tudat
