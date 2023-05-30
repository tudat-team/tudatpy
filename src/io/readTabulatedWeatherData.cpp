/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readTabulatedWeatherData.h"

#include <fstream>

#include <boost/algorithm/string.hpp>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace input_output
{

void DsnWeatherData::readSingleFileWeatherData( const std::string& weatherFile )
{
    // Open file
    std::ifstream stream( weatherFile, std::ios_base::binary );
    if ( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening weather file: " + weatherFile );
    }

     // Define number of characters associated with the weather data
    std::vector< int > weatherDataCharStart = { 10, 19, 28, 39, 54 };
    std::vector< int > weatherDataCharLen = { 5, 5, 6, 6, 3 };

    // Line based parsing
    std::string line;
    int year, month, day;
    while ( stream.good( ) && !stream.eof( ) )
    {
        // Get line from stream
        std::getline( stream, line );

        std::vector< std::string > vectorOfIndividualStrings;
        boost::algorithm::split( vectorOfIndividualStrings, line,
                                 boost::algorithm::is_any_of( " :" ), boost::algorithm::token_compress_on );
        // Check if last string is empty and remove it if so
        if ( vectorOfIndividualStrings.back( ).empty( ) )
        {
            vectorOfIndividualStrings.pop_back( );
        }

        // Skip empty lines
        std::string trimmedLine = line;
        boost::algorithm::trim( trimmedLine );
        if ( trimmedLine.empty( ) )
        {
            continue;
        }

        // Check if first line of day
        if ( line.substr( 0, 4 ) == "DATE" )
        {
            year = std::stoi( line.substr( 6, 2 ) );
            if ( year <= 68 )
            {
                year += 2000;
            }
            else
            {
                year += 1900;
            }
            month = std::stoi( line.substr( 8, 2 ) );
            day = std::stoi( line.substr( 10, 2 ) );

            // Initialize DSN complex number
            if ( dsnStationComplexId_ < 0 )
            {
                dsnStationComplexId_ = std::stoi( line.substr( 26, 3 ) );
            }
            else if ( dsnStationComplexId_ != std::stoi( line.substr( 26, 3 ) ) )
            {
                throw std::runtime_error( "Error when reading tabulated weather data: the DSN complex ID is inconsistent." );
            }

            // Skip following 4 lines
            for ( unsigned int i = 0; i < 4; ++i )
            {
                std::getline( stream, line );
            }
        }
        // If not first line of day
        else
        {
            int hours = std::stoi( line.substr( 1, 2 ) );
            int minutes = std::stoi( line.substr( 3, 2 ) );

            time_.push_back(
                    basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                            year, month, day, hours, minutes, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
                            physical_constants::JULIAN_DAY );

            for ( unsigned int i = 0; i < 5; ++i )
            {
                double value;
                try
                {
                    value = std::stod( line.substr( weatherDataCharStart.at( i ), weatherDataCharLen.at( i ) ) );
                }
                catch( ... )
                {
                    value = TUDAT_NAN;
                }

                switch ( i )
                {
                    // Temperatures converted from celsius to kelvin
                    case 0:
                        dewPoint_.push_back( value + 273.15 );
                        break;
                    case 1:
                        temperature_.push_back( value + 273.15 );
                        break;
                    // Pressures converted from mbar to Pa
                    case 2:
                        pressure_.push_back( value * 1e2 );
                        break;
                    case 3:
                        waterVaporPartialPressure_.push_back( value * 1e2 );
                        break;
                    // Relative humidity converted from percentage to fraction
                    case 4:
                        relativeHumidity_.push_back( value / 1e2 );
                        break;
                    default:
                        throw std::runtime_error("Invalid index when reading weather data.");
                        break;
                }
            }

        }
    }

    // Close file
    stream.close( );
}

std::map< int, std::vector< std::string > > getDefaultDsnStationNamesPerComplex( )
{
    std::map< int, std::vector< std::string > > stationsPerComplex;
    stationsPerComplex[ 10 ] = { "DSS-13", "DSS-14", "DSS-15", "DSS-24", "DSS-25", "DSS-26" };
    stationsPerComplex[ 40 ] = { "DSS-34", "DSS-35", "DSS-36", "DSS-43", "DSS-45" };
    stationsPerComplex[ 60 ] = { "DSS-54", "DSS-55", "DSS-63", "DSS-65" };

    return stationsPerComplex;
}

bool compareDsnWeatherFileStartDate( std::shared_ptr< DsnWeatherData > file1,
                                     std::shared_ptr< DsnWeatherData > file2 )
{

    if ( file1 == nullptr || file2 == nullptr )
    {
        throw std::runtime_error( "Error when comparing DSN weather files: invalid files." );
    }

    if ( file1->time_.front( ) < file2->time_.front( ) )
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::map< int, std::shared_ptr< DsnWeatherData > > readDsnWeatherDataFiles(
        const std::vector< std::string >& weatherFileNames )
{
    // Read all files and store them in vectors. One vector per DSN complex
    std::map< int, std::vector< std::shared_ptr< DsnWeatherData > > > weatherDataVectorPerComplex;
    for ( std::string weatherFileName : weatherFileNames )
    {
        std::shared_ptr< DsnWeatherData > weatherData = std::make_shared< DsnWeatherData >( weatherFileName );
        // Add data to map
        if ( !weatherDataVectorPerComplex.count( weatherData->dsnStationComplexId_ ) )
        {
            weatherDataVectorPerComplex[ weatherData->dsnStationComplexId_ ] =
                    std::vector< std::shared_ptr< DsnWeatherData > >{ weatherData };
        }
        else
        {
            weatherDataVectorPerComplex[ weatherData->dsnStationComplexId_ ].push_back( weatherData );
        }
    }

    // Merge weather files per complex and save them to map
    std::map< int, std::shared_ptr< DsnWeatherData > > weatherDataPerComplex;

    for ( auto complexIdIt = weatherDataVectorPerComplex.begin( ); complexIdIt != weatherDataVectorPerComplex.end( );
            ++complexIdIt )
    {
        // Sort weather files in each vector
        std::stable_sort( complexIdIt->second.begin( ), complexIdIt->second.end( ), &compareDsnWeatherFileStartDate );

        // Merge weather files and save them
        weatherDataPerComplex[ complexIdIt->first ] = std::make_shared< DsnWeatherData >( );
        weatherDataPerComplex[ complexIdIt->first ]->dsnStationComplexId_ = complexIdIt->first;
        for ( unsigned int i = 0; i < complexIdIt->second.size( ); ++i )
        {
            weatherDataPerComplex[ complexIdIt->first ]->fileNames_.insert(
                    weatherDataPerComplex[ complexIdIt->first ]->fileNames_.end( ),
                    complexIdIt->second.at( i )->fileNames_.begin( ),
                    complexIdIt->second.at( i )->fileNames_.end( ) );

            weatherDataPerComplex[ complexIdIt->first ]->time_.insert(
                    weatherDataPerComplex[ complexIdIt->first ]->time_.end( ),
                    complexIdIt->second.at( i )->time_.begin( ),
                    complexIdIt->second.at( i )->time_.end( ) );

            weatherDataPerComplex[ complexIdIt->first ]->dewPoint_.insert(
                    weatherDataPerComplex[ complexIdIt->first ]->dewPoint_.end( ),
                    complexIdIt->second.at( i )->dewPoint_.begin( ),
                    complexIdIt->second.at( i )->dewPoint_.end( ) );

            weatherDataPerComplex[ complexIdIt->first ]->temperature_.insert(
                    weatherDataPerComplex[ complexIdIt->first ]->temperature_.end( ),
                    complexIdIt->second.at( i )->temperature_.begin( ),
                    complexIdIt->second.at( i )->temperature_.end( ) );

            weatherDataPerComplex[ complexIdIt->first ]->pressure_.insert(
                    weatherDataPerComplex[ complexIdIt->first ]->pressure_.end( ),
                    complexIdIt->second.at( i )->pressure_.begin( ),
                    complexIdIt->second.at( i )->pressure_.end( ) );

            weatherDataPerComplex[ complexIdIt->first ]->waterVaporPartialPressure_.insert(
                    weatherDataPerComplex[ complexIdIt->first ]->waterVaporPartialPressure_.end( ),
                    complexIdIt->second.at( i )->waterVaporPartialPressure_.begin( ),
                    complexIdIt->second.at( i )->waterVaporPartialPressure_.end( ) );

            weatherDataPerComplex[ complexIdIt->first ]->relativeHumidity_.insert(
                    weatherDataPerComplex[ complexIdIt->first ]->relativeHumidity_.end( ),
                    complexIdIt->second.at( i )->relativeHumidity_.begin( ),
                    complexIdIt->second.at( i )->relativeHumidity_.end( ) );
        }
    }

    return weatherDataPerComplex;
}

std::function< double ( double ) > createInterpolatingFunction(
        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings,
        const std::vector< double >& keys,
        const std::vector< double >& values )
{
    std::vector< double > validKeys;
    std::vector< double > validValues;

    for ( unsigned int i = 0; i < keys.size( ); ++i )
    {
        if ( !std::isnan( values.at( i ) ) )
        {
            validKeys.push_back( keys.at( i ) );
            validValues.push_back( values.at( i ) );
        }
    }

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > interpolator =
            createOneDimensionalInterpolator(
                    utilities::createMapFromVectors( validKeys, validValues ), interpolatorSettings );
    std::function< double ( double ) > function = [=] ( const double time ) { return interpolator->interpolate( time ); };

    return function;
}

void setDsnWeatherDataInGroundStations(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< int, std::shared_ptr< DsnWeatherData > >& weatherDataPerComplex,
        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings,
        const std::map< int, std::vector< std::string > >& groundStationsPerComplex,
        const std::string& bodyWithGroundStations )
{
    // Loop over DSN complexes
    for ( auto weatherDataIt = weatherDataPerComplex.begin( ); weatherDataIt != weatherDataPerComplex.end( ); ++weatherDataIt )
    {
        int dsnComplex = weatherDataIt->first;
        std::shared_ptr< DsnWeatherData > weatherData = weatherDataIt->second;

        std::vector< std::string > groundStations;
        if ( groundStationsPerComplex.count( dsnComplex ) )
        {
            groundStations = groundStationsPerComplex.at( dsnComplex );
        }
        else
        {
            throw std::runtime_error( "Error when setting weather data in ground station: no ground stations in complex." );
        }

        // Create interpolating functions
        std::function< double ( double ) > dewPointFunction = createInterpolatingFunction(
                interpolatorSettings, weatherData->time_, weatherData->dewPoint_ );
        std::function< double ( double ) > temperatureFunction = createInterpolatingFunction(
                interpolatorSettings, weatherData->time_, weatherData->temperature_ );
        std::function< double ( double ) > pressureFunction = createInterpolatingFunction(
                interpolatorSettings, weatherData->time_, weatherData->pressure_ );
        std::function< double ( double ) > waterVaporPartialPressureFunction = createInterpolatingFunction(
                interpolatorSettings, weatherData->time_, weatherData->waterVaporPartialPressure_ );
        std::function< double ( double ) > relativeHumidityFunction = createInterpolatingFunction(
                interpolatorSettings, weatherData->time_, weatherData->relativeHumidity_ );

        // Set functions in ground stations
        for ( const std::string& groundStation : groundStations )
        {
            bodies.getBody( bodyWithGroundStations )->getGroundStation( groundStation )->setDewPointFunction( dewPointFunction );
            bodies.getBody( bodyWithGroundStations )->getGroundStation( groundStation )->setTemperatureFunction( temperatureFunction );
            bodies.getBody( bodyWithGroundStations )->getGroundStation( groundStation )->setPressureFunction( pressureFunction );
            bodies.getBody( bodyWithGroundStations )->getGroundStation( groundStation )->setWaterVaporPartialPressureFunction( waterVaporPartialPressureFunction );
            bodies.getBody( bodyWithGroundStations )->getGroundStation( groundStation )->setRelativeHumidityFunction( relativeHumidityFunction );
        }
    }
}

} // namespace input_output

} // namespace tudat