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
    std::ifstream stream( weatherFile, std::ios_base::in );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening weather file: " + weatherFile );
    }

    // Define number of characters associated with the weather data
    std::vector< int > weatherDataCharStart = { 10, 19, 28, 39, 54 };
    std::vector< int > weatherDataCharLen = { 5, 5, 6, 6, 3 };

    // Line based parsing
    std::string line;
    int year = -0, month = -1, day = -1;

    double currentTime = TUDAT_NAN;
    Eigen::VectorXd currentMeteoData = Eigen::VectorXd::Zero( 5 );

    while( stream.good( ) && !stream.eof( ) )
    {
        // Get line from stream
        std::getline( stream, line );

        std::vector< std::string > vectorOfIndividualStrings;
        boost::algorithm::split(
                vectorOfIndividualStrings, line, boost::algorithm::is_any_of( " :" ), boost::algorithm::token_compress_on );
        // Check if last string is empty and remove it if so
        if( vectorOfIndividualStrings.back( ).empty( ) )
        {
            vectorOfIndividualStrings.pop_back( );
        }

        // Skip empty lines
        std::string trimmedLine = line;
        boost::algorithm::trim( trimmedLine );
        if( trimmedLine.empty( ) )
        {
            continue;
        }

        // Check if first line of day
        if( line.substr( 0, 4 ) == "DATE" )
        {
            year = std::stoi( line.substr( 6, 2 ) );
            if( year <= 68 )
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
            if( dsnStationComplexId_ < 0 )
            {
                dsnStationComplexId_ = std::stoi( line.substr( 26, 3 ) );
            }
            else if( dsnStationComplexId_ != std::stoi( line.substr( 26, 3 ) ) )
            {
                throw std::runtime_error( "Error when reading tabulated weather data: the DSN complex ID is inconsistent." );
            }

            // Skip following 4 lines
            for( unsigned int i = 0; i < 4; ++i )
            {
                std::getline( stream, line );
            }
        }
        // If not first line of day
        else
        {
            int hours = std::stoi( line.substr( 1, 2 ) );
            int minutes = std::stoi( line.substr( 3, 2 ) );

            currentTime = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                                  year, month, day, hours, minutes, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
                    physical_constants::JULIAN_DAY;

            currentMeteoData.setZero( );
            for( unsigned int i = 0; i < 5; ++i )
            {
                double value;
                // Try to convert the characters associated with the current data field to a double
                try
                {
                    value = std::stod( line.substr( weatherDataCharStart.at( i ), weatherDataCharLen.at( i ) ) );
                }
                // If that fails, set data field to NAN
                catch( ... )
                {
                    value = TUDAT_NAN;
                }

                switch( i )
                {
                    // Temperature converted from celsius to kelvin
                    case 0:
                        currentMeteoData( 0 ) = value + 273.15;
                        break;
                    // Temperature converted from celsius to kelvin
                    case 1:
                        currentMeteoData( 1 ) = value + 273.15;
                        break;
                    // Pressure converted from mbar to Pa
                    case 2:
                        currentMeteoData( 2 ) = value * 1e2;
                        break;
                    // Pressure converted from mbar to Pa
                    case 3:
                        currentMeteoData( 3 ) = value * 1e2;
                        break;
                    // Relative humidity converted from percentage to fraction
                    case 4:
                        currentMeteoData( 4 ) = value / 1e2;
                        break;
                    default:
                        throw std::runtime_error( "Invalid index when reading weather data." );
                        break;
                }
            }
            meteoDataMap_[ currentTime ] = currentMeteoData;
        }
    }

    // Close file
    stream.close( );
}

void EstrackWeatherData::readSingleWeatherDataFile( const std::string& weatherFile )
{
    // Open file
    std::ifstream stream( weatherFile, std::ios_base::in );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening ESTRACK file: " + weatherFile );
    }

    // Line based parsing
    std::string currentLine;
    std::vector< std::string > currentSplitLine;

    std::map< double, Eigen::VectorXd > meteoMap;
    while( stream.good( ) && !stream.eof( ) )
    {
        // Get line from stream
        std::getline( stream, currentLine );

        if( currentLine.size( ) > 0 )
        {
            // Trim input string (removes all leading and trailing whitespaces).
            boost::algorithm::trim( currentLine );
            boost::algorithm::split(
                    currentSplitLine, currentLine, boost::algorithm::is_any_of( " \t" ), boost::algorithm::token_compress_on );
            std::string utcString = currentSplitLine.at( 1 );
            double utc = basic_astrodynamics::dateTimeFromIsoString( utcString ).epoch< double >( );
            double humidity = std::stod( currentSplitLine.at( 4 ) ) / 100.0;
            double pressure = std::stod( currentSplitLine.at( 5 ) ) * 100.0;
            double temperature = std::stod( currentSplitLine.at( 6 ) ) + 273.15;
            meteoMap[ utc ] = ( Eigen::VectorXd( 3 ) << humidity, pressure, temperature ).finished( );
        }
    }
    meteoDataPerFile_.push_back( meteoMap );
}

bool compareDsnWeatherFileStartDate( std::shared_ptr< DsnWeatherData > file1, std::shared_ptr< DsnWeatherData > file2 )
{
    if( file1 == nullptr || file2 == nullptr )
    {
        throw std::runtime_error( "Error when comparing DSN weather files: invalid files." );
    }

    if( file1->meteoDataMap_.begin( )->first < file2->meteoDataMap_.begin( )->first )
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::map< int, std::shared_ptr< DsnWeatherData > > readDsnWeatherDataFiles( const std::vector< std::string >& weatherFileNames )
{
    // Read all files and store them in vectors. One vector per DSN complex
    std::map< int, std::vector< std::shared_ptr< DsnWeatherData > > > weatherDataVectorPerComplex;
    for( std::string weatherFileName: weatherFileNames )
    {
        std::shared_ptr< DsnWeatherData > weatherData = std::make_shared< DsnWeatherData >( weatherFileName );
        // Add data to map
        if( !weatherDataVectorPerComplex.count( weatherData->dsnStationComplexId_ ) )
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

    for( auto complexIdIt = weatherDataVectorPerComplex.begin( ); complexIdIt != weatherDataVectorPerComplex.end( ); ++complexIdIt )
    {
        // Sort weather files in each vector
        std::stable_sort( complexIdIt->second.begin( ), complexIdIt->second.end( ), &compareDsnWeatherFileStartDate );

        // Merge weather files and save them
        weatherDataPerComplex[ complexIdIt->first ] = std::make_shared< DsnWeatherData >( );
        weatherDataPerComplex[ complexIdIt->first ]->dsnStationComplexId_ = complexIdIt->first;
        for( unsigned int i = 0; i < complexIdIt->second.size( ); ++i )
        {
            weatherDataPerComplex[ complexIdIt->first ]->fileNames_.insert( weatherDataPerComplex[ complexIdIt->first ]->fileNames_.end( ),
                                                                            complexIdIt->second.at( i )->fileNames_.begin( ),
                                                                            complexIdIt->second.at( i )->fileNames_.end( ) );

            weatherDataPerComplex[ complexIdIt->first ]->meteoDataMap_.insert( complexIdIt->second.at( i )->meteoDataMap_.begin( ),
                                                                               complexIdIt->second.at( i )->meteoDataMap_.end( ) );
        }
    }

    return weatherDataPerComplex;
}

std::function< double( double ) > createInterpolatingFunction( std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings,
                                                               const std::vector< double >& keys,
                                                               const std::vector< double >& values )
{
    std::vector< double > validKeys;
    std::vector< double > validValues;

    for( unsigned int i = 0; i < keys.size( ); ++i )
    {
        if( !std::isnan( values.at( i ) ) )
        {
            validKeys.push_back( keys.at( i ) );
            validValues.push_back( values.at( i ) );
        }
    }

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > interpolator =
            createOneDimensionalInterpolator( utilities::createMapFromVectors( validKeys, validValues ), interpolatorSettings );
    std::function< double( double ) > function = [ = ]( const double time ) { return interpolator->interpolate( time ); };

    return function;
}

void setDsnWeatherDataInGroundStations( simulation_setup::SystemOfBodies& bodies,
                                        const std::map< int, std::shared_ptr< DsnWeatherData > >& weatherDataPerComplex,
                                        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings,
                                        const std::map< int, std::vector< std::string > >& groundStationsPerComplex,
                                        const std::string& bodyWithGroundStations )
{
    // Loop over DSN complexes
    for( auto weatherDataIt = weatherDataPerComplex.begin( ); weatherDataIt != weatherDataPerComplex.end( ); ++weatherDataIt )
    {
        int dsnComplex = weatherDataIt->first;
        std::shared_ptr< DsnWeatherData > weatherData = weatherDataIt->second;

        std::vector< std::string > groundStations;
        if( groundStationsPerComplex.count( dsnComplex ) )
        {
            groundStations = groundStationsPerComplex.at( dsnComplex );
        }
        else
        {
            throw std::runtime_error( "Error when setting weather data in ground station: no ground stations in complex." );
        }

        std::map< ground_stations::MeteoDataEntries, int > dsnMeteoEntries = { { ground_stations::temperature_meteo_data, 1 },
                                                                               { ground_stations::pressure_meteo_data, 2 },
                                                                               { ground_stations::water_vapor_pressure_meteo_data, 3 },
                                                                               { ground_stations::relative_humidity_meteo_data, 4 },
                                                                               { ground_stations::dew_point_meteo_data, 0 } };

        std::shared_ptr< ground_stations::StationMeteoData > meteoData =
                std::make_shared< ground_stations::ContinuousInterpolatedMeteoData >(
                        interpolators::createOneDimensionalInterpolator( weatherData->meteoDataMap_, interpolatorSettings ),
                        dsnMeteoEntries );
        // Set functions in ground stations
        for( const std::string& groundStation: groundStations )
        {
            bodies.getBody( bodyWithGroundStations )->getGroundStation( groundStation )->setMeteoData( meteoData );
        }
    }
}

void setEstrackWeatherDataInGroundStation( simulation_setup::SystemOfBodies& bodies,
                                           const std::vector< std::string >& weatherFiles,
                                           const std::string groundStation,
                                           std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings,
                                           const std::string& bodyWithGroundStations )
{
    EstrackWeatherData weatherData = EstrackWeatherData( weatherFiles );
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > > interpolators;
    auto meteoDataList = weatherData.getMeteoDataPerFile( );
    for( unsigned int i = 0; i < meteoDataList.size( ); i++ )
    {
        interpolators.push_back( interpolators::createOneDimensionalInterpolator( meteoDataList.at( i ), interpolatorSettings ) );
    }

    std::map< ground_stations::MeteoDataEntries, int > estrackMeteoEntries = { { ground_stations::temperature_meteo_data, 2 },
                                                                               { ground_stations::pressure_meteo_data, 1 },
                                                                               { ground_stations::relative_humidity_meteo_data, 0 } };

    std::shared_ptr< ground_stations::StationMeteoData > meteoData =
            std::make_shared< ground_stations::PiecewiseInterpolatedMeteoData >( interpolators, estrackMeteoEntries );
    bodies.at( bodyWithGroundStations )->getGroundStation( groundStation )->setMeteoData( meteoData );
}

}  // namespace input_output

}  // namespace tudat