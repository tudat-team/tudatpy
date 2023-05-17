/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          820-013 TRK-2-24, Weather Data Interface, Revision A (2006), DSN/JPL
 */

#ifndef TUDAT_READTABULATEDWEATHERDATA_H
#define TUDAT_READTABULATEDWEATHERDATA_H

#include <vector>
#include <string>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/simulation/environment_setup.h"

namespace tudat
{

namespace input_output
{

class DsnWeatherData
{
public:

    DsnWeatherData( const std::string& weatherFile ):
        dsnStationComplexId_( -1 )
    {
        fileNames_.push_back( weatherFile );
        readSingleFileWeatherData( weatherFile );
    }

    DsnWeatherData( ):
        dsnStationComplexId_( -1 )
    { }

    int dsnStationComplexId_;

    std::vector< std::string > fileNames_;

    // Time since J2000 [s]
    std::vector< double > time_;

    // Dew point [K]
    std::vector< double > dewPoint_;

    // Temperature [K]
    std::vector< double > temperature_;

    // Pressure [Pa]
    std::vector< double > pressure_;

    // Water vapor partial pressure [Pa]
    std::vector< double > waterVaporPartialPressure_;

    // Relative humidity [-] (defined in [0,1])
    std::vector< double > relativeHumidity_;

private:

    void readSingleFileWeatherData( const std::string& weatherFile );
};

std::map< int, std::vector< std::string > > getDefaultDsnStationNamesPerComplex( );

bool compareDsnWeatherFileStartDate( std::shared_ptr< DsnWeatherData > file1,
                                     std::shared_ptr< DsnWeatherData > file2 );

std::map< int, std::shared_ptr< DsnWeatherData > > readDsnWeatherDataFiles(
        const std::vector< std::string >& weatherFileNames );

std::function< double ( double ) > createInterpolatingFunction(
        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings,
        const std::vector< double >& keys,
        const std::vector< double >& values );

void setDsnWeatherDataInGroundStations(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< int, std::shared_ptr< DsnWeatherData > >& weatherDataPerComplex,
        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings,
        const std::map< int, std::vector< std::string > >& groundStationsPerComplex,
        const std::string& bodyWithGroundStations );

inline void setDsnWeatherDataInGroundStations(
        simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::string >& weatherFiles,
        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = interpolators::linearInterpolation( ),
        const std::map< int, std::vector< std::string > >& groundStationsPerComplex = getDefaultDsnStationNamesPerComplex( ),
        const std::string& bodyWithGroundStations = "Earth" )
{
    setDsnWeatherDataInGroundStations(
            bodies, readDsnWeatherDataFiles( weatherFiles ), interpolatorSettings, groundStationsPerComplex,
            bodyWithGroundStations );
}

} // namespace input_output

} // namespace tudat

#endif //TUDAT_READTABULATEDWEATHERDATA_H
