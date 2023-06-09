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

// Class containing weather data for a single DSN stations complex. Data is distributed in TRK-2-24 files.
class DsnWeatherData
{
public:

    /*!
     * Constructor. Reads weather file and saves the data.
     *
     * @param weatherFile File name.
     */
    DsnWeatherData( const std::string& weatherFile ):
        dsnStationComplexId_( -1 )
    {
        fileNames_.push_back( weatherFile );
        readSingleFileWeatherData( weatherFile );
    }

    /*!
     * Constructor.
     */
    DsnWeatherData( ):
        dsnStationComplexId_( -1 )
    { }

    // Number of the DSN station complex
    int dsnStationComplexId_;

    // Names of the files from which the data originaets
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

    /*!
     * Extracts the data in the read file and places it in the appropriate vectors. Missing measurements are set to NAN.
     * Data is extracted according to TRK-2-24 (2006).
     *
     * @param weatherFile File name.
     */
    void readSingleFileWeatherData( const std::string& weatherFile );
};

/*!
 * Returns the default DSN station names per DSN station complex id. Stations are named as "DSS-i", following the
 * nomenclature used when retrieving the default DSN ground station settings.
 */
std::map< int, std::vector< std::string > > getDefaultDsnStationNamesPerComplex( );

/*!
 * Checks which file starts first. Used to sort weather files. Returns true if file1 starts first, false otherwise.
 *
 * @param file1 Weather data file.
 * @param file2 Weather data file.
 * @return
 */
bool compareDsnWeatherFileStartDate( std::shared_ptr< DsnWeatherData > file1,
                                     std::shared_ptr< DsnWeatherData > file2 );

/*!
 * Reads multiple DSN weather files. Merges the data associated with each DSN complex.
 *
 * @param weatherFileNames Vector with weather file names.
 * @return Map with a single DsnWeatherData object per DSN complex id.
 */
std::map< int, std::shared_ptr< DsnWeatherData > > readDsnWeatherDataFiles(
        const std::vector< std::string >& weatherFileNames );

/*!
 * Creates interpolation function with the specified settings, keys and values. The only difference with respect to a manual
 * creation of the interpolator is that this function first removes all NAN values and respective keys from the provided
 * vectors.
 *
 * @param interpolatorSettings Interpolator settings.
 * @param keys Vector with interpolation keys (time).
 * @param values Vector with interpolation values (e.g. temperature, pressure, etc.)
 * @return Value as a function of the key.
 */
std::function< double ( double ) > createInterpolatingFunction(
        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings,
        const std::vector< double >& keys,
        const std::vector< double >& values );

/*!
 * Sets the functions to compute the weather data variables (pressure, temperature, etc.) as a function of time in the
 * ground station objects of the system of bodies.
 *
 * @param bodies System of bodies.
 * @param weatherDataPerComplex Map with the weather data per DSN complex.
 * @param interpolatorSettings Interpolator settings to use when interpolating the weather data.
 * @param groundStationsPerComplex Map containing the names of the DSN stations per DSN complex id.
 * @param bodyWithGroundStations Name of the body with the ground stations.
 */
void setDsnWeatherDataInGroundStations(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< int, std::shared_ptr< DsnWeatherData > >& weatherDataPerComplex,
        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings,
        const std::map< int, std::vector< std::string > >& groundStationsPerComplex,
        const std::string& bodyWithGroundStations );

/*!
 * Sets the functions to compute the weather data variables (pressure, temperature, etc.) as a function of time in the
 * ground station objects of the system of bodies.
 *
 * @param bodies System of bodies.
 * @param weatherFiles Vector with the names of the weather files to use. They don't need to have any specific order.
 * @param interpolatorSettings Interpolator settings to use when interpolating the weather data.
 * @param groundStationsPerComplex Map containing the names of the DSN stations per DSN complex id.
 * @param bodyWithGroundStations Name of the body with the ground stations.
 */
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
