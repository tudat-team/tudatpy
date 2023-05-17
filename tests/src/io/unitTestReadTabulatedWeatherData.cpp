/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include "tudat/basics/testMacros.h"

#include "tudat/io/readTabulatedWeatherData.h"

#include "tudat/simulation/environment_setup.h"

namespace tudat
{
namespace unit_tests
{

using namespace input_output;

BOOST_AUTO_TEST_SUITE( test_read_tabulated_weather_data )

BOOST_AUTO_TEST_CASE( readWeatherData )
{

    // Reading data without gaps
    {
        std::shared_ptr< DsnWeatherData > weatherFile = std::make_shared< DsnWeatherData >(
                "/Users/pipas/Documents/mro-data/wea/mromagr20170012017365_10.wea.txt" );

        BOOST_CHECK_EQUAL ( weatherFile->dsnStationComplexId_, 10 );

        // Comparison with values from file
        // Values of time obtained from https://nsidc.org/data/icesat/glas-date-conversion-tool/date_convert/

        double tolerance = 1e-13; // Corresponds to numerical precision for temperature values

        // DATE: 170101 DOY: 001 DSS 10, TIME 0000
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->time_.at( 0 ), 536500800.0, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->dewPoint_.at( 0 ) - 273.15, 2.8, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->temperature_.at( 0 ) - 273.15, 10.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->pressure_.at( 0 ) * 1e-2, 897.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->waterVaporPartialPressure_.at( 0 ) * 1e-2, 7.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->relativeHumidity_.at( 0 ) * 100, 58, tolerance );

        // DATE: 170101 DOY: 001 DSS 10, TIME 0030
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->time_.at( 1 ), 536502600.0, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->dewPoint_.at( 1 ) - 273.15, 1.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->temperature_.at( 1 ) - 273.15, 9.9, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->pressure_.at( 1 ) * 1e-2, 896.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->waterVaporPartialPressure_.at( 1 ) * 1e-2, 6.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->relativeHumidity_.at( 1 ) * 100, 54, tolerance );

        // DATE: 170101 DOY: 001 DSS 10, TIME 2359
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->time_.at( 48 ), 536587140.0, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->dewPoint_.at( 48 ) - 273.15, 1.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->temperature_.at( 48 ) - 273.15, 10.5, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->pressure_.at( 48 ) * 1e-2, 894.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->waterVaporPartialPressure_.at( 48 ) * 1e-2, 6.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->relativeHumidity_.at( 48 ) * 100, 52, tolerance );

        // DATE: 170102 DOY: 002 DSS 10, TIME 0000
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->time_.at( 49 ), 536587200.0, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->dewPoint_.at( 49 ) - 273.15, 0.9, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->temperature_.at( 49 ) - 273.15, 10.5, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->pressure_.at( 49 ) * 1e-2, 894.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->waterVaporPartialPressure_.at( 49 ) * 1e-2, 6.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->relativeHumidity_.at( 49 ) * 100, 51, tolerance );
    }

    // Reading data with gaps
    {
        std::shared_ptr< DsnWeatherData > weatherFile = std::make_shared< DsnWeatherData >(
            "/Users/pipas/Documents/mro-data/wea/mromagr20060012006365_10.wea.txt" );

        BOOST_CHECK_EQUAL ( weatherFile->dsnStationComplexId_, 10 );

        // Comparison with values from file
        // Values of time obtained from https://nsidc.org/data/icesat/glas-date-conversion-tool/date_convert/

        double tolerance = 1e-13; // Corresponds to numerical precision for temperature values

        // DATE: 060109 DOY: 009 DSS 10, TIME 0334
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->time_.at( 344 ), 190049640.0, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->dewPoint_.at( 344 ) - 273.15, -7.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->temperature_.at( 344 ) - 273.15, 11.7, tolerance );
        BOOST_CHECK_EQUAL ( std::isnan( weatherFile->pressure_.at( 344 ) ), true );
        BOOST_CHECK_CLOSE_FRACTION ( weatherFile->waterVaporPartialPressure_.at( 344 ) * 1e-2, 3.6, tolerance );
        BOOST_CHECK_EQUAL ( std::isnan( weatherFile->relativeHumidity_.at( 344 ) ), true );
    }
}

// Check that date is set correctly into the ground stations
// Also checks that files are ordered correctly
BOOST_AUTO_TEST_CASE( setWeatherData )
{
    double tolerance = 1e-13; // Corresponds to numerical precision for temperature values

    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodiesToCreate = { "Earth" };
    simulation_setup::BodyListSettings bodySettings = simulation_setup::getDefaultBodySettings( bodiesToCreate );
    bodySettings.at( "Earth" )->groundStationSettings = simulation_setup::getDsnStationSettings( );
    simulation_setup::SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    setDsnWeatherDataInGroundStations(
            bodies,
            std::vector< std::string >
                    { "/Users/pipas/Documents/mro-data/wea/mromagr20180012018365_10.wea.txt",
                      "/Users/pipas/Documents/mro-data/wea/mromagr20170012017365_60.wea.txt",
                      "/Users/pipas/Documents/mro-data/wea/mromagr20170012017365_10.wea.txt",
                      "/Users/pipas/Documents/mro-data/wea/mromagr20170012017365_40.wea.txt",
                      "/Users/pipas/Documents/mro-data/wea/mromagr20160012016366_10.wea.txt" } );

    std::vector< std::string > complex10GroundStations = { "DSS-13", "DSS-14", "DSS-15", "DSS-24", "DSS-25", "DSS-26" };
    for ( std::string groundStation : complex10GroundStations )
    {
        std::shared_ptr< ground_stations::GroundStation > gs = bodies.getBody( "Earth" )->getGroundStation( groundStation );
        // DATE: 170101 DOY: 001 DSS 10, TIME 0000
        double time = 536500800.0;
        BOOST_CHECK_CLOSE_FRACTION ( gs->getDewPointFunction( )( time ) - 273.15, 2.8, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getTemperatureFunction( )( time ) - 273.15, 10.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getPressureFunction( )( time ) * 1e-2, 897.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 7.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getRelativeHumidityFunction( )( time ) * 100, 58, tolerance );

        // DATE: 170101 DOY: 001 DSS 10, TIME 0030
        time = 536502600.0;
        BOOST_CHECK_CLOSE_FRACTION ( gs->getDewPointFunction( )( time ) - 273.15, 1.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getTemperatureFunction( )( time ) - 273.15, 9.9, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getPressureFunction( )( time ) * 1e-2, 896.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 6.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getRelativeHumidityFunction( )( time ) * 100, 54, tolerance );

        // DATE: 170101 DOY: 001 DSS 10, TIME 2359
        // Why do some tolerances need larger values here? Loss of significant digits in the linear interpolation?
        time = 536587140.0;
        BOOST_CHECK_CLOSE_FRACTION ( gs->getDewPointFunction( )( time ) - 273.15, 1.1, 1e-9 );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getTemperatureFunction( )( time ) - 273.15, 10.5, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getPressureFunction( )( time ) * 1e-2, 894.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 6.7, 1e-10 );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getRelativeHumidityFunction( )( time ) * 100, 52, 1e-10 );

        // DATE: 170102 DOY: 002 DSS 10, TIME 0000
        time = 536587200.0;
        BOOST_CHECK_CLOSE_FRACTION ( gs->getDewPointFunction( )( time ) - 273.15, 0.9, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getTemperatureFunction( )( time ) - 273.15, 10.5, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getPressureFunction( )( time ) * 1e-2, 894.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 6.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION ( gs->getRelativeHumidityFunction( )( time ) * 100, 51, tolerance );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat