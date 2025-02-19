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
#include "tudat/simulation/estimation_setup/createLightTimeCorrection.h"
namespace tudat
{
namespace unit_tests
{

using namespace tudat;

using namespace input_output;
using namespace observation_models;

BOOST_AUTO_TEST_SUITE( test_read_tabulated_weather_data )

BOOST_AUTO_TEST_CASE( readWeatherData )
{
    // Reading data without gaps
    {
        std::shared_ptr< DsnWeatherData > weatherFile =
                std::make_shared< DsnWeatherData >( tudat::paths::getTudatTestDataPath( ) + "mromagr20170012017365_10.wea.txt" );

        BOOST_CHECK_EQUAL( weatherFile->dsnStationComplexId_, 10 );

        // Comparison with values from file
        // Values of time obtained from https://nsidc.org/data/icesat/glas-date-conversion-tool/date_convert/

        double tolerance = 1e-13;  // Corresponds to numerical precision for temperature values

        // DATE: 170101 DOY: 001 DSS 10, TIME 0000
        double currentTime = weatherFile->meteoDataMap_.begin( )->first;
        BOOST_CHECK_CLOSE_FRACTION( currentTime, 536500800.0, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 0 ) - 273.15, 2.8, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 1 ) - 273.15, 10.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 2 ) * 1e-2, 897.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 3 ) * 1e-2, 7.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 4 ) * 100, 58, tolerance );

        // DATE: 170101 DOY: 001 DSS 10, TIME 0030
        currentTime = tudat::utilities::getKeyByIndex( weatherFile->meteoDataMap_, 1 );
        BOOST_CHECK_CLOSE_FRACTION( currentTime, 536502600.0, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 0 ) - 273.15, 1.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 1 ) - 273.15, 9.9, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 2 ) * 1e-2, 896.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 3 ) * 1e-2, 6.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 4 ) * 100, 54, tolerance );

        // DATE: 170101 DOY: 001 DSS 10, TIME 2359
        currentTime = tudat::utilities::getKeyByIndex( weatherFile->meteoDataMap_, 48 );
        BOOST_CHECK_CLOSE_FRACTION( currentTime, 536587140.0, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 0 ) - 273.15, 1.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 1 ) - 273.15, 10.5, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 2 ) * 1e-2, 894.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 3 ) * 1e-2, 6.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 4 ) * 100, 52, tolerance );

        // DATE: 170102 DOY: 002 DSS 10, TIME 0000
        currentTime = tudat::utilities::getKeyByIndex( weatherFile->meteoDataMap_, 49);
        BOOST_CHECK_CLOSE_FRACTION( currentTime, 536587200.0, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 0 ) - 273.15, 0.9, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 1 ) - 273.15, 10.5, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 2 ) * 1e-2, 894.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 3 ) * 1e-2, 6.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 4 ) * 100, 51, tolerance );
    }

    // Reading data with gaps
    {
        std::shared_ptr< DsnWeatherData > weatherFile =
                std::make_shared< DsnWeatherData >( tudat::paths::getTudatTestDataPath( ) + "mromagr20060012006365_10.wea.txt" );

        BOOST_CHECK_EQUAL( weatherFile->dsnStationComplexId_, 10 );

        // Comparison with values from file
        // Values of time obtained from https://nsidc.org/data/icesat/glas-date-conversion-tool/date_convert/

        double tolerance = 1e-13;  // Corresponds to numerical precision for temperature values

        // DATE: 060109 DOY: 009 DSS 10, TIME 0334
        double currentTime = tudat::utilities::getKeyByIndex( weatherFile->meteoDataMap_, 344);
        BOOST_CHECK_CLOSE_FRACTION( currentTime, 190049640.0, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 0 ) - 273.15, -7.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 1 ) - 273.15, 11.7, tolerance );
        BOOST_CHECK_EQUAL( std::isnan( weatherFile->meteoDataMap_.at( currentTime )( 2 ) ), true );
        BOOST_CHECK_CLOSE_FRACTION( weatherFile->meteoDataMap_.at( currentTime )( 3 ) * 1e-2, 3.6, tolerance );
        BOOST_CHECK_EQUAL( std::isnan( weatherFile->meteoDataMap_.at( currentTime )( 4 ) ), true );
    }
}

// Check that date is set correctly into the ground stations
// Also checks that files are ordered correctly
BOOST_AUTO_TEST_CASE( setWeatherData )
{
    using namespace tudat;
    using namespace tudat::input_output;

    double tolerance = 1e-13;  // Corresponds to numerical precision for temperature values

    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodiesToCreate = { "Earth" };
    simulation_setup::BodyListSettings bodySettings = simulation_setup::getDefaultBodySettings( bodiesToCreate );
    bodySettings.at( "Earth" )->groundStationSettings = simulation_setup::getDsnStationSettings( );
    simulation_setup::SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    setDsnWeatherDataInGroundStations(
            bodies,
            std::vector< std::string >{ tudat::paths::getTudatTestDataPath( ) + "mromagr20170012017365_60.wea.txt",
                                        tudat::paths::getTudatTestDataPath( ) + "mromagr20170012017365_10.wea.txt",
                                        tudat::paths::getTudatTestDataPath( ) + "mromagr20170012017365_40.wea.txt" } );

    simulation_setup::SystemOfBodies bodies2 = createSystemOfBodies( bodySettings );
    observation_models::setVmfTroposphereCorrections(
        std::vector< std::string >( { paths::getTudatTestDataPath( ) + "y2017.vmf3_r.txt" } ),
            true, false,  bodies2, false, true );


    std::vector< std::string > complex10GroundStations = { "DSS-13", "DSS-14", "DSS-15", "DSS-24", "DSS-25", "DSS-26" };
    for( std::string groundStation: complex10GroundStations )
    {
        std::shared_ptr< ground_stations::GroundStation > gs = bodies.getBody( "Earth" )->getGroundStation( groundStation );
        std::shared_ptr< ground_stations::GroundStation > gs2 = bodies2.getBody( "Earth" )->getGroundStation( groundStation );

        // DATE: 170101 DOY: 001 DSS 10, TIME 0000
        double time = 536500800.0;
        if( groundStation == "DSS-13" )
        {
            BOOST_CHECK_SMALL( gs->getDewPointFunction( )( time ) - gs2->getDewPointFunction( )( time ), 1.0 );
            BOOST_CHECK_SMALL( gs->getTemperatureFunction( )( time ) - gs2->getTemperatureFunction( )( time ), 1.0E-3 );
            BOOST_CHECK_SMALL( gs->getPressureFunction( )( time ) - gs2->getPressureFunction( )( time ), 0.1 );
            BOOST_CHECK_SMALL( gs->getWaterVaporPartialPressureFunction( )( time ) - gs2->getWaterVaporPartialPressureFunction( )( time ), 1.0E-2 );
            BOOST_CHECK_SMALL( gs->getRelativeHumidityFunction( )( time ) - gs2->getRelativeHumidityFunction( )( time ), 2.0E-2 );
        }
        else
        {
            BOOST_CHECK_SMALL( gs->getTemperatureFunction( )( time ) - gs2->getTemperatureFunction( )( time ), 5.0 );
            BOOST_CHECK_SMALL( gs->getPressureFunction( )( time ) - gs2->getPressureFunction( )( time ), 200.0 * ( groundStation == "DSS-14"  ? 2.5 : 1.0 ) );
            BOOST_CHECK_SMALL( gs->getWaterVaporPartialPressureFunction( )( time ) - gs2->getWaterVaporPartialPressureFunction( )( time ), 100.0 );


        }

        BOOST_CHECK_CLOSE_FRACTION( gs->getDewPointFunction( )( time ) - 273.15, 2.8, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getTemperatureFunction( )( time ) - 273.15, 10.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getPressureFunction( )( time ) * 1e-2, 897.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 7.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getRelativeHumidityFunction( )( time ) * 100, 58, tolerance );

        // DATE: 170101 DOY: 001 DSS 10, TIME 0030
        time = 536502600.0;
        BOOST_CHECK_CLOSE_FRACTION( gs->getDewPointFunction( )( time ) - 273.15, 1.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getTemperatureFunction( )( time ) - 273.15, 9.9, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getPressureFunction( )( time ) * 1e-2, 896.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 6.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getRelativeHumidityFunction( )( time ) * 100, 54, tolerance );

        // DATE: 170101 DOY: 001 DSS 10, TIME 2359
        // Why do some tolerances need larger values here? Loss of significant digits in the linear interpolation?
        time = 536587140.0;
        BOOST_CHECK_CLOSE_FRACTION( gs->getDewPointFunction( )( time ) - 273.15, 1.1, 1e-9 );
        BOOST_CHECK_CLOSE_FRACTION( gs->getTemperatureFunction( )( time ) - 273.15, 10.5, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getPressureFunction( )( time ) * 1e-2, 894.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 6.7, 1e-10 );
        BOOST_CHECK_CLOSE_FRACTION( gs->getRelativeHumidityFunction( )( time ) * 100, 52, 1e-10 );

        // DATE: 170102 DOY: 002 DSS 10, TIME 0000
        time = 536587200.0;
        BOOST_CHECK_CLOSE_FRACTION( gs->getDewPointFunction( )( time ) - 273.15, 0.9, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getTemperatureFunction( )( time ) - 273.15, 10.5, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getPressureFunction( )( time ) * 1e-2, 894.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 6.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getRelativeHumidityFunction( )( time ) * 100, 51, tolerance );
    }
}

BOOST_AUTO_TEST_CASE( compareDsnAndWmfWeatherData )
{
    double tolerance = 1e-13;  // Corresponds to numerical precision for temperature values

    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodiesToCreate = { "Earth" };
    simulation_setup::BodyListSettings bodySettings = simulation_setup::getDefaultBodySettings( bodiesToCreate );
    bodySettings.at( "Earth" )->groundStationSettings = simulation_setup::getDsnStationSettings( );
    simulation_setup::SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    setDsnWeatherDataInGroundStations(
        bodies,
        std::vector< std::string >{ tudat::paths::getTudatTestDataPath( ) + "mromagr20180012018365_10.wea.txt",
                                    tudat::paths::getTudatTestDataPath( ) + "mromagr20170012017365_60.wea.txt",
                                    tudat::paths::getTudatTestDataPath( ) + "mromagr20170012017365_10.wea.txt",
                                    tudat::paths::getTudatTestDataPath( ) + "mromagr20170012017365_40.wea.txt",
                                    tudat::paths::getTudatTestDataPath( ) + "mromagr20160012016366_10.wea.txt" } );

    simulation_setup::SystemOfBodies bodies2 = createSystemOfBodies( bodySettings );


    std::vector< std::string > complex10GroundStations = { "DSS-13", "DSS-14", "DSS-15", "DSS-24", "DSS-25", "DSS-26" };
    for( std::string groundStation: complex10GroundStations )
    {
        std::shared_ptr< ground_stations::GroundStation > gs = bodies.getBody( "Earth" )->getGroundStation( groundStation );
        // DATE: 170101 DOY: 001 DSS 10, TIME 0000
        double time = 536500800.0;

        BOOST_CHECK_CLOSE_FRACTION( gs->getDewPointFunction( )( time ) - 273.15, 2.8, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getTemperatureFunction( )( time ) - 273.15, 10.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getPressureFunction( )( time ) * 1e-2, 897.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 7.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getRelativeHumidityFunction( )( time ) * 100, 58, tolerance );

        // DATE: 170101 DOY: 001 DSS 10, TIME 0030
        time = 536502600.0;
        BOOST_CHECK_CLOSE_FRACTION( gs->getDewPointFunction( )( time ) - 273.15, 1.1, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getTemperatureFunction( )( time ) - 273.15, 9.9, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getPressureFunction( )( time ) * 1e-2, 896.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 6.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getRelativeHumidityFunction( )( time ) * 100, 54, tolerance );

        // DATE: 170101 DOY: 001 DSS 10, TIME 2359
        // Why do some tolerances need larger values here? Loss of significant digits in the linear interpolation?
        time = 536587140.0;
        BOOST_CHECK_CLOSE_FRACTION( gs->getDewPointFunction( )( time ) - 273.15, 1.1, 1e-9 );
        BOOST_CHECK_CLOSE_FRACTION( gs->getTemperatureFunction( )( time ) - 273.15, 10.5, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getPressureFunction( )( time ) * 1e-2, 894.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 6.7, 1e-10 );
        BOOST_CHECK_CLOSE_FRACTION( gs->getRelativeHumidityFunction( )( time ) * 100, 52, 1e-10 );

        // DATE: 170102 DOY: 002 DSS 10, TIME 0000
        time = 536587200.0;
        BOOST_CHECK_CLOSE_FRACTION( gs->getDewPointFunction( )( time ) - 273.15, 0.9, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getTemperatureFunction( )( time ) - 273.15, 10.5, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getPressureFunction( )( time ) * 1e-2, 894.7, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getWaterVaporPartialPressureFunction( )( time ) * 1e-2, 6.6, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( gs->getRelativeHumidityFunction( )( time ) * 100, 51, tolerance );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat