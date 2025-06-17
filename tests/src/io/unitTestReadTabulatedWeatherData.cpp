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
        currentTime = tudat::utilities::getKeyByIndex( weatherFile->meteoDataMap_, 49 );
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
        double currentTime = tudat::utilities::getKeyByIndex( weatherFile->meteoDataMap_, 344 );
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
BOOST_AUTO_TEST_CASE( compareDsnAndWmfWeatherData )
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
                                        tudat::paths::getTudatTestDataPath( ) + "mromagr20170012017365_40.wea.txt" },
            interpolators::linearInterpolation( ) );

    simulation_setup::SystemOfBodies bodies2 = createSystemOfBodies( bodySettings );
    observation_models::setVmfTroposphereCorrections(
            std::vector< std::string >( { paths::getTudatTestDataPath( ) + "y2017.vmf3_r.txt" } ), true, false, bodies2, false, true );

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
            BOOST_CHECK_SMALL( gs->getWaterVaporPartialPressureFunction( )( time ) - gs2->getWaterVaporPartialPressureFunction( )( time ),
                               1.0E-2 );
            BOOST_CHECK_SMALL( gs->getRelativeHumidityFunction( )( time ) - gs2->getRelativeHumidityFunction( )( time ), 2.0E-2 );
        }
        else
        {
            BOOST_CHECK_SMALL( gs->getTemperatureFunction( )( time ) - gs2->getTemperatureFunction( )( time ), 5.0 );
            BOOST_CHECK_SMALL( gs->getPressureFunction( )( time ) - gs2->getPressureFunction( )( time ),
                               200.0 * ( groundStation == "DSS-14" ? 2.5 : 1.0 ) );
            BOOST_CHECK_SMALL( gs->getWaterVaporPartialPressureFunction( )( time ) - gs2->getWaterVaporPartialPressureFunction( )( time ),
                               100.0 );
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

BOOST_AUTO_TEST_CASE( setDsnWeatherData )
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
                                        tudat::paths::getTudatTestDataPath( ) + "mromagr20160012016366_10.wea.txt" },
            interpolators::linearInterpolation( ) );

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

BOOST_AUTO_TEST_CASE( setEstrackWeatherData )
{
    double tolerance = 1e-13;  // Corresponds to numerical precision for temperature values

    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodiesToCreate = { "Earth" };
    simulation_setup::BodyListSettings bodySettings = simulation_setup::getDefaultBodySettings( bodiesToCreate );
    bodySettings.at( "Earth" )->groundStationSettings = simulation_setup::getRadioTelescopeStationSettings( );
    simulation_setup::SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    std::vector< std::string > fullEstrackDataFiles =
            std::vector< std::string >{ tudat::paths::getTudatTestDataPath( ) + "M32ICL2L1B_MET_133621727_00.TAB",
                                        tudat::paths::getTudatTestDataPath( ) + "M32ICL2L1B_MET_133621727_01.TAB",
                                        tudat::paths::getTudatTestDataPath( ) + "M32ICL2L1B_MET_133621727_02.TAB",
                                        tudat::paths::getTudatTestDataPath( ) + "M32ICL2L1B_MET_133621727_03.TAB",
                                        tudat::paths::getTudatTestDataPath( ) + "M32ICL2L1B_MET_133621727_04.TAB",
                                        tudat::paths::getTudatTestDataPath( ) + "M32ICL2L1B_MET_133621727_05.TAB" };

    std::vector< double > fullTimes;
    std::vector< Eigen::VectorXd > fullMeteo;

    for( int test = 0; test < 4; test++ )
    {
        std::vector< std::string > estrackDataFiles;
        if( test == 0 )
        {
            estrackDataFiles = fullEstrackDataFiles;
        }
        else if( test == 1 )
        {
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 1 ) );
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 0 ) );
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 4 ) );
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 5 ) );
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 3 ) );
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 2 ) );
        }
        else if( test == 2 )
        {
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 0 ) );
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 1 ) );
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 3 ) );
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 4 ) );
        }
        else if( test == 3 )
        {
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 0 ) );
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 2 ) );
            estrackDataFiles.push_back( fullEstrackDataFiles.at( 4 ) );
        }
        setEstrackWeatherDataInGroundStation( bodies, estrackDataFiles, "NWNORCIA", interpolators::linearInterpolation( ) );

        std::shared_ptr< ground_stations::GroundStation > gs = bodies.getBody( "Earth" )->getGroundStation( "NWNORCIA" );
        std::shared_ptr< ground_stations::PiecewiseInterpolatedMeteoData > meteoDataObject =
                std::dynamic_pointer_cast< ground_stations::PiecewiseInterpolatedMeteoData >( gs->getMeteoData( ) );

        if( test == 0 )
        {
            BOOST_CHECK_EQUAL( meteoDataObject->getMeteoDataInterpolators( ).size( ), 1 );

            fullTimes = meteoDataObject->getMeteoDataInterpolators( ).at( 0 )->getIndependentValues( );
            fullMeteo = meteoDataObject->getMeteoDataInterpolators( ).at( 0 )->getDependentValues( );
        }
        else if( test == 1 )
        {
            BOOST_CHECK_EQUAL( meteoDataObject->getMeteoDataInterpolators( ).size( ), 1 );

            std::vector< double > currentTimes = meteoDataObject->getMeteoDataInterpolators( ).at( 0 )->getIndependentValues( );
            std::vector< Eigen::VectorXd > currentMeteo = meteoDataObject->getMeteoDataInterpolators( ).at( 0 )->getDependentValues( );

            BOOST_CHECK_EQUAL( currentTimes.size( ), fullTimes.size( ) );

            for( int i = 0; i < currentTimes.size( ); i++ )
            {
                BOOST_CHECK_EQUAL( currentTimes.at( i ), fullTimes.at( i ) );
                BOOST_CHECK_EQUAL( currentMeteo.at( i )( 0 ), fullMeteo.at( i )( 0 ) );
                BOOST_CHECK_EQUAL( currentMeteo.at( i )( 1 ), fullMeteo.at( i )( 1 ) );
                BOOST_CHECK_EQUAL( currentMeteo.at( i )( 2 ), fullMeteo.at( i )( 2 ) );
            }
        }
        else if( test == 2 || test == 3 )
        {
            BOOST_CHECK_EQUAL( meteoDataObject->getMeteoDataInterpolators( ).size( ), test );

            std::vector< int > shiftIndex;
            int blockSize;
            if( test == 2 )
            {
                shiftIndex = { 0, 180 };
                blockSize = 120;
            }
            else
            {
                shiftIndex = { 0, 120, 240 };
                blockSize = 60;
            }
            for( int j = 0; j < test; j++ )
            {
                std::vector< double > currentTimes = meteoDataObject->getMeteoDataInterpolators( ).at( j )->getIndependentValues( );
                std::vector< Eigen::VectorXd > currentMeteo = meteoDataObject->getMeteoDataInterpolators( ).at( j )->getDependentValues( );

                BOOST_CHECK_EQUAL( currentTimes.size( ), blockSize );

                for( int i = 0; i < blockSize; i++ )
                {
                    BOOST_CHECK_EQUAL( currentTimes.at( i ), fullTimes.at( i + shiftIndex.at( j ) ) );
                    BOOST_CHECK_EQUAL( currentMeteo.at( i )( 0 ), fullMeteo.at( i + shiftIndex.at( j ) )( 0 ) );
                    BOOST_CHECK_EQUAL( currentMeteo.at( i )( 1 ), fullMeteo.at( i + shiftIndex.at( j ) )( 1 ) );
                    BOOST_CHECK_EQUAL( currentMeteo.at( i )( 2 ), fullMeteo.at( i + shiftIndex.at( j ) )( 2 ) );
                }
            }
        }

        {
            double testEpoch = basic_astrodynamics::DateTime::fromIsoString( "2013-12-28T17:43:20.000" ).epoch< double >( );

            double testTemperature = 19.6 + 273.15;
            double testHumidity = 0.549;
            double testPressure = 97950.0;

            BOOST_CHECK_CLOSE_FRACTION(
                    gs->getDewPointFunction( )( testEpoch ), ground_stations::computeDewPoint( testHumidity, testTemperature ), tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getTemperatureFunction( )( testEpoch ), testTemperature, tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getPressureFunction( )( testEpoch ), testPressure, tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getWaterVaporPartialPressureFunction( )( testEpoch ),
                                        testHumidity * ground_stations::computeSaturationWaterVaporPressure( testTemperature ),
                                        tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getRelativeHumidityFunction( )( testEpoch ), testHumidity, tolerance );
        }

        {
            double testEpoch = basic_astrodynamics::DateTime::fromIsoString( "2013-12-28T22:21:20.000" ).epoch< double >( );

            double testTemperature = 21.9 + 273.15;
            double testHumidity = 0.385;
            double testPressure = 98020.0;

            BOOST_CHECK_CLOSE_FRACTION(
                    gs->getDewPointFunction( )( testEpoch ), ground_stations::computeDewPoint( testHumidity, testTemperature ), tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getTemperatureFunction( )( testEpoch ), testTemperature, tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getPressureFunction( )( testEpoch ), testPressure, tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getWaterVaporPartialPressureFunction( )( testEpoch ),
                                        testHumidity * ground_stations::computeSaturationWaterVaporPressure( testTemperature ),
                                        tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getRelativeHumidityFunction( )( testEpoch ), testHumidity, tolerance );
        }

        if( test < 3 )
        {
            double testEpoch = basic_astrodynamics::DateTime::fromIsoString( "2013-12-28T18:26:50.000" ).epoch< double >( );

            double testTemperature = ( 19.8 + 19.7 ) / 2.0 + 273.15;
            double testHumidity = ( 0.529 + 0.533 ) / 2.0;
            double testPressure = ( 97920.0 + 97980.0 ) / 2.0;

            BOOST_CHECK_CLOSE_FRACTION(
                    gs->getDewPointFunction( )( testEpoch ), ground_stations::computeDewPoint( testHumidity, testTemperature ), tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getTemperatureFunction( )( testEpoch ), testTemperature, tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getPressureFunction( )( testEpoch ), testPressure, tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getWaterVaporPartialPressureFunction( )( testEpoch ),
                                        testHumidity * ground_stations::computeSaturationWaterVaporPressure( testTemperature ),
                                        tolerance );
            BOOST_CHECK_CLOSE_FRACTION( gs->getRelativeHumidityFunction( )( testEpoch ), testHumidity, tolerance );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat