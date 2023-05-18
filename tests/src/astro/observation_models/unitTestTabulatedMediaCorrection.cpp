/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include "tudat/basics/testMacros.h"

#include <numeric>

#include "tudat/astro/observation_models.h"
#include "tudat/simulation/estimation_setup.h"
#include "tudat/simulation/environment_setup.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/io/readTabulatedMediaCorrections.h"

namespace tudat
{
namespace unit_tests
{

using namespace observation_models;

BOOST_AUTO_TEST_SUITE( test_tabulated_media_corrections )

// Only testing wet simplified Chao model. Consistency of dry Chao model with Niell checked in following test
// Anyway, wet and dry Chao are identical but with different coefficients.
BOOST_AUTO_TEST_CASE( testSimplifiedChaoMappingFunction )
{

    // Estefan and Sovers (1994), A Comparative Survey of Current and Proposed Tropospheric Refraction-Delay Models for
    // DSN Radio Metric Data Calibration, JPL/NASA, 94-24
    // Appendix
    std::vector< double > chaoTabulatedWetMapping = {
            35.3955, 23.4935, 17.2402, 13.4500, 10.9967, 9.2827, 8.0235, 7.0621, 6.3054, 5.6951, 5.1929,
            4.7728, 4.4164, 4.1104, 3.8449, 3.6125, 3.4075, 3.2253, 3.0625, 2.9160, 2.7838, 2.6637, 2.5543
        };

    std::vector< double > chaoTabulatedMappingElevation( chaoTabulatedWetMapping.size( ) );
    for ( unsigned int i = 0; i < chaoTabulatedMappingElevation.size( ); ++i )
    {
        chaoTabulatedMappingElevation.at( i ) = ( i + 1.0 ) * mathematical_constants::PI / 180.0;
    }

    for ( unsigned int i = 0; i < chaoTabulatedWetMapping.size( ); ++i )
    {
        std::function< double ( Eigen::Vector3d, double ) > elevationFunction =
                [=]( Eigen::Vector3d inertialVectorAwayFromStation, double time ){ return chaoTabulatedMappingElevation.at( i ); };

        // isUplinkCorrection set to true. Irrelevant if true or false, since the elevation function returns a  fixed value
        SimplifiedChaoTroposphericMapping simplifiedChaoModel = SimplifiedChaoTroposphericMapping( elevationFunction, true );

        // Moyer (2000) says that error should be smaller than 1% for elevation larger than 1 deg. Not sure
        // why it needs higher tolerance below 3.5 deg
        double tolerance;
        if ( chaoTabulatedMappingElevation.at( i ) * 180.0 / mathematical_constants::PI > 3.5 )
        {
            tolerance = 0.01;
        }
        else
        {
            tolerance = 0.025;
        }
        // Arguments to computeWetTroposphericMapping are irrelevant
        BOOST_CHECK_CLOSE_FRACTION(
                chaoTabulatedWetMapping.at( i ),
                simplifiedChaoModel.computeWetTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), 0.0, 0.0 ), tolerance );
    }

}

// Check consistency between Niell and Chao mapping
BOOST_AUTO_TEST_CASE( testNiellChaoMappingFunctionConsistency )
{

    std::vector< double > elevation( 17 );
    for ( unsigned int i = 0; i < elevation.size( ); ++i )
    {
        elevation.at( i ) = ( i + 2 ) * 5.0 * mathematical_constants::PI / 180.0;
    }

    // Geodetic position: [altitude, latitude, longitude]
    std::function< Eigen::Vector3d ( double ) > groundStationGeodeticPositionFunction =
            [=]( double time ){ return ( Eigen::Vector3d( ) << 800.0, 35.0 * mathematical_constants::PI / 180.0, TUDAT_NAN ).finished( ); };

    // Comparing wet mapping
    for ( unsigned int i = 0; i < elevation.size( ); ++i )
    {
        std::function< double ( Eigen::Vector3d, double ) > elevationFunction =
                [=]( Eigen::Vector3d inertialVectorAwayFromStation, double time ){ return elevation.at( i ); };


        // isUplinkCorrection set to true. Irrelevant if true or false, since the elevation function returns a  fixed value
        SimplifiedChaoTroposphericMapping simplifiedChaoModel = SimplifiedChaoTroposphericMapping( elevationFunction, true );
        NiellTroposphericMapping niellModel = NiellTroposphericMapping(
                elevationFunction, groundStationGeodeticPositionFunction, true );

        double chaoMapping = simplifiedChaoModel.computeWetTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), 0.0, 0.0 );
        double niellMapping = niellModel.computeWetTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), 0.0, 0.0 );
        BOOST_CHECK_CLOSE_FRACTION( chaoMapping, niellMapping, 0.01 );
    }

    // Comparing dry mapping
    for ( unsigned int i = 0; i < elevation.size( ); ++i )
    {
        std::function< double ( Eigen::Vector3d, double ) > elevationFunction =
                [=]( Eigen::Vector3d inertialVectorAwayFromStation, double time ){ return elevation.at( i ); };

        // isUplinkCorrection set to true. Irrelevant if true or false, since the elevation function returns a  fixed value
        SimplifiedChaoTroposphericMapping simplifiedChaoModel = SimplifiedChaoTroposphericMapping( elevationFunction, true );
        NiellTroposphericMapping niellModel = NiellTroposphericMapping(
                elevationFunction, groundStationGeodeticPositionFunction, true );

        double chaoMapping = simplifiedChaoModel.computeDryTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), 0.0, 0.0 );
        double niellMapping = niellModel.computeDryTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), 0.0, 0.0 );

        BOOST_CHECK_CLOSE_FRACTION( chaoMapping, niellMapping, 0.01 );
    }
}

// Compare values of Niell mapping function with
// Niell (1996), Global mapping functions for the atmosphere delay at radio wavelengths, Journal of Geophysics Research
BOOST_AUTO_TEST_CASE( testNiellMappingFunctionNiellPaper )
{

    double elevation = 5.0 * mathematical_constants::PI / 180.0;
    double geodeticLatitude = 64.92 * mathematical_constants::PI / 180.0;
    double altitude = 132.0;

    double time1987 = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
            1987, 1, 1, 0, 0, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    // Figure 3 of Niell (1996)
    std::vector< double > dryMapping = { 10.1875, 10.18, 10.1475, 10.1375, 10.185 };
    std::vector< double > dryTolerance = { 0.0025, 0.0025, 0.0025, 0.0025, 0.005 };
    std::vector< double > dryTime { time1987 + 0.1 * physical_constants::JULIAN_YEAR,
                                    time1987 + 0.2 * physical_constants::JULIAN_YEAR,
                                    time1987 + 0.4 * physical_constants::JULIAN_YEAR,
                                    time1987 + 0.7 * physical_constants::JULIAN_YEAR,
                                    time1987 + 1.0 * physical_constants::JULIAN_YEAR };

    // Figure 5 of Niell (1996)
    double wetMapping = 10.73;
    double wetTolerance = 0.01;
    std::vector< double > wetTime = dryTime;

    // Create mapping function
    std::function< double ( Eigen::Vector3d, double ) > elevationFunction =
                [=]( Eigen::Vector3d inertialVectorAwayFromStation, double time ){ return elevation; };
    // Geodetic position: [altitude, latitude, longitude]
    std::function< Eigen::Vector3d ( double ) > groundStationGeodeticPositionFunction =
            [=]( double time ){ return ( Eigen::Vector3d( ) << altitude, geodeticLatitude, TUDAT_NAN ).finished( ); };

    NiellTroposphericMapping niellModel = NiellTroposphericMapping(
            elevationFunction, groundStationGeodeticPositionFunction, true );

    // Test dry mapping
    for ( unsigned int i = 0; i < dryMapping.size( ); ++i )
    {
        double calculatedMapping = niellModel.computeDryTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), dryTime.at( i ), dryTime.at( i ) );

        BOOST_CHECK_SMALL( std::abs( calculatedMapping - dryMapping.at( i ) ), dryTolerance.at( i ) );
    }

    // Test wet mapping
    for ( unsigned int i = 0; i < wetTime.size( ); ++i )
    {
        double calculatedMapping = niellModel.computeWetTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), wetTime.at( i ), wetTime.at( i ) );

        BOOST_CHECK_SMALL( std::abs( calculatedMapping - wetMapping ), wetTolerance );
    }

}

// Compare Niell mapping values with GODOT
// Reference values extracted from GODOT for XMM-Newton satellite
BOOST_AUTO_TEST_CASE( testNiellMappingFunctionGodot )
{
    double toleranceWet = 1e-13;
    // Dry tolerance needs to be much larger because GODOT has two errors:
    // - Instead of using 28 days as reference time, they use 27 days
    // - When computing the dry coefficients they are missing the 1/2 term inside the cosine
    // See https://gitlab.space-codev.org/godot/godot/-/blob/f0e574576361b2b227d60672d6f48edc564e9f14/godot/model/obs/TroposphereMappingFunctions.cpp#L74-75
    double toleranceDry = 1e-5;

    for ( unsigned int gs = 0; gs < 3; ++gs )
    {
        std::vector< double > times = { 517557600.0 };
        Eigen::Vector3d groundStationGeodeticPosition;
        std::vector< double > elevations;
        std::vector< double > expectedWetMapping;
        std::vector< double > expectedDryMapping;
        if ( gs == 0 ) // Kourou
        {
            groundStationGeodeticPosition = ( Eigen::Vector3d( ) <<
                    -0.0142714925195833e3, 0.0916549210572268, TUDAT_NAN ).finished( );
            elevations = { 0.279499332478995 };
            expectedWetMapping = { 3.59993033305367 };
            expectedDryMapping = { 3.57132955130791 };
        }
        else if ( gs == 1 ) // Kiruna
        {
            groundStationGeodeticPosition = ( Eigen::Vector3d( ) <<
                    0.402633377266284e3, 1.18433033971699, TUDAT_NAN ).finished( );
            elevations = { -1.03692130895857 };
            expectedWetMapping = { -1.16140864286966 };
            expectedDryMapping = { -1.16115998824861 };
        }
        else // Canberra_34
        {
            groundStationGeodeticPosition = ( Eigen::Vector3d( ) <<
                    0.692418697207358e3, -0.61781997683548, TUDAT_NAN ).finished( );
            elevations = { 0.187111816282427 };
            expectedWetMapping = { 5.29410018021229 };
            expectedDryMapping = { 5.20743179879896 };
        }

        for ( unsigned int j = 0; j < times.size( ); ++j )
        {
            // Create mapping function
            std::function< double ( Eigen::Vector3d, double ) > elevationFunction =
                    [=]( Eigen::Vector3d inertialVectorAwayFromStation, double time ){ return elevations.at( j ); };
            // Geodetic position: [altitude, latitude, longitude]
            std::function< Eigen::Vector3d ( double ) > groundStationGeodeticPositionFunction =
                    [=]( double time ){ return groundStationGeodeticPosition; };

            NiellTroposphericMapping niellModel = NiellTroposphericMapping(
                    elevationFunction, groundStationGeodeticPositionFunction, true );

            BOOST_CHECK_CLOSE_FRACTION(
                    niellModel.computeWetTroposphericMapping(
                            Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), times.at( j ), times.at( j ) ),
                    expectedWetMapping.at( j ), toleranceWet );
            BOOST_CHECK_CLOSE_FRACTION(
                    niellModel.computeDryTroposphericMapping(
                            Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), times.at( j ), times.at( j ) ),
                    expectedDryMapping.at( j ), toleranceDry );
        }
    }
}

// Compare Saastamoinen correction values with GODOT
// Reference values extracted from GODOT using godot::model::obs::saastamoinenDryCorrection and
// godot::model::obs::saastamoinenWetCorrection functions
BOOST_AUTO_TEST_CASE( testSaastamoinenTroposphericCorrectionGodot )
{
    double tolerance = 1e-14;

    std::vector< double > temperature = { 10.7 + 273.15, 5.9 + 273.15, 3.3 + 273.15, 4.7 + 273.15 }; // K
    std::vector< double > pressure = { 897.1e2, 897.8e2, 896.3e2, 897.8e2 }; // Pa
    std::vector< double > relativeHumidity = { 0.58, 0.95, 1.0, 1.0 }; // [-]

    std::vector< double > expectedDryCorrection = { 6.81310428429791e-09, 6.81842049542154e-09, 6.80702861444233e-09, 6.81842049542154e-09 };
    std::vector< double > expectedWetCorrection = { 2.53513266881891e-10, 3.04776408738257e-10, 2.69896652188657e-10, 2.96360948712566e-10 };
    // GODOT uses 0.0022768 instead of 0.002277 for the wet correction
    double factor = 0.002277 / 0.0022768;
    for ( unsigned int i = 0; i < expectedWetCorrection.size( ); ++i )
    {
        expectedWetCorrection.at( i ) *= factor;
    }

    // Geodetic position: selected so that position correction is zero (because GODOT doesn't correct the Saastamoinen
    // correction for the ground station position)
    std::function< Eigen::Vector3d ( double ) > groundStationGeodeticPositionFunction =
            [=]( double time ){ return ( Eigen::Vector3d( ) <<
                    0.0, mathematical_constants::PI / 4, TUDAT_NAN ).finished( ); };

    for ( unsigned int i = 0; i < temperature.size( ); ++i )
    {
        std::function< double ( double ) > temperatureFunction = [=]( double time ){ return temperature.at( i ); };
        std::function< double ( double ) > pressureFunction = [=]( double time ){ return pressure.at( i ); };
        std::function< double ( double ) > humidityFunction = [=]( double time ){ return relativeHumidity.at( i ); };

        std::function< double ( double ) > waterVaporPartialPressureFunction = getBeanAndDuttonWaterVaporPartialPressureFunction(
                humidityFunction, temperatureFunction );

        SaastamoinenTroposphericCorrection troposphericCorrection = SaastamoinenTroposphericCorrection(
                groundStationGeodeticPositionFunction, pressureFunction, temperatureFunction, waterVaporPartialPressureFunction,
                nullptr, true );

        BOOST_CHECK_CLOSE_FRACTION(
                troposphericCorrection.getDryZenithRangeCorrectionFunction( )( TUDAT_NAN ) / physical_constants::SPEED_OF_LIGHT,
                expectedDryCorrection.at( i ),
                tolerance );
        BOOST_CHECK_CLOSE_FRACTION(
                troposphericCorrection.getWetZenithRangeCorrectionFunction( )( TUDAT_NAN ) / physical_constants::SPEED_OF_LIGHT,
                expectedWetCorrection.at( i ),
                tolerance );
    }

}

// Check consistency between tabulated and Saastamoinen tropospheric corrections
BOOST_AUTO_TEST_CASE( testTabulatedAndSaastamoinenTroposphericCorrectionsConsistency )
{
    // Create bodies
    spice_interface::loadStandardSpiceKernels( { "/Users/pipas/Documents/mro-spice/mro_psp43.bsp" } );

    std::vector< std::string > bodiesToCreate = { "Earth", "Mars" };
    simulation_setup::BodyListSettings bodySettings = simulation_setup::getDefaultBodySettings( bodiesToCreate );
    bodySettings.at( "Earth" )->groundStationSettings = simulation_setup::getDsnStationSettings( );

    std::string spacecraftName = "MRO";
    bodySettings.addSettings( spacecraftName );
    bodySettings.at( spacecraftName )->ephemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >(  );

    simulation_setup::SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Load weather data
    input_output::setDsnWeatherDataInGroundStations(
            bodies,
            std::vector< std::string >
                    { "/Users/pipas/Documents/mro-data/wea/mromagr20170012017365_10.wea.txt",
                      "/Users/pipas/Documents/mro-data/wea/mromagr20170012017365_40.wea.txt",
                      "/Users/pipas/Documents/mro-data/wea/mromagr20170012017365_60.wea.txt" } );

    // Load tabulated corrections data
    std::shared_ptr< input_output::CspRawFile > troposphericCspFile = std::make_shared< input_output::CspRawFile >(
            "/Users/pipas/Documents/mro-data/tro/mromagr2017_091_2017_121.tro.txt" );

    // Create link ends
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "MRO" );
    linkEnds[ receiver ] = LinkEndId( "Earth", "DSS-26" );

    // Create Saastamoinen corrections
    std::shared_ptr< LightTimeCorrectionSettings > saastamoinenCorrectionSettings =
            std::make_shared< SaastamoinenTroposphericCorrectionSettings >( );
    std::shared_ptr< LightTimeCorrection > saastamoinenCorrectionBase = createLightTimeCorrections(
            saastamoinenCorrectionSettings, bodies, linkEnds, transmitter, receiver,
            observation_models::dsn_n_way_averaged_doppler );
    std::shared_ptr< MappedTroposphericCorrection > saastamoinenCorrection =
            std::dynamic_pointer_cast< MappedTroposphericCorrection >( saastamoinenCorrectionBase );

    // Create tabulated corrections
    std::shared_ptr< LightTimeCorrectionSettings > tabulatedCorrectionSettings =
            std::make_shared< observation_models::TabulatedTroposphericCorrectionSettings >(
                    input_output::createTroposphericDryCorrectionAdjustment( { troposphericCspFile } ),
                    input_output::createTroposphericWetCorrectionAdjustment( { troposphericCspFile } ) );
    std::shared_ptr< LightTimeCorrection > tabulatedCorrectionBase = createLightTimeCorrections(
            tabulatedCorrectionSettings, bodies, linkEnds, transmitter, receiver,
            observation_models::dsn_n_way_averaged_doppler );
    std::shared_ptr< MappedTroposphericCorrection > tabulatedCorrection =
            std::dynamic_pointer_cast< MappedTroposphericCorrection >( tabulatedCorrectionBase );

    double initialTime = 544795200.0;
    double timeStep = 3600.0;
    unsigned int numberOfPoints = 144;

    for ( unsigned int i = 0; i < numberOfPoints; ++i )
    {
        double time = initialTime + i * timeStep;

        std::cout << "SAS: " << std::setprecision(15) <<
            saastamoinenCorrection->getDryZenithRangeCorrectionFunction( )( time ) << "\t" <<
            saastamoinenCorrection->getWetZenithRangeCorrectionFunction( )( time ) << std::endl;
        std::cout << "TAB: " <<
            tabulatedCorrection->getDryZenithRangeCorrectionFunction( )( time ) << "\t" <<
            tabulatedCorrection->getWetZenithRangeCorrectionFunction( )( time ) << std::endl << std::endl;
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}