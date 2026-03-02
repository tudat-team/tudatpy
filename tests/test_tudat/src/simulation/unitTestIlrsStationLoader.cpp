/*    Copyright (c) 2010-2026, Delft University of Technology
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

#include <cstdio>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/io/readSinexFile.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

namespace tudat
{
namespace unit_tests
{

namespace
{

std::string createTempPath( const std::string& suffix )
{
    const boost::filesystem::path tempPath =
            boost::filesystem::temp_directory_path( ) / boost::filesystem::unique_path( "tudat-ilrs-%%%%%%" + suffix );
    return tempPath.string( );
}

void writeSyntheticSinexFile( const std::string& fileName )
{
    std::ofstream file( fileName );
    file << "+SITE/ID\n";
    file << "*CODE PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_     _SOD#___\n";
    file << " 7110  A 12345M001 L ARECIBO    FIXED      66 45  0.0  18 20  0.0   400.0     71100901\n";
    file << "-SITE/ID\n";
    file << "+SOLUTION/ESTIMATE\n";
    file << "*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S ____ESTIMATED_VALUE ____STD_DEV___\n";
    file << "1 STAX 7110 A 1 25:150:00000 m 2 2390490.0 0.001\n";
    file << "2 STAY 7110 A 1 25:150:00000 m 2 -5564763.0 0.001\n";
    file << "3 STAZ 7110 A 1 25:150:00000 m 2 1994727.0 0.001\n";
    file << "4 VELX 7110 A 1 25:150:00000 m/y 2 31557600.0 0.001\n";
    file << "5 VELY 7110 A 1 25:150:00000 m/y 2 0.0 0.001\n";
    file << "6 VELZ 7110 A 1 25:150:00000 m/y 2 -31557600.0 0.001\n";
    file << "-SOLUTION/ESTIMATE\n";
}

void writeSyntheticSinexEccentricityFile( const std::string& fileName )
{
    std::ofstream file( fileName );
    file << "+SITE/ID\n";
    file << "*CODE PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_     _SOD#___\n";
    file << " 7110  A 12345M001 L ARECIBO    FIXED      66 45  0.0  18 20  0.0   400.0     71100901\n";
    file << "-SITE/ID\n";
    file << "+SITE/ECCENTRICITY\n";
    file << "*SITE PT SOLN T DATA_START__ DATA_END____ XYZ X_______ Y_______ Z_______        CDP-SOD_\n";
    file << " 7110  A    1 L 24:001:00000 24:365:86399 XYZ   0.1000   0.2000   0.3000        71100901\n";
    file << " 7110  A    1 L 25:001:00000 00:000:00000 XYZ   0.4000   0.5000   0.6000        71100901\n";
    file << "-SITE/ECCENTRICITY\n";
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_ilrs_station_loader )

BOOST_AUTO_TEST_CASE( testLoadIlrsStationFromDomesAndSinex )
{
    const std::string sinexFile = createTempPath( ".snx" );
    const std::string eccentricityFile = createTempPath( ".snx" );
    writeSyntheticSinexFile( sinexFile );
    writeSyntheticSinexEccentricityFile( eccentricityFile );

    const std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > stationSettings =
            simulation_setup::getIlrsStationSettingsFromSinexDomes(
                    { "12345M001" }, sinexFile, eccentricityFile, TUDAT_NAN, true );

    BOOST_CHECK_EQUAL( stationSettings.size( ), 1 );
    BOOST_CHECK_EQUAL( stationSettings.at( 0 )->getStationName( ), "12345M001" );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 0 ), 2390490.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 1 ), -5564763.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 2 ), 1994727.0, 1.0E-15 );

    bool foundLinearMotion = false;
    bool foundPiecewiseEccentricityMotion = false;
    for( const auto& motionSettings: stationSettings.at( 0 )->getStationMotionSettings( ) )
    {
        const std::shared_ptr< simulation_setup::LinearGroundStationMotionSettings > linearSettings =
                std::dynamic_pointer_cast< simulation_setup::LinearGroundStationMotionSettings >( motionSettings );
        if( linearSettings != nullptr )
        {
            foundLinearMotion = true;
            BOOST_CHECK_CLOSE_FRACTION( linearSettings->linearVelocity_( 0 ), 1.0, 1.0E-15 );
            BOOST_CHECK_CLOSE_FRACTION( linearSettings->linearVelocity_( 1 ), 0.0, 1.0E-15 );
            BOOST_CHECK_CLOSE_FRACTION( linearSettings->linearVelocity_( 2 ), -1.0, 1.0E-15 );
        }

        const std::shared_ptr< simulation_setup::PiecewiseConstantGroundStationMotionSettings > piecewiseSettings =
                std::dynamic_pointer_cast< simulation_setup::PiecewiseConstantGroundStationMotionSettings >( motionSettings );
        if( piecewiseSettings != nullptr )
        {
            foundPiecewiseEccentricityMotion = true;
            BOOST_CHECK_EQUAL( piecewiseSettings->displacementList_.size( ), 3 );

            const double firstArcEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:001:00000" );
            const double secondArcEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "25:001:00000" );

            BOOST_CHECK_CLOSE_FRACTION( piecewiseSettings->displacementList_.at( firstArcEpoch )( 0 ), 0.1, 1.0E-15 );
            BOOST_CHECK_CLOSE_FRACTION( piecewiseSettings->displacementList_.at( firstArcEpoch )( 1 ), 0.2, 1.0E-15 );
            BOOST_CHECK_CLOSE_FRACTION( piecewiseSettings->displacementList_.at( firstArcEpoch )( 2 ), 0.3, 1.0E-15 );
            BOOST_CHECK_CLOSE_FRACTION( piecewiseSettings->displacementList_.at( secondArcEpoch )( 0 ), 0.4, 1.0E-15 );
            BOOST_CHECK_CLOSE_FRACTION( piecewiseSettings->displacementList_.at( secondArcEpoch )( 1 ), 0.5, 1.0E-15 );
            BOOST_CHECK_CLOSE_FRACTION( piecewiseSettings->displacementList_.at( secondArcEpoch )( 2 ), 0.6, 1.0E-15 );
        }
    }
    BOOST_CHECK_EQUAL( foundLinearMotion, true );
    BOOST_CHECK_EQUAL( foundPiecewiseEccentricityMotion, true );

    std::remove( sinexFile.c_str( ) );
    std::remove( eccentricityFile.c_str( ) );
}

BOOST_AUTO_TEST_CASE( testEccentricitySelectionByEpoch )
{
    const std::string sinexFile = createTempPath( ".snx" );
    const std::string eccentricityFile = createTempPath( ".snx" );
    writeSyntheticSinexFile( sinexFile );
    writeSyntheticSinexEccentricityFile( eccentricityFile );

    const double evaluationEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:200:00000" );
    const std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > stationSettings =
            simulation_setup::getIlrsStationSettingsFromSinexDomes(
                    { "12345M001" }, sinexFile, eccentricityFile, evaluationEpoch, true );

    BOOST_CHECK_EQUAL( stationSettings.size( ), 1 );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 0 ), 2390490.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 1 ), -5564763.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 2 ), 1994727.0, 1.0E-15 );

    bool foundPiecewiseEccentricityMotion = false;
    for( const auto& motionSettings: stationSettings.at( 0 )->getStationMotionSettings( ) )
    {
        const std::shared_ptr< simulation_setup::PiecewiseConstantGroundStationMotionSettings > piecewiseSettings =
                std::dynamic_pointer_cast< simulation_setup::PiecewiseConstantGroundStationMotionSettings >( motionSettings );
        if( piecewiseSettings != nullptr )
        {
            foundPiecewiseEccentricityMotion = true;
            BOOST_CHECK_EQUAL( piecewiseSettings->displacementList_.size( ), 3 );
        }
    }
    BOOST_CHECK_EQUAL( foundPiecewiseEccentricityMotion, true );

    std::remove( sinexFile.c_str( ) );
    std::remove( eccentricityFile.c_str( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
