/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <cstdio>
#include <fstream>
#include <limits>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/io/readSinexFile.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/support/testFileUtilities.h"

namespace tudat
{
namespace unit_tests
{

namespace
{

std::string createTempPath( const std::string& suffix )
{
    return createTemporaryFilePath( "tudat-ilrs", suffix );
}

//! Write a complete fixed-width solution for the ILRS loader integration cases.
void writeSyntheticSinexFile( const std::string& fileName )
{
    std::ofstream file( fileName );
    file << "%=SNX 2.02 TUD 25:150:00000 TUD 25:001:00000 25:365:86399 L 00006 2 X V\n";
    file << "+SITE/ID\n";
    file << " 7110  A 12345M001 L ARECIBO               \n";
    file << "-SITE/ID\n";
    file << "+SOLUTION/ESTIMATE\n";
    file << "     1 STAX   7110  A    1 25:150:00000 m    2  2.39049000000000E+06 1.00000E-03\n";
    file << "     2 STAY   7110  A    1 25:150:00000 m    2 -5.56476300000000E+06 1.00000E-03\n";
    file << "     3 STAZ   7110  A    1 25:150:00000 m    2  1.99472700000000E+06 1.00000E-03\n";
    file << "     4 VELX   7110  A    1 25:150:00000 m/y  2  3.15576000000000E+07 1.00000E-03\n";
    file << "     5 VELY   7110  A    1 25:150:00000 m/y  2  0.00000000000000E+00 1.00000E-03\n";
    file << "     6 VELZ   7110  A    1 25:150:00000 m/y  2 -3.15576000000000E+07 1.00000E-03\n";
    file << "-SOLUTION/ESTIMATE\n";
    file << "%ENDSNX\n";
}

struct SyntheticEccentricityArc {
    std::string startEpoch_;
    std::string endEpoch_;
    Eigen::Vector3d eccentricity_ = Eigen::Vector3d::Zero( );
};

void writeSyntheticSinexEccentricityFile( const std::string& fileName, const std::vector< SyntheticEccentricityArc >& arcs )
{
    std::ofstream file( fileName );
    file << "+SITE/ID\n";
    file << "*CODE PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_     _SOD#___\n";
    file << " 7110  A 12345M001 L ARECIBO    FIXED      66 45  0.0  18 20  0.0   400.0     71100901\n";
    file << "-SITE/ID\n";
    file << "+SITE/ECCENTRICITY\n";
    file << "*SITE PT SOLN T DATA_START__ DATA_END____ XYZ X_______ Y_______ Z_______        CDP-SOD_\n";
    for( const SyntheticEccentricityArc& arc : arcs )
    {
        file << " 7110  A    1 L " << arc.startEpoch_ << " " << arc.endEpoch_ << " XYZ   " << arc.eccentricity_( 0 ) << "   "
             << arc.eccentricity_( 1 ) << "   " << arc.eccentricity_( 2 ) << "        71100901\n";
    }
    file << "-SITE/ECCENTRICITY\n";
}

void writeSyntheticSinexEccentricityFile( const std::string& fileName )
{
    writeSyntheticSinexEccentricityFile( fileName,
                                         { { "24:001:00000", "24:365:86399", Eigen::Vector3d( 0.1, 0.2, 0.3 ) },
                                           { "25:001:00000", "00:000:00000", Eigen::Vector3d( 0.4, 0.5, 0.6 ) } } );
}

std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > loadIlrsStationSettingsFromSyntheticFiles(
        const std::vector< SyntheticEccentricityArc >& arcs,
        const bool includeEccentricityFile = true )
{
    const std::string sinexFile = createTempPath( ".snx" );
    writeSyntheticSinexFile( sinexFile );

    std::string eccentricityFile = "";
    if( includeEccentricityFile )
    {
        eccentricityFile = createTempPath( ".snx" );
        writeSyntheticSinexEccentricityFile( eccentricityFile, arcs );
    }

    const std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > stationSettings =
            simulation_setup::getIlrsStationSettingsFromSinexDomes( { "12345M001" }, sinexFile, eccentricityFile, true );

    std::remove( sinexFile.c_str( ) );
    if( includeEccentricityFile )
    {
        std::remove( eccentricityFile.c_str( ) );
    }

    return stationSettings;
}

std::shared_ptr< simulation_setup::PiecewiseConstantGroundStationMotionSettings > getPiecewiseMotionSettings(
        const std::shared_ptr< simulation_setup::GroundStationSettings >& stationSettings )
{
    for( const auto& motionSettings : stationSettings->getStationMotionSettings( ) )
    {
        const std::shared_ptr< simulation_setup::PiecewiseConstantGroundStationMotionSettings > piecewiseSettings =
                std::dynamic_pointer_cast< simulation_setup::PiecewiseConstantGroundStationMotionSettings >( motionSettings );
        if( piecewiseSettings != nullptr )
        {
            return piecewiseSettings;
        }
    }

    return nullptr;
}

std::shared_ptr< simulation_setup::LinearGroundStationMotionSettings > getLinearMotionSettings(
        const std::shared_ptr< simulation_setup::GroundStationSettings >& stationSettings )
{
    for( const auto& motionSettings : stationSettings->getStationMotionSettings( ) )
    {
        const std::shared_ptr< simulation_setup::LinearGroundStationMotionSettings > linearSettings =
                std::dynamic_pointer_cast< simulation_setup::LinearGroundStationMotionSettings >( motionSettings );
        if( linearSettings != nullptr )
        {
            return linearSettings;
        }
    }

    return nullptr;
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_ilrs_station_loader )

//! Build an ILRS station from its registry identity, Cartesian state and eccentricity.
BOOST_AUTO_TEST_CASE( testLoadIlrsStationFromDomesAndSinex )
{
    const std::string sinexFile = createTempPath( ".snx" );
    const std::string eccentricityFile = createTempPath( ".snx" );
    writeSyntheticSinexFile( sinexFile );
    writeSyntheticSinexEccentricityFile( eccentricityFile );

    const std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > stationSettings =
            simulation_setup::getIlrsStationSettingsFromSinexDomes( { "12345M001" }, sinexFile, eccentricityFile, true );

    BOOST_CHECK_EQUAL( stationSettings.size( ), 1 );
    BOOST_CHECK_EQUAL( stationSettings.at( 0 )->getStationName( ), "12345M001" );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 0 ), 2390490.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 1 ), -5564763.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 2 ), 1994727.0, 1.0E-15 );

    bool foundLinearMotion = false;
    bool foundPiecewiseEccentricityMotion = false;
    for( const auto& motionSettings : stationSettings.at( 0 )->getStationMotionSettings( ) )
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
            const double firstArcEndEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:365:86399" );
            const double secondArcEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "25:001:00000" );

            BOOST_CHECK_CLOSE_FRACTION( piecewiseSettings->displacementList_.at( firstArcEpoch )( 0 ), 0.1, 1.0E-15 );
            BOOST_CHECK_CLOSE_FRACTION( piecewiseSettings->displacementList_.at( firstArcEpoch )( 1 ), 0.2, 1.0E-15 );
            BOOST_CHECK_CLOSE_FRACTION( piecewiseSettings->displacementList_.at( firstArcEpoch )( 2 ), 0.3, 1.0E-15 );
            BOOST_CHECK_SMALL( piecewiseSettings->displacementList_.at( firstArcEndEpoch ).norm( ), 1.0E-15 );
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

//! Select each eccentricity vector in its declared validity interval.
BOOST_AUTO_TEST_CASE( testEccentricitySelectionByEpoch )
{
    const std::string sinexFile = createTempPath( ".snx" );
    const std::string eccentricityFile = createTempPath( ".snx" );
    writeSyntheticSinexFile( sinexFile );
    writeSyntheticSinexEccentricityFile( eccentricityFile );

    const std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > stationSettings =
            simulation_setup::getIlrsStationSettingsFromSinexDomes( { "12345M001" }, sinexFile, eccentricityFile, true );

    BOOST_CHECK_EQUAL( stationSettings.size( ), 1 );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 0 ), 2390490.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 1 ), -5564763.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( stationSettings.at( 0 )->getGroundStationPosition( )( 2 ), 1994727.0, 1.0E-15 );

    bool foundPiecewiseEccentricityMotion = false;
    for( const auto& motionSettings : stationSettings.at( 0 )->getStationMotionSettings( ) )
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

//! Apply an initially unbounded eccentricity before the first dated change.
BOOST_AUTO_TEST_CASE( testFirstOpenEndedEccentricityUsesMinimumEpoch )
{
    const std::vector< SyntheticEccentricityArc > arcs = { { "24:001:00000", "00:000:00000", Eigen::Vector3d( 0.1, 0.2, 0.3 ) },
                                                           { "25:001:00000", "25:200:00000", Eigen::Vector3d( 0.4, 0.5, 0.6 ) } };
    const std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > stationSettings =
            loadIlrsStationSettingsFromSyntheticFiles( arcs, true );

    BOOST_REQUIRE_EQUAL( stationSettings.size( ), 1 );
    const std::shared_ptr< simulation_setup::PiecewiseConstantGroundStationMotionSettings > piecewiseSettings =
            getPiecewiseMotionSettings( stationSettings.at( 0 ) );
    BOOST_REQUIRE( piecewiseSettings != nullptr );

    const std::map< double, Eigen::Vector3d >& displacementList = piecewiseSettings->displacementList_;
    const double minimumEpoch = -std::numeric_limits< double >::max( );
    const double firstArcEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:001:00000" );
    const double secondArcEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "25:001:00000" );
    const double secondArcEndEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "25:200:00000" );

    BOOST_CHECK_EQUAL( displacementList.count( minimumEpoch ), 1 );
    BOOST_CHECK_CLOSE_FRACTION( displacementList.at( minimumEpoch )( 0 ), 0.1, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( displacementList.at( minimumEpoch )( 1 ), 0.2, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( displacementList.at( minimumEpoch )( 2 ), 0.3, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( displacementList.at( firstArcEpoch )( 0 ), 0.1, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( displacementList.at( secondArcEpoch )( 0 ), 0.4, 1.0E-15 );
    BOOST_CHECK_SMALL( displacementList.at( secondArcEndEpoch ).norm( ), 1.0E-15 );
}

//! Leave eccentricity corrections zero between disjoint closed intervals.
BOOST_AUTO_TEST_CASE( testClosedEccentricityArcsCreateZeroGaps )
{
    const std::vector< SyntheticEccentricityArc > arcs = { { "24:001:00000", "24:100:00000", Eigen::Vector3d( 0.1, 0.2, 0.3 ) },
                                                           { "24:200:00000", "24:300:00000", Eigen::Vector3d( 0.4, 0.5, 0.6 ) } };
    const std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > stationSettings =
            loadIlrsStationSettingsFromSyntheticFiles( arcs, true );

    BOOST_REQUIRE_EQUAL( stationSettings.size( ), 1 );
    const std::shared_ptr< simulation_setup::PiecewiseConstantGroundStationMotionSettings > piecewiseSettings =
            getPiecewiseMotionSettings( stationSettings.at( 0 ) );
    BOOST_REQUIRE( piecewiseSettings != nullptr );

    const std::map< double, Eigen::Vector3d >& displacementList = piecewiseSettings->displacementList_;
    const double minimumEpoch = -std::numeric_limits< double >::max( );
    const double firstArcEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:001:00000" );
    const double firstArcEndEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:100:00000" );
    const double secondArcEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:200:00000" );
    const double secondArcEndEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:300:00000" );

    BOOST_CHECK_EQUAL( displacementList.count( minimumEpoch ), 0 );
    BOOST_CHECK_EQUAL( displacementList.size( ), 4 );
    BOOST_CHECK_CLOSE_FRACTION( displacementList.at( firstArcEpoch )( 0 ), 0.1, 1.0E-15 );
    BOOST_CHECK_SMALL( displacementList.at( firstArcEndEpoch ).norm( ), 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( displacementList.at( secondArcEpoch )( 0 ), 0.4, 1.0E-15 );
    BOOST_CHECK_SMALL( displacementList.at( secondArcEndEpoch ).norm( ), 1.0E-15 );
}

//! Apply chronological eccentricity changes even when source records are unordered.
BOOST_AUTO_TEST_CASE( testUnsortedEccentricityArcsAreSorted )
{
    const std::vector< SyntheticEccentricityArc > arcs = { { "24:200:00000", "24:300:00000", Eigen::Vector3d( 0.4, 0.5, 0.6 ) },
                                                           { "24:001:00000", "24:100:00000", Eigen::Vector3d( 0.1, 0.2, 0.3 ) } };
    const std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > stationSettings =
            loadIlrsStationSettingsFromSyntheticFiles( arcs, true );

    BOOST_REQUIRE_EQUAL( stationSettings.size( ), 1 );
    const std::shared_ptr< simulation_setup::PiecewiseConstantGroundStationMotionSettings > piecewiseSettings =
            getPiecewiseMotionSettings( stationSettings.at( 0 ) );
    BOOST_REQUIRE( piecewiseSettings != nullptr );

    const std::map< double, Eigen::Vector3d >& displacementList = piecewiseSettings->displacementList_;
    const double firstArcEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:001:00000" );
    const double firstArcEndEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:100:00000" );
    const double secondArcEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:200:00000" );

    BOOST_CHECK_EQUAL( displacementList.begin( )->first, firstArcEpoch );
    BOOST_CHECK_CLOSE_FRACTION( displacementList.at( firstArcEpoch )( 0 ), 0.1, 1.0E-15 );
    BOOST_CHECK_SMALL( displacementList.at( firstArcEndEpoch ).norm( ), 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( displacementList.at( secondArcEpoch )( 0 ), 0.4, 1.0E-15 );
}

//! Respect both finite and unbounded validity for a single eccentricity record.
BOOST_AUTO_TEST_CASE( testSingleEntryClosedAndOpenEndedEccentricity )
{
    {
        const std::vector< SyntheticEccentricityArc > arcs = { { "24:001:00000", "24:100:00000", Eigen::Vector3d( 0.1, 0.2, 0.3 ) } };
        const std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > stationSettings =
                loadIlrsStationSettingsFromSyntheticFiles( arcs, true );

        BOOST_REQUIRE_EQUAL( stationSettings.size( ), 1 );
        const std::shared_ptr< simulation_setup::PiecewiseConstantGroundStationMotionSettings > piecewiseSettings =
                getPiecewiseMotionSettings( stationSettings.at( 0 ) );
        BOOST_REQUIRE( piecewiseSettings != nullptr );

        const std::map< double, Eigen::Vector3d >& displacementList = piecewiseSettings->displacementList_;
        const double minimumEpoch = -std::numeric_limits< double >::max( );
        const double firstArcEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:001:00000" );
        const double firstArcEndEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:100:00000" );

        BOOST_CHECK_EQUAL( displacementList.count( minimumEpoch ), 0 );
        BOOST_CHECK_EQUAL( displacementList.size( ), 2 );
        BOOST_CHECK_CLOSE_FRACTION( displacementList.at( firstArcEpoch )( 0 ), 0.1, 1.0E-15 );
        BOOST_CHECK_SMALL( displacementList.at( firstArcEndEpoch ).norm( ), 1.0E-15 );
    }

    {
        const std::vector< SyntheticEccentricityArc > arcs = { { "24:001:00000", "00:000:00000", Eigen::Vector3d( 0.4, 0.5, 0.6 ) } };
        const std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > stationSettings =
                loadIlrsStationSettingsFromSyntheticFiles( arcs, true );

        BOOST_REQUIRE_EQUAL( stationSettings.size( ), 1 );
        const std::shared_ptr< simulation_setup::PiecewiseConstantGroundStationMotionSettings > piecewiseSettings =
                getPiecewiseMotionSettings( stationSettings.at( 0 ) );
        BOOST_REQUIRE( piecewiseSettings != nullptr );

        const std::map< double, Eigen::Vector3d >& displacementList = piecewiseSettings->displacementList_;
        const double minimumEpoch = -std::numeric_limits< double >::max( );
        const double firstArcEpoch = input_output::convertSinexDateTimeToSecondsSinceEpoch( "24:001:00000" );

        BOOST_CHECK_EQUAL( displacementList.count( minimumEpoch ), 1 );
        BOOST_CHECK_EQUAL( displacementList.size( ), 2 );
        BOOST_CHECK_CLOSE_FRACTION( displacementList.at( minimumEpoch )( 0 ), 0.4, 1.0E-15 );
        BOOST_CHECK_CLOSE_FRACTION( displacementList.at( firstArcEpoch )( 0 ), 0.4, 1.0E-15 );
    }
}

//! Load the Cartesian station without adding a missing eccentricity correction.
BOOST_AUTO_TEST_CASE( testNoEccentricityDataSkipsPiecewiseMotion )
{
    const std::vector< std::shared_ptr< simulation_setup::GroundStationSettings > > stationSettings =
            loadIlrsStationSettingsFromSyntheticFiles( {}, false );

    BOOST_REQUIRE_EQUAL( stationSettings.size( ), 1 );
    BOOST_CHECK_EQUAL( getPiecewiseMotionSettings( stationSettings.at( 0 ) ), nullptr );
    BOOST_REQUIRE( getLinearMotionSettings( stationSettings.at( 0 ) ) != nullptr );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
