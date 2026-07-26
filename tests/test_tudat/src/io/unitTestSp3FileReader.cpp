/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/io/readSp3File.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"

namespace tudat
{
namespace unit_tests
{

namespace
{

const std::string velocitySp3Contents = R"sp3(#dV2001  8  8  0  0  0.00000000       2 ORBIT IGS97 HLM MGEX
## 1126 259200.00000000   900.00000000 52129 0.0000000000000
+    1   G01
%c G  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
/* Adapted from IGS SP3-d example 2, as distributed with the Orekit tests.
/* https://gitlab.orekit.org/orekit/orekit/-/raw/develop/src/test/resources/sp3/example-d-2.sp3
*  2001  8  8  0  0  0.00000000
PG01 -11044.805800 -10475.672350  21929.418200    189.163300 18 18 18 219
EP    55   55   55     222  1234567 -1234567  5999999      -30       21 -1230000
VG01  20298.880364 -18462.044804   1381.387685     -4.534317 14 14 14 191
EV    22   22   22     111  1234567  1234567  1234567  1234567  1234567  1234567
*  2001  8  9 23 45  0.00000000
PG01 -11044.805800 -10475.672350  21929.418200    189.163300 18 18 18 219  P   P
EP    55   55   55     222  1234567 -1234567  5999999      -30       21 -1230000
VG01  20298.880364 -18462.044804   1381.387685     -4.534317 14 14 14 191
EV    22   22   22     111  1234567  1234567  1234567  1234567  1234567  1234567
EOF
)sp3";

const std::string positionOnlySp3Contents = R"sp3(#cP2001  8  8  0  0  0.00000000       3 ORBIT WGS84 FIT TEST
## 1126 259200.00000000    10.00000000 52129 0.0000000000000
+    1   E11
%c M  cc GAL ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
*  2001  8  8  0  0  0.00000000
PE11      1.000000     -2.000000      3.000000 999999.999999
*  2001  8  8  0  0 10.00000000
PE11     71.000000     23.000000    -37.000000 999999.999999
*  2001  8  8  0  0 30.00000000
PE11    511.000000    223.000000   -117.000000 999999.999999
EOF
)sp3";

boost::filesystem::path writeTemporarySp3File( const std::string& contents )
{
    const boost::filesystem::path path =
            boost::filesystem::temp_directory_path( ) / boost::filesystem::unique_path( "tudat-sp3-%%%%%%.sp3" );
    std::ofstream stream( path.string( ) );
    stream << contents;
    return path;
}

std::string makeStateRecord( const char recordType, const std::string& satelliteId, const double x, const double y, const double z )
{
    std::ostringstream stream;
    stream << recordType << std::setw( 3 ) << satelliteId << std::fixed << std::setprecision( 6 ) << std::setw( 14 ) << x << std::setw( 14 )
           << y << std::setw( 14 ) << z << std::setw( 14 ) << 0.0 << '\n';
    return stream.str( );
}

std::string makeMultiConstellationSp3Contents( )
{
    const std::vector< std::string > satelliteIds = { "G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09",
                                                      "E01", "E02", "E03", "E04", "E05", "E06", "E07", "E08", "E09" };

    std::ostringstream stream;
    stream << "#dV2001  8  8  0  0  0.00000000       1 ORBIT IGb14 FIT TEST\n"
           << "## 1126 259200.00000000   900.00000000 52129 0.0000000000000\n"
           << "+   18   ";
    for( unsigned int i = 0; i < 17; ++i )
    {
        stream << satelliteIds.at( i );
    }
    stream << "\n+        " << satelliteIds.back( ) << "\n%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n"
           << "*  2001  8  8  0  0  0.00000000\n";

    for( unsigned int i = 0; i < satelliteIds.size( ); ++i )
    {
        stream << makeStateRecord( 'P', satelliteIds.at( i ), 10000.0 + i, 20000.0 + i, 30000.0 + i );
        stream << makeStateRecord( 'V', satelliteIds.at( i ), 10.0 + i, 20.0 + i, 30.0 + i );
    }
    stream << "EOF\n";
    return stream.str( );
}

void checkReadThrows( const std::string& contents )
{
    const boost::filesystem::path path = writeTemporarySp3File( contents );
    BOOST_CHECK_THROW( input_output::readSp3File( path.string( ) ), std::runtime_error );
    boost::filesystem::remove( path );
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_sp3_file_reader )

BOOST_AUTO_TEST_CASE( testVelocityFileMetadataStateAndLegacyAlias )
{
    const boost::filesystem::path path = writeTemporarySp3File( velocitySp3Contents );
    const std::shared_ptr< input_output::Sp3FileContents > contents = input_output::readSp3File( path.string( ) );

    BOOST_CHECK_EQUAL( contents->formatVersion, 'd' );
    BOOST_CHECK( contents->hasVelocityRecords );
    BOOST_CHECK( !contents->velocitiesWereDerived );
    BOOST_CHECK_EQUAL( contents->declaredNumberOfEpochs, 2 );
    BOOST_CHECK_EQUAL( contents->declaredNumberOfSatellites, 1 );
    BOOST_CHECK_CLOSE_FRACTION( contents->declaredEpochInterval, 900.0, 1.0e-15 );
    BOOST_CHECK_EQUAL( contents->frameName, "IGS97" );
    BOOST_CHECK_EQUAL( contents->analysisCenter, "MGEX" );
    BOOST_CHECK_EQUAL( contents->timeScale, "GPS" );
    BOOST_REQUIRE_EQUAL( contents->satelliteStates.count( "G01" ), 1 );
    BOOST_REQUIRE_EQUAL( contents->satelliteStates.at( "G01" ).size( ), 2 );

    const double expectedStartEpoch = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch< double >(
                                              2001, 8, 8, 0, 0, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;
    BOOST_CHECK_SMALL( contents->startEpoch - expectedStartEpoch, 1.0e-9 );

    const Eigen::VectorXd& firstState = contents->satelliteStates.at( "G01" ).begin( )->second;
    BOOST_CHECK_CLOSE_FRACTION( firstState( 0 ), -11044805.8, 1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( firstState( 1 ), -10475672.35, 1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( firstState( 2 ), 21929418.2, 1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( firstState( 3 ), 2029.8880364, 1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( firstState( 4 ), -1846.2044804, 1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( firstState( 5 ), 138.1387685, 1.0e-15 );

    const std::shared_ptr< input_output::SP3cFileContents > legacyContents = input_output::readSp3cFile( path.string( ) );
    BOOST_REQUIRE_EQUAL( legacyContents->satelliteStates.at( "G01" ).size( ), contents->satelliteStates.at( "G01" ).size( ) );
    BOOST_CHECK( legacyContents->satelliteStates.at( "G01" ).begin( )->second.isApprox( firstState ) );
    boost::filesystem::remove( path );
}

BOOST_AUTO_TEST_CASE( testPositionOnlyFileUsesSecondOrderFiniteDifferences )
{
    const boost::filesystem::path path = writeTemporarySp3File( positionOnlySp3Contents );
    const std::shared_ptr< input_output::Sp3FileContents > contents = input_output::readSp3File( path.string( ) );

    BOOST_CHECK_EQUAL( contents->formatVersion, 'c' );
    BOOST_CHECK( !contents->hasVelocityRecords );
    BOOST_CHECK( contents->velocitiesWereDerived );
    BOOST_CHECK_EQUAL( contents->timeScale, "GAL" );
    BOOST_REQUIRE_EQUAL( contents->satelliteStates.at( "E11" ).size( ), 3 );

    const std::vector< Eigen::Vector3d > expectedVelocities = { Eigen::Vector3d( 2000.0, 0.0, -4000.0 ),
                                                                Eigen::Vector3d( 12000.0, 5000.0, -4000.0 ),
                                                                Eigen::Vector3d( 32000.0, 15000.0, -4000.0 ) };
    unsigned int stateIndex = 0;
    for( const auto& state : contents->satelliteStates.at( "E11" ) )
    {
        BOOST_CHECK_SMALL( ( state.second.segment< 3 >( 3 ) - expectedVelocities.at( stateIndex ) ).norm( ), 1.0e-10 );
        stateIndex++;
    }

    const auto settings = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >(
            simulation_setup::sp3EphemerisSettings( contents, "E11" ) );
    BOOST_REQUIRE( settings != nullptr );
    BOOST_CHECK_EQUAL( settings->getFrameOrientation( ), "WGS84" );
    boost::filesystem::remove( path );
}

BOOST_AUTO_TEST_CASE( testEphemerisFactoryFrameTransformations )
{
    const boost::filesystem::path path = writeTemporarySp3File( velocitySp3Contents );
    const std::shared_ptr< input_output::Sp3FileContents > contents = input_output::readSp3File( path.string( ) );
    const Eigen::Vector6d sourceState = contents->satelliteStates.at( "G01" ).begin( )->second;

    const auto itrf2014Settings = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >(
            simulation_setup::sp3EphemerisSettings( contents, "G01", "Earth", "ITRF2014" ) );
    BOOST_REQUIRE( itrf2014Settings != nullptr );
    BOOST_CHECK_EQUAL( itrf2014Settings->getFrameOrientation( ), "ITRF2014" );
    const Eigen::Vector6d itrf2014State = itrf2014Settings->getBodyStateHistory( ).begin( )->second;
    BOOST_CHECK_GT( ( itrf2014State - sourceState ).norm( ), 1.0e-3 );
    BOOST_CHECK_LT( ( itrf2014State.segment< 3 >( 0 ) - sourceState.segment< 3 >( 0 ) ).norm( ), 1.0 );

    const auto j2000Settings = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >(
            simulation_setup::sp3EphemerisSettings( contents, "G01", "Earth", "J2000" ) );
    BOOST_REQUIRE( j2000Settings != nullptr );
    const Eigen::Vector6d j2000State = j2000Settings->getBodyStateHistory( ).begin( )->second;
    BOOST_CHECK_CLOSE_FRACTION( j2000State.segment< 3 >( 0 ).norm( ), sourceState.segment< 3 >( 0 ).norm( ), 1.0e-14 );
    BOOST_CHECK_GT( ( j2000State.segment< 3 >( 3 ) - sourceState.segment< 3 >( 3 ) ).norm( ), 100.0 );

    const std::shared_ptr< input_output::Sp3FileContents > j2000Contents = std::make_shared< input_output::Sp3FileContents >( *contents );
    j2000Contents->frameName = "J2000";
    for( const auto& state : j2000Settings->getBodyStateHistory( ) )
    {
        j2000Contents->satelliteStates.at( "G01" ).at( state.first ) = state.second;
    }
    const auto roundTripSettings = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >(
            simulation_setup::sp3EphemerisSettings( j2000Contents, "G01", "Earth", "IGS97" ) );
    BOOST_REQUIRE( roundTripSettings != nullptr );
    BOOST_CHECK_SMALL( ( roundTripSettings->getBodyStateHistory( ).begin( )->second - sourceState ).norm( ), 1.0e-3 );

    const std::shared_ptr< input_output::Sp3FileContents > aliasContents = std::make_shared< input_output::Sp3FileContents >( *contents );
    aliasContents->frameName = "IGb14";
    const auto aliasSettings = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >(
            simulation_setup::sp3EphemerisSettings( aliasContents, "G01", "Earth", "ITRF2014" ) );
    BOOST_REQUIRE( aliasSettings != nullptr );
    BOOST_CHECK( aliasSettings->getBodyStateHistory( ).begin( )->second.isApprox( sourceState, 0.0 ) );

    const std::shared_ptr< input_output::Sp3FileContents > wgs84Contents = std::make_shared< input_output::Sp3FileContents >( *contents );
    wgs84Contents->frameName = "WGS84";
    BOOST_CHECK_THROW( simulation_setup::sp3EphemerisSettings( wgs84Contents, "G01", "Earth", "ITRF2014" ), std::invalid_argument );

    const std::shared_ptr< input_output::Sp3FileContents > unknownFrameContents =
            std::make_shared< input_output::Sp3FileContents >( *contents );
    unknownFrameContents->frameName = "LOCAL";
    BOOST_CHECK_THROW( simulation_setup::sp3EphemerisSettings( unknownFrameContents, "G01", "Earth", "J2000" ), std::invalid_argument );
    boost::filesystem::remove( path );
}

BOOST_AUTO_TEST_CASE( testTudatEphemerisFactories )
{
    const boost::filesystem::path path = writeTemporarySp3File( velocitySp3Contents );
    const std::shared_ptr< input_output::Sp3FileContents > contents = input_output::readSp3File( path.string( ) );

    const std::shared_ptr< simulation_setup::EphemerisSettings > fromContents = simulation_setup::sp3EphemerisSettings( contents, "G01" );
    const std::shared_ptr< simulation_setup::EphemerisSettings > fromFile = simulation_setup::sp3EphemerisSettings( path.string( ), "G01" );

    const auto tabulatedFromContents = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >( fromContents );
    const auto tabulatedFromFile = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >( fromFile );
    BOOST_REQUIRE( tabulatedFromContents != nullptr );
    BOOST_REQUIRE( tabulatedFromFile != nullptr );
    BOOST_CHECK_EQUAL( tabulatedFromContents->getFrameOrigin( ), "Earth" );
    BOOST_CHECK_EQUAL( tabulatedFromContents->getFrameOrientation( ), "IGS97" );
    BOOST_REQUIRE_EQUAL( tabulatedFromContents->getBodyStateHistory( ).size( ), tabulatedFromFile->getBodyStateHistory( ).size( ) );
    BOOST_CHECK( tabulatedFromContents->getBodyStateHistory( ).begin( )->second.isApprox(
            tabulatedFromFile->getBodyStateHistory( ).begin( )->second ) );

    BOOST_CHECK_THROW( simulation_setup::sp3EphemerisSettings( contents, "G99" ), std::invalid_argument );

    contents->satelliteStates.at( "G01" ).begin( )->second( 3 ) = std::numeric_limits< double >::quiet_NaN( );
    BOOST_CHECK_THROW( simulation_setup::sp3EphemerisSettings( contents, "G01" ), std::runtime_error );
    boost::filesystem::remove( path );
}

BOOST_AUTO_TEST_CASE( testSupportedVersionsTimeSystemsAndFrameMetadata )
{
    for( const char version : { 'a', 'b', 'c', 'd' } )
    {
        std::string contents = velocitySp3Contents;
        contents.at( 1 ) = version;
        if( version == 'a' || version == 'b' )
        {
            contents.replace( contents.find( "G01" ), 3, "  1" );
            std::string::size_type stateIdPosition = contents.find( "G01" );
            while( stateIdPosition != std::string::npos )
            {
                contents.replace( stateIdPosition, 3, "  1" );
                stateIdPosition = contents.find( "G01", stateIdPosition + 3 );
            }
            contents.replace( contents.find( "GPS" ), 3, "ccc" );
        }

        const boost::filesystem::path path = writeTemporarySp3File( contents );
        const auto parsed = input_output::readSp3File( path.string( ) );
        BOOST_CHECK_EQUAL( parsed->formatVersion, version );
        BOOST_CHECK_EQUAL( parsed->timeScale, "GPS" );
        BOOST_CHECK_EQUAL( parsed->satelliteStates.count( version < 'c' ? "1" : "G01" ), 1 );
        boost::filesystem::remove( path );
    }

    const std::vector< std::string > timeSystems = { "GPS", "GLO", "GAL", "BDT", "TAI", "UTC", "IRN", "QZS" };
    for( const std::string& timeSystem : timeSystems )
    {
        std::string contents = velocitySp3Contents;
        contents.replace( contents.find( "GPS" ), 3, timeSystem );
        const boost::filesystem::path path = writeTemporarySp3File( contents );
        BOOST_CHECK_EQUAL( input_output::readSp3File( path.string( ) )->timeScale, timeSystem );
        boost::filesystem::remove( path );
    }

    for( const std::string frameName : { "IGS97", "IGb14", "WGS84", "GCRF " } )
    {
        std::string contents = velocitySp3Contents;
        contents.replace( contents.find( "IGS97" ), 5, frameName );
        const boost::filesystem::path path = writeTemporarySp3File( contents );
        BOOST_CHECK_EQUAL( input_output::readSp3File( path.string( ) )->frameName, boost::algorithm::trim_copy( frameName ) );
        boost::filesystem::remove( path );
    }
}

BOOST_AUTO_TEST_CASE( testMultipleSatelliteListLinesAndConstellations )
{
    const boost::filesystem::path path = writeTemporarySp3File( makeMultiConstellationSp3Contents( ) );
    const auto contents = input_output::readSp3File( path.string( ) );
    BOOST_CHECK_EQUAL( contents->declaredNumberOfSatellites, 18 );
    BOOST_CHECK_EQUAL( contents->satelliteStates.size( ), 18 );
    BOOST_CHECK_EQUAL( contents->satelliteStates.count( "G01" ), 1 );
    BOOST_CHECK_EQUAL( contents->satelliteStates.count( "E09" ), 1 );
    boost::filesystem::remove( path );
}

BOOST_AUTO_TEST_CASE( testBadOrAbsentStateValuesAreMarkedNonFinite )
{
    for( const char missingRecordType : { 'P', 'V' } )
    {
        std::string contents = velocitySp3Contents;
        const std::string recordPrefix = std::string( 1, missingRecordType ) + "G01";
        const std::string::size_type recordStart = contents.find( recordPrefix );
        const std::string::size_type recordEnd = contents.find( '\n', recordStart );
        contents.replace( recordStart, recordEnd - recordStart + 1, makeStateRecord( missingRecordType, "G01", 0.0, 0.0, 0.0 ) );

        const boost::filesystem::path path = writeTemporarySp3File( contents );
        const auto parsed = input_output::readSp3File( path.string( ) );
        const Eigen::VectorXd& firstState = parsed->satelliteStates.at( "G01" ).begin( )->second;
        const int missingSegmentStart = missingRecordType == 'P' ? 0 : 3;
        BOOST_CHECK( !firstState.segment< 3 >( missingSegmentStart ).allFinite( ) );
        BOOST_CHECK_THROW( simulation_setup::sp3EphemerisSettings( parsed, "G01" ), std::runtime_error );
        boost::filesystem::remove( path );
    }
}

BOOST_AUTO_TEST_CASE( testMalformedAndIncompleteFilesAreRejected )
{
    std::string contents = velocitySp3Contents;
    contents.replace( contents.find( "       2 ORBIT" ), 8, "       3" );
    checkReadThrows( contents );

    contents = velocitySp3Contents;
    contents.at( 1 ) = 'e';
    checkReadThrows( contents );

    contents = velocitySp3Contents;
    contents.at( 2 ) = 'X';
    checkReadThrows( contents );

    contents = velocitySp3Contents;
    contents.replace( contents.find( "GPS" ), 3, "UNK" );
    checkReadThrows( contents );

    contents = velocitySp3Contents;
    contents.erase( contents.find( "VG01" ), contents.find( '\n', contents.find( "VG01" ) ) - contents.find( "VG01" ) + 1 );
    checkReadThrows( contents );

    contents = positionOnlySp3Contents;
    const std::string::size_type lastEpoch = contents.rfind( "*  " );
    contents.erase( lastEpoch, contents.find( "EOF", lastEpoch ) - lastEpoch );
    contents.replace( contents.find( "       3 ORBIT" ), 8, "       2" );
    checkReadThrows( contents );

    contents = positionOnlySp3Contents;
    contents.insert( contents.find( "EOF" ), "VE11      1.000000      2.000000      3.000000      0.000000\n" );
    checkReadThrows( contents );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
