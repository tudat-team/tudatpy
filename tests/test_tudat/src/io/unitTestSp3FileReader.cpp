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
#include "tudat/simulation/environment_setup/createSp3Ephemeris.h"

namespace tudat
{
namespace unit_tests
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

//! Write SP3 text to a uniquely named temporary file.
boost::filesystem::path writeTemporarySp3File( const std::string& contents )
{
    // Give each fixture a unique path so test cases remain independent when the suite runs concurrently.
    const boost::filesystem::path path =
            boost::filesystem::temp_directory_path( ) / boost::filesystem::unique_path( "tudat-sp3-%%%%%%.sp3" );
    // Write exactly the supplied fixed-width records; the individual tests remove the file after inspection.
    std::ofstream stream( path.string( ) );
    stream << contents;
    return path;
}

//! Format one fixed-width SP3 position or velocity record.
std::string makeStateRecord( const char recordType, const std::string& satelliteId, const double x, const double y, const double z )
{
    // Match the SP3 P/V column widths so tests can replace a complete state record without disturbing adjacent records.
    std::ostringstream stream;
    stream << recordType << std::setw( 3 ) << satelliteId << std::fixed << std::setprecision( 6 ) << std::setw( 14 ) << x << std::setw( 14 )
           << y << std::setw( 14 ) << z << std::setw( 14 ) << 0.0 << '\n';
    return stream.str( );
}

//! Create an SP3 file that exercises continuation lines in the satellite list.
std::string makeMultiConstellationSp3Contents( )
{
    // Eighteen identifiers force one continuation line because a '+' record contains at most seventeen entries.
    const std::vector< std::string > satelliteIds = { "G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09",
                                                      "E01", "E02", "E03", "E04", "E05", "E06", "E07", "E08", "E09" };

    // Place the first seventeen identifiers on the initial list record and the last identifier on its continuation.
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

    // Supply a complete P/V pair for every declared satellite so this fixture isolates satellite-list parsing.
    for( unsigned int i = 0; i < satelliteIds.size( ); ++i )
    {
        stream << makeStateRecord( 'P', satelliteIds.at( i ), 10000.0 + i, 20000.0 + i, 30000.0 + i );
        stream << makeStateRecord( 'V', satelliteIds.at( i ), 10.0 + i, 20.0 + i, 30.0 + i );
    }
    stream << "EOF\n";
    return stream.str( );
}

//! Verify that parsing the supplied SP3 text raises a runtime error.
void checkReadThrows( const std::string& contents )
{
    const boost::filesystem::path path = writeTemporarySp3File( contents );
    // Confirm the malformed fixture is rejected by the reader rather than producing partial contents.
    BOOST_CHECK_THROW( input_output::readSp3File( path.string( ) ), std::runtime_error );
    boost::filesystem::remove( path );
}

BOOST_AUTO_TEST_SUITE( test_sp3_file_reader )

BOOST_AUTO_TEST_CASE( testVelocityFileMetadataAndState )
{
    // This test reads an SP3-d velocity product and verifies its header metadata, epoch conversion, and SI-unit Cartesian state.
    const boost::filesystem::path path = writeTemporarySp3File( velocitySp3Contents );
    const std::shared_ptr< input_output::Sp3FileContents > contents = input_output::readSp3File( path.string( ) );

    // Confirm the reader retains the format version that determines the applicable SP3 grammar.
    BOOST_CHECK_EQUAL( contents->formatVersion, 'd' );
    // Confirm a V-mode header is exposed so callers know velocities came from the file.
    BOOST_CHECK( contents->hasVelocityRecords );
    // Confirm file-provided velocities are not mislabeled as numerically reconstructed values.
    BOOST_CHECK( !contents->velocitiesWereDerived );
    // Confirm the declared epoch count is parsed for later completeness validation.
    BOOST_CHECK_EQUAL( contents->declaredNumberOfEpochs, 2 );
    // Confirm the declared satellite count is parsed from the satellite-list header.
    BOOST_CHECK_EQUAL( contents->declaredNumberOfSatellites, 1 );
    // Confirm the nominal sample interval is parsed without a unit conversion.
    BOOST_CHECK_CLOSE_FRACTION( contents->declaredEpochInterval, 900.0, 1.0e-15 );
    // Confirm the terrestrial realization tag remains available for environment-level frame conversion.
    BOOST_CHECK_EQUAL( contents->frameName, "IGS97" );
    // Confirm the analysis-center metadata is retained for product provenance.
    BOOST_CHECK_EQUAL( contents->analysisCenter, "MGEX" );
    // Confirm the time-system tag remains available for environment-level Earth-rotation setup.
    BOOST_CHECK_EQUAL( contents->timeScale, "GPS" );
    // Require the declared satellite to exist before accessing its state history below.
    BOOST_REQUIRE_EQUAL( contents->satelliteStates.count( "G01" ), 1 );
    // Require both declared epochs to be present before selecting the first state.
    BOOST_REQUIRE_EQUAL( contents->satelliteStates.at( "G01" ).size( ), 2 );

    const double expectedStartEpoch = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch< double >(
                                              2001, 8, 8, 0, 0, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;
    // Confirm calendar fields are converted to seconds since J2000, which is the reader's default reference.
    BOOST_CHECK_SMALL( contents->startEpoch - expectedStartEpoch, 1.0e-9 );

    const Eigen::Vector6d& firstState = contents->satelliteStates.at( "G01" ).begin( )->second;
    // Confirm the x position is converted from the SP3 kilometre field to metres.
    BOOST_CHECK_CLOSE_FRACTION( firstState( 0 ), -11044805.8, 1.0e-15 );
    // Confirm the y position is converted from the SP3 kilometre field to metres.
    BOOST_CHECK_CLOSE_FRACTION( firstState( 1 ), -10475672.35, 1.0e-15 );
    // Confirm the z position is converted from the SP3 kilometre field to metres.
    BOOST_CHECK_CLOSE_FRACTION( firstState( 2 ), 21929418.2, 1.0e-15 );
    // Confirm the x velocity is converted from the SP3 decimetre-per-second field to metres per second.
    BOOST_CHECK_CLOSE_FRACTION( firstState( 3 ), 2029.8880364, 1.0e-15 );
    // Confirm the y velocity is converted from the SP3 decimetre-per-second field to metres per second.
    BOOST_CHECK_CLOSE_FRACTION( firstState( 4 ), -1846.2044804, 1.0e-15 );
    // Confirm the z velocity is converted from the SP3 decimetre-per-second field to metres per second.
    BOOST_CHECK_CLOSE_FRACTION( firstState( 5 ), 138.1387685, 1.0e-15 );

    boost::filesystem::remove( path );
}

BOOST_AUTO_TEST_CASE( testPositionOnlyFileUsesSecondOrderFiniteDifferences )
{
    // This test reads an SP3-c position-only product and verifies that complete states are formed with quadratic finite-difference
    // velocities.
    const boost::filesystem::path path = writeTemporarySp3File( positionOnlySp3Contents );
    const std::shared_ptr< input_output::Sp3FileContents > contents = input_output::readSp3File( path.string( ) );

    // Confirm the older SP3-c version is accepted and reported correctly.
    BOOST_CHECK_EQUAL( contents->formatVersion, 'c' );
    // Confirm the reader distinguishes position-only input from a file containing V records.
    BOOST_CHECK( !contents->hasVelocityRecords );
    // Confirm the completed velocity components are explicitly marked as derived.
    BOOST_CHECK( contents->velocitiesWereDerived );
    // Confirm non-GPS time-system metadata is retained without conversion in data_input.
    BOOST_CHECK_EQUAL( contents->timeScale, "GAL" );
    // Require all three epochs because each second-order derivative uses a three-point stencil.
    BOOST_REQUIRE_EQUAL( contents->satelliteStates.at( "E11" ).size( ), 3 );

    const std::vector< Eigen::Vector3d > expectedVelocities = { Eigen::Vector3d( 2000.0, 0.0, -4000.0 ),
                                                                Eigen::Vector3d( 12000.0, 5000.0, -4000.0 ),
                                                                Eigen::Vector3d( 32000.0, 15000.0, -4000.0 ) };
    unsigned int stateIndex = 0;
    for( const auto& state : contents->satelliteStates.at( "E11" ) )
    {
        // Confirm the forward, centered, or backward stencil reproduces the analytic derivative at this epoch.
        BOOST_CHECK_SMALL( ( state.second.segment< 3 >( 3 ) - expectedVelocities.at( stateIndex ) ).norm( ), 1.0e-10 );
        stateIndex++;
    }

    const auto settings = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >(
            simulation_setup::sp3EphemerisSettings( contents, "E11" ) );
    // Require the SP3 factory to return tabulated settings before inspecting their frame metadata.
    BOOST_REQUIRE( settings != nullptr );
    // Confirm an omitted target frame preserves the frame declared by the SP3 file.
    BOOST_CHECK_EQUAL( settings->getFrameOrientation( ), "WGS84" );
    boost::filesystem::remove( path );
}

BOOST_AUTO_TEST_CASE( testEphemerisFactoryFrameTransformations )
{
    // This test exercises terrestrial-realization and terrestrial/inertial transformations performed during SP3 ephemeris creation.
    const boost::filesystem::path path = writeTemporarySp3File( velocitySp3Contents );
    const std::shared_ptr< input_output::Sp3FileContents > contents = input_output::readSp3File( path.string( ) );
    const Eigen::Vector6d sourceState = contents->satelliteStates.at( "G01" ).begin( )->second;

    const auto itrf2014Settings = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >(
            simulation_setup::sp3EphemerisSettings( contents, "G01", "Earth", "ITRF2014" ) );
    // Require tabulated output so the transformed state history can be inspected.
    BOOST_REQUIRE( itrf2014Settings != nullptr );
    // Confirm the factory labels the result with the requested ITRF realization.
    BOOST_CHECK_EQUAL( itrf2014Settings->getFrameOrientation( ), "ITRF2014" );
    const Eigen::Vector6d itrf2014State = itrf2014Settings->getBodyStateHistory( ).begin( )->second;
    // Confirm the IGS97-to-ITRF2014 Helmert parameters actually alter the Cartesian state.
    BOOST_CHECK_GT( ( itrf2014State - sourceState ).norm( ), 1.0e-3 );
    // Confirm the realization change remains at the expected sub-metre position scale.
    BOOST_CHECK_LT( ( itrf2014State.segment< 3 >( 0 ) - sourceState.segment< 3 >( 0 ) ).norm( ), 1.0 );

    const auto j2000Settings = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >(
            simulation_setup::sp3EphemerisSettings( contents, "G01", "Earth", "J2000" ) );
    // Require tabulated output before checking the Earth-orientation state transformation.
    BOOST_REQUIRE( j2000Settings != nullptr );
    const Eigen::Vector6d j2000State = j2000Settings->getBodyStateHistory( ).begin( )->second;
    // Confirm a pure frame rotation preserves position magnitude.
    BOOST_CHECK_CLOSE_FRACTION( j2000State.segment< 3 >( 0 ).norm( ), itrf2014State.segment< 3 >( 0 ).norm( ), 1.0e-14 );
    // Confirm the rotating-frame velocity term is included, rather than rotating position and velocity independently.
    BOOST_CHECK_GT( ( j2000State.segment< 3 >( 3 ) - sourceState.segment< 3 >( 3 ) ).norm( ), 100.0 );

    const std::shared_ptr< input_output::Sp3FileContents > j2000Contents = std::make_shared< input_output::Sp3FileContents >( *contents );
    j2000Contents->frameName = "J2000";
    for( const auto& state : j2000Settings->getBodyStateHistory( ) )
    {
        j2000Contents->satelliteStates.at( "G01" ).at( state.first ) = state.second;
    }
    const auto roundTripSettings = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >(
            simulation_setup::sp3EphemerisSettings( j2000Contents, "G01", "Earth", "IGS97" ) );
    // Require the reverse conversion to produce inspectable tabulated settings.
    BOOST_REQUIRE( roundTripSettings != nullptr );
    // Confirm terrestrial-to-inertial-to-terrestrial conversion recovers the original state within numerical tolerance.
    BOOST_CHECK_SMALL( ( roundTripSettings->getBodyStateHistory( ).begin( )->second - sourceState ).norm( ), 1.0e-3 );

    const std::shared_ptr< input_output::Sp3FileContents > aliasContents = std::make_shared< input_output::Sp3FileContents >( *contents );
    aliasContents->frameName = "IGb14";
    const auto aliasSettings = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >(
            simulation_setup::sp3EphemerisSettings( aliasContents, "G01", "Earth", "ITRF2014" ) );
    // Require the IGS alias conversion to return tabulated output.
    BOOST_REQUIRE( aliasSettings != nullptr );
    // Confirm IGb14 is recognized as the ITRF2014 realization and therefore needs no numerical transformation.
    BOOST_CHECK( aliasSettings->getBodyStateHistory( ).begin( )->second.isApprox( sourceState, 0.0 ) );

    const std::shared_ptr< input_output::Sp3FileContents > genericItrsContents =
            std::make_shared< input_output::Sp3FileContents >( *contents );
    genericItrsContents->frameName = "ITRS";
    // Confirm generic ITRS cannot be converted to a specific realization without realization metadata.
    BOOST_CHECK_THROW( simulation_setup::sp3EphemerisSettings( genericItrsContents, "G01", "Earth", "ITRF2014" ), std::invalid_argument );

    const std::shared_ptr< input_output::Sp3FileContents > wgs84Contents = std::make_shared< input_output::Sp3FileContents >( *contents );
    wgs84Contents->frameName = "WGS84";
    // Confirm WGS84 is not silently treated as an ITRF realization when no unambiguous Helmert model exists.
    BOOST_CHECK_THROW( simulation_setup::sp3EphemerisSettings( wgs84Contents, "G01", "Earth", "ITRF2014" ), std::invalid_argument );

    const std::shared_ptr< input_output::Sp3FileContents > unknownFrameContents =
            std::make_shared< input_output::Sp3FileContents >( *contents );
    unknownFrameContents->frameName = "LOCAL";
    // Confirm an unknown frame tag is rejected instead of applying an assumed Earth-orientation rotation.
    BOOST_CHECK_THROW( simulation_setup::sp3EphemerisSettings( unknownFrameContents, "G01", "Earth", "J2000" ), std::invalid_argument );
    boost::filesystem::remove( path );
}

BOOST_AUTO_TEST_CASE( testTudatEphemerisFactories )
{
    // This test compares the parsed-contents and filename factory overloads and verifies their input validation.
    const boost::filesystem::path path = writeTemporarySp3File( velocitySp3Contents );
    const std::shared_ptr< input_output::Sp3FileContents > contents = input_output::readSp3File( path.string( ) );

    const std::shared_ptr< simulation_setup::EphemerisSettings > fromContents = simulation_setup::sp3EphemerisSettings( contents, "G01" );
    const std::shared_ptr< simulation_setup::EphemerisSettings > fromFile = simulation_setup::sp3EphemerisSettings( path.string( ), "G01" );

    const auto tabulatedFromContents = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >( fromContents );
    const auto tabulatedFromFile = std::dynamic_pointer_cast< simulation_setup::TabulatedEphemerisSettings >( fromFile );
    // Require the parsed-contents overload to produce tabulated ephemeris settings.
    BOOST_REQUIRE( tabulatedFromContents != nullptr );
    // Require the filename overload to produce the same settings type.
    BOOST_REQUIRE( tabulatedFromFile != nullptr );
    // Confirm the factory default origin is Earth, matching the geocentric SP3 state definition.
    BOOST_CHECK_EQUAL( tabulatedFromContents->getFrameOrigin( ), "Earth" );
    // Confirm the factory default orientation is taken from the SP3 header.
    BOOST_CHECK_EQUAL( tabulatedFromContents->getFrameOrientation( ), "IGS97" );
    // Require both overloads to preserve the same number of epochs before comparing their states.
    BOOST_REQUIRE_EQUAL( tabulatedFromContents->getBodyStateHistory( ).size( ), tabulatedFromFile->getBodyStateHistory( ).size( ) );
    // Confirm parsing before the call and parsing inside the filename overload yield identical state data.
    BOOST_CHECK( tabulatedFromContents->getBodyStateHistory( ).begin( )->second.isApprox(
            tabulatedFromFile->getBodyStateHistory( ).begin( )->second ) );

    // Confirm requesting an absent satellite fails before constructing invalid ephemeris settings.
    BOOST_CHECK_THROW( simulation_setup::sp3EphemerisSettings( contents, "G99" ), std::invalid_argument );

    contents->satelliteStates.at( "G01" ).begin( )->second( 3 ) = std::numeric_limits< double >::quiet_NaN( );
    // Confirm missing state components are rejected because tabulated ephemerides require finite states.
    BOOST_CHECK_THROW( simulation_setup::sp3EphemerisSettings( contents, "G01" ), std::runtime_error );
    boost::filesystem::remove( path );
}

BOOST_AUTO_TEST_CASE( testSupportedVersionsTimeSystemsAndFrameMetadata )
{
    // This test varies independent header fields to verify all supported versions/time systems and opaque frame metadata are parsed.
    for( const char version : { 'a', 'b', 'c', 'd' } )
    {
        // Adapt identifiers and the legacy placeholder time tag to the grammar used by SP3-a/b fixtures.
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
        // Confirm each supported version reaches the parsed metadata unchanged.
        BOOST_CHECK_EQUAL( parsed->formatVersion, version );
        // Confirm the legacy "ccc" placeholder is interpreted as GPS while modern files retain GPS directly.
        BOOST_CHECK_EQUAL( parsed->timeScale, "GPS" );
        // Confirm legacy numeric and modern constellation-prefixed satellite identifiers are both retained.
        BOOST_CHECK_EQUAL( parsed->satelliteStates.count( version < 'c' ? "1" : "G01" ), 1 );
        boost::filesystem::remove( path );
    }

    const std::vector< std::string > timeSystems = { "GPS", "GLO", "GAL", "BDT", "TAI", "UTC", "IRN", "QZS" };
    for( const std::string& timeSystem : timeSystems )
    {
        std::string contents = velocitySp3Contents;
        contents.replace( contents.find( "GPS" ), 3, timeSystem );
        const boost::filesystem::path path = writeTemporarySp3File( contents );
        // Confirm every time-system code accepted by the reader is exposed unchanged for later environment setup.
        BOOST_CHECK_EQUAL( input_output::readSp3File( path.string( ) )->timeScale, timeSystem );
        boost::filesystem::remove( path );
    }

    for( const std::string frameName : { "IGS97", "IGb14", "WGS84", "GCRF " } )
    {
        std::string contents = velocitySp3Contents;
        contents.replace( contents.find( "IGS97" ), 5, frameName );
        const boost::filesystem::path path = writeTemporarySp3File( contents );
        // Confirm data_input trims fixed-width padding but otherwise preserves the frame tag for environment-level interpretation.
        BOOST_CHECK_EQUAL( input_output::readSp3File( path.string( ) )->frameName, boost::algorithm::trim_copy( frameName ) );
        boost::filesystem::remove( path );
    }
}

BOOST_AUTO_TEST_CASE( testMultipleSatelliteListLinesAndConstellations )
{
    // This test uses an 18-satellite mixed GPS/Galileo product to exercise '+' continuation records and constellation identifiers.
    const boost::filesystem::path path = writeTemporarySp3File( makeMultiConstellationSp3Contents( ) );
    const auto contents = input_output::readSp3File( path.string( ) );
    // Confirm the declared count is read from the first satellite-list record.
    BOOST_CHECK_EQUAL( contents->declaredNumberOfSatellites, 18 );
    // Confirm identifiers from both the first and continuation records receive state histories.
    BOOST_CHECK_EQUAL( contents->satelliteStates.size( ), 18 );
    // Confirm the first GPS identifier is not lost while parsing the initial '+' record.
    BOOST_CHECK_EQUAL( contents->satelliteStates.count( "G01" ), 1 );
    // Confirm the final Galileo identifier is read from the continuation '+' record.
    BOOST_CHECK_EQUAL( contents->satelliteStates.count( "E09" ), 1 );
    boost::filesystem::remove( path );
}

BOOST_AUTO_TEST_CASE( testBadOrAbsentStateValuesAreMarkedNonFinite )
{
    // This test replaces either a P or V record with SP3's all-zero missing-value marker and checks downstream rejection.
    for( const char missingRecordType : { 'P', 'V' } )
    {
        std::string contents = velocitySp3Contents;
        const std::string recordPrefix = std::string( 1, missingRecordType ) + "G01";
        const std::string::size_type recordStart = contents.find( recordPrefix );
        const std::string::size_type recordEnd = contents.find( '\n', recordStart );
        contents.replace( recordStart, recordEnd - recordStart + 1, makeStateRecord( missingRecordType, "G01", 0.0, 0.0, 0.0 ) );

        const boost::filesystem::path path = writeTemporarySp3File( contents );
        const auto parsed = input_output::readSp3File( path.string( ) );
        const Eigen::Vector6d& firstState = parsed->satelliteStates.at( "G01" ).begin( )->second;
        const int missingSegmentStart = missingRecordType == 'P' ? 0 : 3;
        // Confirm the all-zero marker becomes a non-finite position or velocity segment instead of a physical zero vector.
        BOOST_CHECK( !firstState.segment< 3 >( missingSegmentStart ).allFinite( ) );
        // Confirm ephemeris construction rejects the incomplete state history before interpolation.
        BOOST_CHECK_THROW( simulation_setup::sp3EphemerisSettings( parsed, "G01" ), std::runtime_error );
        boost::filesystem::remove( path );
    }
}

BOOST_AUTO_TEST_CASE( testMalformedAndIncompleteFilesAreRejected )
{
    // This test corrupts one structural or numeric requirement at a time and verifies strict whole-file validation.
    std::string contents = velocitySp3Contents;
    contents.replace( contents.find( "       2 ORBIT" ), 8, "       3" );
    // Confirm a declared epoch count that exceeds the parsed count is rejected.
    checkReadThrows( contents );

    contents = velocitySp3Contents;
    contents.at( 1 ) = 'e';
    // Confirm unsupported versions are rejected before version-dependent records are interpreted.
    checkReadThrows( contents );

    contents = velocitySp3Contents;
    contents.at( 2 ) = 'X';
    // Confirm the first header must declare either position-only or position-and-velocity mode.
    checkReadThrows( contents );

    contents = velocitySp3Contents;
    contents.replace( contents.find( "       2 ORBIT" ), 8, "      2X" );
    // Confirm integer fields with trailing non-numeric characters are not accepted as valid prefixes.
    checkReadThrows( contents );

    contents = velocitySp3Contents;
    contents.replace( contents.find( " -11044.805800" ), 14, " -11044.80580X" );
    // Confirm floating-point state fields with trailing non-numeric characters are rejected.
    checkReadThrows( contents );

    contents = velocitySp3Contents;
    contents.replace( contents.find( "GPS" ), 3, "UNK" );
    // Confirm an unknown time-system code is rejected because later frame conversion could not assign epochs safely.
    checkReadThrows( contents );

    contents = velocitySp3Contents;
    contents.erase( contents.find( "VG01" ), contents.find( '\n', contents.find( "VG01" ) ) - contents.find( "VG01" ) + 1 );
    // Confirm every epoch in a V-mode file contains a velocity record for each declared satellite.
    checkReadThrows( contents );

    contents = positionOnlySp3Contents;
    const std::string::size_type lastEpoch = contents.rfind( "*  " );
    contents.erase( lastEpoch, contents.find( "EOF", lastEpoch ) - lastEpoch );
    contents.replace( contents.find( "       3 ORBIT" ), 8, "       2" );
    // Confirm position-only files with fewer than three epochs cannot claim second-order derived velocities.
    checkReadThrows( contents );

    contents = positionOnlySp3Contents;
    contents.insert( contents.find( "EOF" ), "VE11      1.000000      2.000000      3.000000      0.000000\n" );
    // Confirm a velocity record is rejected when the header declares a position-only product.
    checkReadThrows( contents );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
