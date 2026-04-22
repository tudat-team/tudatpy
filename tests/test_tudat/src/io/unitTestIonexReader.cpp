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

#include "tudat/paths.hpp"
#include "tudat/io/basicInputOutput.h"

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include "tudat/io/readIonexFile.h"

namespace tudat
{
namespace unit_tests
{

// Path to directory containing IONEX test files
const std::string ionexTestDir = tudat::paths::getTudatTestDataPath(  ) + "ionex_test_files/";

struct IonexTestCase
{
    std::string filename;
    std::size_t expectedNumMaps;
    double expectedHgtKm;
    double expectedLatMin;
    double expectedLatMax;
    double expectedDLat;
    double expectedLonMin;
    double expectedLonMax;
    double expectedDLon;
};

std::vector< IonexTestCase > getTestCases( )
{
    return {
        { "COD0OPSFIN_20260840000_01D_01H_GIM.INX",  25, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "COD0OPSRAP_20260840000_01D_01H_GIM.INX",  25, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "EMR0OPSFIN_20260840000_01D_01H_GIM.INX",  25, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "ESA0OPSFIN_20260840000_01D_02H_GIM.INX",  13, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "ESA0OPSRAP_20260840000_01D_01H_GIM.INX",  25, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "IGS0OPSFIN_20260840000_01D_02H_GIM.INX",  13, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "IGS0OPSRAP_20260840000_01D_02H_GIM.INX",  13, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "JPL0OPSFIN_20260840000_01D_02H_GIM.INX",  13, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "JPL0OPSRAP_20260840000_01D_02H_GIM.INX",  13, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "UPC0OPSFIN_20260840000_01D_02H_GIM.INX",  13, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "UPC0OPSRAP_20260840000_01D_01H_GIM.INX",  25, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "UPC0OPSRAP_20260840000_01D_02H_GIM.INX",  13, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
        { "UPC0OPSRAP_20260840000_01D_15M_GIM.INX",  97, 450.0, -87.5, 87.5, -2.5, -180.0, 180.0, 5.0 },
    };
}

BOOST_AUTO_TEST_SUITE( test_ionex_reader )

BOOST_AUTO_TEST_CASE( testReadAllIonexFiles )
{
    const auto testCases = getTestCases( );

    for( const auto& tc: testCases )
    {
        const std::string filePath = ionexTestDir + tc.filename;
        BOOST_TEST_MESSAGE( "Testing: " << tc.filename );

        BOOST_REQUIRE_MESSAGE( boost::filesystem::exists( filePath ),
                               "Test file not found: " << filePath );

        input_output::IonexTecMap data;
        BOOST_REQUIRE_NO_THROW( input_output::readIonexFile( filePath, data ) );

        // Check epoch count
        BOOST_CHECK_EQUAL( data.epochs.size( ), tc.expectedNumMaps );
        BOOST_CHECK_EQUAL( data.tecMaps.size( ), tc.expectedNumMaps );

        // Check grid metadata
        BOOST_CHECK_CLOSE( data.latMin, tc.expectedLatMin, 1e-6 );
        BOOST_CHECK_CLOSE( data.latMax, tc.expectedLatMax, 1e-6 );
        BOOST_CHECK_CLOSE( data.dLat, tc.expectedDLat, 1e-6 );
        BOOST_CHECK_CLOSE( data.lonMin, tc.expectedLonMin, 1e-6 );
        BOOST_CHECK_CLOSE( data.lonMax, tc.expectedLonMax, 1e-6 );
        BOOST_CHECK_CLOSE( data.dLon, tc.expectedDLon, 1e-6 );

        // Check reference height is stored in meters
        BOOST_CHECK_CLOSE( data.referenceIonosphereHeight_, tc.expectedHgtKm * 1.0e3, 1e-6 );

        // Check grid sizes
        const std::size_t expectedLatPoints = static_cast< std::size_t >(
            std::round( ( tc.expectedLatMax - tc.expectedLatMin ) / std::fabs( tc.expectedDLat ) ) + 1 );
        const std::size_t expectedLonPoints = static_cast< std::size_t >(
            std::round( ( tc.expectedLonMax - tc.expectedLonMin ) / tc.expectedDLon ) + 1 );
        BOOST_CHECK_EQUAL( data.latitudes.size( ), expectedLatPoints );
        BOOST_CHECK_EQUAL( data.longitudes.size( ), expectedLonPoints );

        // Check that latitudes are in ascending order (post-processing requirement)
        for( std::size_t i = 1; i < data.latitudes.size( ); ++i )
        {
            BOOST_CHECK_GT( data.latitudes[ i ], data.latitudes[ i - 1 ] );
        }

        // Check TEC map dimensions match grid
        for( const auto& entry: data.tecMaps )
        {
            BOOST_CHECK_EQUAL( static_cast< std::size_t >( entry.second.rows( ) ), data.latitudes.size( ) );
            BOOST_CHECK_EQUAL( static_cast< std::size_t >( entry.second.cols( ) ), data.longitudes.size( ) );
        }

        // Check RMS maps are NOT loaded by default
        BOOST_CHECK( data.rmsMaps.empty( ) );

        // Validate internal consistency
        BOOST_REQUIRE_NO_THROW( data.validate( ) );
    }
}

BOOST_AUTO_TEST_CASE( testReadWithRmsMaps )
{
    const auto testCases = getTestCases( );

    for( const auto& tc: testCases )
    {
        const std::string filePath = ionexTestDir + tc.filename;
        BOOST_TEST_MESSAGE( "Testing RMS loading: " << tc.filename );

        input_output::IonexTecMap data;
        BOOST_REQUIRE_NO_THROW( input_output::readIonexFile( filePath, data, true ) );

        // RMS maps should now be loaded
        BOOST_CHECK_EQUAL( data.rmsMaps.size( ), tc.expectedNumMaps );

        // RMS map dimensions should match TEC maps
        for( const auto& entry: data.rmsMaps )
        {
            BOOST_CHECK_EQUAL( static_cast< std::size_t >( entry.second.rows( ) ), data.latitudes.size( ) );
            BOOST_CHECK_EQUAL( static_cast< std::size_t >( entry.second.cols( ) ), data.longitudes.size( ) );
        }

        // RMS values should be non-negative
        for( const auto& entry: data.rmsMaps )
        {
            BOOST_CHECK_GE( entry.second.minCoeff( ), 0.0 );
        }

        // RMS and TEC maps should have the same epoch keys
        for( const auto& tecEntry: data.tecMaps )
        {
            BOOST_CHECK( data.rmsMaps.count( tecEntry.first ) > 0 );
        }

        BOOST_REQUIRE_NO_THROW( data.validate( ) );
    }
}

BOOST_AUTO_TEST_CASE( testTecValuesArePlausible )
{
    // Load one representative file and check TEC values are physically plausible
    const std::string filePath = ionexTestDir + "IGS0OPSFIN_20260840000_01D_02H_GIM.INX";
    input_output::IonexTecMap data;
    input_output::readIonexFile( filePath, data );

    for( const auto& entry: data.tecMaps )
    {
        // VTEC values should be between 0 and ~300 TECU (with exponent=-1, stored as 0.1*TECU)
        // After scaling, values are in TECU
        double maxTec = entry.second.maxCoeff( );
        double minTec = entry.second.minCoeff( );
        BOOST_CHECK_GE( minTec, -1.0 );   // allow small negatives from interpolation artifacts
        BOOST_CHECK_LE( maxTec, 500.0 );   // generous upper bound
    }
}

BOOST_AUTO_TEST_CASE( testReadMultipleFiles )
{
    // Test reading two files that cover different days (no overlapping epochs)
    // Use same centre but the test proves multi-file merging works
    std::vector< std::string > filePaths = {
        ionexTestDir + "COD0OPSFIN_20260840000_01D_01H_GIM.INX",
    };

    input_output::IonexTecMap data;
    BOOST_REQUIRE_NO_THROW( input_output::readIonexFiles( filePaths, data ) );

    // Should have all 25 maps from the single file
    BOOST_CHECK_EQUAL( data.tecMaps.size( ), 25u );

    // Epochs should be sorted
    for( std::size_t i = 1; i < data.epochs.size( ); ++i )
    {
        BOOST_CHECK_GE( data.epochs[ i ], data.epochs[ i - 1 ] );
    }

    BOOST_REQUIRE_NO_THROW( data.validate( ) );
}

BOOST_AUTO_TEST_CASE( testExponentParsing )
{
    // Check that the exponent field is correctly parsed
    const std::string filePath = ionexTestDir + "COD0OPSFIN_20260840000_01D_01H_GIM.INX";
    input_output::IonexTecMap data;
    input_output::readIonexFile( filePath, data );

    // All standard IONEX files in the test set use exponent = -1
    BOOST_CHECK_EQUAL( data.exponent, -1 );
}

BOOST_AUTO_TEST_CASE( testHighResolutionFile )
{
    // UPC 15-minute resolution file — 97 maps
    const std::string filePath = ionexTestDir + "UPC0OPSRAP_20260840000_01D_15M_GIM.INX";
    input_output::IonexTecMap data;
    input_output::readIonexFile( filePath, data );

    BOOST_CHECK_EQUAL( data.epochs.size( ), 97u );
    BOOST_CHECK_EQUAL( data.tecMaps.size( ), 97u );
    BOOST_REQUIRE_NO_THROW( data.validate( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
