/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved.
 *    This file is part of Tudat and is distributed under the Modified BSD license.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "tudat/io/readDsnAmcFile.h"

namespace tudat
{
namespace unit_tests
{
namespace
{

class TemporaryAmcFile
{
public:
    TemporaryAmcFile( const std::string& contents )
    {
        path_ = std::filesystem::temp_directory_path( ) /
                ( "tudat_amc_" + std::to_string( std::chrono::steady_clock::now( ).time_since_epoch( ).count( ) ) + ".txt" );
        std::ofstream file( path_ );
        file << contents;
    }

    ~TemporaryAmcFile( )
    {
        std::error_code error;
        std::filesystem::remove( path_, error );
    }

    std::string name( ) const
    {
        return path_.string( );
    }

private:
    std::filesystem::path path_;
};

// Actual PJ6 rows, retaining an intentional gap and the undocumented surface flag 3.
const std::string pj6FirstRow = "P 1495152060 -076.742 0168.835 0004.186 0205.390 -00005 010 0 0 1 0 0 2 007.6 295\n";
const std::string pj6LaterRow = "P 1495173540 0024.985 0132.403 0004.261 0206.209 -00002 010 0 0 0 3 0 0 009.6 023\n";

}  // namespace

BOOST_AUTO_TEST_SUITE( test_read_dsn_amc_file )

BOOST_AUTO_TEST_CASE( readJunoRecordsWithoutQualityFiltering )
{
    TemporaryAmcFile file(
            "# Path delay data file\r\n# wet pdzen in cm\r\n# dry pdzen in cm\r\n"
            "# time az el wetdly drydly cloud alg f f f f f f ws wd\r\n\n" +
            pj6FirstRow + pj6LaterRow );
    const auto data = input_output::readDsnAmcFile( file.name( ) );
    BOOST_CHECK_EQUAL( data->getFileName( ), file.name( ) );
    BOOST_REQUIRE_EQUAL( data->getRecords( ).size( ), 2 );
    const auto& first = data->getRecords( ).front( );
    BOOST_CHECK_EQUAL( first.utcTime, 548424060.0 );
    BOOST_CHECK_CLOSE_FRACTION( first.zenithWetDelay, 0.04186, 1.0e-14 );
    BOOST_CHECK_CLOSE_FRACTION( first.zenithDryDelay, 2.05390, 1.0e-14 );
    BOOST_CHECK_CLOSE_FRACTION( first.azimuth, -76.742, 1.0e-14 );
    BOOST_CHECK_CLOSE_FRACTION( first.elevation, 168.835, 1.0e-14 );
    BOOST_CHECK_EQUAL( first.cloudLiquid, -5.0 );
    BOOST_CHECK_EQUAL( first.algorithmId, "010" );
    BOOST_CHECK_EQUAL( first.mtpFlag, 1 );
    BOOST_CHECK_EQUAL( first.pathDelayFlag, 2 );
    BOOST_CHECK_EQUAL( first.cloudFlag, 0 );
    BOOST_CHECK_EQUAL( first.wvrFlag, 0 );
    BOOST_CHECK_EQUAL( first.windFlag, 0 );
    BOOST_CHECK_EQUAL( first.windSpeed, 7.6 );
    BOOST_CHECK_EQUAL( first.windDirection, 295.0 );
    const auto& last = data->getRecords( ).back( );
    BOOST_CHECK_EQUAL( last.utcTime, 548445540.0 );
    BOOST_CHECK_EQUAL( last.surfaceFlag, 3 );
    BOOST_CHECK_EQUAL( last.utcTime - first.utcTime, 21480.0 );
}

BOOST_AUTO_TEST_CASE( preserveSentinelsAndUnknownFlags )
{
    TemporaryAmcFile file( "P 1481454570 -70.873 164.935 999 206.224 3 002 3 1 2 2 9 1 9.3 257\n" );
    const auto data = input_output::readDsnAmcFile( file.name( ) );
    const auto& retainedRecord = data->getRecords( ).front( );
    BOOST_CHECK_EQUAL( retainedRecord.utcTime, 534726570.0 );
    BOOST_CHECK_EQUAL( retainedRecord.algorithmId, "002" );
    BOOST_CHECK_CLOSE_FRACTION( retainedRecord.zenithWetDelay, 9.99, 1.0e-14 );
    BOOST_CHECK_EQUAL( retainedRecord.surfaceFlag, 2 );
    BOOST_CHECK_EQUAL( retainedRecord.windFlag, 9 );
    BOOST_CHECK_EQUAL( retainedRecord.cloudFlag, 3 );
    BOOST_CHECK_EQUAL( retainedRecord.wvrFlag, 1 );
    BOOST_CHECK_EQUAL( retainedRecord.mtpFlag, 2 );
    BOOST_CHECK_EQUAL( retainedRecord.pathDelayFlag, 1 );
}

BOOST_AUTO_TEST_CASE( rejectMalformedAndUnorderedFiles )
{
    for( const auto& contents : std::vector< std::string >{ "# Empty file\n",
                                                            "Q 1495152060\n",
                                                            "P not-a-time\n",
                                                            pj6FirstRow + "P 1495152090 0 0\n",
                                                            pj6FirstRow + "unexpected record\n",
                                                            pj6FirstRow + pj6FirstRow,
                                                            pj6LaterRow + pj6FirstRow,
                                                            "P 1495152060 0 0 4 205 0 010 0 0 0 0 0 0 2 90 extra\n",
                                                            "P 1495152060 0 0 nan 205 0 010 0 0 0 0 0 0 2 90\n" } )
    {
        TemporaryAmcFile file( contents );
        BOOST_CHECK_THROW( input_output::readDsnAmcFile( file.name( ) ), std::runtime_error );
    }
    TemporaryAmcFile missing( "" );
    std::filesystem::remove( missing.name( ) );
    BOOST_CHECK_THROW( input_output::readDsnAmcFile( missing.name( ) ), std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
