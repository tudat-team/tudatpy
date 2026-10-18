/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include <filesystem>
#include <fstream>

#include "tudat/io/readDsnAmcFile.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_read_dsn_amc_file )

BOOST_AUTO_TEST_CASE( testReadDsnAmcFile )
{
    // Juno PJ6 records (DSS-25), with a nominal record, a record with non-zero quality flag, and a missing retrieval
    std::string fileName = ( std::filesystem::temp_directory_path( ) / "tudat_test_amc.txt" ).string( );
    {
        std::ofstream file( fileName );
        file << "# Path delay data file\n"
             << "# Header lines start with \"#\"\n"
             << "P 1495173540 0024.985 0132.403 0004.261 0206.209 -00002 010 0 0 0 0 0 0 009.6 023\n"
             << "P 1495173570 0024.985 0132.403 0004.263 0206.209 -00002 010 0 0 0 3 0 0 009.6 023\n"
             << "P 1495173600 0024.985 0132.403 0999.000 0206.209 -00002 010 0 0 0 0 0 0 009.6 023\n";
    }
    std::vector< input_output::DsnAmcRecord > records = input_output::readDsnAmcFile( fileName );
    std::filesystem::remove( fileName );

    BOOST_REQUIRE_EQUAL( records.size( ), 3 );
    // 2017-05-19 06:39:00 UTC, in seconds since J2000
    BOOST_CHECK_EQUAL( records.at( 0 ).utcTime, 548445540.0 );
    BOOST_CHECK_CLOSE_FRACTION( records.at( 0 ).zenithWetDelay, 0.04261, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( records.at( 0 ).zenithDryDelay, 2.06209, 1.0E-14 );
    BOOST_CHECK_EQUAL( records.at( 0 ).algorithmId, "010" );
    BOOST_CHECK_EQUAL( records.at( 0 ).cloudLiquid, -2.0 );
    BOOST_CHECK_EQUAL( records.at( 1 ).surfaceFlag, 3 );

    // Only the first record is valid
    std::map< double, Eigen::Vector2d > zenithDelays = input_output::getValidDsnAmcZenithDelays( records );
    BOOST_REQUIRE_EQUAL( zenithDelays.size( ), 1 );
    BOOST_CHECK_EQUAL( zenithDelays.begin( )->first, 548445540.0 );
    BOOST_CHECK_EQUAL( zenithDelays.begin( )->second( 0 ), records.at( 0 ).zenithDryDelay );
    BOOST_CHECK_EQUAL( zenithDelays.begin( )->second( 1 ), records.at( 0 ).zenithWetDelay );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
