/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>
#include "tudat/basics/testMacros.h"
#include "tudat/io/readJplRangeFile.h"


namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_jpl_range_file_reader );

//! Checks parsed binary file (odf07155.odf) against values in txt file (odf07155.txt)
//! Files available at https://pds-geosciences.wustl.edu/dataserv/radio_science.htm (see radio science documentation bundle)

BOOST_AUTO_TEST_CASE( JplRangeDataBlock )
{

    // The test string below is arbitrary data.
    std::string testString = " 3 15 16 1997 JUL  6 04:32:13   1287957826.329   .067";
    input_output::JplRangeDataBlock dataBlock = input_output::JplRangeDataBlock( testString );

    BOOST_CHECK_EQUAL( dataBlock.spacecraftId_, 3 );
    BOOST_CHECK_EQUAL( dataBlock.transmittingStationId_, 15 );
    BOOST_CHECK_EQUAL( dataBlock.receivingStationId_, 16 );
    BOOST_CHECK_EQUAL( dataBlock.utcYear_, 1997 );
    BOOST_CHECK_EQUAL( dataBlock.utcMonth_, 7 );
    BOOST_CHECK_EQUAL( dataBlock.utcHour_, 4 );
    BOOST_CHECK_EQUAL( dataBlock.utcMinute_, 32 );
    BOOST_CHECK_EQUAL( dataBlock.utcSecond_, 13 );
    BOOST_CHECK_EQUAL( dataBlock.roundTripLightTime_, 1287957826.329 );
    BOOST_CHECK_EQUAL( dataBlock.measurementAccuracy_, 0.067 );
}


BOOST_AUTO_TEST_CASE( JplRangeFile )
{
    std::string file = tudat::paths::getTudatTestDataPath( ) + "mars-pathfinder-range.txt";
    std::shared_ptr< input_output::JplRangeFileContents > jplRangeContents = input_output::readJplRangeFile( file );

    BOOST_CHECK_EQUAL( file, jplRangeContents->fileName_ );

    // Check the fourth line in the file.
    input_output::JplRangeDataBlock dataBlock3 = *jplRangeContents->dataBlocks_.at( 3 );
    BOOST_CHECK_EQUAL( dataBlock3.spacecraftId_, 3 );
    BOOST_CHECK_EQUAL( dataBlock3.transmittingStationId_, 65 );
    BOOST_CHECK_EQUAL( dataBlock3.receivingStationId_, 65 );
    BOOST_CHECK_EQUAL( dataBlock3.utcYear_, 1997 );
    BOOST_CHECK_EQUAL( dataBlock3.utcMonth_, 7 );
    BOOST_CHECK_EQUAL( dataBlock3.utcDay_, 25 );
    BOOST_CHECK_EQUAL( dataBlock3.utcHour_, 18 );
    BOOST_CHECK_EQUAL( dataBlock3.utcMinute_, 03 );
    BOOST_CHECK_EQUAL( dataBlock3.utcSecond_, 31 );
    BOOST_CHECK_EQUAL( dataBlock3.roundTripLightTime_, 1420414209.718 - 420 ); // Transponder delay
    BOOST_CHECK_EQUAL( dataBlock3.measurementAccuracy_, .067 );

    // Check the last line in the file.
    input_output::JplRangeDataBlock dataBlockLast = *jplRangeContents->dataBlocks_.back( );
    BOOST_CHECK_EQUAL( dataBlockLast.spacecraftId_, 3 );
    BOOST_CHECK_EQUAL( dataBlockLast.transmittingStationId_, 34 );
    BOOST_CHECK_EQUAL( dataBlockLast.receivingStationId_, 34 );
    BOOST_CHECK_EQUAL( dataBlockLast.utcYear_, 1997 );
    BOOST_CHECK_EQUAL( dataBlockLast.utcMonth_, 9 );
    BOOST_CHECK_EQUAL( dataBlockLast.utcDay_, 25 );
    BOOST_CHECK_EQUAL( dataBlockLast.utcHour_, 05 );
    BOOST_CHECK_EQUAL( dataBlockLast.utcMinute_, 35 );
    BOOST_CHECK_EQUAL( dataBlockLast.utcSecond_, 54 );
    BOOST_CHECK_EQUAL( dataBlockLast.roundTripLightTime_, 1762272411.593 - 420); // Transponder delay
    BOOST_CHECK_EQUAL( dataBlockLast.measurementAccuracy_, .133 );
}



BOOST_AUTO_TEST_SUITE_END( );
}
}
