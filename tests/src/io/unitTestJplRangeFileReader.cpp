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

BOOST_AUTO_TEST_CASE( Early )
{
    std::string file = tudat::paths::getTudatTestDataPath( ) + "mars-pathfinder-range.txt";
    std::shared_ptr< input_output::JplRangeFileContents > jplRangeContents = std::make_shared< input_output::JplRangeFileContents >(
            file );

    std::cout << "Test case with file " << jplRangeContents->getFileName( ) << "\n";

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

BOOST_AUTO_TEST_SUITE_END( );
}
}
