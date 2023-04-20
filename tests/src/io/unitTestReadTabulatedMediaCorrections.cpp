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

#include "tudat/basics/testMacros.h"

#include "tudat/io/readTabulatedMediaCorrections.h"

namespace tudat
{
namespace unit_tests
{

using namespace input_output;

BOOST_AUTO_TEST_SUITE( test_read_tabulated_media_corrections )

BOOST_AUTO_TEST_CASE( testName )
{

//    std::vector< std::string > cspCommands = readCspCommandsFile( "/Users/pipas/Documents/mro-data/tro/mromagr2017_091_2017_121.tro.txt" );
//
//    std::cout << cspCommands.at(0) << std::endl << std::endl;
//
//    std::vector< std::string > vectorOfIndividualStrings;
//    boost::algorithm::split( vectorOfIndividualStrings,
//                             cspCommands.at(0),
//                             boost::algorithm::is_any_of( ",()" ),
//                             boost::algorithm::token_compress_on );
//
//    for ( int i = 0; i < vectorOfIndividualStrings.size(); ++i )
//    {
//        std::cout << vectorOfIndividualStrings.at(i) << std::endl;
//    }
//
//    tdbName( vectorOfIndividualStrings );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat