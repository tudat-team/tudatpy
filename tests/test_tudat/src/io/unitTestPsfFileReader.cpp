/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

// #define BOOST_TEST_DYN_LINK
// #define BOOST_TEST_MAIN

#include <vector>

#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readPsfFile.h"

// namespace tudat
// {
// namespace unit_tests
// {

using namespace tudat;
using namespace tudat::input_output;
using namespace tudat::input_output::psf;
//
// BOOST_AUTO_TEST_SUITE( test_psf_file_reader )
//
// BOOST_AUTO_TEST_CASE( testSinglePsfFileReader )
int main( )
{
    std::cout << "A" << std::endl;
    std::string file = "/home/dominic/Downloads/psf_vgr2_neptune.txt";
    RawPsfFileContents psfFile = readPsfFile( file );
}

// BOOST_AUTO_TEST_SUITE_END( )

// }  // namespace unit_tests
//
// }  // namespace tudat
