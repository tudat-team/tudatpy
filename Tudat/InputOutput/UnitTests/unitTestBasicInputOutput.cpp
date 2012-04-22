/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120207    K. Kumar          File created.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 */

#define BOOST_TEST_MAIN

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_basic_inputoutput )

//! Test if listing all files in specified directory works correctly.
BOOST_AUTO_TEST_CASE( testListAllFilesInDirectory )
{
    // Set path to new directory.
    boost::filesystem::path pathToNewDirectory(
                tudat::input_output::getTudatRootPath( ) + "InputOutput/UnitTests/TestDirectory" );

    // Set number of files in directory.
    unsigned int numberOfFiles = 10;

    // Create new directory.
    boost::filesystem::create_directory( pathToNewDirectory );

    // List all files in directory and check that there are none.
    std::vector< boost::filesystem::path > emptyListOfFilenames =
            tudat::input_output::listAllFilesInDirectory( pathToNewDirectory.string( ) );

    BOOST_CHECK_EQUAL( emptyListOfFilenames.size( ), 0 );

    // Create test files.
    for ( unsigned int i = 0; i < numberOfFiles; i++ )
    {
        // Create stream new filename.
        std::stringstream newFile;
        newFile << pathToNewDirectory.string( ) << "/testFile" << i << ".txt" << std::endl;

        // Create test file and fill with random contents.
        std::ofstream testFile( newFile.str( ).c_str( ) );
        testFile << "tastes good!\n";
        testFile.close( );
    }

    // List all files in directory and check that they are as expected.
    std::vector< boost::filesystem::path > listOfFilenames =
            tudat::input_output::listAllFilesInDirectory( pathToNewDirectory.string( ) );
    std::sort( listOfFilenames.begin(), listOfFilenames.end() );

    for ( unsigned int i = 0; i < listOfFilenames.size( ); i++ )
    {
        std::stringstream newFile;
        newFile << "testFile" << i << ".txt" << std::endl;
        BOOST_CHECK_EQUAL( newFile.str( ), listOfFilenames.at( i ).string( ) );
    }

    // Remove new directory.
    boost::filesystem::remove_all( pathToNewDirectory );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
