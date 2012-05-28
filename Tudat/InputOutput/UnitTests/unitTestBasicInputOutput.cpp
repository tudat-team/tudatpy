/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
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
