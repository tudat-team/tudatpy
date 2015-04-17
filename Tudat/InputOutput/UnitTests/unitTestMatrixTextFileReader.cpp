/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      110530    F.M. Engelen      Code created.
 *      120519    K. Kumar          Boostified unit test.
 *      130109    K. Kumar          Ported from Tudat; added note that Test 1 also tests for spaces
 *                                  at start and end of line.
 *      130111    K. Kumar          Added unit tests to test that errors are thrown when matrix is
 *                                  incomplete or input file doesn't exist.
 *
 *    References
 *
 *    Notes
 *      If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>
#include <stdexcept>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"

namespace tudat
{
namespace unit_tests
{

// Test if matrix text file reader is working correctly.
BOOST_AUTO_TEST_CASE( testMatrixTextFileReader )
{
    // Set expected matrix.
    const Eigen::MatrixXd expectedMatrix = ( Eigen::MatrixXd( 4, 3 )
                                             << 1.0,  2.0,  3.0,
                                                4.0,  5.0,  6.0,
                                                7.0,  8.0,  9.0,
                                                10.0, 11.0, 12.0 ).finished( );

    // Test 1: Test for semi-colon-separated files. This test also ensures that lines containing
    //         spaces at the start and end are trimmed correctly.
    {
        // Read input file and store data in matrix.
        const Eigen::MatrixXd inputFileMatrix = input_output::readMatrixFromFile(
                    input_output::getTudatRootPath( )
                    + "/InputOutput/UnitTests/testMatrix.txt", ";"  );

        // Check if data input file matrix matches expected matrix.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMatrix, inputFileMatrix,
                                           std::numeric_limits< double >::epsilon( ) );

    }

    // Test 2: Test for space-separated files.
    {
        // Read input file and store data in matrix.
        const Eigen::MatrixXd inputFileMatrix = input_output::readMatrixFromFile(
                    input_output::getTudatRootPath( )
                    + "/InputOutput/UnitTests/testMatrix2.txt", " \t", "#" );

        // Check if data input file matrix matches expected matrix.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMatrix, inputFileMatrix,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test 3: Test that runtime error is thrown if matrix is incomplete.
    {
        // Set flag indicating if runtime error is thrown to false.
        bool isRuntimeError = false;

        // Try to read input file.
        try
        {
            // Read input file and store data in matrix.
            const Eigen::MatrixXd inputFileMatrix = input_output::readMatrixFromFile(
                        input_output::getTudatRootPath( )
                        + "/InputOutput/UnitTests/testMatrix3.txt", " \t", "#" );
        }

        // Catch expected runtime error.
        catch( std::runtime_error& matrixIncomplete )
        {
            // Set runtime error flag to true.
            isRuntimeError = true;
        }

        // Check that runtime error was thrown because matrix is incomplete.
        BOOST_CHECK( isRuntimeError );
    }

    // Test 4: Test that runtime error is thrown if input file doesn't exist.
    {
        // Set flag indicating if runtime error is thrown to false.
        bool isRuntimeError = false;

        // Try to read input file.
        try
        {
            // Read fake input file.
            const Eigen::MatrixXd inputFileMatrix = input_output::readMatrixFromFile(
                        input_output::getTudatRootPath( )
                        + "/InputOutput/UnitTests/fakeFile.txt", " \t", "#" );
        }

        // Catch expected runtime error.
        catch( std::runtime_error& fileDoesNotExist )
        {
            // Set runtime error flag to true.
            isRuntimeError = true;
        }

        // Check that runtime error was thrown because file does not exist
        BOOST_CHECK( isRuntimeError );
    }
}

} // namespace unit_tests
} // namespace tudat
