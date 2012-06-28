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
 *      110530    F.M. Engelen      First creation of code.
 *      120519    K. Kumar          Boostified unit test.
 *
 *    References
 *
 *    If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_matrix_text_file_reader )

// Test if matrix text file reader is working correctly.
BOOST_AUTO_TEST_CASE( testMatrixTextFileReader )
{
    // Set expected matrix.
    const Eigen::MatrixXd expectedMatrix = ( Eigen::MatrixXd( 4, 3 )
                                             << 1.0,  2.0,  3.0,
                                                4.0,  5.0,  6.0,
                                                7.0,  8.0,  9.0,
                                                10.0, 11.0, 12.0 ).finished( );

    // Test 1: test for semi-colon-separated files.
    {
        // Read input file and store data in matrix.
        const Eigen::MatrixXd inputFileMatrix = tudat::input_output::readMatrixFromFile(
                    tudat::input_output::getTudatRootPath( )
                    + "/InputOutput/UnitTests/testMatrix.txt", ";"  );

        // Check if data input file matrix matches expected matrix.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMatrix, inputFileMatrix,
                                           std::numeric_limits< double >::epsilon( ) );

    }

    // Test 2: test for space-separated files.
    {
        // Read input file and store data in matrix.
        const Eigen::MatrixXd inputFileMatrix = tudat::input_output::readMatrixFromFile(
                    tudat::input_output::getTudatRootPath( )
                    + "/InputOutput/UnitTests/testMatrix2.txt", " \t", "#" );

        // Check if data input file matrix matches expected matrix.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMatrix, inputFileMatrix,
                                           std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
