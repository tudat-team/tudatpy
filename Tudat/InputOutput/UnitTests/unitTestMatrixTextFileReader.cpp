/*    Copyright (c) 2010-2018, Delft University of Technology
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

        // Read input file and store data in matrix (with empty lines).
        const Eigen::MatrixXd inputFileMatrix2 = input_output::readMatrixFromFile(
                    input_output::getTudatRootPath( )
                    + "/InputOutput/UnitTests/testMatrix4.txt", ";"  );

        // Check if data input file matrix matches expected matrix.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMatrix, inputFileMatrix2,
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
