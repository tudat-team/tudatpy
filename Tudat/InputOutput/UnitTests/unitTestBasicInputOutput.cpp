/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *
 */

#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/streamFilters.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

//! Reads a text file and tokenizes it into lines.
std::vector< std::string > readLinesFromFile(
        const std::string& dataFileAbsolutePath, const std::string& delimiter = " ",
        const char commentCharacter = '#' )
{
    // Read in data file.
    std::ifstream inputFile;
    inputFile.open( dataFileAbsolutePath.c_str( ) );

    // Declare filtered data string.
    std::string filteredData;

    // If input file could be opened, read contents into data buffer, filter out comment lines,
    // and close file.
    if ( inputFile )
    {
        // Create unfiltered data buffer.
        std::stringstream unfilteredDataBuffer;

        // Read input file contents into data buffer.
        unfilteredDataBuffer << inputFile.rdbuf( );

        // Close input file.
        inputFile.close( );

        // Create a filter chain.
        boost::iostreams::filtering_ostream filterProcessor;

        // Add remove comment lines filter.
        filterProcessor.push( input_output::stream_filters::RemoveComment(
                                  commentCharacter, true ) );

        // Last step in the chain; store the resulting string in filteredData.
        filterProcessor.push( boost::iostreams::back_inserter( filteredData ) );

        // Push in unfiltered data.
        filterProcessor << unfilteredDataBuffer.str( );

        // Make sure the data is processed by the chain.
        filterProcessor.flush( );
    }

    else
    {
        throw std::runtime_error( "Data file could not be opened: " + dataFileAbsolutePath );
    }

    // Trim all stray whitespaces.
    boost::trim( filteredData );

    // Erase all non-whitespace delimiter characters.
    if ( delimiter.compare( " " ) != 0 )
    {
        boost::erase_all( filteredData, delimiter );
    }

    // Tokenize the row data by splitting for any of a list of characters.
    std::vector< std::string > inputFileRowTokens;
    boost::split( inputFileRowTokens, filteredData, boost::is_any_of( "\n\r\t" ),
                  boost::token_compress_on );

    // Return tokenized input file (tokenized per row).
    return inputFileRowTokens;
}

BOOST_AUTO_TEST_SUITE( test_basic_inputoutput )

BOOST_AUTO_TEST_CASE( testPrintingNumberInFormattedScientificNotation )
{
    using input_output::printInFormattedScientificNotation;

    // Set floating-point numbers used for testing.
    const double number1 = 1.23e-0045;
    const double number2 = 9.87e65;
    const double number3 = -3.240112e201;

    // Test 1: format numbers with "digits10" precision for the base and two digits for exponent if
    //         possible (if exponent has more significant digits, the minimum number required
    //         should be printed).
    {
        // Set expected output for first number.
        const std::string expectedOutputForNumber1 = "1.230000000000000E-45";

        // Check that output generated for first number matches expected output.
        BOOST_CHECK_EQUAL( printInFormattedScientificNotation( number1 ),
                           expectedOutputForNumber1 );

        // Set expected output for second number.
        const std::string expectedOutputForNumber2 = "9.870000000000000E+65";

        // Check that output generated for first number matches expected output.
        BOOST_CHECK_EQUAL( printInFormattedScientificNotation( number2 ),
                           expectedOutputForNumber2 );

        // Set expected output for third number.
        const std::string expectedOutputForNumber3 = "-3.240112000000000E+201";

        // Check that output generated for first number matches expected output.
        BOOST_CHECK_EQUAL( printInFormattedScientificNotation( number3 ),
                           expectedOutputForNumber3 );
    }

    // Test 2: format numbers with five-digit precision for the base and four digits for exponent
    //         if possible (if exponent has more significant digits, the minimum number required
    //         should be printed).
    {
        // Set expected output for first number.
        const std::string expectedOutputForNumber1 = "1.23000E-0045";

        // Check that output generated for first number matches expected output.
        BOOST_CHECK_EQUAL( printInFormattedScientificNotation( number1, 5, 4 ),
                           expectedOutputForNumber1 );

        // Set expected output for second number.
        const std::string expectedOutputForNumber2 = "9.87000E+0065";

        // Check that output generated for first number matches expected output.
        BOOST_CHECK_EQUAL( printInFormattedScientificNotation( number2, 5, 4 ),
                           expectedOutputForNumber2 );

        // Set expected output for third number.
        const std::string expectedOutputForNumber3 = "-3.24011E+0201";

        // Check that output generated for first number matches expected output.
        BOOST_CHECK_EQUAL( printInFormattedScientificNotation( number3, 5, 4 ),
                           expectedOutputForNumber3 );
    }
}

// Test if listing all files in specified directory works correctly.
BOOST_AUTO_TEST_CASE( testListAllFilesInDirectory )
{
    // Set path to new directory.
    const boost::filesystem::path pathToNewDirectory(
                input_output::getTudatRootPath( ) + "InputOutput/UnitTests/ListAllFiles" );

    // Set number of files in directory.
    const unsigned int numberOfFiles = 10;

    // Create new directory.
    boost::filesystem::create_directory( pathToNewDirectory );

    // List all files in directory and check that there are none.
    const std::vector< boost::filesystem::path > emptyListOfFilenames =
            input_output::listAllFilesInDirectory( pathToNewDirectory.string( ) );

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
            input_output::listAllFilesInDirectory( pathToNewDirectory.string( ) );
    std::sort( listOfFilenames.begin( ), listOfFilenames.end( ) );

    for ( unsigned int i = 0; i < listOfFilenames.size( ); i++ )
    {
        std::stringstream newFile;
        newFile << "testFile" << i << ".txt" << std::endl;
        BOOST_CHECK_EQUAL( newFile.str( ), listOfFilenames.at( i ).string( ) );
    }

    // Remove new directory.
    boost::filesystem::remove_all( pathToNewDirectory );
}

// Test if writing data map to text file works correctly.
BOOST_AUTO_TEST_CASE( testWriteDataMapToTextFile )
{
    // Set path to output directory.
    const boost::filesystem::path pathToOutputDirectory(
                input_output::getTudatRootPath( ) + "InputOutput/UnitTests/WriteDataMap" );

    // Case 1: write key=double, value=double map to file.
    {
        // Set map of data to write to file.
        using input_output::DoubleKeyTypeDoubleValueTypeMap;
        DoubleKeyTypeDoubleValueTypeMap keyDoubleValueDoubleMap;
        keyDoubleValueDoubleMap[ std::sqrt( 3.0 ) ] = 1.0 / std::sqrt( 2.0 );
        keyDoubleValueDoubleMap[ 4.5 ] = 56.89;
        keyDoubleValueDoubleMap[ 12.65 ] = 1.0 / 3.0;

        // Set delimiter.
        std::string delimiter = ",";

        // Set file header.
        std::string fileHeader = "# This is a comment line.\n";

        // Write map of data to file.
        input_output::writeDataMapToTextFile(
                    keyDoubleValueDoubleMap,
                    "keyDoubleValueDoubleMapDataFile",
                    pathToOutputDirectory,
                    fileHeader,
                    std::numeric_limits< double >::digits10,
                    std::numeric_limits< double >::digits10,
                    delimiter );

        // Set absolute path to data file.
        std::string dataFileAbsolutePath = pathToOutputDirectory.string( )
                + "/keyDoubleValueDoubleMapDataFile";

        // Read data file and tokenize per row.
        std::vector< std::string > inputFileRowTokens
                = readLinesFromFile( dataFileAbsolutePath, delimiter );

        // Declare row counter.
        unsigned int rowCounter = 0;

        // Loop over map data.
        for ( DoubleKeyTypeDoubleValueTypeMap::iterator iteratorDataMap
              = keyDoubleValueDoubleMap.begin( );
              iteratorDataMap != keyDoubleValueDoubleMap.end( );
              iteratorDataMap++ )
        {
            // Tokenize the row data by splitting for whitespaces.
            std::vector< std::string > rowTokens;
            boost::split( rowTokens, inputFileRowTokens.at( rowCounter ), boost::is_any_of( " " ),
                          boost::token_compress_on );

            // Check if map key written to file is as expected.
            BOOST_CHECK_CLOSE_FRACTION( iteratorDataMap->first,
                                        std::stod( rowTokens.at( 0 ) ),
                                        1.0e-14 );

            // Check if map value written to file is as expected.
            BOOST_CHECK_CLOSE_FRACTION( iteratorDataMap->second,
                                        std::stod( rowTokens.at( 1 ) ),
                                        1.0e-15 );

            // Increment row counter.
            rowCounter++;
        }

        // Remove output directory.
        boost::filesystem::remove_all( pathToOutputDirectory );
    }

    // Case 2: write key=double, value=double map to file with default options.
    {
        // Set map of data to write to file.
        using input_output::DoubleKeyTypeDoubleValueTypeMap;
        DoubleKeyTypeDoubleValueTypeMap keyDoubleValueDoubleMap;
        keyDoubleValueDoubleMap[ std::sqrt( 3.0 ) ] = 1.0 / std::sqrt( 2.0 );
        keyDoubleValueDoubleMap[ 4.5 ] = 56.89;
        keyDoubleValueDoubleMap[ 12.65 ] = 1.0 / 3.0;

        // Write map of data to file.
        input_output::writeDataMapToTextFile(
                    keyDoubleValueDoubleMap, "keyDoubleValueDoubleMapDataFileWithDefaults" );

        // Set absolute path to data file.
        std::string dataFileAbsolutePath = input_output::getTudatRootPath( )
                + "/keyDoubleValueDoubleMapDataFileWithDefaults";

        // Read data file and tokenize per row.
        std::vector< std::string > inputFileRowTokens = readLinesFromFile( dataFileAbsolutePath );

        // Declare row counter.
        unsigned int rowCounter = 0;

        // Loop over map data.
        for ( DoubleKeyTypeDoubleValueTypeMap::iterator iteratorDataMap
              = keyDoubleValueDoubleMap.begin( );
              iteratorDataMap != keyDoubleValueDoubleMap.end( );
              iteratorDataMap++ )
        {
            // Tokenize the row data by splitting for whitespaces.
            std::vector< std::string > rowTokens;
            boost::split( rowTokens, inputFileRowTokens.at( rowCounter ), boost::is_any_of( " " ),
                          boost::token_compress_on );

            // Check if map key written to file is as expected.
            BOOST_CHECK_CLOSE_FRACTION( iteratorDataMap->first,
                                        std::stod( rowTokens.at( 0 ) ),
                                        1.0e-14 );

            // Check if map value written to file is as expected.
            BOOST_CHECK_CLOSE_FRACTION( iteratorDataMap->second,
                                        std::stod( rowTokens.at( 1 ) ),
                                        1.0e-15 );

            // Increment row counter.
            rowCounter++;
        }

        // Remove output file.
        boost::filesystem::remove( input_output::getTudatRootPath( )
                                   + "/keyDoubleValueDoubleMapDataFileWithDefaults" );
    }

    // Case 3: write key=int, value=double map to file.
    {
        // Set map of data to write to file.
        using input_output::IntKeyTypeDoubleValueTypeMap;
        IntKeyTypeDoubleValueTypeMap keyIntValueDoubleMap;
        keyIntValueDoubleMap[ 1 ] = 1.0 / std::sqrt( 2.0 );
        keyIntValueDoubleMap[ 7 ] = 56.89;
        keyIntValueDoubleMap[ -9 ] = 1.0 / 3.0;

        // Set delimiter.
        std::string delimiter = ",";

        // Set file header.
        std::string fileHeader = "# This is a comment line.\n";

        // Write map of data to file.
        input_output::writeDataMapToTextFile(
                    keyIntValueDoubleMap,
                    "keyIntValueDoubleMapDataFile",
                    pathToOutputDirectory,
                    fileHeader,
                    std::numeric_limits< int >::digits10,
                    std::numeric_limits< double >::digits10,
                    delimiter );

        // Set absolute path to data file.
        std::string dataFileAbsolutePath = pathToOutputDirectory.string( )
                + "/keyIntValueDoubleMapDataFile";

        // Read data file and tokenize per row.
        std::vector< std::string > inputFileRowTokens
                = readLinesFromFile( dataFileAbsolutePath, delimiter );

        // Declare row counter.
        unsigned int rowCounter = 0;

        // Loop over map data.
        for ( IntKeyTypeDoubleValueTypeMap::iterator iteratorDataMap
              = keyIntValueDoubleMap.begin( );
              iteratorDataMap != keyIntValueDoubleMap.end( );
              iteratorDataMap++ )
        {
            // Tokenize the row data by splitting for whitespaces.
            std::vector< std::string > rowTokens;
            boost::split( rowTokens, inputFileRowTokens.at( rowCounter ), boost::is_any_of( " " ),
                          boost::token_compress_on );

            // Check if map key written to file is as expected.
            BOOST_CHECK_EQUAL( iteratorDataMap->first,
                               std::stoi( rowTokens.at( 0 ) ) );

            // Check if map value written to file is as expected.
            BOOST_CHECK_CLOSE_FRACTION( iteratorDataMap->second,
                                        std::stod( rowTokens.at( 1 ) ),
                                        1.0e-15 );

            // Increment row counter.
            rowCounter++;
        }

        // Remove output directory.
        boost::filesystem::remove_all( pathToOutputDirectory );
    }

    // Case 4: write key=double, value=Vector3d map to file with default options.
    {
        // Set map of data to write to file.
        using input_output::DoubleKeyTypeVector3dValueTypeMap;
        DoubleKeyTypeVector3dValueTypeMap keyDoubleValueVector3dMap;
        keyDoubleValueVector3dMap[ 1.1 ] = Eigen::Vector3d( 0.0, 1.3, -6.54 );
        keyDoubleValueVector3dMap[ 6.5 ] = Eigen::Vector3d( -4.56, 1.23, -9.98 );
        keyDoubleValueVector3dMap[ 10.9 ] = Eigen::Vector3d( -46.13, 1.0 / 3.0, std::sqrt( 2.0 ) );

        // Write map of data to file.
        input_output::writeDataMapToTextFile(
                    keyDoubleValueVector3dMap, "keyDoubleValueVector3dMapDataFile" );

        // Set absolute path to data file.
        std::string dataFileAbsolutePath = input_output::getTudatRootPath( )
                + "/keyDoubleValueVector3dMapDataFile";

        // Read data file and tokenize per row.
        std::vector< std::string > inputFileRowTokens = readLinesFromFile( dataFileAbsolutePath );

        // Declare row counter.
        unsigned int rowCounter = 0;

        // Loop over map data.
        for ( DoubleKeyTypeVector3dValueTypeMap::iterator iteratorDataMap
              = keyDoubleValueVector3dMap.begin( );
              iteratorDataMap != keyDoubleValueVector3dMap.end( );
              iteratorDataMap++ )
        {
            // Tokenize the row data by splitting for whitespaces.
            std::vector< std::string > rowTokens;
            boost::split( rowTokens, inputFileRowTokens.at( rowCounter ), boost::is_any_of( " " ),
                          boost::token_compress_on );

            // Check if map key written to file is as expected.
            BOOST_CHECK_CLOSE_FRACTION( iteratorDataMap->first,
                                        std::stod( rowTokens.at( 0 ) ),
                                        1.0e-15 );

            // Create Vector3d from tokens.
            Eigen::Vector3d valueVector3d( std::stod( rowTokens.at( 1 ) ),
                                           std::stod( rowTokens.at( 2 ) ),
                                           std::stod( rowTokens.at( 3 ) ) );

            // Check if map value written to file is as expected.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( iteratorDataMap->second,
                                               valueVector3d,
                                               1.0e-14 );

            // Increment row counter.
            rowCounter++;
        }

        // Remove output file.
        boost::filesystem::remove( input_output::getTudatRootPath( )
                                   + "/keyDoubleValueVector3dMapDataFile" );
    }

    // Case 4: write key=double, value=Matrix3d map to file.
    {
        // Set map of data to write to file.
        using input_output::DoubleKeyTypeMatrix3dValueTypeMap;
        DoubleKeyTypeMatrix3dValueTypeMap keyDoubleValueMatrix3dMap;

        keyDoubleValueMatrix3dMap[ 1.1 ] << ( Eigen::Matrix3d( )
                                              << 0.0,       1.3,   -6.54,
                                              12.65,     10.23,  0.61,
                                              3.0 / 4.7, 12.345, 70.908 ).finished( );
        keyDoubleValueMatrix3dMap[ 6.5 ] << ( Eigen::Matrix3d( )
                                              << -4.56,     1.23,   -9.98,
                                              5.0 / 7.0, 4.65,   std::sqrt( 7.0 ),
                                              64.65,     -7.645, -1001.2908 ).finished( );
        keyDoubleValueMatrix3dMap[ -10.9 ] << ( Eigen::Matrix3d( )
                                                << -46.13, 1.0 / 3.0, std::sqrt( 2.0 ),
                                                -34.65, 987.1025,  1.0 / 3.0,
                                                12.65,  0.999,     6.544 ).finished( );

        // Set delimiter.
        std::string delimiter = "|";

        // Set file header.
        std::string fileHeader = "# This is a comment line.\n";

        // Write map of data to file.
        input_output::writeDataMapToTextFile(
                    keyDoubleValueMatrix3dMap,
                    "keyDoubleValueMatrix3dMapDataFile",
                    pathToOutputDirectory,
                    fileHeader,
                    std::numeric_limits< double >::digits10,
                    std::numeric_limits< double >::digits10,
                    delimiter );

        // Set absolute path to data file.
        std::string dataFileAbsolutePath = pathToOutputDirectory.string( )
                + "/keyDoubleValueMatrix3dMapDataFile";

        // Read data file and tokenize per row.
        std::vector< std::string > inputFileRowTokens
                = readLinesFromFile( dataFileAbsolutePath, delimiter );

        // Declare row counter.
        unsigned int rowCounter = 0;

        // Loop over map data.
        for ( DoubleKeyTypeMatrix3dValueTypeMap::iterator iteratorDataMap
              = keyDoubleValueMatrix3dMap.begin( );
              iteratorDataMap != keyDoubleValueMatrix3dMap.end( );
              iteratorDataMap++ )
        {
            // Tokenize the row data by splitting for whitespaces.
            std::vector< std::string > rowTokens;
            boost::split( rowTokens, inputFileRowTokens.at( rowCounter ), boost::is_any_of( " " ),
                          boost::token_compress_on );

            // Check if map key written to file is as expected.
            BOOST_CHECK_CLOSE_FRACTION( iteratorDataMap->first,
                                        std::stod( rowTokens.at( 0 ) ),
                                        1.0e-15 );

            // Create Matrix3d from tokens.
            Eigen::Matrix3d valueMatrix3d;
            valueMatrix3d << ( Eigen::Matrix3d( )
                               << std::stod( rowTokens.at( 1 ) ),
                               std::stod( rowTokens.at( 2 ) ),
                               std::stod( rowTokens.at( 3 ) ),
                               std::stod( rowTokens.at( 4 ) ),
                               std::stod( rowTokens.at( 5 ) ),
                               std::stod( rowTokens.at( 6 ) ),
                               std::stod( rowTokens.at( 7 ) ),
                               std::stod( rowTokens.at( 8 ) ),
                               std::stod( rowTokens.at( 9 ) ) ).finished( );

            // Check if map value written to file is as expected.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( iteratorDataMap->second,
                                               valueMatrix3d,
                                               1.0e-14 );

            // Increment row counter.
            rowCounter++;
        }

        // Remove output directory.
        boost::filesystem::remove_all( pathToOutputDirectory );
    }
}

//! Test if matrix is correctly written to a file
BOOST_AUTO_TEST_CASE( testMatrixFileWriting )
{
    // Create random matrix
    Eigen::MatrixXd randomMatrix = Eigen::MatrixXd::Random( 20, 10 );

    // Write matrix to file
    std::string fileName = "randomTestMatrix.save";
    std::string outputPath = input_output::getTudatRootPath( ) + "InputOutput/UnitTests/";
    input_output::writeMatrixToFile( randomMatrix, fileName, 16, outputPath, "\t" );

    //  Read matrix to file and check if it is equal to original
    Eigen::MatrixXd retrievedRandomMatrix = input_output::readMatrixFromFile( outputPath + fileName, "\t" );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( randomMatrix,
                                       retrievedRandomMatrix,
                                       ( 2.0 * std::numeric_limits< double >::epsilon( ) ) );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
