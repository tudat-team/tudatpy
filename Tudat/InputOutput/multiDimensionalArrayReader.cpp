/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <map>
#include <fstream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "Tudat/InputOutput/multiDimensionalArrayReader.h"

namespace tudat
{

namespace input_output
{

//! Function to parse a block of values read from a file into a multi-array of size 1.
boost::multi_array< double, 1 > parseRawOneDimensionalCoefficientsFromFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock )
{
    boost::multi_array< double, 1 > coefficientMultiarray;

    // Check input consistency
    if( independentVariableSize.size( ) == 1 )
    {
        // Define size of multi-array
        coefficientMultiarray.resize( boost::extents[ independentVariableSize.at( 0 ) ] );

        // Parse data
        for( int i = 0; i < independentVariableSize.at( 0 ); i++ )
        {
            coefficientMultiarray[ i ] = coefficientsBlock( i, 0 );
        }
    }
    else
    {
        throw std::runtime_error( "Error, expected size of 1 dimension when parsing 1-dimensional data into multi-array" );
    }
    return coefficientMultiarray;
}

//! Function to parse a block of values read from a file into a multi-array of size 2.
boost::multi_array< double, 2 > parseRawTwoDimensionalCoefficientsFromFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock )
{
    boost::multi_array< double, 2 > coefficientMultiarray;

    // Check input consistency
    if( independentVariableSize.size( ) == 2 )
    {
        // Define size of multi-array
        coefficientMultiarray.resize( boost::extents[ independentVariableSize.at( 0 ) ][ independentVariableSize.at( 1 ) ]  );

        // Parse data
        for( int i = 0; i < independentVariableSize.at( 0 ); i++ )
        {
            for( int j = 0; j < independentVariableSize.at( 1 ); j++ )
            {
                coefficientMultiarray[ i ][ j ] = coefficientsBlock( i, j );
            }
        }
    }
    else
    {
        throw std::runtime_error( "Error, expected size of 2 dimensions when parsing 1-dimensional data into multi-array" );
    }
    return coefficientMultiarray;
}

//! Function to parse a block of values read from a file into a multi-array of size 3.
boost::multi_array< double, 3 > parseRawThreeDimensionalCoefficientsFromFile(
        const std::vector< int > independentVariableSize,
        const Eigen::MatrixXd& coefficientsBlock )
{
    boost::multi_array< double, 3 > coefficientMultiarray;

    // Check input consistency
    if( independentVariableSize.size( ) == 3 )
    {
        // Define size of multi-array
        coefficientMultiarray.resize( boost::extents[ independentVariableSize.at( 0 ) ]
                [ independentVariableSize.at( 1 ) ][ independentVariableSize.at( 2 ) ] );

        // Parse data
        int currentStartRow = 0;
        for( int k = 0; k < independentVariableSize.at( 2 ); k++ )
        {
            for( int i = 0; i < independentVariableSize.at( 0 ); i++ )
            {
                for( int j = 0; j < independentVariableSize.at( 1 ); j++ )
                {
                    coefficientMultiarray[ i ][ j ][ k ] = coefficientsBlock( i + currentStartRow, j );
                }
            }
            currentStartRow += independentVariableSize.at( 0 );
        }
    }

    else
    {
        throw std::runtime_error( "Error, expected size of 3 dimensions when parsing 1-dimensional data into multi-array" );
    }
    return coefficientMultiarray;
}

//! Function to read a coefficient file (data on a structured grid as a function of N independent variables)
void readCoefficientsFile(
        const std::string fileName,
        std::vector< std::vector< double > >& independentVariables,
        Eigen::MatrixXd& coefficientBlock )
{
    // Open file and create file stream.
    std::fstream stream( fileName.c_str( ), std::ios::in );

    // Check if file opened correctly.
    if ( stream.fail( ) )
    {
        throw std::runtime_error( "Data file could not be opened: " + fileName );
    }

    // Initialize boolean that gets set to true once the file header is passed.
    bool isHeaderPassed = 0;
    bool isFirstLinePassed = 0;


    // Line based parsing
    std::string line;
    std::vector< std::string > vectorOfIndividualStrings;

    // Read file line-by-line
    int numberOfDataLinesParsed = 0;
    int numberOfIndependentVariables = -1;
    while ( !stream.fail( ) && !stream.eof( ) )
    {
        // Get line from stream
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        // Skip empty and comment lines
        if( line.size( ) > 0 && !( line.at( 0 ) == '#' ) )
        {
            // Split string into multiple strings, each containing one element from a line from the data file.
            boost::algorithm::split( vectorOfIndividualStrings,
                                     line,
                                     boost::algorithm::is_any_of( "\t ;, " ),
                                     boost::algorithm::token_compress_on );

            // If this is the first line that is read, it should contain the number of independent variables
            if( !isFirstLinePassed )
            {
                if( vectorOfIndividualStrings.size( ) != 1 )
                {
                    throw std::runtime_error( "Error when reading multi-array, expected number of independent variables" );
                }
                numberOfIndependentVariables = std::stoi( vectorOfIndividualStrings.at( 0 ) );
                isFirstLinePassed = true;
            }            
            // If the file header is not passed, this should contain the independent variables
            else if( !isHeaderPassed )
            {
                std::vector< double > currentDataPoints;
                for( unsigned int i = 0; i < vectorOfIndividualStrings.size( ); i++ )
                {
                    currentDataPoints.push_back( std::stod( vectorOfIndividualStrings.at( i ) ) );
                }
                independentVariables.push_back( currentDataPoints );
            }
            else if( isHeaderPassed )
            {
                // Check line consistency
                if( vectorOfIndividualStrings.size( ) != static_cast< unsigned int >( coefficientBlock.cols( ) ) )
                {
                    throw std::runtime_error(
                                "Error on data line " + std::to_string( numberOfDataLinesParsed ) +
                                " found " + std::to_string( vectorOfIndividualStrings.size( ) ) +
                                " columns, but expected " + std::to_string( coefficientBlock.cols( ) ) );
                }
                else if( numberOfDataLinesParsed > coefficientBlock.rows( ) )
                {
                    throw std::runtime_error(
                                "Error on data line " + std::to_string( numberOfDataLinesParsed ) +
                                " expected " +  std::to_string( coefficientBlock.rows( ) ) + "rows" );
                }
                else
                {
                    // Parse data from current line into output matrix.
                    for( unsigned int i = 0; i < vectorOfIndividualStrings.size( ); i++ )
                    {
                        coefficientBlock( numberOfDataLinesParsed, i ) =
                                std::stod( vectorOfIndividualStrings.at( i ) );
                    }
                    numberOfDataLinesParsed++;
                }
            }

            // If the number of independent variables read is the same as the number of independent variables that
            // should be defined, allocate memory for output matrix.
            if( ( static_cast< int >( independentVariables.size( ) ) == numberOfIndependentVariables ) && !isHeaderPassed )
            {
                // Check input consistency
                if( independentVariables.size( ) == 0 )
                {
                    throw std::runtime_error( "Error when reading multi-array, no header found" );
                }
                else
                {
                    // Define size of output matrix, and allocate memory
                    int numberOfRows = independentVariables.at( 0 ).size( );
                    int numberOfColumns = 1;
                    if( independentVariables.size( ) > 1 )
                    {
                        numberOfColumns = independentVariables.at( 1 ).size( );
                        for( unsigned int i = 2; i < independentVariables.size( ); i++ )
                        {
                            numberOfRows *= independentVariables.at( i ).size( );
                        }
                    }
                    coefficientBlock.setZero( numberOfRows, numberOfColumns );
                }
                isHeaderPassed = true;
            }

        }
    }

    if( numberOfDataLinesParsed != coefficientBlock.rows( ) )
    {
        throw std::runtime_error(
                    "Error at end of coefficient file reader, found " +
                    std::to_string( numberOfDataLinesParsed ) +
                    " lines, but expected " +  std::to_string( coefficientBlock.rows( ) ) + "rows" );
    }
}

//! Function to retrieve the number of independent variables that the coefficients in a file are given for.
int getNumberOfIndependentVariablesInCoefficientFile( const std::string& fileName )
{
    // Open file and create file stream.
    std::fstream stream( fileName.c_str( ), std::ios::in );

    // Check if file opened correctly.
    if ( stream.fail( ) )
    {
       throw std::runtime_error( "Data file could not be opened: " + fileName );
    }

    std::string line;
    std::vector< std::string > vectorOfIndividualStrings;

    // Retrieve number of independent variables from file.
    int numberOfIndependentVariables = -1;
    while ( !stream.fail( ) && !stream.eof( ) && ( numberOfIndependentVariables < 0 ) )
    {
        // Get line from stream
        std::getline( stream, line );
        boost::algorithm::trim( line );

        if( line.size( ) > 0 && !( line.at( 0 ) == '#' ) )
        {
            boost::algorithm::split( vectorOfIndividualStrings,
                                     line,
                                     boost::algorithm::is_any_of( "\t ;, " ),
                                     boost::algorithm::token_compress_on );
            try
            {
                numberOfIndependentVariables = std::stod( vectorOfIndividualStrings.at( 0 ) );
            }
            catch( std::runtime_error )
            {
                throw std::runtime_error( "Error when reading coefficicent file size, input is inconsistent" );
            }
        }
    }
    return numberOfIndependentVariables;
}

}

}
