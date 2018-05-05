/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPARTA_DATA_READER_H
#define TUDAT_SPARTA_DATA_READER_H

#include <string>
#include <vector>

#include <Eigen/Core>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "Tudat/InputOutput/spartaInputOutput.h"

namespace tudat
{

namespace input_output
{

//! Function to read the geometry file for a SPARTA rarefied flow simulation.
/*!
 *  Function to read the geometry file for a SPARTA rarefied flow simulation.
 *  \param geometryFile File name for the geometry.
 *  \return Pair: first entry containing list of points, second containing list of triangles.
 */
std::pair< Eigen::Matrix< double, Eigen::Dynamic, 3 >, Eigen::Matrix< int, Eigen::Dynamic, 3 > >
readSpartaGeometryFile( const std::string& geometryFile )
{
    // Initialize output vectors
    Eigen::Matrix< double, Eigen::Dynamic, 3 > shapePoints;
    Eigen::Matrix< int, Eigen::Dynamic, 3 > shapeTriangles;
    int numberOfPoints = 0;
    int numberOfTriangles = 0;

    // Open file and create file stream.
    std::fstream stream( geometryFile.c_str( ), std::ios::in );

    // Check if file opened correctly.
    if ( stream.fail( ) )
    {
        throw std::runtime_error( "Data file could not be opened: " + geometryFile );
    }

    // Initialize booleans that specifies once parts of file have been passed.
    bool isNumberOfPointsPassed = false;
    bool isNumberOfTrianglesPassed = false;
    bool isListOfPointsPassed = false;

    // Line based parsing
    std::string line;
    std::vector< std::string > vectorOfIndividualStrings;

    // Read file line-by-line
    int numberOfPointsParsed = 0;
    int numberOfTrianglesParsed = 0;
    while ( !stream.fail( ) && !stream.eof( ) )
    {
        // Get line from stream
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        // Skip empty and comment lines
        if ( line.size( ) > 0 && !( line.at( 0 ) == '#' ) )
        {
            // Split string into multiple strings, each containing one element from a line from the data file.
            boost::algorithm::split( vectorOfIndividualStrings,
                                     line,
                                     boost::algorithm::is_any_of( "\t ;, " ),
                                     boost::algorithm::token_compress_on );

            // If this is the first line that is read, it should contain the number of points
            if ( !isNumberOfPointsPassed )
            {
                if ( vectorOfIndividualStrings.size( ) != 2 )
                {
                    throw std::runtime_error( "Error when reading SPARTA geometry file, expected number of points." );
                }
                numberOfPoints = std::stoi( vectorOfIndividualStrings.at( 0 ) );
                isNumberOfPointsPassed = true;
            }
            // If this is the first line that is read, it should contain the number of triangles
            else if ( !isNumberOfTrianglesPassed )
            {
                if ( vectorOfIndividualStrings.size( ) != 2 )
                {
                    throw std::runtime_error( "Error when reading SPARTA geometry file, expected number of triangles." );
                }
                numberOfTriangles = std::stoi( vectorOfIndividualStrings.at( 0 ) );
                isNumberOfTrianglesPassed = true;
            }
            else if ( !isListOfPointsPassed )
            {
                if ( vectorOfIndividualStrings.at( 0 ) != "Points" )
                {
                    // Check line consistency
                    if ( vectorOfIndividualStrings.size( ) != static_cast< unsigned int >( 4 ) )
                    {
                        throw std::runtime_error(
                                    "Error on data line " + std::to_string( numberOfPointsParsed ) +
                                    " found " + std::to_string( vectorOfIndividualStrings.size( ) ) +
                                    " columns, but expected " + std::to_string( 4 ) );
                    }
                    else
                    {
                        // Parse data from current line into output matrix.
                        for ( unsigned int i = 0; i < ( vectorOfIndividualStrings.size( ) - 1 ); i++ )
                        {
                            shapePoints( numberOfPointsParsed, i ) = std::stod( vectorOfIndividualStrings.at( i + 1 ) );
                        }
                        numberOfPointsParsed++;
                    }
                }

                if ( numberOfPointsParsed == numberOfPoints )
                {
                    isListOfPointsPassed = true;
                }
            }
            else if ( isListOfPointsPassed )
            {
                if ( vectorOfIndividualStrings.at( 0 ) != "Triangles" )
                {
                    // Check line consistency
                    if ( vectorOfIndividualStrings.size( ) != static_cast< unsigned int >( 4 ) )
                    {
                        throw std::runtime_error(
                                    "Error on data line " + std::to_string( numberOfTrianglesParsed ) +
                                    " found " + std::to_string( vectorOfIndividualStrings.size( ) ) +
                                    " columns, but expected " + std::to_string( 4 ) );
                    }
                    else
                    {
                        // Parse data from current line into output matrix.
                        for ( unsigned int i = 0; i < ( vectorOfIndividualStrings.size( ) - 1 ); i++ )
                        {
                            shapeTriangles( numberOfTrianglesParsed, i ) = std::stod( vectorOfIndividualStrings.at( i + 1 ) );
                        }
                        numberOfTrianglesParsed++;
                    }
                }

                if ( numberOfTrianglesParsed > numberOfTriangles )
                {
                    throw std::runtime_error( "Number of triangles in file does not match file header." );
                }
            }

            // Allocate memory for point and triangle matrices.
            if ( isNumberOfPointsPassed && isNumberOfTrianglesPassed )
            {
                // Check input consistency
                if ( ( numberOfPoints == 0 ) || ( numberOfTriangles == 0 ) )
                {
                    throw std::runtime_error( "Error when reading shape file, expected to find a non-zero number of points and triangles." );
                }
                else
                {
                    // Define size of output vectors
                    shapePoints.resize( numberOfPoints, 3 );
                    shapeTriangles.resize( numberOfTriangles, 3 );
                }
            }
        }
    }

    return std::make_pair( shapePoints, shapeTriangles );
}

//! Function to read the input file template format for a SPARTA rarefied flow simulation.
/*!
 *  Function to read the input file template format for a SPARTA rarefied flow simulation.
 *  \param geometryFile File name for the template.
 *  \return String of input template format.
 */
std::string readSpartaInputFileTemplate( )
{
    // Initialize output variable
    std::string inputFile = getSpartaInputFileTemplate( );
    std::string inputTemplate;

    // Open file and create file stream.
    std::fstream stream( inputFile.c_str( ), std::ios::in );

    // Check if file opened correctly.
    if ( stream.fail( ) )
    {
        throw std::runtime_error( "Data file could not be opened: " + inputFile );
    }

    // Line based parsing
    std::string line;

    // Read file line-by-line
    while ( !stream.fail( ) && !stream.eof( ) )
    {
        // Get line from stream
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        // Skip empty and comment lines
        if ( line.size( ) > 0 && !( line.at( 0 ) == '#' ) )
        {
            // Append line to string
            inputTemplate.append( line + "\n" );
        }
    }
    return inputTemplate;
}

} // namespace input_output

} // namespace tudat

#endif // TUDAT_SPARTA_DATA_READER_H
