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
 *      If tabs are used as spaces, it doesn't work. The separator should also be tabs then.
 *
 */

#include <vector>
#include <fstream>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/InputOutput/streamFilters.h"

namespace tudat
{
namespace input_output
{

//! Read the file and return the data matrix.
Eigen::MatrixXd readMatrixFromFile( const std::string& relativePath, const std::string& separators,
                                    const std::string& skipLinesCharacter )
{
    // Open input and output.
    std::fstream file( relativePath.c_str( ), std::ios::in );
    if ( file.fail( ) )
    {
        throw std::runtime_error( "Data file could not be opened: " + relativePath );
    }

    std::stringstream filteredStream( std::ios::in | std::ios::out );
    {
        // Filter the file stream. This needs to be in its own scope, because filtering_stream::
        // flush( ) does not work if the underlying end point is a stringstream, so the flush has
        // to be forced by letting the filtering_stream go out of scope.
        boost::iostreams::filtering_ostream filterProcessor;
        for ( unsigned int i = 0; i < skipLinesCharacter.size( ); i++ )
        {
            // Remove all comments from the stream.
            filterProcessor.push( input_output::stream_filters::RemoveComment(
                                      skipLinesCharacter[ i ], true ) );
        }

        // Add the output to the filter chain.
        filterProcessor.push( filteredStream );

        // Copy the input to the filter.
        boost::iostreams::copy( file, filterProcessor );
    }

    // Seek stream back to start.
    filteredStream.seekg( 0, std::ios::beg );

    // Read the filtered stream into lines.
    std::vector< std::string > lines_;
    while ( !filteredStream.eof( ) )
    {
        std::string line_;
        getline( filteredStream, line_ );
        if ( !line_.empty( ) )
        {
            boost::trim_all( line_ );
            lines_.push_back( line_ );
        }
    }

    // If there are no lines, return an empty matrix.
    if ( lines_.empty( ) )
    {
        return Eigen::MatrixXd( );
    }

    const std::string realSeparators = std::string( separators ) + " ";

    // Determine the number of columns from.
    std::vector< std::string > lineSplit_;
    boost::algorithm::split( lineSplit_, lines_[ 0 ], boost::is_any_of( realSeparators ),
            boost::algorithm::token_compress_on );
    const unsigned int numberOfColumns = lineSplit_.size( );

    // Initialize the matrix with sizes obtained from the number of lines and the entries in the
    // first line.
    Eigen::MatrixXd dataMatrix_( lines_.size( ), numberOfColumns );
    for ( int rowIndex = 0; rowIndex < dataMatrix_.rows( ); rowIndex++ )
    {
        lineSplit_.clear( );

        // Read current line and split into separate entries.
        boost::algorithm::split( lineSplit_, lines_[ rowIndex ],
                                 boost::is_any_of( realSeparators ),
                                 boost::algorithm::token_compress_on );

        // Check if number of column entries in line matches the number of columns in the matrix.
        // If not, throw a runtime error.
        if ( lineSplit_.size( ) != numberOfColumns )
        {
            throw std::runtime_error(
                        "Number of columns in row " + std::to_string( rowIndex ) + " is " +
                        std::to_string( lineSplit_.size( ) ) + " should be " +
                        std::to_string( numberOfColumns ) );
        }

        // Put single line entries into matrix as doubles.
        for ( int columnIndex = 0; columnIndex < dataMatrix_.cols( ); columnIndex++ )
        {
            boost::trim( lineSplit_.at( columnIndex ) );
            dataMatrix_( rowIndex, columnIndex ) =
                    std::stod( lineSplit_.at( columnIndex ) );
        }
    }

    return dataMatrix_;
}

} // namespace input_output
} // namespace tudat
