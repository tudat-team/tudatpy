/*    Copyright (c) 2010 Delft University of Technology.
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
 *      110530    F.M. Engelen      First creation of code.
 *      120111    F.M. Engelen      Replaced part of the code with boost alternatives.
 *      120206    K. Kumar          Added Boost::trim( ) function to trim output of filter
 *                                  before casting it to doubles.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
// 

#include <boost/format.hpp>
#include <boost/throw_exception.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <fstream>
#include <sstream>
#include <TudatCore/InputOutput/streamFilters.h>
#include "Tudat/InputOutput/matrixTextFileReader.h"

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
        boost::throw_exception( std::runtime_error( boost::str(
                boost::format( "Data file '%s' could not be opened." ) % relativePath.c_str( ) ) ) );
    }

    std::stringstream filteredStream( std::ios::in | std::ios::out );
    {
        // Filter the file stream. This needs to be in its own scope, because filtering_stream::
        // flush( ) does not work if the underlying end point is a stringstream, so the flush has
        // to be forced by letting the filtering_stream go out of scope
        boost::iostreams::filtering_ostream filterProcessor;
        for ( unsigned int i = 0; i < skipLinesCharacter.size( ); i++ )
        {
            // Remove all comments from the stream
            filterProcessor.push( tudat::input_output::stream_filters::RemoveComment(
                                      skipLinesCharacter[i], true ) );
        }
        // Add the output to the filter chain
        filterProcessor.push( filteredStream );

        // Copy the input to the filter
        boost::iostreams::copy( file, filterProcessor );
    }
    // Seek stream back to start
    filteredStream.seekg( 0, std::ios::beg );

    // Read the filtered stream into lines
    std::vector< std::string > lines_;
    while ( !filteredStream.eof( ) )
    {
        std::string line_;
        getline( filteredStream, line_ );
        if ( !line_.empty( ) )
        {
            lines_.push_back( line_ );
        }
    }

    // If there are no lines, return an empty matrix
    if ( lines_.empty( ) )
    {
        return Eigen::MatrixXd( );
    }

    const std::string realSeparators = std::string( separators ) + " ";

    // Determine the number of columns from
    std::vector< std::string > lineSplit_;
    boost::algorithm::split( lineSplit_, lines_[0], boost::is_any_of( realSeparators ),
                             boost::algorithm::token_compress_on );

    // Initialize the matrix with sizes obtained from the number of lines and the entries in the
    // first line.
    Eigen::MatrixXd dataMatrix_( lines_.size( ), lineSplit_.size( ) );
    for ( int rowIndex = 0; rowIndex < dataMatrix_.rows( ); rowIndex++ )
    {
        lineSplit_.clear( );

        boost::algorithm::split( lineSplit_, lines_[rowIndex], boost::is_any_of( realSeparators ),
                                 boost::algorithm::token_compress_on );
        for ( int columnIndex = 0; columnIndex < dataMatrix_.cols( ); columnIndex++ )
        {
            boost::trim( lineSplit_.at( columnIndex ) );
            dataMatrix_( rowIndex, columnIndex ) =
                    boost::lexical_cast< double >( lineSplit_.at( columnIndex ) );
        }
    }

    return dataMatrix_;
}

} // namespace input_output
} // namespace tudat
