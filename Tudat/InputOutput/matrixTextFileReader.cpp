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
 *      120111    F.M. Engelen      Replaced part of the code with boost alternatives.
 *      120206    K. Kumar          Added Boost::trim( ) function to trim output of filter
 *                                  before casting it to doubles.
 *      130109    K. Kumar          Ported from Tudat, added trim_all function-call to eliminate
 *                                  spaces at start and end of line.
 *      130111    K. Kumar          Added runtime error to check for incomplete matrix.
 *
 *    References
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
#include <boost/format.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>

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
        boost::throw_exception(
                    std::runtime_error(
                                    boost::str(
                            boost::format( "Data file '%s' could not be opened." )
                            % relativePath.c_str( ) ) ) );
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
            boost::throw_exception(
                        std::runtime_error(
                            boost::str(
                                boost::format(
                                    "Number of columns in row %1% is %2%; should be %3%." )
                                % rowIndex % lineSplit_.size( ) % numberOfColumns ) ) );
        }

        // Put single line entries into matrix as doubles.
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
