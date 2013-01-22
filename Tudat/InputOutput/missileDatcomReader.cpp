/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      110530    F.M. Engelen      File created.
 *      120326    D. Dirkx          Modified code to be consistent with latest Tudat/TudatCore;
 *                                  moved relevant functionality of (text)FileReader to this class.
 *
 *    References
 *      Blake, W.B. Missile Datcom User's Manual - 1997 Fortran 90 Version, AFRL-VA-WP-TR-1998-3009
 *          Air Force Research Laboratory, 1998.
 *
 *    Notes
 *
 */

#include <algorithm>
#include <functional>
#include <sstream>
#include <stdexcept>

#include <boost/format.hpp>
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>

#include "Tudat/InputOutput/missileDatcomReader.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace input_output
{

using std::string;
using std::stringstream;
using std::vector;

//! Default constructor.
MissileDatcomReader::MissileDatcomReader ( const std::string& fileNameAndPath ) : dataFile_( )
{
    readFor004( fileNameAndPath );
}

//! Function to read the for004.dat file and return one long vector.
void MissileDatcomReader::readFor004( const std::string& fileNameAndPath )
{
    // Set path and open test text file, skip lines with "CONTENTS".
    openFile( fileNameAndPath );
    readAndStoreData( "CONTENTS" );

    // Retreive data from the textFileReader
    missileDatcomDataContainer_ = containerOfDataFromFile_;

    // Loop through all the strings stored in the container.
    for ( iteratorContainerOfData_ = missileDatcomDataContainer_.begin( ); iteratorContainerOfData_
          != missileDatcomDataContainer_.end( ); iteratorContainerOfData_++ )
    {
        // Empty the dataVector such that it can be reused.
        dataVector_.clear( );

        // Split the long string with muliple entries into a vector of small strings.
        split( iteratorContainerOfData_->second,' ' ,dataVector_ );

        // Convert Strings to doubles and place in the result vector.
        for ( unsigned int iteratorInt_ = 0; iteratorInt_ < dataVector_.size( ); iteratorInt_++ )
        {
            missileDatcomData_.push_back( stringToDouble( dataVector_[ iteratorInt_ ] ) );
        }
    }

    dataFile_.close( );
}

//! Open data file.
void MissileDatcomReader::openFile( const string& fileNameAndPath )
{
    // Open data file.
    dataFile_.open( fileNameAndPath.c_str( ), std::ios::binary );

    // Check if file could be opened. Throw exception with error message if file could not be
    // opened.
    if ( !dataFile_ )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            boost::str( boost::format( "Data file '%s' could not be opened." )
                                        % fileNameAndPath.c_str( ) ) ) )
                    << boost::errinfo_file_name( fileNameAndPath.c_str( ) )
                    << boost::errinfo_file_open_mode( "std::ios::binary" )
                    << boost::errinfo_api_function( "std::ifstream::open" ) );
    }
}

//! Function to split a string consisting of line of entries in smaller strings.
void MissileDatcomReader::split( const string& dataString, char separator,
                                 vector< string >& dataVector )
{
    string::size_type i = 0;
    string::size_type j = dataString.find( separator );

    // Loop over the whole string and place characters into new string until the separator is
    // found. Then start with a new string at the next entry in the vector. Until the end of the
    // original string is reached.
    while ( j != string::npos )
    {
        // to prevent placements of empty entries if multiple seperators after eachother are placed.
        if ( !dataString.substr( i, j-i ).empty() )
        {
            dataVector.push_back( dataString.substr( i, j-i ) );
        }

        i = ++j;
        j = dataString.find( separator, j );

        if ( j == string::npos )
        {
            dataVector.push_back( dataString.substr(i, dataString.length( ) ) );
        }
    }
}

//! Read and store data.
void MissileDatcomReader::readAndStoreData( const string& skipKeyword )
{
    // Local variable for reading a single line
    std::string stringOfData;

    int lineCounter = 1;

    // Whilst the end of the data file has not been reached, continue reading
    // from lines from data file.
    while( !dataFile_.eof( ) )
    {
        // Get next line of data from data file and store in a string.
        getline( dataFile_, stringOfData );

        // Check if string doesn't start with set starting character, if string
        // is not empty, and if the skip keyword is not in the string.
        if ( ( ( !skipKeyword.empty( ) && stringOfData.find( skipKeyword ) == string::npos )
               || ( skipKeyword.empty( ) ) ) && !stringOfData.empty( ) )
        {
            // Store string in container.
            containerOfDataFromFile_[ lineCounter ] = stringOfData;
        }

        // Increment line counter.
        lineCounter++;
    }
}

//! Function to Convert a string to a double.
double MissileDatcomReader::stringToDouble( std::string const& inputString )
{
    // Create a local input string stream object.
    std::istringstream iStringStream( inputString );

    // Create a local double.
    double value;

    // Convert string to double and trow exception if this is not possible.
    if ( !( iStringStream >> value ) ) throw std::runtime_error( "invalid double" );

    return value;
}

} // namespace input_output
} // namespace tudat
