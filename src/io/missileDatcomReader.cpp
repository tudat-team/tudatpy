/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Blake, W.B. Missile Datcom User's Manual - 1997 Fortran 90 Version, AFRL-VA-WP-TR-1998-3009
 *          Air Force Research Laboratory, 1998.
 *
 */

#include <algorithm>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <boost/algorithm/string.hpp>

#include "tudat/io/missileDatcomReader.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace input_output
{

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
void MissileDatcomReader::openFile( const std::string& fileNameAndPath )
{
    // Open data file.
    dataFile_.open( fileNameAndPath.c_str( ), std::ios::binary );

    // Check if file could be opened. Throw exception with error message if file could not be
    // opened.
    if ( !dataFile_ )
    {
        throw std::runtime_error(  "Data file could not be opened." + fileNameAndPath );
    }
}

//! Function to split a string consisting of line of entries in smaller strings.
void MissileDatcomReader::split( const std::string& dataString, char separator,
                                 std::vector< std::string >& dataVector )
{
    std::string::size_type i = 0;
    std::string::size_type j = dataString.find( separator );

    // Loop over the whole string and place characters into new string until the separator is
    // found. Then start with a new string at the next entry in the vector. Until the end of the
    // original string is reached.
    while ( j != std::string::npos )
    {
        // to prevent placements of empty entries if multiple seperators after eachother are placed.
        if ( !dataString.substr( i, j-i ).empty() )
        {
            dataVector.push_back( dataString.substr( i, j-i ) );
        }

        i = ++j;
        j = dataString.find( separator, j );

        if ( j == std::string::npos )
        {
            dataVector.push_back( dataString.substr(i, dataString.length( ) ) );
        }
    }
}

//! Read and store data.
void MissileDatcomReader::readAndStoreData( const std::string& skipKeyword )
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
        if ( ( ( !skipKeyword.empty( ) && stringOfData.find( skipKeyword ) == std::string::npos )
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

    // Remove leading or trailing spaces
    std::string inputStringNoSpaces = inputString;
    boost::algorithm::trim(inputStringNoSpaces);

    // Return NaN if the string represent Not a Number
    if (inputStringNoSpaces == "NaN")
    {
        return TUDAT_NAN;
    }

    // Create a local double.
    double value;

    // Convert string to double and trow exception if this is not possible.
    if ( !( iStringStream >> value ) ) throw std::runtime_error( "invalid double" );

    return value;
}

} // namespace input_output
} // namespace tudat
