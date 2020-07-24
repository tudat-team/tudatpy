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
 *      Leloux, J. Filtering Techniques for Orbital Debris Conjunction Analysis
 *          - applied to SSN TLE catalog data and including astrodynamics and
 *          collision probability theory, MSc Literature Research, Delft
 *          University of Technology, 2010.
 *      Celestrak (a). Space Track TLE Retriever Help,
 *          http://celestrak.com/SpaceTrack/TLERetrieverHelp.asp, 2011. Last
 *          accessed: 5 August, 2011.
 *      Space Track. TLE Format, http://www.space-track.org/tle_format.html,
 *          2004. Last accessed: 5 August, 2011.
 *      Celestrak (b). FAQs: Two-Line Element Set Format,
 *          http://celestrak.com/columns/v04n03/, 2006. Last accessed:
 *          5 August, 2011.
 *      Celestrak (c). NORAD Two-Line Element Set Format,
 *          http://celestrak.com/NORAD/documentation/tle-fmt.asp, 2004. Last
 *          accessed: 5 August, 2011.
 *
 *    Notes
 *      Raw TLE data can be obtained from (Celestrak (a), 2011). Explanations of the TLE data
 *      format can be viewed in (Space Track, 2004), (Celestrak (b), 2006), and
 *      (Celestrak (c), 2004).
 *
 */ 

#include <cmath>
#include <map>
#include <string>
#include <utility>

#include <boost/algorithm/string/trim.hpp>
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/twoLineElementsFileReader.h"

namespace tudat
{
namespace input_output
{

// Using declarations.
using mathematical_constants::PI;

//! Open data file.
void TwoLineElementsFileReader::openFile( )
{
    if ( absoluteDirectoryPath_.empty( ) )
    {
        absoluteFilePath_ = getTudatRootPath( ) + relativeDirectoryPath_ + fileName_;
    }

    else
    {
        absoluteFilePath_ = absoluteDirectoryPath_ + fileName_;
    }

    // Open data file.
    dataFile_.open( absoluteFilePath_.c_str( ) );

    // Check if file could be opened. Throw exception with error message if file could not be
    // opened.
    if ( !dataFile_.is_open( ) )
    {
       throw std::runtime_error( "Data file could not be opened. Check the file path and if no other"
								 "programs are using the file.\nCurrent filepath: " + absoluteFilePath_ );
    }
}

//! Read and store data.
void TwoLineElementsFileReader::readAndStoreData( )
{
	// Initialize the satellite name temporary variable
	std::string satelliteName;

    // Initialize line counter
    unsigned int lineCounter = 1;

    // Placeholder for current line in the file
    std::string currentLine;

    // Whilst the end of the data file has not been reached, continue reading
    // from lines from data file.
    while( !dataFile_.eof( ) )
    {
        // Get next line of data from data file and store in a string
        std::getline( dataFile_, currentLine );

        // Check if line of data is header line
        if ( lineCounter <= numberOfHeaderLines_ )
        {
            // Store header line data.
            headerLines_.push_back( currentLine );
        }

        // Else process non-header data line
        else
        {
        	// Temporary storage for first TLE line
        	std::string firstLine;

            // Check if first character is 1 or 2
            if( currentLine.substr( 0, 2 ) == "1 " )
			{
				firstLine = currentLine;

				// Get second TLE line
				std::string secondLine;
				std::getline( dataFile_, secondLine );
				lineCounter++;

				// Remove trailing newline character
				secondLine.erase( std::remove( secondLine.begin( ), secondLine.end( ), '\r' ) );

				// Concatenate the two lines and save in map
				tlePairs_[ satelliteName ] = firstLine + secondLine;
			}
            // Line contains the satellite name
            else
			{
            	// Trim satellite name (remove whitespaces and end-of-line characters
            	satelliteName = currentLine;
            	unsigned int end = satelliteName.find_last_not_of( " \n\r\t\f\v" );
            	satelliteName = satelliteName.substr( 0, end + 1 );
			}
        }

        // Increment line counter
        lineCounter++;
    }

    numberOfObjects_ = tlePairs_.size( );
}

} // namespace input_output
} // namespace tudat
