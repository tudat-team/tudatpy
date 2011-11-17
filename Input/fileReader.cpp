/*! \file fileReader.cpp
 *    This source file contains the definition of a base class for all file
 *    readers included in Tudat.
 *
 *    Path              : /Input/
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Author/Checker    : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author/Checker    : S. Billemont
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : simon@angelcorp.be
 *
 *    Checker           : J. Leloux
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.Leloux@student.tudelft.nl
 *
 *    Date created      : 24 February, 2011
 *    Last modified     : 17 November, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110223    K. Kumar          First creation of code.
 *      110224    K. Kumar          Changed vector container to map container.
 *      110224    J. Leloux         Checked code and fixed typo.
 *      110627    K. Kumar          Moved skipLinesWithKeyword() from TextFileReader.
 *      111115    S. Billemont      Added exception handling for invalid input file.
 *      111117    K. Kumar          Added header-line functionality.
 */

// Include statements.
#include <boost/format.hpp>
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>
#include "Input/fileReader.h"

//! Tudat library namespace.
namespace tudat
{

// Using declarations.
using std::ifstream;
using std::string;
using std::map;
using std::cerr;
using std::endl;

//! Open data file.
void FileReader::openFile( )
{
    if ( absoluteDirectoryPath_.compare( "" ) == 0 )
    {
        absoluteFilePath_ = basic_functions::getRootPath( ) + relativeDirectoryPath_ + fileName_;
    }

    else
    {
        absoluteFilePath_ = absoluteDirectoryPath_ + fileName_;
    }

    // Open data file.
    dataFile_.open( absoluteFilePath_.c_str( ), std::ios::binary );

    // Check if file could be opened. Throw exception with error message if file could not be
    // opened.
    if ( !dataFile_ )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            boost::str( boost::format( "Data file '%s' could not be opened." )
                                 % absoluteFilePath_.c_str( ) ) ) )
            << boost::errinfo_file_name( absoluteFilePath_.c_str( ) )
            << boost::errinfo_file_open_mode( "std::ios::binary" )
            << boost::errinfo_api_function( "std::ifstream::open" ) );
    }
}

//! Skip lines.
void FileReader::skipLines( unsigned int numberOfLines )
{
    // Call getline( ) function for set number of lines to be skipped.
    for ( unsigned int i = 0; i < numberOfLines; i++ )
    {
        // Get next line of data from file and don't do anything.
        std::getline( dataFile_, stringOfData_ );

        // Increment line counter.
        lineCounter_++;
    }
}

}

// End of file.
