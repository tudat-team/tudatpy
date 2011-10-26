/*! \file fileReader.cpp
 *    This source file contains the definition of a base class for all file
 *    readers included in Tudat.
 *
 *    Path              : /Input/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Leloux
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.Leloux@student.tudelft.nl
 *
 *    Date created      : 24 February, 2011
 *    Last modified     : 27 June, 2011
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
 */

// Include statements.
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
    // Set absolute file path.
    absolutePath_ = basic_functions::getRootPath( ) + relativePath_ + fileName_;

    // Open data file.
    dataFile_.open( absolutePath_.c_str( ), std::ios::binary );

    // Check if file could be opened.
    if ( !dataFile_ )
    {
        cerr << "Error: Data file could not be opened." << endl;
        cerr << "File location is stored as: " <<  absolutePath_.c_str( ) << endl;
    }
}

//! Skip lines.
void FileReader::skipLines( unsigned int numberOfLines )
{
    // Call getline( ) function for set number of lines to be skipped.
    for ( unsigned int i = 0; i < numberOfLines; i++ )
    {
        // Get next line of data from file and don't do anything.
        getline( dataFile_, stringOfData_ );

        // Increment line counter.
        lineCounter_++;
    }
}

}

// End of file.
