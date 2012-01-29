/*! \file matrixTextFileReader.cpp
 *    This file contains the source file of the matrix text file reader included in Tudat.
 *
 *    Path              : /Input/
 *    Version           : 2
 *    Check status      : Unchecked
 *
 *    Author            : F. M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : S. Billemont
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : simon@angelcorp.be
 *
 *    Date created      : 30 May, 2011
 *    Last modified     : 11 January, 2011
 *
 *    References
 *
 *    Notes
 *      If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 */

// Include statements.
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <map>
#include <vector>
#include <stdexcept>
#include "Input/matrixTextFileReader.h"
#include "Input/textFileReader.h"

//! Tudat library namespace.
namespace tudat
{

// Using statements.
using std::vector;
using std::map;
using std::iterator;

//! Read the file and return the data matrix.
Eigen::MatrixXd MatrixTextFileReader::getMatrixFromFile( const string& filename,
                                                         const string& seperators,
                                                         const string& skipLinesCharacter,
                                                         const string& relativePath )
{
    TextFileReader myTextFileReader_;

    myTextFileReader_.setRelativeDirectoryPath( relativePath );
    myTextFileReader_.setFileName( filename );
    myTextFileReader_.openFile( );
    myTextFileReader_.skipLinesStartingWithCharacter( skipLinesCharacter );
    myTextFileReader_.readAndStoreData( );
    myTextFileReader_.stripEndOfLineCharacters( myTextFileReader_.getContainerOfData( ) );
    myTextFileReader_.closeFile( );

    // initialize the map and its itterator and store the file data in it.
    map< unsigned int, string > fileData = myTextFileReader_.getContainerOfData( );
    map< unsigned int, string >::iterator fileDataIterator;

    // read the first string to determine the amount of columns.
    string line_ = fileData.begin( )->second;
    vector< string > splitStringsVector_;

    boost::algorithm::split( splitStringsVector_,
                             line_,
                             boost::is_any_of( seperators ),
                             boost::algorithm::token_compress_on );

    int numberOfColumns_ = splitStringsVector_.size( );

    // Go to the end of the map to determine the number of rows.
    unsigned int numberOfLines = fileData.size( );

    dataMatrix_ = Eigen::MatrixXd( numberOfLines, numberOfColumns_ );

    int rowIndex = 0;

    // Loop over the whole map and convert to a matrix of doubles.
    for ( fileDataIterator = fileData.begin( ); fileDataIterator != fileData.end( );
          fileDataIterator++ )
    {
        splitStringsVector_.clear( );
        line_ = fileDataIterator->second;

        boost::algorithm::split( splitStringsVector_,
                                 line_,
                                 boost::is_any_of(seperators + " "),
                                 boost::algorithm::token_compress_on );


        for ( int columIndex = 0; columIndex < numberOfColumns_; columIndex++ )
        {
            dataMatrix_( rowIndex, columIndex ) =
                    boost::lexical_cast<double>( splitStringsVector_[ columIndex ] );
        }

        rowIndex++;
    }

    return dataMatrix_;
}

} // Namespace tudat.

// End of file.
