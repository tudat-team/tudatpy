/*! \file textFileReader.cpp
 *    This source file contains the definition of a text file reader class.
 *
 *    Path              : /Input/
 *    Version           : 7
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
 *    Last modified     : 10 August, 2011
 *
 *    References
 *      ASCII Table, http://www.asciitable.com/, last accessed: 21st May, 2011.
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
 *      110315    J. Leloux         Checked code.
 *      110316    K. Kumar          Added ostream >> operator overload.
 *      110521    K. Kumar          Modified stripEndOfLineCharacters().
 *      110607    F.M. Engelen      Added skipKeyWord Feature.
 *      110627    K. Kumar          Moved skipLinesWithKeyword() to FileReader.
 *      110810    J. Leloux         Changed if-statement of readAndStoreData function.
 */

// Include statements.
#include <string>
#include "Input/textFileReader.h"

//! Tudat library namespace.
namespace tudat
{

// Using declarations.
using std::string;
using std::endl;

//! Read and store data.
void TextFileReader::readAndStoreData( )
{
    // Whilst the end of the data file has not been reached, continue reading
    // from lines from data file.
    while( !dataFile_.eof( ) )
    {
        // Get next line of data from data file and store in a string.
        getline( dataFile_, stringOfData_ );

        // Check if string doesn't start with set starting character, if string
        // is not empty, and if the skip keyword is not in the string.
        if ( ( ( !startingCharacter_.empty( ) && stringOfData_.substr( 0, 1 )
                 .compare( startingCharacter_ ) != 0 )
               || ( !skipKeyword_.empty( ) && stringOfData_.find( skipKeyword_ )
                    == string::npos )
               || ( startingCharacter_.empty( ) && skipKeyword_.empty( ) ) )
             && !stringOfData_.empty( ) )
        {
            // Store string in container.
            containerOfDataFromFile_[ lineCounter_ ] = stringOfData_;
        }

        // Increment line counter.
        lineCounter_++;
    }
}

//! Read and store data.
void TextFileReader::readAndStoreData( unsigned int numberOfLines )
{
    // Loop over number of lines of data to read and stored from data file.
    for ( unsigned int i = 0; i < numberOfLines; i++ )
    {
        // Get next line of data from data file and store in a string.
        getline( dataFile_, stringOfData_ );

        // Check string is not empty.
        if ( !stringOfData_.empty( ) )
        {
            // Store string in container.
            containerOfDataFromFile_[ lineCounter_ ] = stringOfData_;
        }

        // Increment line counter.
        lineCounter_++;
    }
}

//! Strip End-Of-Line characters.
void TextFileReader::stripEndOfLineCharacters( )
{
    // Declare local variables.
    // Declare string iterator.
    string::iterator iteratorString_;

    // Loop through all the strings stored in the container.
    for ( LineBasedStringDataMap::iterator iteratorContainerOfDataFromFile_
          = containerOfDataFromFile_.begin( );
          iteratorContainerOfDataFromFile_ != containerOfDataFromFile_.end( );
          iteratorContainerOfDataFromFile_++ )
    {
        // Loop through all the characters in the string.
        for ( iteratorString_ = iteratorContainerOfDataFromFile_
                                ->second.begin( );
              iteratorString_ != iteratorContainerOfDataFromFile_
                                 ->second.end( );
              iteratorString_++ )
        {
            // Check if end-of-line characters are present in string.
            // The end-of-line characters are checked for their integer
            // equivalents. See: http://www.asciitable.com/.
            if ( static_cast< int >( *iteratorString_ ) == 10
                 || static_cast< int >( *iteratorString_ ) == 13 )
            {
                // Strip end-of-line character from string.
                iteratorContainerOfDataFromFile_
                        ->second.erase( iteratorString_ );

                // Decrement string iterator since character was erased from
                // string.
                iteratorString_--;
            }
        }
    }
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          TextFileReader& pointerToTextFileReader )
{
    stream << "This is a TextFileReader object." << endl;
    stream << "The input data file name is: " << pointerToTextFileReader.fileName_ << endl;
    stream << "The absolute path to the input data file is: "
           << pointerToTextFileReader.absolutePath_ << endl;

    if ( pointerToTextFileReader.startingCharacter_.empty( ) )
    {
        stream << "The starting character to skip when reading the file "
               << "contents has been defined as: "
               << pointerToTextFileReader.startingCharacter_ << endl;
    }

    stream << "The line counter is set currently to: "
           << pointerToTextFileReader.lineCounter_ << endl;

    // Return stream.
    return stream;
}

}

// End of file.
