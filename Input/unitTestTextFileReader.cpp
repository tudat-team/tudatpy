/*! \file unitTestTextFileReader.h
 *    This header file contains the definition of a unit test for the text file
 *    reader class.
 *
 *    Path              : /Input/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : J. Leloux
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.Leloux@tudelft.nl, J.Leloux@student.tudelft.nl,
 *                        jleloux@gmail.com
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 24 March, 2011
 *    Last modified     : 21 May, 2011
 *
 *    References
 *
 *    Notes
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
 *      110324    J. Leloux         Creation of unit test.
 *      110330    K. Kumar          Minor corrections; test unsuccessful.
 *      110407    J. Leloux         Changed setup/order of file and tests ->
 *                                  successful.
 *      110408    K. Kumar          Added carriage return to string ( temp
 *                                  solution ).
 *      110521    K. Kumar          Updated to reflect changes to
 *                                  stripEndOfLineCharacters()
 */

// Include statements.
#include "unitTestTextFileReader.h"

// Using directives.
using std::cerr;
using std::endl;
using std::cout;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test implementation of textFileReader class.
bool testTextFileReader( )
{
    // Test result initialised to false.
    bool isTextFileReaderErroneous = false;

    // Definition of text file strings for check process.
    string line1 = "% Starting character % or skip line 1.";
    string line2 = "% Starting character % or skip line 2.";
    string line3 = "Test line 1.";
    string line4 = "Test line 2.";
    string line5 = "Test line 3.";
    string line6 = "Test line 4.";
    string line7 = "% Starting character % or skip line 3.";
    string line8 = "Test line 5.";
    string line9 = "Test line 6.";

    // Declare test TextFileReader object.
    TextFileReader testTextFileReader;

    // Set path and open test text file, skip lines with
    // starting character %, read rest of the file, and strip End-Of-Line
    // characters from data.
    testTextFileReader.setRelativePath( "Input/" );
    testTextFileReader.setFileName( "testTextFile.txt" );
    testTextFileReader.openFile( );
    testTextFileReader.skipLinesStartingWithCharacter( "%" );
    testTextFileReader.readAndStoreData( );

    // Create a map variable with the read input file string data.
    map< unsigned int, string > testContainerOfData =
            testTextFileReader.getContainerOfData( );

    // Strip End-Of-Line characters from string data.
    testTextFileReader.stripEndOfLineCharacters( );

    // Create a map variable with the read input file string data.
    testContainerOfData =
            testTextFileReader.getContainerOfData( );

    // Check if the defined lines match the stored string data, if it doesn't,
    // set the test boolean to true and give an error output message.
    if ( testContainerOfData[ 3 ] != line3 ||
         testContainerOfData[ 4 ] != line4 ||
         testContainerOfData[ 5 ] != line5 ||
         testContainerOfData[ 6 ] != line6 ||
         testContainerOfData[ 8 ] != line8 ||
         testContainerOfData[ 9 ] != line9    )
    {
        isTextFileReaderErroneous = true;
        cerr << "Test text file reading with lines with "
                "starting character % and string storage was unsuccessful."
             << endl;
    }

    // Close the file.
    testTextFileReader.closeFile( );

    // Declare second test TextFileReader object.
    TextFileReader testTextFileReader2;

    // Set path and open test text file, skip first two lines,
    // read four lines, skip one line, and read two lines, and strip
    // End-Of-Line characters from data.
    testTextFileReader2.setRelativePath( "Input/" );
    testTextFileReader2.setFileName( "testTextFile.txt" );
    testTextFileReader2.openFile( );
    testTextFileReader2.skipLines( 2 );
    testTextFileReader2.readAndStoreData( 4 );
    testTextFileReader2.skipLines( 1 );
    testTextFileReader2.readAndStoreData( 2 );

    // Set map variable to the read input file string data.
    map< unsigned int, string > testContainerOfData2 =
                testTextFileReader2.getContainerOfData( );

    // Strip End-Of-Line characters from string data.
    testTextFileReader2.stripEndOfLineCharacters( );

    // Set map variable to the read input file string data.
    testContainerOfData2 =
            testTextFileReader2.getContainerOfData( );

    // Check if the defined lines match the stored string data, if it doesn't,
    // set the test boolean to true and give an error output message.
    if ( testContainerOfData2[ 3 ] != line3 ||
         testContainerOfData2[ 4 ] != line4 ||
         testContainerOfData2[ 5 ] != line5 ||
         testContainerOfData2[ 6 ] != line6 ||
         testContainerOfData2[ 8 ] != line8 ||
         testContainerOfData2[ 9 ] != line9    )
    {
        isTextFileReaderErroneous = true;
        cerr << "Test text file reading with skipping of lines and reading "
                "integer number of lines and string storage was unsuccessful."
             << endl;
    }

    // Close the file.
    testTextFileReader2.closeFile( );

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isTextFileReaderErroneous;
}

}

// End of file.
