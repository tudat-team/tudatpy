/*! \file unitTestTextFileReader.h
 *    This header file contains the definition of a unit test for the text file
 *    reader class.
 *
 *    Path              : /Input/
 *    Version           : 7
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
 *    Last modified     : 27 June, 2011
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
 *                                  stripEndOfLineCharacters().
 *      110607    F.M. Engelen      Updated with skipKeyword() test.
 *      110627    K. Kumar          Fixed skipKeyword() test.
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
    string line6 = "Test line 4 with KEYWORD in the line.";
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

    // Strip End-Of-Line characters from string data.
    testTextFileReader.stripEndOfLineCharacters( );

    // Create a map variable with the read input file string data.
    map< unsigned int, string > testContainerOfData =
            testTextFileReader.getContainerOfData( );

    // Check if the defined lines match the stored string data, if it doesn't,
    // set the test boolean to true and give an error output message.
    if ( testContainerOfData[ 3 ] != line3 ||
         testContainerOfData[ 4 ] != line4 ||
         testContainerOfData[ 5 ] != line5 ||
         testContainerOfData[ 6 ] != line6 ||
         testContainerOfData[ 8 ] != line8 ||
         testContainerOfData[ 9 ] != line9 )
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

    // Strip End-Of-Line characters from string data.
    testTextFileReader2.stripEndOfLineCharacters( );

    // Set map variable to the read input file string data.
    map< unsigned int, string > testContainerOfData2 =
                testTextFileReader2.getContainerOfData( );

    // Check if the defined lines match the stored string data, if it doesn't,
    // set the test boolean to true and give an error output message.
    if ( testContainerOfData2[ 3 ] != line3 ||
         testContainerOfData2[ 4 ] != line4 ||
         testContainerOfData2[ 5 ] != line5 ||
         testContainerOfData2[ 6 ] != line6 ||
         testContainerOfData2[ 8 ] != line8 ||
         testContainerOfData2[ 9 ] != line9 )
    {
        isTextFileReaderErroneous = true;
        cerr << "Test text file reading with skipping of lines and reading "
                "integer number of lines and string storage was unsuccessful."
             << endl;
    }

    // Close the file.
    testTextFileReader2.closeFile( );

    // Declare third test TextFileReader object.
    TextFileReader testTextFileReader3;

    // Declare a skip keyword string;
    string skipKeyword = "KEYWORD";

    // Set path and open test text file, Skip lines with keyword.
    testTextFileReader3.setRelativePath( "Input/" );
    testTextFileReader3.setFileName( "testTextFile.txt" );
    testTextFileReader3.openFile( );
    testTextFileReader3.skipLinesWithKeyword( skipKeyword );
    testTextFileReader3.readAndStoreData( );

    // Strip End-Of-Line characters from string data.
    testTextFileReader3.stripEndOfLineCharacters( );

    // Set map variable to the read input file string data.
    map< unsigned int, string > testContainerOfData3
            = testTextFileReader3.getContainerOfData( );

    // Check if the defined lines match the stored string data, if it doesn't,
    // set the test boolean to true and give an error output message.
    if ( testContainerOfData3[ 1 ] != line1 ||
         testContainerOfData3[ 2 ] != line2 ||
         testContainerOfData3[ 3 ] != line3 ||
         testContainerOfData3[ 4 ] != line4 ||
         testContainerOfData3[ 5 ] != line5 ||
         testContainerOfData3[ 7 ] != line7 ||
         testContainerOfData3[ 8 ] != line8 )
    {
        isTextFileReaderErroneous = true;
        cerr << "Test text file reading with skipping of lines with a keyword "
                "was unsuccessful."
             << endl;
    }

    // Close the file.
    testTextFileReader3.closeFile( );

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isTextFileReaderErroneous;
}

}

// End of file.
