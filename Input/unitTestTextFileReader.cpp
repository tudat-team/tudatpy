/*! \file unitTestTextFileReader.cpp
 *    This source file contains the definition of a unit test for the text file
 *    reader class.
 *
 *    Path              : /Input/
 *    Version           : 8
 *    Check status      : Checked
 *
 *    Author            : J. Leloux
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.Leloux@tudelft.nl, J.Leloux@student.tudelft.nl,
 *                        jleloux@gmail.com
 *
 *    Author/Checker    : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : S. Billemont
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : simon@angelcorp.be
 *
 *    Date created      : 24 March, 2011
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
 *      110324    J. Leloux         Creation of unit test.
 *      110330    K. Kumar          Minor corrections; test unsuccessful.
 *      110407    J. Leloux         Changed setup/order of file and tests -> successful.
 *      110408    K. Kumar          Added carriage return to string ( temp solution ).
 *      110521    K. Kumar          Updated to reflect changes to stripEndOfLineCharacters().
 *      110607    F.M. Engelen      Updated with skipKeyword() test.
 *      110627    K. Kumar          Fixed skipKeyword() test.
 *      111117    K. Kumar          Added invalid input file exception test based on suggestion by
 *                                  S. Billemont.
 */

// Include statements.
#include <boost/exception/all.hpp>
#include <cmath>
#include "Input/textFileReader.h"

//! Test implementation of textFileReader class.
int main( )
{
    // Using declarations.
    using std::string;
    using std::map;
    using std::cerr;
    using std::endl;
    using namespace tudat;

    // Summary of tests.
    // Test 1: Test file reader, skipping lines starting with %.
    // Test 2: Test file reader, skipping set number of lines and reading and storing data.
    // Test 3: Test file reader, skipping words containing KEYWORD.
    // Test 4: Test file reader, setting number of lines of header data.
    // Test 5: Test file reader, setting an invalid input file.

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

    // Test 1: Test file reader, skipping lines starting with %.

    // Declare test TextFileReader object.
    TextFileReader testTextFileReader;

    // Set path and open test text file, skip lines with starting character %, read rest of the
    // file, and strip End-Of-Line characters from data.
    testTextFileReader.setRelativeDirectoryPath( "Input/" );
    testTextFileReader.setFileName( "testTextFile.txt" );
    testTextFileReader.openFile( );
    testTextFileReader.skipLinesStartingWithCharacter( "%" );
    testTextFileReader.readAndStoreData( );

    // Create a map variable with the read input file string data.
    map< unsigned int, string > testContainerOfData = testTextFileReader.getContainerOfData( );

    // Strip End-Of-Line characters from string data.
    testTextFileReader.stripEndOfLineCharacters( testContainerOfData );

    // Check if the defined lines match the stored string data, if it doesn't,
    // set the test boolean to true and give an error output message.
    if ( testContainerOfData.at( 3 ) != line3 ||
         testContainerOfData.at( 4 ) != line4 ||
         testContainerOfData.at( 5 ) != line5 ||
         testContainerOfData.at( 6 ) != line6 ||
         testContainerOfData.at( 8 ) != line8 ||
         testContainerOfData.at( 9 ) != line9 )
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

    // Test 2: Test file reader, skipping set number of lines and reading and storing data.

    // Set path and open test text file, skip first two lines, read four lines, skip one line, and
    // read two lines, and strip End-Of-Line characters from data.
    testTextFileReader2.setRelativeDirectoryPath( "Input/" );
    testTextFileReader2.setFileName( "testTextFile.txt" );
    testTextFileReader2.openFile( );
    testTextFileReader2.skipLines( 2 );
    testTextFileReader2.readAndStoreData( 4 );
    testTextFileReader2.skipLines( 1 );
    testTextFileReader2.readAndStoreData( 2 );

    // Set map variable to the read input file string data.
    map< unsigned int, string > testContainerOfData2 = testTextFileReader2.getContainerOfData( );

    // Strip End-Of-Line characters from string data.
    testTextFileReader2.stripEndOfLineCharacters( testContainerOfData2 );

    // Check if the defined lines match the stored string data, if it doesn't,
    // set the test boolean to true and give an error output message.
    if ( testContainerOfData2.at( 3 ) != line3 ||
         testContainerOfData2.at( 4 ) != line4 ||
         testContainerOfData2.at( 5 ) != line5 ||
         testContainerOfData2.at( 6 ) != line6 ||
         testContainerOfData2.at( 8 ) != line8 ||
         testContainerOfData2.at( 9 ) != line9 )
    {
        isTextFileReaderErroneous = true;
        cerr << "Test text file reading with skipping of lines and reading "
                "integer number of lines and string storage was unsuccessful."
             << endl;
    }

    // Close the file.
    testTextFileReader2.closeFile( );

    // Test 3: Test file reader, skipping words containing KEYWORD.

    // Declare third test TextFileReader object.
    TextFileReader testTextFileReader3;

    // Declare a skip keyword string;
    string skipKeyword = "KEYWORD";

    // Set path and open test text file, Skip lines with keyword.
    testTextFileReader3.setRelativeDirectoryPath( "Input/" );
    testTextFileReader3.setFileName( "testTextFile.txt" );
    testTextFileReader3.openFile( );
    testTextFileReader3.skipLinesWithKeyword( skipKeyword );
    testTextFileReader3.readAndStoreData( );

    // Set map variable to the read input file string data.
    map< unsigned int, string > testContainerOfData3 = testTextFileReader3.getContainerOfData( );

    // Strip End-Of-Line characters from string data.
    testTextFileReader3.stripEndOfLineCharacters( testContainerOfData3 );

    // Check if the defined lines match the stored string data, if it doesn't,
    // set the test boolean to true and give an error output message.
    if ( testContainerOfData3.at( 1 ) != line1 ||
         testContainerOfData3.at( 2 ) != line2 ||
         testContainerOfData3.at( 3 ) != line3 ||
         testContainerOfData3.at( 4 ) != line4 ||
         testContainerOfData3.at( 5 ) != line5 ||
         testContainerOfData3.at( 7 ) != line7 ||
         testContainerOfData3.at( 8 ) != line8 )
    {
        isTextFileReaderErroneous = true;
        cerr << "Test text file reading with skipping of lines with a keyword was unsuccessful."
             << endl;
    }

    // Test 4: Test file reader, setting number of lines of header data.

    // Declare fourth test TextFileReader object.
    TextFileReader testTextFileReader4;

    // Set path and open test text file, Skip lines with keyword.
    testTextFileReader4.setRelativeDirectoryPath( "Input/" );
    testTextFileReader4.setFileName( "testTextFile.txt" );
    testTextFileReader4.setNumberOfHeaderLines( 2 );
    testTextFileReader4.openFile( );
    testTextFileReader4.readAndStoreData( );

    // Set map variable to the read input file string data.
    map< unsigned int, string > testContainerOfData4 = testTextFileReader4.getContainerOfData( );

    // Strip End-Of-Line characters from string data.
    testTextFileReader4.stripEndOfLineCharacters( testContainerOfData4 );

    // Set map variable to the read input file header string data.
    map< unsigned int, string > testContainerOfHeaderData4
            = testTextFileReader4.getContainerOfHeaderData( );

    // Strip End-Of-Line characters from string header data.
    testTextFileReader4.stripEndOfLineCharacters( testContainerOfHeaderData4 );

    // Check if the defined lines match the stored string data, if it doesn't,
    // set the test boolean to true and give an error output message.
    if ( testContainerOfHeaderData4.at( 1 ) != line1 ||
         testContainerOfHeaderData4.at( 2 ) != line2 ||
         testContainerOfData4.at( 3 ) != line3 ||
         testContainerOfData4.at( 4 ) != line4 ||
         testContainerOfData4.at( 5 ) != line5 ||
         testContainerOfData4.at( 6 ) != line6 ||
         testContainerOfData4.at( 7 ) != line7 ||
         testContainerOfData4.at( 8 ) != line8 ||
         testContainerOfData4.at( 9 ) != line9 )
    {
        isTextFileReaderErroneous = true;
        cerr << "Test text file reading with set number of header lines was unsuccessful." << endl;
    }

    // Close the file.
    testTextFileReader4.closeFile( );

    // Test 5: Test file reader, setting an invalid input file.
    try
    {
        // Declare fourth test TextFileReader object.
        TextFileReader testTextFileReader5;

        // Set relative directory path and invalid file name.
        testTextFileReader5.setRelativeDirectoryPath( "Input/" );
        testTextFileReader5.setFileName( "invalidInputFile.txt" );

        // Open invalid input file.
        testTextFileReader5.openFile( );

        // Set test boolean to true if the test reaches this point.
        isTextFileReaderErroneous = true;

        // Output error statement.
        cerr << "Exception handling for invalid input handling failed." << endl;

        // Close the file.
        testTextFileReader5.closeFile( );
    }

    catch ( boost::exception& invalidInputFileException ) { }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isTextFileReaderErroneous )
    { cerr << "testTextFileReader failed!" << std::endl; }
}

// End of file.
