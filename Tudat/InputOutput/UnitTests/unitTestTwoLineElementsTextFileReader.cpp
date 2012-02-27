/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      110304    J. Leloux         First setup of unit test.
 *      110805    J. Leloux         Unit test made Tudat-ready for codecheck.
 *      110807    K. Kumar          Minor layout and comment corrections.
 *      110810    J. Leloux         Minor comment corrections.
 *      110826    J. Leloux         Updated test for 2-line and 3-line cases.
 *      111027    K. Kumar          Removed dynamic memory allocation.
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
 *          http://celestrak.com/columns/v04n03/, 2006. Last accessed: 5 August, 2011.
 *      Celestrak (c). NORAD Two-Line Element Set Format,
 *          http://celestrak.com/NORAD/documentation/tle-fmt.asp, 2004. Last
 *          accessed: 5 August, 2011.
 *
 */

// Temporary notes (move to class/function doxygen):
// Raw TLE data can be obtained from (Celestrak (a), 2011). Explanations of the TLE data
// format can be viewed in (Space Track, 2004), (Celestrak (b), 2006), and
// (Celestrak (c), 2004).
// 

#include <ctime>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include "Tudat/InputOutput/twoLineElementsTextFileReader.h"

//! Test implementation of TLE text file reader class.
int main( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;
    using std::multimap;
    using std::vector;
    using std::string;
    using namespace tudat;
    using namespace tudat::input_output;


    // Test result initialized to false.
    bool isTwoLineElementsTextFileReaderErroneous = false;

    // Declare multimap of corrupted TLE data errors.
    multimap< int, string > corruptedTwoLineElementDataErrors_;

    // Initialization of two TLE text file reader, one two-line and one three-line.
    TwoLineElementsTextFileReader twoLineElementsTextFileReaderForTwoLine;
    TwoLineElementsTextFileReader twoLineElementsTextFileReaderForThreeLine;
    twoLineElementsTextFileReaderForTwoLine.setLineNumberTypeForTwoLineElementInputData(
                TwoLineElementsTextFileReader::twoLineType );
    twoLineElementsTextFileReaderForThreeLine.setLineNumberTypeForTwoLineElementInputData(
                TwoLineElementsTextFileReader::threeLineType );

    // Create vector to store TLE text file readers and add the file readers.
    vector< TwoLineElementsTextFileReader* > twoLineElementsTextFileReaders;
    twoLineElementsTextFileReaders.push_back( &twoLineElementsTextFileReaderForTwoLine );
    twoLineElementsTextFileReaders.push_back( &twoLineElementsTextFileReaderForThreeLine );

    for ( unsigned int i = 0; i < 2; i++ )
    {
        // Set input directory, file name, then open file and store data strings
        // with Textfilereader class, and finally close file.
        twoLineElementsTextFileReaders[ i ]->setRelativeDirectoryPath( "InputOutput/UnitTests/" );

        if ( i == 0 )
        {
            twoLineElementsTextFileReaders[ i ]->setFileName(
                        "testTwoLineElementsTextFile2Line.txt" );
        }

        else if ( i == 1 )
        {
            twoLineElementsTextFileReaders[ i ]->setFileName(
                        "testTwoLineElementsTextFile3Line.txt" );
        }

        // Set input directory, file name, then open file and store data strings
        // with Textfilereader class, and finally close file.
        twoLineElementsTextFileReaders[ i ]->openFile( );
        twoLineElementsTextFileReaders[ i ]->readAndStoreData( );
        twoLineElementsTextFileReaders[ i ]->closeFile( );

        // Store TLE data variables from input file using TLE data class container.
        twoLineElementsTextFileReaders[ i ]->setCurrentYear( 2011 );
        twoLineElementsTextFileReaders[ i ]->storeTwoLineElementData( );

        // Get TLE data object.
        vector< TwoLineElementData > twoLineElementData =
                twoLineElementsTextFileReaders[ i ]->getTwoLineElementData( );

        // Some random checks of TLE data variables of (corrupt) object 9 to see if
        // the data is stored correctly.
        if ( twoLineElementData.at( 8 ).fourDigitEpochYear != 2011 ||
             twoLineElementData.at( 8 ).exponentOfBStar != -2 ||
             twoLineElementData.at( 8 ).tleNumber != 26 ||
             twoLineElementData.at( 8 ).objectIdentificationNumberLine2 != 30303 ||
             twoLineElementData.at( 8 ).TLEKeplerianElements.getInclination( ) != 24.6237 )
        {
            isTwoLineElementsTextFileReaderErroneous = true;

            cerr << "TLE input data of (corrupt) object 9 vs stored variables mismatch." << endl;
            cerr << "TLE text file reader unit test failed!" << endl;
        }

        // Stored TLE data is checked for integrity and number of corrupted TLEs is saved,
        // while the corrupted TLEs have been erased from the data
        // and the reason and the corrupt objects involved are written to output.
        // First object is stored as #0!
        corruptedTwoLineElementDataErrors_
                = twoLineElementsTextFileReaders[ i ]->checkTwoLineElementsFileIntegrity( );

        // Get updated TLE data object.
        vector< TwoLineElementData > twoLineElementDataAfterIntegrityCheck =
                twoLineElementsTextFileReaders[ i ]->getTwoLineElementData( );


        // Test input file has been setup to contain 7 corrupted TLEs,
        // if this is not detected the test fails.
        unsigned int numberOfObjects = twoLineElementsTextFileReaders[ i ]->getNumberOfObjects( );

        if ( numberOfObjects != 3 )
        {
            isTwoLineElementsTextFileReaderErroneous = true;
            cerr << "Amount of corrupt TLEs or number of valid TLEs does not match the "
                 << "predefined 7." << endl;
        }

        // Some random checks of TLE data variables to see if the corrupted TLEs were in fact
        // deleted and the valid one are still stored correctly.
        if ( twoLineElementDataAfterIntegrityCheck.at( 0 ).launchNumber != 2 ||
             twoLineElementDataAfterIntegrityCheck.at( 0 ).firstDerivativeOfMeanMotionDividedByTwo
             != 0.00000290 ||
             twoLineElementDataAfterIntegrityCheck.at( 1 ).launchPart != "B  " ||
             twoLineElementDataAfterIntegrityCheck.at( 1 ).TLEKeplerianElements.getEccentricity( )
             != 0.0024687 ||
             twoLineElementDataAfterIntegrityCheck.at( 2 ).revolutionNumber != 57038 )
        {
            isTwoLineElementsTextFileReaderErroneous = true;

            cerr << "Valid TLE input data vs stored variables mismatch." << endl;
            cerr << "TLE text file reader (file integrity) unit test failed!" << endl;
        }
    }

    return isTwoLineElementsTextFileReaderErroneous;
}
