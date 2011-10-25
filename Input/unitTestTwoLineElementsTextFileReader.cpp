/*! \file unitTestTwoLineElementsTextFileReader.cpp
 *    Source file that implements the unit test for the Two-Line Elements (TLE)
 *    text file reader. A seperate test input TLE catalog file has been made,
 *    which contains 7 corrupt objects (objects 2 and 4-9), with the other 3
 *    objects being valid.
 *
 *    Path              : /Input/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : J. Leloux
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : j.leloux@tudelft.nl, j.leloux@student.tudelft.nl,
 *                        j.leloux@gmail.com
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 4 March, 2011
 *    Last modified     : 10 August, 2011
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
 *      Raw TLE data can be obtained from (Celestrak (a), 2011). Explanations
 *      of the TLE data format can be viewed in (Space Track, 2004),
 *      (Celestrak (b), 2006), and (Celestrak (c), 2004).
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
 *      110304    J. Leloux         First setup of unit test.
 *      110805    J. Leloux         Unit test made Tudat-ready for codecheck.
 *      110807    K. Kumar          Minor layout and comment corrections.
 *      110810    J. Leloux         Minor comment corrections.
 */

// Include statements.
#include <ctime>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include "Input/twoLineElementsTextFileReader.h"
#include "Input/unitTestTwoLineElementsTextFileReader.h"

// Using declarations.
using std::cerr;
using std::endl;
using std::multimap;
using std::vector;
using std::string;
using basic_functions::outputCurrentRunningTime;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test implementation of TLE text file reader class.
bool testTwoLineElementsTextFileReader( )
{
    // Starting clock initialised.
    clock_t start_clock = clock( );

    // Test result initialized to false.
    bool isTwoLineElementsTextFileReaderErroneous = false;

    // Declare multimap of corrupted TLE data errors.
    multimap< int, string > corruptedTwoLineElementDataErrors_;

    // Initialization of TLE text file reader.
    TwoLineElementsTextFileReader twoLineElementsTextFileReader;

    // Set input directory, file name, then open file and store data strings
    // with Textfilereader class, and finally close file.
    twoLineElementsTextFileReader.setRelativePath( "Input/" );
    twoLineElementsTextFileReader.setFileName( "testTwoLineElementsTextFile.txt" );
    twoLineElementsTextFileReader.openFile( );
    twoLineElementsTextFileReader.readAndStoreData( );
    twoLineElementsTextFileReader.stripEndOfLineCharacters( );
    twoLineElementsTextFileReader.closeFile( );

    // Store TLE data variables from input file using TLE data class container.
    twoLineElementsTextFileReader.setCurrentYear( 2011 );
    twoLineElementsTextFileReader.storeTwoLineElementData( );

    // Get TLE data object.
    vector< TwoLineElementData > twoLineElementData =
            twoLineElementsTextFileReader.getTwoLineElementData( );

    // Some random checks of TLE data variables of (corrupt) object 9 to see if
    // the data is stored correctly.
    if ( twoLineElementData.at( 8 ).objectName.at( 1 ) != "R/B" ||
         twoLineElementData.at( 8 ).fourDigitEpochYear != 2011 ||
         twoLineElementData.at( 8 ).exponentOfBStar != -2 ||
         twoLineElementData.at( 8 ).objectIdentificationNumberLine2 != 30303 ||
         twoLineElementData.at( 8 ).TLEKeplerianElements.getInclination( ) != 24.6237 )
    {
        isTwoLineElementsTextFileReaderErroneous = true;

        cerr << "TLE input data of (corrupt) object 9 vs stored variables mismatch." << endl;
        cerr << "TLE text file reader unit test failed!" << endl;
    }

    // Current running time and status written to vector.
    vector< string > runningTimeAndStatusContainer;
    runningTimeAndStatusContainer = outputCurrentRunningTime(
                start_clock, "Storing of TwoLineElement data is done." );

    // Stored TLE data is checked for integrity and number of corrupted TLEs is saved,
    // while the corrupted TLEs have been erased from the data
    // and the reason and the corrupt objects involved are written to output.
    // First object is stored as #0!
    corruptedTwoLineElementDataErrors_
            = twoLineElementsTextFileReader.checkTwoLineElementsFileIntegrity( );

    // Get updated TLE data object.
    vector< TwoLineElementData > twoLineElementDataAfterIntegrityCheck =
            twoLineElementsTextFileReader.getTwoLineElementData( );

    // Current running time and status written to vector.
    runningTimeAndStatusContainer = outputCurrentRunningTime(
                start_clock, "Checking of TwoLineElement file integrity is done." );

    // Test input file has been setup to contain 7 corrupted TLEs,
    // if this is not detected the test fails.
    unsigned int numberOfObjects = twoLineElementsTextFileReader.getNumberOfObjects( );
    if ( numberOfObjects != 3 )
    {
        isTwoLineElementsTextFileReaderErroneous = true;
        cerr << "Amount of corrupt TLEs or number of valid TLEs does not match the predefined 7."
             << endl;
    }

    // Some random checks of TLE data variables to see if the corrupted TLEs were in fact deleted
    // and the valid one are still stored correctly.
    if ( twoLineElementDataAfterIntegrityCheck.at( 0 ).objectNameString
         != "VANGUARD 1              " ||
         twoLineElementDataAfterIntegrityCheck.at( 0 ).launchNumber != 2 ||
         twoLineElementDataAfterIntegrityCheck.at( 1 ).launchPart != "B  " ||
         twoLineElementDataAfterIntegrityCheck.at( 1 ).TLEKeplerianElements.getEccentricity( )
         != 0.0024687 ||
         twoLineElementDataAfterIntegrityCheck.at( 2 ).revolutionNumber != 57038 )
    {
        isTwoLineElementsTextFileReaderErroneous = true;

        cerr << "Valid TLE input data vs stored variables mismatch." << endl;
        cerr << "TLE text file reader (file integrity) unit test failed!" << endl;
    }

    // Current running time and status written to vector.
    runningTimeAndStatusContainer = outputCurrentRunningTime(
                start_clock, "Unit test for TwoLineElementsTextFileReader is done.");

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isTwoLineElementsTextFileReaderErroneous;
}

}

// End of file.
