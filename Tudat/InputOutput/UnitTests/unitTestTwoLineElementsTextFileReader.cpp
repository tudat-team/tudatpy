/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110304    J. Leloux         First setup of unit test.
 *      110805    J. Leloux         Unit test made Tudat-ready for codecheck.
 *      110807    K. Kumar          Minor layout and comment corrections.
 *      110810    J. Leloux         Minor comment corrections.
 *      110826    J. Leloux         Updated test for 2-line and 3-line cases.
 *      111027    K. Kumar          Removed dynamic memory allocation.
 *      120627    A. Ronse          Boostified unit test.
 *
 *    References
 *      Leloux, J. Filtering Techniques for Orbital Debris Conjunction Analysis - applied to SSN
 *          TLE catalog data and including astrodynamics and collision probability theory, MSc
 *          Literature Research, Delft University of Technology, 2010.
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
 *    Raw TLE data can be obtained from (Celestrak (a), 2011). Explanations of the TLE data
 *    format can be viewed in (Space Track, 2004), (Celestrak (b), 2006), and
 *    (Celestrak (c), 2004).
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <map>
#include <string>
#include <vector>

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/InputOutput/twoLineElementsTextFileReader.h"

namespace tudat
{
namespace unit_tests
{

//! Test implementation of TLE text file reader class.
BOOST_AUTO_TEST_SUITE( test_Two_Line_Elements_Text_File_Reader )

//! Test two-line TLE catalog.
BOOST_AUTO_TEST_CASE( testTwoLineElementsTextFileReaderForTwoLine )
{
    // Declare multimap of corrupted TLE data errors.
    std::multimap< int, std::string > corruptedTwoLineElementDataErrors_;

    // Initialization of two TLE text file reader, one two-line and one three-line.
    using tudat::input_output::TwoLineElementsTextFileReader;
    TwoLineElementsTextFileReader twoLineElementsTextFileReaderForTwoLine;
    twoLineElementsTextFileReaderForTwoLine.setLineNumberTypeForTwoLineElementInputData(
                TwoLineElementsTextFileReader::twoLineType );

    // Set input directory, file name, then open file and store data strings
    // with Textfilereader class, and finally close file.
    twoLineElementsTextFileReaderForTwoLine.setRelativeDirectoryPath( "InputOutput/UnitTests/" );

    twoLineElementsTextFileReaderForTwoLine.setFileName(
                "testTwoLineElementsTextFile2Line.txt" );

    // Set input directory, file name, then open file and store data strings
    // with Textfilereader class, and finally close file.
    twoLineElementsTextFileReaderForTwoLine.openFile( );
    twoLineElementsTextFileReaderForTwoLine.readAndStoreData( );
    twoLineElementsTextFileReaderForTwoLine.closeFile( );

    // Store TLE data variables from input file using TLE data class container.
    twoLineElementsTextFileReaderForTwoLine.setCurrentYear( 2011 );
    twoLineElementsTextFileReaderForTwoLine.storeTwoLineElementData( );

    // Get TLE data object.
    std::vector< tudat::input_output::TwoLineElementData > twoLineElementData =
            twoLineElementsTextFileReaderForTwoLine.getTwoLineElementData( );

    // Some random checks of TLE data variables of (corrupt) object 9 to see if
    // the data is stored correctly.
    BOOST_CHECK_EQUAL( twoLineElementData.at( 8 ).fourDigitEpochYear, 2011 );
    BOOST_CHECK_EQUAL( twoLineElementData.at( 8 ).exponentOfBStar, -2 );
    BOOST_CHECK_EQUAL( twoLineElementData.at( 8 ).tleNumber, 26 );
    BOOST_CHECK_EQUAL( twoLineElementData.at( 8 ).objectIdentificationNumberLine2, 30303 );
    BOOST_CHECK_EQUAL( twoLineElementData.at( 8 ).TLEKeplerianElements(
                           basic_astrodynamics::inclinationIndex ), 24.6237 );

    // Stored TLE data is checked for integrity and number of corrupted TLEs is saved,
    // while the corrupted TLEs have been erased from the data
    // and the reason and the corrupt objects involved are written to output.
    // First object is stored as #0!
    corruptedTwoLineElementDataErrors_
            = twoLineElementsTextFileReaderForTwoLine.checkTwoLineElementsFileIntegrity( );

    // Get updated TLE data object.
    std::vector< tudat::input_output::TwoLineElementData > twoLineElementDataAfterIntegrityCheck
            = twoLineElementsTextFileReaderForTwoLine.getTwoLineElementData( );

    // Test input file has been setup to contain 7 corrupted TLEs,
    // if this is not detected the test fails.
    BOOST_CHECK_EQUAL( twoLineElementsTextFileReaderForTwoLine.getNumberOfObjects( ), 3 );

    // Some random checks of TLE data variables to see if the corrupted TLEs were in fact
    // deleted and the valid one are still stored correctly.
    BOOST_CHECK_EQUAL( twoLineElementDataAfterIntegrityCheck.at( 0 ).launchNumber, 2 );
    BOOST_CHECK_EQUAL( twoLineElementDataAfterIntegrityCheck.at( 0 ).
                       firstDerivativeOfMeanMotionDividedByTwo, 0.00000290 );
    BOOST_CHECK_EQUAL( twoLineElementDataAfterIntegrityCheck.at( 1 ).launchPart, "B  " );
    BOOST_CHECK_EQUAL( twoLineElementDataAfterIntegrityCheck.at( 1 ).TLEKeplerianElements(
                           basic_astrodynamics::eccentricityIndex ), 0.0024687 );
    BOOST_CHECK_EQUAL( twoLineElementDataAfterIntegrityCheck.at( 2 ).revolutionNumber, 57038 );
}

//! Test three-line TLE catalog.
BOOST_AUTO_TEST_CASE( testTwoLineElementsTextFileReaderForThreeLine )
{
    // Declare multimap of corrupted TLE data errors.
    std::multimap< int, std::string > corruptedTwoLineElementDataErrors_;

    // Initialization of two TLE text file reader, one two-line and one three-line.
    using tudat::input_output::TwoLineElementsTextFileReader;
    TwoLineElementsTextFileReader twoLineElementsTextFileReaderForThreeLine;
    twoLineElementsTextFileReaderForThreeLine.setLineNumberTypeForTwoLineElementInputData(
                TwoLineElementsTextFileReader::threeLineType );

    // Set input directory, file name, then open file and store data strings
    // with Textfilereader class, and finally close file.
    twoLineElementsTextFileReaderForThreeLine.setRelativeDirectoryPath(
                "InputOutput/UnitTests/" );

    twoLineElementsTextFileReaderForThreeLine.setFileName("testTwoLineElementsTextFile3Line.txt" );

    // Set input directory, file name, then open file and store data strings
    // with Textfilereader class, and finally close file.
    twoLineElementsTextFileReaderForThreeLine.openFile( );
    twoLineElementsTextFileReaderForThreeLine.readAndStoreData( );
    twoLineElementsTextFileReaderForThreeLine.closeFile( );

    // Store TLE data variables from input file using TLE data class container.
    twoLineElementsTextFileReaderForThreeLine.setCurrentYear( 2011 );
    twoLineElementsTextFileReaderForThreeLine.storeTwoLineElementData( );

    // Get TLE data object.
    std::vector< tudat::input_output::TwoLineElementData > twoLineElementData =
            twoLineElementsTextFileReaderForThreeLine.getTwoLineElementData( );

    // Some random checks of TLE data variables of (corrupt) object 9 to see if
    // the data is stored correctly.
    BOOST_CHECK_EQUAL( twoLineElementData.at( 8 ).fourDigitEpochYear, 2011 );
    BOOST_CHECK_EQUAL( twoLineElementData.at( 8 ).exponentOfBStar, -2 );
    BOOST_CHECK_EQUAL( twoLineElementData.at( 8 ).tleNumber, 26 );
    BOOST_CHECK_EQUAL( twoLineElementData.at( 8 ).objectIdentificationNumberLine2, 30303 );
    BOOST_CHECK_EQUAL( twoLineElementData.at( 8 ).TLEKeplerianElements(
                           basic_astrodynamics::inclinationIndex ), 24.6237 );

    // Stored TLE data is checked for integrity and number of corrupted TLEs is saved,
    // while the corrupted TLEs have been erased from the data
    // and the reason and the corrupt objects involved are written to output.
    // First object is stored as #0!
    corruptedTwoLineElementDataErrors_
            = twoLineElementsTextFileReaderForThreeLine.checkTwoLineElementsFileIntegrity( );

    // Get updated TLE data object.
    std::vector< tudat::input_output::TwoLineElementData > twoLineElementDataAfterIntegrityCheck =
            twoLineElementsTextFileReaderForThreeLine.getTwoLineElementData( );

    // Test input file has been setup to contain 7 corrupted TLEs,
    // if this is not detected the test fails.
    BOOST_CHECK_EQUAL( twoLineElementsTextFileReaderForThreeLine.getNumberOfObjects( ), 3 );

    // Some random checks of TLE data variables to see if the corrupted TLEs were in fact
    // deleted and the valid one are still stored correctly.
    BOOST_CHECK_EQUAL( twoLineElementDataAfterIntegrityCheck.at( 0 ).launchNumber, 2 );
    BOOST_CHECK_EQUAL( twoLineElementDataAfterIntegrityCheck.at( 0 ).
                       firstDerivativeOfMeanMotionDividedByTwo, 0.00000290 );
    BOOST_CHECK_EQUAL( twoLineElementDataAfterIntegrityCheck.at( 1 ).launchPart, "B  " );
    BOOST_CHECK_EQUAL( twoLineElementDataAfterIntegrityCheck.at( 1 ).TLEKeplerianElements(
                           basic_astrodynamics::eccentricityIndex ), 0.0024687 );
    BOOST_CHECK_EQUAL( twoLineElementDataAfterIntegrityCheck.at( 2 ).revolutionNumber, 57038 );
}

BOOST_AUTO_TEST_SUITE_END( )

}   // namespace unit_tests
}   // namespace tudat
