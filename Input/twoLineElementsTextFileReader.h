/*! \file twoLineElementsTextFileReader.h
 *    This header file defines a class which can:
 *     - read a 3-line Two-Line Element (TLE) catalog file and store its data,
 *     - check if the data file is valid.
 *    Note that 2-line TLE handling still has to be added.
 *
 *    Path              : /Input/
 *    Version           : 7
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
 *    Date created      : 21 February, 2011
 *    Last modified     : 27 October, 2011
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
 *      110221    J. Leloux         Startup of TLE header file and class.
 *      110301    J. Leloux         Adjusting header file to parent classes and Tudat rules.
 *      110803    J. Leloux         First setup for codecheck.
 *      110805    K. Kumar          Layout and comment corrections; added
 *                                  get-function for std::vector container of TLE data.
 *      110810    J. Leloux         Tested new setup and changed descriptions.
 *      110826    J. Leloux         Added functionality for 2-line and 3-line data.
 *      111027    K. Kumar          Modified 2-line and 3-line options using enum.
 */

#ifndef TWOLINEELEMENTSTEXTFILEREADER_H
#define TWOLINEELEMENTSTEXTFILEREADER_H

// Include statements.
#include "Input/textFileReader.h"
#include "Input/twoLineElementData.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! TLE catalog text file reader class.
/*!
 * Definition of TLE catalog text file reader class.
 */
class TwoLineElementsTextFileReader : public TextFileReader
{
public:

    //! Line-number types for TLE input data.
    /*!
     * Line-number types for TLE input data.
     */
    enum LineNumberTypesForTwoLineElementInputData { twoLineType, threeLineType };

    //! Default constructor.
    /*!
     * Default constructor.
     */
    TwoLineElementsTextFileReader( ) : currentYear_( -0 ), numberOfObjects_( -0 ),
        numberOfLinesPerTwoLineElementDatum_( 3 ), twoLineElementData_( ) { }

    //! Set current year.
    /*!
     * Sets the current year.
     * \param currentYear Current year.
     */
    void setCurrentYear( const unsigned int& currentYear ) { currentYear_ = currentYear; }

    //! Get TLE data.
    /*!
     * Returns container of TwoLineElementData objects, containing all the TLE
     * data retrieved from the catalog file and stored in objects.
     * \return TLE data stored in TwoLineElementData objects.
     */
    std::vector< TwoLineElementData >& getTwoLineElementData( ) { return twoLineElementData_; }

    //! Get number of objects.
    /*!
     * Returns number of objects in TLE data catalog file.
     * \return Number of objects in TLE data catalog file.
     */
    unsigned int& getNumberOfObjects( ) { return numberOfObjects_; }

    //! Convert and store TLE data.
    /*!
     * Converts strings read by TextFileReader to the variables contained in
     * their TLE format and stores their data according to their variable-type
     * in the variables contained in this class.
     */
    void storeTwoLineElementData( );

    //! Checks the integrity of the TLE input file.
    /*!
     * Checks the integrity of the TLE input file. For instance:
     * line number 1 = 1 and 2 = 2, object identification number of
     * line 1 = line 2, U for unclassified data, modulo-10 checksum check, etc.
     * Deletes objects whose TLE data are corrupted.
     * \return 0 for success, integer with number of corrupted objects for
     *          failure (with error output).
     */
    std::multimap< int, std::string > checkTwoLineElementsFileIntegrity( );

    //! Set line number type for TLE input data.
    /*!
      * Sets the line number type for TLE input data. This can be either 2-line or 3-line.
      */
    void setLineNumberTypeForTwoLineElementInputData(
        LineNumberTypesForTwoLineElementInputData lineNumberType )
    {
        if ( lineNumberType == twoLineType ) { numberOfLinesPerTwoLineElementDatum_ = 2; }
        else if ( lineNumberType == threeLineType ) { numberOfLinesPerTwoLineElementDatum_ = 3; }
    }

protected:

private:

    //! Current year.
    /*!
     * Current year.
     */
    unsigned int currentYear_;

    //! Number of object in input file.
    /*!
     * Number of objects in input file, linecounter divided by 3.
     */
    unsigned int numberOfObjects_;

    //! Number of lines per TLE datum.
    /*!
     * Number of lines per TLE datum. (Can be 2 or 3).
     */
    unsigned int numberOfLinesPerTwoLineElementDatum_;

    //! Vector of TwoLineElementData objects.
    /*!
     * Vector of TwoLineElementData objects, used to store multiple objects
     * with their TLE data from catalog file.
     */
    std::vector< TwoLineElementData > twoLineElementData_;
};

}

#endif // TWOLINEELEMENTSTEXTFILEREADER_H

// End of file.
