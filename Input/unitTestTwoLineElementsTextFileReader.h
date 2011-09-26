/*! \file unitTestTwoLineElementsTextFileReader.h
 *    Header file that defined the unit test for the Two-Line Elements (TLE)
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

#ifndef UNITTESTTWOLINEELEMENTSTEXTFILEREADER_H
#define UNITTESTTWOLINEELEMENTSTEXTFILEREADER_H

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test implementation of TLE text file reader class.
/*!
 * Tests implementation of the TwoLineElementsTextFileReader class.
 * A test input file has been written and is located in the Input folder.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testTwoLineElementsTextFileReader( );

}

#endif // UNITTESTTWOLINEELEMENTSTEXTFILEREADER_H

// End of file.
