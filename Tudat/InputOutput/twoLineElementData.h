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
 *      110304    J. Leloux         First setup of TLE data class.
 *      110415    K. Kumar          Minor corrections.
 *      110722    J. Leloux         Added variables and changed variables into the
 *                                  KeplerianElements class.
 *      110802    K. Kumar          Minor changes.
 *      110802    J. Leloux         Minor correction and checked for compliance with TLEreader in
 *                                  rev 138.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      110826    J. Leloux         Added TLE string vector container.
 *
 *    References
 *      Literature research including information on TLE's:
 *      J. Leloux (2010), Filtering Techniques for Orbital Debris Conjunction
 *                        Analysis - applied to SSN TLE catalog data and
 *                        including astrodynamics and collision probability
 *                        theory.
 *      Program to download clean TLE's:
 *                        http://celestrak.com/SpaceTrack/TLERetrieverHelp.asp.
 *      TLE format explanation:
 *                        http://www.space-track.org/tle_format.html
 *                        http://celestrak.com/columns/v04n03/
 *                        http://celestrak.com/NORAD/documentation/tle-fmt.asp.
 *
 */

#ifndef TUDAT_TWO_LINE_ELEMENT_DATA_H
#define TUDAT_TWO_LINE_ELEMENT_DATA_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "Tudat/Astrodynamics/States/keplerianElements.h"

namespace tudat
{
namespace input_output
{

//! TLE data class.
/*!
 * Class containing all variables of one TLE for one space debris
 * object, according to the above format definition.
 * See reference for explanation of the variables.
 */
struct TwoLineElementData
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    TwoLineElementData( );

    //! TLE strings vector.
    /*!
     * Vector containing all strings of the TLE.
     */
    std::vector< std::string > twoLineElementStrings;

    // Line 0 strings.

    //! Object name string vector.
    /*!
     * Vector containing the separate words of the name of the object as strings,
     * range could go as far as 10 strings containing parts of the name.
     */
    std::vector< std::string > objectName;

    //! Object name string.
    /*!
     * String containing the name of the object,
     */
    std::string objectNameString;

    // Line 1 variables.

    //! Line number of line 1.
    /*!
     * Line number of line 1, always '1'.
     */
    unsigned int lineNumberLine1;

    //! Object identfication number.
    /*!
     * The object identfication number, range 0 to 99999.
     */
    unsigned int objectIdentificationNumber;

    //! TLE classification.
    /*!
     * TLE classification, U is unclassified.
     */
    char tleClassification;

    //! Launch year two digits.
    /*!
     * Launch year, two digits, 57-99 (1900s) and 00-56 (2000s).
     */
    unsigned int launchYear;

    //! Launch year four digits.
    /*!
     * Launch year, four digits, converted from launchYear,
     * range 1957 to present year.
     */
    unsigned int fourDigitlaunchYear;

    //! Launch number.
    /*!
     * Launch number of above launch year, range 0 to 999.
     */
    unsigned int launchNumber;

    //! Part of the launch.
    /*!
     * Part/piece/fragment of the launch, range AAA - ZZZ.
     */
    std::string launchPart;

    //! TLE epoch year two digits.
    /*!
     * TLE epoch year, two digits, 57-99 (1900s) and 00-56 (2000s).
     */
    unsigned int epochYear;

    //! TLE epoch year four digits.
    /*!
     * TLE epoch year, four digits, converted from epochYear, range 1957 to
     * present year.
     */
    unsigned int fourDigitEpochYear;

    //! Epoch day of the year.
    /*!
     * Epoch day of the year as a double, range 0.0 to about 367.0.
     */
    double epochDay;

    //! First derivative of the mean motion divided by two.
    /*!
     * First derivative of the mean motion divided by two, in units of
     * revolutions per day squared.
     */
    double firstDerivativeOfMeanMotionDividedByTwo;

    //! Coefficient of second derivative of the mean motion divided by six.
    /*!
     * Coefficient of scientific notation of the second derivative of the mean
     * motion divided by six.
     */
    double coefficientOfSecondDerivativeOfMeanMotionDividedBySix;

    //! Exponent of second derivative of the mean motion divided by six.
    /*!
     * Exponent of scientific notation of the second derivative of the mean
     * motion divided by six.
     */
    double exponentOfSecondDerivativeOfMeanMotionDividedBySix;

    //! Second derivative of the mean motion divided by six.
    /*!
     * The second derivative of the mean motion divided by six, in units of
     * revolutions per day cubed.
     */
    double secondDerivativeOfMeanMotionDividedBySix;

    //! Coefficient of B* (bStar) drag term.
    /*!
     * Coefficient of scientific notation of B* (bStar) drag term.
     */
    double coefficientOfBStar;

    //! Exponent of B* (bStar) drag term.
    /*!
     * Exponent of scientific notation of B* (bStar) drag term.
     */
    int exponentOfBStar;

    //! B* (bStar) drag term.
    /*!
     * B* (bStar) drag term, in units of inverse Earth radius.
     */
    double bStar;

    //! Orbital model.
    /*!
     * Orbital model, always '0' nowadays (SGP4/SDP4).
     */
    unsigned int orbitalModel;

    //! TLE number.
    /*!
     * TLE number, range 0 to 9999, loops back to 0 when 9999 is reached.
     */
    unsigned int tleNumber;

    //! Modulo-10 checksum of line 1.
    /*!
     * Modulo-10 checksum of TLE line 1, all integers except for the
     * modulo-10 checksum integer are added, minus signs count as one,
     * rest is ignored.
     */
    unsigned int modulo10CheckSumLine1;

    // Line 2 variables.

    //! Line number of line 2.
    /*!
     * Line number of line 2, always '2'.
     */
    unsigned int lineNumberLine2;

    //! Object identification number of line 2.
    /*!
     * Object identification number of line 2, same as the one of line 1.
     */
    unsigned int objectIdentificationNumberLine2;

    //! TLE Keplerian Elements.
    /*!
     * Including:
     * Inclination of the object, ranging between 0 and 180 degrees.
     * Right ascension of ascending node, ranging between 0 and 360 degrees.
     * Eccentricity, ranging between 0 and 180 degrees,
     * Argument of Perigee, ranging between 0 and 360 degrees.
     * Semi-Major Axis, calculated from the other TLE variables.
     */
    astrodynamics::states::KeplerianElements TLEKeplerianElements;

    //! Mean anomaly.
    /*!
     *  Mean anomaly, range between 0 and 360 degrees
     */
    double meanAnomaly;

    //! Mean motion.
    /*!
     *  Mean motion, in revolutions per day.
     */
    double meanMotionInRevolutionsPerDay;

    //! Revolution number.
    /*!
     * Revolution number of the object, range 0 to 99999,
     * loops from 99999 back to 0.
     */
    unsigned int revolutionNumber;

    //! Total revolution number.
    /*!
     * Total revolution number calculated with mean motion.
     */
    int totalRevolutionNumber;

    //! Modulo-10 checksum of line 2.
    /*!
     * Modulo-10 checksum of TLE line 2, all integers except for the
     * modulo-10 checksum integer are added, rest is ignored.
     */
    unsigned int modulo10CheckSumLine2;

    //! Perigee of object.
    /*!
     * Perigee of the object is calculated from the other TLE variables.
     */
    double perigee;

    //! Apogee of object.
    /*!
     * Apogee of the object is calculated from the other TLE variables.
     */
    double apogee;

    //! Original line numbers of the TLE.
    /*!
     * The line numbers of this TLE in the original input file are saved here.
     */
    std::vector< unsigned int > lineNumbers;

    //! Overloaded ostream to print class information.
    /*!
     *  Overloaded ostream to print class information; prints all
     *  converted TLE variables obtained from TLE.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     TwoLineElementData& twoLineElementData );

protected:

private:
};

} // namespace input_output
} // namespace tudat

#endif // TUDAT_TWO_LINE_ELEMENT_DATA_H
