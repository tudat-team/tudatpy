/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"

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
    Eigen::Vector6d TLEKeplerianElements;

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

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information; prints all converted TLE variables obtained
     * from TLE.
     * \param stream Stream to print class data to.
     * \param twoLineElementData TLE data to print.
     * \return Stream handler.
     */
    friend std::ostream& operator << ( std::ostream& stream,
                                     TwoLineElementData& twoLineElementData );

protected:

private:
};

//! Typedef for shared-pointer to TwoLineElementData object.
typedef boost::shared_ptr< TwoLineElementData > TwoLineElementDataPointer;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_TWO_LINE_ELEMENT_DATA_H
