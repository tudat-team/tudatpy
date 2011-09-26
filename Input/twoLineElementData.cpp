/*! \file twoLineElementData.cpp
 *    This source file defines a class which stores (Two-Line Element) TLE data
 *    read by the TwoLineElementsTextFileReader.
 *
 *    Path              : /Input/
 *    Version           : 6
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
 *      110304    J. Leloux         First setup of TLE data class.
 *      110415    K. Kumar          Minor corrections.
 *      110722    J. Leloux         Added variables and changed variables into
 *                                  the KeplerianElements class.
 *      110802    K. Kumar          Minor changes.
 *      110802    J. Leloux         Minor correction and checked for compliance
 *                                  with TLEreader in rev 138.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

// Include statements.
#include "twoLineElementData.h"

// Using declarations.
using std::endl;

//! Default constructor.
TwoLineElementData::TwoLineElementData( ) : objectNameString( "" ), lineNumberLine1( -0 ),
    objectIdentificationNumber( -0 ), tleClassification( 'x' ), launchYear( -0 ),
    fourDigitlaunchYear( -0 ), launchNumber( -0 ), launchPart( "" ), epochYear( -0 ),
    fourDigitEpochYear( -0 ), epochDay( -0.0 ), firstDerivativeOfMeanMotionDividedByTwo( -0.0 ),
    coefficientOfSecondDerivativeOfMeanMotionDividedBySix( -0.0 ),
    exponentOfSecondDerivativeOfMeanMotionDividedBySix( -0 ),
    secondDerivativeOfMeanMotionDividedBySix( -0.0 ), coefficientOfBStar( -0.0 ),
    exponentOfBStar( -0 ), bStar(-0.0  ), orbitalModel( 1 ), tleNumber( -0 ),
    modulo10CheckSumLine1( 10 ), lineNumberLine2( 0 ), objectIdentificationNumberLine2( -0 ),
    meanAnomaly( -0.0 ), meanMotionInRevolutionsPerDay( -0.0 ), revolutionNumber( -0 ),
    totalRevolutionNumber( -0 ), modulo10CheckSumLine2( 10 ), perigee( -0.0 ), apogee( -0.0 ) { }

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          TwoLineElementData& twoLineElementData )
{
    stream << "This is a TLE data object." << endl;
    stream << "The converted TLE information is stored as: " << endl;

    stream << "TLE line 0: " << endl;
    for ( unsigned int i = 0; i < twoLineElementData.objectName.size( ); i++ )
    {
        stream << twoLineElementData.objectName.at( i ) << " ";
    }
    stream << endl;
    stream << twoLineElementData.objectNameString << endl;

    stream << "TLE line 1: " << endl;
    stream << twoLineElementData.lineNumberLine1 << " ";
    stream << twoLineElementData.objectIdentificationNumber << " ";
    stream << twoLineElementData.tleClassification << " ";
    stream << twoLineElementData.launchYear << " ";
    stream << twoLineElementData.fourDigitlaunchYear << " ";
    stream << twoLineElementData.launchNumber << " ";
    stream << twoLineElementData.launchPart << " ";
    stream << twoLineElementData.epochYear << " ";
    stream << twoLineElementData.fourDigitEpochYear << " ";
    stream << twoLineElementData.epochDay << " ";
    stream << twoLineElementData.firstDerivativeOfMeanMotionDividedByTwo << " ";
    stream << twoLineElementData.coefficientOfSecondDerivativeOfMeanMotionDividedBySix << " ";
    stream << twoLineElementData.exponentOfSecondDerivativeOfMeanMotionDividedBySix << " ";
    stream << twoLineElementData.secondDerivativeOfMeanMotionDividedBySix << " ";
    stream << twoLineElementData.coefficientOfBStar << " ";
    stream << twoLineElementData.exponentOfBStar << " ";
    stream << twoLineElementData.bStar << " ";
    stream << twoLineElementData.orbitalModel << " ";
    stream << twoLineElementData.tleNumber << " ";
    stream << twoLineElementData.modulo10CheckSumLine1 << endl;

    stream << "TLE line 2: " << endl;
    stream << twoLineElementData.lineNumberLine2 << " ";
    stream << twoLineElementData.objectIdentificationNumberLine2 << " ";
    stream << twoLineElementData.TLEKeplerianElements.getInclination( ) << " ";
    stream << twoLineElementData.TLEKeplerianElements.getLongitudeOfAscendingNode( ) << " ";
    stream << twoLineElementData.TLEKeplerianElements.getEccentricity( ) << " ";
    stream << twoLineElementData.TLEKeplerianElements.getArgumentOfPeriapsis( ) << " ";
    stream << twoLineElementData.meanAnomaly << " ";
    stream << twoLineElementData.meanMotionInRevolutionsPerDay << " ";
    stream << twoLineElementData.revolutionNumber << " ";
    stream << twoLineElementData.totalRevolutionNumber << " ";
    stream << twoLineElementData.modulo10CheckSumLine2 << endl;

    stream << "Calculated values: " << endl;
    stream << "Semi-Major Axis = " << twoLineElementData.TLEKeplerianElements.getSemiMajorAxis( )
           << endl;
    stream << "Perigee = " << twoLineElementData.perigee << endl;
    stream << "Apogee = " << twoLineElementData.apogee << endl;

    // Return stream.
    return stream;
}

// End of file.
