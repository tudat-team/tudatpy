/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      110722    J. Leloux         Added variables and changed variables into
 *                                  the KeplerianElements class.
 *      110802    K. Kumar          Minor changes.
 *      110802    J. Leloux         Minor correction and checked for compliance
 *                                  with TLEreader in rev 138.
 *      110810    J. Leloux         Corrected doxygen documentation.
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
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/InputOutput/twoLineElementData.h"

namespace tudat
{
namespace input_output
{

//! Default constructor.
TwoLineElementData::TwoLineElementData( )
    : objectNameString( "" ), lineNumberLine1( -0 ),
      objectIdentificationNumber( -0 ), tleClassification( 'x' ), launchYear( -0 ),
      fourDigitlaunchYear( -0 ), launchNumber( -0 ), launchPart( "" ), epochYear( -0 ),
      fourDigitEpochYear( -0 ), epochDay( -0.0 ), firstDerivativeOfMeanMotionDividedByTwo( -0.0 ),
      coefficientOfSecondDerivativeOfMeanMotionDividedBySix( -0.0 ),
      exponentOfSecondDerivativeOfMeanMotionDividedBySix( -0 ),
      secondDerivativeOfMeanMotionDividedBySix( -0.0 ), coefficientOfBStar( -0.0 ),
      exponentOfBStar( -0 ), bStar( -0.0  ), orbitalModel( 1 ), tleNumber( -0 ),
      modulo10CheckSumLine1( 10 ), lineNumberLine2( 0 ), objectIdentificationNumberLine2( -0 ),
      meanAnomaly( -0.0 ), meanMotionInRevolutionsPerDay( -0.0 ), revolutionNumber( -0 ),
      totalRevolutionNumber( -0 ), modulo10CheckSumLine2( 10 ), perigee( -0.0 ), apogee( -0.0 )
{ }

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          TwoLineElementData& twoLineElementData )
{
    // Using declarations.
    using std::endl;

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
    stream << twoLineElementData.TLEKeplerianElements(
                  basic_astrodynamics::inclinationIndex ) << " ";
    stream << twoLineElementData.TLEKeplerianElements(
                  basic_astrodynamics::longitudeOfAscendingNodeIndex ) << " ";
    stream << twoLineElementData.TLEKeplerianElements(
                  basic_astrodynamics::eccentricityIndex ) << " ";
    stream << twoLineElementData.TLEKeplerianElements(
                  basic_astrodynamics::argumentOfPeriapsisIndex ) << " ";
    stream << twoLineElementData.meanAnomaly << " ";
    stream << twoLineElementData.meanMotionInRevolutionsPerDay << " ";
    stream << twoLineElementData.revolutionNumber << " ";
    stream << twoLineElementData.totalRevolutionNumber << " ";
    stream << twoLineElementData.modulo10CheckSumLine2 << endl;

    stream << "Calculated values: " << endl;
    stream << "Semi-Major Axis = " << twoLineElementData.TLEKeplerianElements(
                  basic_astrodynamics::semiMajorAxisIndex )
           << endl;
    stream << "Perigee = " << twoLineElementData.perigee << endl;
    stream << "Apogee = " << twoLineElementData.apogee << endl;

    // Return stream.
    return stream;
}

} // namespace input_output
} // namespace tudat
