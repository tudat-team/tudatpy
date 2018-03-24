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
std::ostream& operator << ( std::ostream& stream,
                          TwoLineElementData& twoLineElementData )
{

    stream << "This is a TLE data object." << std::endl;
    stream << "The converted TLE information is stored as: " << std::endl;

    stream << "TLE line 0: " << std::endl;
    for ( unsigned int i = 0; i < twoLineElementData.objectName.size( ); i++ )
    {
        stream << twoLineElementData.objectName.at( i ) << " ";
    }
    stream << std::endl;
    stream << twoLineElementData.objectNameString << std::endl;

    stream << "TLE line 1: " << std::endl;
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
    stream << twoLineElementData.modulo10CheckSumLine1 << std::endl;

    stream << "TLE line 2: " << std::endl;
    stream << twoLineElementData.lineNumberLine2 << " ";
    stream << twoLineElementData.objectIdentificationNumberLine2 << " ";
    stream << twoLineElementData.TLEKeplerianElements(
                  orbital_element_conversions::inclinationIndex ) << " ";
    stream << twoLineElementData.TLEKeplerianElements(
                  orbital_element_conversions::longitudeOfAscendingNodeIndex ) << " ";
    stream << twoLineElementData.TLEKeplerianElements(
                  orbital_element_conversions::eccentricityIndex ) << " ";
    stream << twoLineElementData.TLEKeplerianElements(
                  orbital_element_conversions::argumentOfPeriapsisIndex ) << " ";
    stream << twoLineElementData.meanAnomaly << " ";
    stream << twoLineElementData.meanMotionInRevolutionsPerDay << " ";
    stream << twoLineElementData.revolutionNumber << " ";
    stream << twoLineElementData.totalRevolutionNumber << " ";
    stream << twoLineElementData.modulo10CheckSumLine2 << std::endl;

    stream << "Calculated values: " << std::endl;
    stream << "Semi-Major Axis = " << twoLineElementData.TLEKeplerianElements(
                  orbital_element_conversions::semiMajorAxisIndex )
           << std::endl;
    stream << "Perigee = " << twoLineElementData.perigee << std::endl;
    stream << "Apogee = " << twoLineElementData.apogee << std::endl;

    // Return stream.
    return stream;
}

} // namespace input_output
} // namespace tudat
