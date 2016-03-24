/*    Copyright (c) 2010-2016, Delft University of Technology
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
 *      120701    A. Ronse          Creation of code.
 *      160324    R. Hoogendoorn    Update for use in current version of Tudat
 *
 *    References
 *      Data source for validation figures:
 *                        http://celestrak.com/SpaceData/sw19571001.txt
 *
 */

#define BOOST_TEST_MAIN

#include <istream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/Astrodynamics/Radiation/parseSolarActivityData.h"
#include "Tudat/Astrodynamics/Radiation/extractSolarActivityData.h"
#include "Tudat/Astrodynamics/Radiation/solarActivityData.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include <tudat/Basics/testMacros.h>

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_Solar_Activity_Parser_Extractor )

//! Test parsing and extraction process of solar dummy activity file
BOOST_AUTO_TEST_CASE( test_parsing_and_extraction )
{
tudat::input_output::parsed_data_vector_utilities::ParsedDataVectorPtr parsedDataVectorPtr;

tudat::radiation::solar_activity::ParseSolarActivityData solarActivityParser;
tudat::radiation::solar_activity::ExtractSolarActivityData solarActivityExtractor;

// Parse file
// save path of cpp file
std::string cppPath( __FILE__ );

// Strip filename from temporary string and return root-path string.
std::string folder = cppPath.substr( 0, cppPath.find_last_of("/\\")+1);
std::string filePath = folder + "testSolarActivity.txt" ;
//std::string fileName = tudat::input_output::getTudatRootPath( ) +
//                       "Astrodynamics/Radiation/UnitTests/testSolarActivity.txt";

// Open dataFile
std::ifstream dataFile;
//dataFile.open(fileName.c_str(), std::ifstream::in);
dataFile.open(filePath.c_str(), std::ifstream::in);
parsedDataVectorPtr = solarActivityParser.parse( dataFile );
dataFile.close();

// Extract data to object of solarActivityData class
std::vector< boost::shared_ptr<tudat::radiation::solar_activity::SolarActivityData> >
        solarActivityData( 6 );
solarActivityData[0] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 0 ) );
solarActivityData[1] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 3 ) );
solarActivityData[2] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 8 ) );
solarActivityData[3] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 9 ) );
solarActivityData[4] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 12 ) );
solarActivityData[5] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 15 ) );

// Check single parameters
BOOST_CHECK_EQUAL(  solarActivityData[0]->dataType, 1 );
BOOST_CHECK_EQUAL(  solarActivityData[1]->planetaryDailyCharacterFigure, 0.0 );
BOOST_CHECK_EQUAL(  solarActivityData[1]->solarRadioFlux107Adjusted, 103.0 );
BOOST_CHECK_EQUAL(  solarActivityData[2]->month, 5 );
BOOST_CHECK_EQUAL(  solarActivityData[2]->internationalSunspotNumber, 0 );
BOOST_CHECK_EQUAL(  solarActivityData[3]->planetaryRangeIndexSum, 176 );
BOOST_CHECK_EQUAL(  solarActivityData[3]->planetaryDailyCharacterFigureConverted, 0 );
BOOST_CHECK_EQUAL(  solarActivityData[4]->bartelsSolarRotationNumber, 2441 );
BOOST_CHECK_EQUAL(  solarActivityData[5]->last81DaySolarRadioFlux107Observed, 187.8 );

// Check planetaryEquivalentAmplitudeVectors
Eigen::VectorXd correctPlanetaryEquivalentAmplitudeVector1( 8 );
correctPlanetaryEquivalentAmplitudeVector1 << 32, 27, 15, 7, 22, 9, 32, 22;
Eigen::VectorXd correctPlanetaryEquivalentAmplitudeVector2 = Eigen::VectorXd::Zero( 8 );

TUDAT_CHECK_MATRIX_CLOSE(  solarActivityData[0]->planetaryEquivalentAmplitudeVector,
                    correctPlanetaryEquivalentAmplitudeVector1, 0 );
TUDAT_CHECK_MATRIX_CLOSE(  solarActivityData[4]->planetaryEquivalentAmplitudeVector,
                    correctPlanetaryEquivalentAmplitudeVector2 , 0);
TUDAT_CHECK_MATRIX_CLOSE(  solarActivityData[5]->planetaryEquivalentAmplitudeVector,
                    correctPlanetaryEquivalentAmplitudeVector2, 0 );
}

BOOST_AUTO_TEST_SUITE_END( )

}   // unit_tests
}   // tudat
