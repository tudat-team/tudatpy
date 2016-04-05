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
 *      120420    A. Ronse          First setup of solar activity data class.
 *
 *    References
 *      Data file:
 *                        http://celestrak.com/SpaceData/sw19571001.txt
 *      Data format explanation:
 *                        http://celestrak.com/SpaceData/SpaceWx-format.asp
 *
 */

#include <istream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/InputOutput/solarActivityData.h"
#include "Tudat/InputOutput/solarActivityData.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/InputOutput/parseSolarActivityData.h"
#include "Tudat/InputOutput/extractSolarActivityData.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include <tudat/Basics/testMacros.h>

namespace tudat
{
namespace input_output
{
namespace solar_activity
{

 //! Default constructor.
 SolarActivityData::SolarActivityData( ) : year( 0 ), month( 0 ), day( 0 ),
 bartelsSolarRotationNumber( 0 ), dayOfBartelsCycle( 0 ), planetaryRangeIndexSum( 0 ),
 planetaryEquivalentAmplitudeAverage( 0 ), planetaryDailyCharacterFigure( -0.0 ),
 planetaryDailyCharacterFigureConverted( 0 ), internationalSunspotNumber( 0 ),
 solarRadioFlux107Adjusted( -0.0 ), fluxQualifier( 0 ),
 centered81DaySolarRadioFlux107Adjusted( -0.0 ), last81DaySolarRadioFlux107Adjusted ( -0.0 ),
solarRadioFlux107Observed( -0.0 ), centered81DaySolarRadioFlux107Observed( -0.0 ),
last81DaySolarRadioFlux107Observed ( -0.0 ), planetaryRangeIndexVector( Eigen::VectorXd::Zero(8) ),
planetaryEquivalentAmplitudeVector( Eigen::VectorXd::Zero(8) ), dataType( 0 ) { }

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                              SolarActivityData& solarActivityData )
{
    stream << "This is a Solar Activity data object." << std::endl;
    stream << "The solar activity information is stored as: " << std::endl;

    stream << "Year: " << solarActivityData.year << std::endl;
    stream << "Month: " << solarActivityData.month << std::endl;
    stream << "Day: " << solarActivityData.day << std::endl;
    stream << "BSRN: " << solarActivityData.bartelsSolarRotationNumber << std::endl;
    stream << "ND: " << solarActivityData.dayOfBartelsCycle << std::endl;
    stream << "Kp (0000-0300 UT): " << solarActivityData.planetaryRangeIndexVector(0)<< std::endl;
    stream << "Kp (0300-0600 UT): " << solarActivityData.planetaryRangeIndexVector(1) << std::endl;
    stream << "Kp (0600-0900 UT): " << solarActivityData.planetaryRangeIndexVector(2) << std::endl;
    stream << "Kp (0900-1200 UT): " << solarActivityData.planetaryRangeIndexVector(3) << std::endl;
    stream << "Kp (1200-1500 UT): " << solarActivityData.planetaryRangeIndexVector(4) << std::endl;
    stream << "Kp (1500-1800 UT): " << solarActivityData.planetaryRangeIndexVector(5) << std::endl;
    stream << "Kp (1800-2100 UT): " << solarActivityData.planetaryRangeIndexVector(6) << std::endl;
    stream << "Kp (2100-0000 UT): " << solarActivityData.planetaryRangeIndexVector(7) << std::endl;
    stream << "Kp-sum: " << solarActivityData.planetaryRangeIndexSum << std::endl;
    stream << "Ap (0000-0300 UT): " << solarActivityData.planetaryEquivalentAmplitudeVector(0) << std::endl;
    stream << "Ap (0300-0600 UT): " << solarActivityData.planetaryEquivalentAmplitudeVector(1) << std::endl;
    stream << "Ap (0600-0900 UT): " << solarActivityData.planetaryEquivalentAmplitudeVector(2) << std::endl;
    stream << "Ap (0900-1200 UT): " << solarActivityData.planetaryEquivalentAmplitudeVector(3) << std::endl;
    stream << "Ap (1200-1500 UT): " << solarActivityData.planetaryEquivalentAmplitudeVector(4) << std::endl;
    stream << "Ap (1500-1800 UT): " << solarActivityData.planetaryEquivalentAmplitudeVector(5) << std::endl;
    stream << "Ap (1800-2100 UT): " << solarActivityData.planetaryEquivalentAmplitudeVector(6) << std::endl;
    stream << "Ap (2100-0000 UT): " << solarActivityData.planetaryEquivalentAmplitudeVector(7) << std::endl;
    stream << "Ap-average: " << solarActivityData.planetaryEquivalentAmplitudeAverage << std::endl;
    stream << "Cp: " << solarActivityData.planetaryDailyCharacterFigure << std::endl;
    stream << "C9: " << solarActivityData.planetaryDailyCharacterFigureConverted << std::endl;
    stream << "ISN: "<< solarActivityData.internationalSunspotNumber << std::endl;
    stream << "F10.7 (Adjusted): " << solarActivityData.solarRadioFlux107Adjusted << std::endl;
    stream << "Q: " << solarActivityData.fluxQualifier << std::endl;
    stream << "Ctr81 (Adjusted): " << solarActivityData.centered81DaySolarRadioFlux107Adjusted << std::endl;
    stream << "Lst81 (Adjusted): " << solarActivityData.last81DaySolarRadioFlux107Adjusted << std::endl;
    stream << "F10.7 (Observed): " << solarActivityData.solarRadioFlux107Observed << std::endl;
    stream << "Ctr81 (Observed): " << solarActivityData.centered81DaySolarRadioFlux107Observed << std::endl;
    stream << "Lst81 (Observed): " << solarActivityData.last81DaySolarRadioFlux107Observed << std::endl;
    stream << "Datatype: " << solarActivityData.dataType << std::endl;

    // Return stream.
    return stream;
}

//! This function reads a SpaceWeather data file and returns a map with SolarActivityData
SolarActivityDataMap readSolarActivityData(std::string filePath){
    // Data Vector container
    tudat::input_output::parsed_data_vector_utilities::ParsedDataVectorPtr parsedDataVectorPtr;

    // datafile Parser and Extractor
    tudat::input_output::solar_activity::ParseSolarActivityData solarActivityParser;
    tudat::input_output::solar_activity::ExtractSolarActivityData solarActivityExtractor;


    // Open dataFile and Parse
    std::ifstream dataFile;
    dataFile.open(filePath.c_str(), std::ifstream::in);
    parsedDataVectorPtr = solarActivityParser.parse( dataFile );
    dataFile.close();

    int numberOfLines = parsedDataVectorPtr->size() ;
    SolarActivityDataMap DataMap ;
    double JulianDate ;

    // Save each line to datamap
    for(int i = 0 ; i < numberOfLines ; i++ ){
        JulianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
                    solarActivityExtractor.extract( parsedDataVectorPtr->at( i ) )->year,
                    solarActivityExtractor.extract( parsedDataVectorPtr->at( i ) )->month,
                    solarActivityExtractor.extract( parsedDataVectorPtr->at( i ) )->day,
                    0, 0, 0.0) ;
        DataMap[ JulianDate ] = solarActivityExtractor.extract( parsedDataVectorPtr->at( i ) ) ;
    }

    return DataMap;

}

} // solar_activity
} // input_output
} // tudat
