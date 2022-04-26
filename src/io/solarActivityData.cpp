/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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
#include <stdexcept>

#include "tudat/basics/testMacros.h"
#include "tudat/io/solarActivityData.h"
#include "tudat/io/parsedDataVectorUtilities.h"
#include "tudat/io/parseSolarActivityData.h"
#include "tudat/io/extractSolarActivityData.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/math/basic/mathematicalConstants.h"

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
    last81DaySolarRadioFlux107Observed ( -0.0 ), planetaryRangeIndexVector( Eigen::VectorXd::Zero( 8 ) ),
    planetaryEquivalentAmplitudeVector( Eigen::VectorXd::Zero( 8 ) ), dataType( 0 ) { }

//! Overload ostream to print class information.
std::ostream& operator << ( std::ostream& stream,
                          SolarActivityData& solarActivityData )
{
    stream << "This is a Solar Activity data object." << std::endl;
    stream << "The solar activity information is stored as: " << std::endl;

    stream << "Year: " << solarActivityData.year << std::endl;
    stream << "Month: " << solarActivityData.month << std::endl;
    stream << "Day: " << solarActivityData.day << std::endl;
    stream << "BSRN: " << solarActivityData.bartelsSolarRotationNumber << std::endl;
    stream << "ND: " << solarActivityData.dayOfBartelsCycle << std::endl;
    stream << "Kp (0000-0300 UT): " << solarActivityData.planetaryRangeIndexVector(0) << std::endl;
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
    stream << "ISN: " << solarActivityData.internationalSunspotNumber << std::endl;
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
SolarActivityDataMap readSolarActivityData( std::string filePath )
{
    // Data Vector container
    tudat::input_output::parsed_data_vector_utilities::ParsedDataVectorPtr parsedDataVector;

    // datafile Parser and Extractor
    tudat::input_output::solar_activity::ParseSolarActivityData solarActivityParser;
    tudat::input_output::solar_activity::ExtractSolarActivityData solarActivityExtractor;


    // Open dataFile and Parse
    std::ifstream dataFile;
    dataFile.open( filePath.c_str( ), std::ifstream::in );
    if( !dataFile.is_open( ) )
    {
        throw std::runtime_error( "Error when reading space weather data file. Requested file: <" + filePath + "> cannot be opened" );
    }

    parsedDataVector = solarActivityParser.parse( dataFile );
    dataFile.close( );

    int numberOfLines = parsedDataVector->size( );
    SolarActivityDataMap dataMap;
    double julianDate = TUDAT_NAN;

    // Save each line to datamap
    for(int i = 0 ; i < numberOfLines ; i++ ){
        julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
                    solarActivityExtractor.extract( parsedDataVector->at( i ) )->year,
                    solarActivityExtractor.extract( parsedDataVector->at( i ) )->month,
                    solarActivityExtractor.extract( parsedDataVector->at( i ) )->day,
                    0, 0, 0.0 ) ;
        dataMap[ julianDate ] = solarActivityExtractor.extract( parsedDataVector->at( i ) ) ;
    }

    return dataMap;

}

} // solar_activity
} // input_output
} // tudat
