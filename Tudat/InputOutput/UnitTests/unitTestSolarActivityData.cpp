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
 *      Data source for validation figures:
 *                        http://celestrak.com/SpaceData/sw19571001.txt
 *
 */

#define BOOST_TEST_MAIN

#include <istream>
#include <string>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/InputOutput/parseSolarActivityData.h"
#include "Tudat/InputOutput/extractSolarActivityData.h"
#include "Tudat/InputOutput/solarActivityData.h"
#include "Tudat/InputOutput/basicInputOutput.h"


namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_Solar_Activity_Parser_Extractor )

//! Test parsing and extraction process of solar dummy activity file
BOOST_AUTO_TEST_CASE( test_parsing_and_extraction )
{
    tudat::input_output::parsed_data_vector_utilities::ParsedDataVectorPtr parsedDataVectorPtr;
    tudat::input_output::solar_activity::ParseSolarActivityData solarActivityParser;
    tudat::input_output::solar_activity::ExtractSolarActivityData solarActivityExtractor;

    // Parse file
    // save path of cpp file
    std::string cppPath( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string folder = cppPath.substr( 0, cppPath.find_last_of( "/\\" ) + 1 );
    std::string filePath = folder + "testSolarActivity.txt";

    // Open dataFile
    std::ifstream dataFile;
    dataFile.open( filePath.c_str( ), std::ifstream::in );
    parsedDataVectorPtr = solarActivityParser.parse( dataFile );
    dataFile.close( );

    // Extract data to object of solarActivityData class
    std::vector< boost::shared_ptr< tudat::input_output::solar_activity::SolarActivityData > >
            solarActivityData( 6 );
    solarActivityData[ 0 ] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 0 ) );
    solarActivityData[ 1 ] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 3 ) );
    solarActivityData[ 2 ] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 8 ) );
    solarActivityData[ 3 ] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 9 ) );
    solarActivityData[ 4 ] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 12 ) );
    solarActivityData[ 5 ] = solarActivityExtractor.extract( parsedDataVectorPtr->at( 15 ) );

    // Check single parameters
    BOOST_CHECK_EQUAL(  solarActivityData[ 0 ]->dataType, 1 );
    BOOST_CHECK_EQUAL(  solarActivityData[ 1 ]->planetaryDailyCharacterFigure, 0.0 );
    BOOST_CHECK_EQUAL(  solarActivityData[ 1 ]->solarRadioFlux107Adjusted, 103.0 );
    BOOST_CHECK_EQUAL(  solarActivityData[ 2 ]->month, 5 );
    BOOST_CHECK_EQUAL(  solarActivityData[ 2 ]->internationalSunspotNumber, 0 );
    BOOST_CHECK_EQUAL(  solarActivityData[ 3 ]->planetaryRangeIndexSum, 176 );
    BOOST_CHECK_EQUAL(  solarActivityData[ 3 ]->planetaryDailyCharacterFigureConverted, 0 );
    BOOST_CHECK_EQUAL(  solarActivityData[ 4 ]->bartelsSolarRotationNumber, 2441 );
    BOOST_CHECK_EQUAL(  solarActivityData[ 5 ]->last81DaySolarRadioFlux107Observed, 187.8 );

    // Check planetaryEquivalentAmplitudeVectors
    Eigen::VectorXd correctPlanetaryEquivalentAmplitudeVector1( 8 );
    correctPlanetaryEquivalentAmplitudeVector1 << 32, 27, 15, 7, 22, 9, 32, 22;
    Eigen::VectorXd correctPlanetaryEquivalentAmplitudeVector2 = Eigen::VectorXd::Zero( 8 );

    TUDAT_CHECK_MATRIX_CLOSE( solarActivityData[ 0 ]->planetaryEquivalentAmplitudeVector,
            correctPlanetaryEquivalentAmplitudeVector1, 0 );
    TUDAT_CHECK_MATRIX_CLOSE( solarActivityData[ 4 ]->planetaryEquivalentAmplitudeVector,
            correctPlanetaryEquivalentAmplitudeVector2 , 0);
    TUDAT_CHECK_MATRIX_CLOSE( solarActivityData[ 5 ]->planetaryEquivalentAmplitudeVector,
            correctPlanetaryEquivalentAmplitudeVector2, 0 );
}

BOOST_AUTO_TEST_CASE( test_function_readSolarActivityData )
{
    using tudat::input_output::solar_activity::SolarActivityDataMap;
    using tudat::input_output::solar_activity::SolarActivityData;

    // Parse file
    // save path of cpp file
    std::string cppPath( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string folder = cppPath.substr( 0, cppPath.find_last_of( "/\\" ) + 1 );
    std::string filePath = folder + "testSolarActivity.txt";

    // Read file.
    SolarActivityDataMap solarActivity = tudat::input_output::solar_activity::readSolarActivityData( filePath );

    // Select date and test selected entries.
    double julianDate;
    julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
                1957, 10, 3, 0, 0, 0.0 );

    SolarActivityDataMap::iterator solarActivityIterator;
    solarActivityIterator = solarActivity.find( julianDate );

    BOOST_CHECK_EQUAL( solarActivity.size( ), 18 );
    BOOST_CHECK_EQUAL( solarActivity.count( julianDate ), 1 );

    BOOST_CHECK_EQUAL( solarActivityIterator->second->day , 3 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->centered81DaySolarRadioFlux107Observed , 268.1 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->planetaryEquivalentAmplitudeAverage , 19 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->planetaryRangeIndexSum , 250 );

    // Select date and test selected entries.
    julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
                2023, 1, 1, 0, 0, 0.0 );

    solarActivityIterator = solarActivity.find( julianDate );

    BOOST_CHECK_EQUAL( solarActivity.count( julianDate ), 1 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->day , 1 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->bartelsSolarRotationNumber , 2583 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->centered81DaySolarRadioFlux107Observed , 0.0 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->last81DaySolarRadioFlux107Adjusted , 0.0 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->planetaryEquivalentAmplitudeAverage , 0.0 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->planetaryRangeIndexSum , 0.0 );

    // Test full line from full activity file
    std::string filePath2 = folder + "sw19571001.txt";

    // Read file.
    SolarActivityDataMap solarActivity2 = tudat::input_output::solar_activity::readSolarActivityData( filePath2 );

    // Define test time.
    julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
                1993, 12, 11, 0, 0, 0.0 );
    solarActivityIterator = solarActivity2.find( julianDate );
    BOOST_CHECK_EQUAL( solarActivity2.count( julianDate ), 1 );

    // Check contents at given julian day.
    BOOST_CHECK_EQUAL( solarActivityIterator->second->day , 11 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->bartelsSolarRotationNumber , 2190 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->dayOfBartelsCycle , 9 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->planetaryRangeIndexSum , 143 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->planetaryEquivalentAmplitudeAverage , 7 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->planetaryDailyCharacterFigure , 0.4 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->planetaryDailyCharacterFigureConverted , 2 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->internationalSunspotNumber , 31 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->solarRadioFlux107Adjusted , 89.7 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->fluxQualifier , 0 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->centered81DaySolarRadioFlux107Adjusted , 101.0 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->last81DaySolarRadioFlux107Adjusted , 97.6 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->solarRadioFlux107Observed , 92.5 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->centered81DaySolarRadioFlux107Observed , 104.0 );
    BOOST_CHECK_EQUAL( solarActivityIterator->second->last81DaySolarRadioFlux107Observed , 99.0 );

    Eigen::VectorXd expectedPlanetaryRangeIndices = ( Eigen::VectorXd( 8 ) << 33, 23, 20, 13, 17, 13, 17, 7 ).finished( );
    Eigen::VectorXd expectedPlanetaryEquivalentAmplitudes = ( Eigen::VectorXd( 8 ) << 18, 9, 7, 5, 6, 5, 6, 3 ).finished( );

    for( unsigned int i = 0; i < 8; i++ )
    {
        BOOST_CHECK_EQUAL( expectedPlanetaryRangeIndices( i ),
                           solarActivityIterator->second->planetaryRangeIndexVector( i ) );
        BOOST_CHECK_EQUAL( expectedPlanetaryEquivalentAmplitudes( i ),
                           solarActivityIterator->second->planetaryEquivalentAmplitudeVector( i ) );

    }
}

BOOST_AUTO_TEST_SUITE_END( )

}   // unit_tests
}   // tudat
