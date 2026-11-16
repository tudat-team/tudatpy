/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iterator>

#include "tudat/basics/testMacros.h"
#include "tudat/io/readCrdFile.h"
#include "tudat/io/readSinexFile.h"
#include "tudat/support/testFileUtilities.h"

namespace tudat
{
namespace unit_tests
{

namespace
{

std::string createTempPath( const std::string& suffix )
{
    return createTemporaryFilePath( "tudat-ilrs", suffix );
}

std::string makeStationInfoLine( const std::string& stationName, const std::string& domesId, const int monumentId )
{
    std::string line( 64, ' ' );
    const std::string monumentString = std::to_string( monumentId );
    line.replace( 0, std::min< unsigned int >( 24, stationName.size( ) ), stationName.substr( 0, 24 ) );
    line.replace( 33, std::min< unsigned int >( 9, domesId.size( ) ), domesId.substr( 0, 9 ) );
    line.replace( 48, std::min< unsigned int >( 4, monumentString.size( ) ), monumentString.substr( 0, 4 ) );
    return line;
}

void writeSyntheticStationInfoFile( const std::string& fileName )
{
    std::ofstream file( fileName );
    file << makeStationInfoLine( "ARECIBO", "12345M001", 7110 ) << std::endl;
}

void writeSyntheticSinexFile( const std::string& fileName )
{
    std::ofstream file( fileName );
    file << "+SITE/ID\n";
    file << "*CODE PT DOMES____ T STATION DESCRIPTION         APPROX_LON_ APPROX_LAT_ APPROX_H_\n";
    file << "ARCB A 12345M001 P ARECIBO                       0 0 0\n";
    file << "-SITE/ID\n";
    file << "+SOLUTION/ESTIMATE\n";
    file << "*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S ____ESTIMATED_VALUE ____STD_DEV___\n";
    file << "1 STAX ARCB A 1 25:150:00000 m 2 2390490.0 0.001\n";
    file << "2 STAY ARCB A 1 25:150:00000 m 2 -5564763.0 0.001\n";
    file << "3 STAZ ARCB A 1 25:150:00000 m 2 1994727.0 0.001\n";
    file << "4 VELX ARCB A 1 25:150:00000 m/y 2 31557600.0 0.001\n";
    file << "5 VELY ARCB A 1 25:150:00000 m/y 2 0.0 0.001\n";
    file << "6 VELZ ARCB A 1 25:150:00000 m/y 2 -31557600.0 0.001\n";
    file << "-SOLUTION/ESTIMATE\n";
}

void writeSyntheticSinexSiteIdOnlyFile( const std::string& fileName )
{
    std::ofstream file( fileName );
    file << "+SITE/ID\n";
    file << "*CODE PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_     _SOD#___\n";
    file << " 7110  A 12345M001 L ARECIBO    FIXED      66 45  0.0  18 20  0.0   400.0     71100901\n";
    file << "-SITE/ID\n";
}

void writeSyntheticCrdFile( const std::string& fileName )
{
    std::ofstream file( fileName );
    file << "H2 ARECIBO 7110 0 0 7\n";
    file << "H3 LAGEOS2 0 0 0 0 0\n";
    file << "H4 0 2025 5 30 23 59 0 2025 5 31 0 1 0 0 0 0 0 0 0 0 0\n";
    file << "C0 0 532.0 SYS1\n";
    file << "20 86398.0 1013.2 293.0 55.0\n";
    file << "10 86398.5 0.001 SYS1 2 0 0 0\n";
    file << "11 86399.0 0.002 SYS1 2 5.0 12 39.0\n";
    file << "11 1.0 0.004 SYS1 2 5.0 11 41.0\n";
    file << "H8\n";
    file << "H9\n";
}

void writeSyntheticSinexEccentricityFile( const std::string& fileName )
{
    std::ofstream file( fileName );
    file << "+SITE/ID\n";
    file << "*CODE PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_     _SOD#___\n";
    file << " 7110  A 12345M001 L ARECIBO    FIXED      66 45  0.0  18 20  0.0   400.0     71100901\n";
    file << "-SITE/ID\n";
    file << "+SITE/ECCENTRICITY\n";
    file << "*SITE PT SOLN T DATA_START__ DATA_END____ XYZ X_______ Y_______ Z_______        CDP-SOD_\n";
    file << " 7110  A    1 L 24:001:00000 24:365:86399 XYZ   0.1000   0.2000   0.3000        71100901\n";
    file << " 7110  A    1 L 25:001:00000 00:000:00000 XYZ   0.4000   0.5000   0.6000        71100901\n";
    file << "-SITE/ECCENTRICITY\n";
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_crd_and_sinex_readers )

BOOST_AUTO_TEST_CASE( testCrdParsingAndConversion )
{
    const std::string stationInfoFile = createTempPath( ".txt" );
    const std::string crdFile = createTempPath( ".npt" );

    writeSyntheticStationInfoFile( stationInfoFile );
    writeSyntheticCrdFile( crdFile );

    std::vector< input_output::CrdPass > crdPasses = input_output::readCrdFile( crdFile );
    BOOST_CHECK_EQUAL( crdPasses.size( ), 1 );
    BOOST_CHECK_EQUAL( crdPasses.at( 0 ).configuration_.cdpPadId_, 7110 );
    BOOST_CHECK_EQUAL( crdPasses.at( 0 ).data_.fullRateData_.size( ), 1 );
    BOOST_CHECK_EQUAL( crdPasses.at( 0 ).data_.normalPointData_.size( ), 2 );
    BOOST_CHECK_EQUAL( crdPasses.at( 0 ).data_.meteorologicalData_.size( ), 1 );
    BOOST_CHECK_CLOSE_FRACTION( crdPasses.at( 0 ).data_.normalPointData_.at( 0 ).normalPointWindowLength_, 5.0, 1.0E-15 );
    BOOST_CHECK_EQUAL( crdPasses.at( 0 ).data_.normalPointData_.at( 0 ).numberOfReturns_, 12 );
    BOOST_CHECK_CLOSE_FRACTION( crdPasses.at( 0 ).data_.normalPointData_.at( 0 ).binRms_, 39.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( crdPasses.at( 0 ).data_.fullRateData_.at( 0 ).oneWayRange_,
                                input_output::convertCrdTwoWayTimeOfFlightToSlrRange( 0.001 ),
                                1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( crdPasses.at( 0 ).data_.meteorologicalData_.at( 0 ).pressure_, 1013.2, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( crdPasses.at( 0 ).data_.meteorologicalData_.at( 0 ).temperature_, 293.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( crdPasses.at( 0 ).data_.meteorologicalData_.at( 0 ).humidity_, 55.0, 1.0E-15 );

    std::map< double, double > parsedNormalPoints = input_output::extractNormalPointMeasurements( crdPasses.at( 0 ) );
    BOOST_CHECK_EQUAL( parsedNormalPoints.size( ), 2 );
    BOOST_CHECK_CLOSE_FRACTION(
            parsedNormalPoints.begin( )->second, input_output::convertCrdTwoWayTimeOfFlightToSlrRange( 0.002 ), 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( std::next( parsedNormalPoints.begin( ) )->first - parsedNormalPoints.begin( )->first, 2.0, 1.0E-15 );
    std::map< double, double > parsedFullRate = input_output::extractFullRateMeasurements( crdPasses.at( 0 ) );
    BOOST_CHECK_EQUAL( parsedFullRate.size( ), 1 );

    std::map< int, std::string > monumentToStationMap = input_output::readMonumentNumbers( stationInfoFile );
    std::map< std::string, std::vector< input_output::CrdPass > > groupedData =
            input_output::groupCrdDataPerStation( crdPasses, monumentToStationMap );

    BOOST_CHECK_EQUAL( groupedData.size( ), 1 );
    BOOST_CHECK_EQUAL( groupedData.begin( )->first, "ARECIBO" );

    std::map< std::string, double > stationWavelengths = input_output::getStationWavelengths( groupedData );
    BOOST_CHECK_CLOSE_FRACTION( stationWavelengths.at( "ARECIBO" ), 532.0E-9, 1.0E-15 );

    std::remove( stationInfoFile.c_str( ) );
    std::remove( crdFile.c_str( ) );
}

BOOST_AUTO_TEST_CASE( testSinexStateParsing )
{
    const std::string stationInfoFile = createTempPath( ".txt" );
    const std::string sinexFile = createTempPath( ".snx" );

    writeSyntheticStationInfoFile( stationInfoFile );
    writeSyntheticSinexFile( sinexFile );

    std::map< std::string, std::string > stationToDomes = input_output::readDomesIdNumbers( stationInfoFile );
    std::map< std::string, std::string > domesToStation = input_output::readGroundStationNames( stationInfoFile );

    BOOST_CHECK_EQUAL( stationToDomes.at( "ARECIBO" ), "12345M001" );
    BOOST_CHECK_EQUAL( domesToStation.at( "12345M001" ), "ARECIBO" );

    std::map< std::string, input_output::SinexStationState > sinexData = input_output::readSinexStationData( sinexFile );
    BOOST_CHECK_EQUAL( sinexData.size( ), 1 );
    BOOST_CHECK( sinexData.count( "12345M001" ) > 0 );

    const input_output::SinexStationState state = sinexData.at( "12345M001" );
    BOOST_CHECK_EQUAL( state.siteCode_, "ARCB" );
    BOOST_CHECK_CLOSE_FRACTION( state.position_( 0 ), 2390490.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( state.position_( 1 ), -5564763.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( state.position_( 2 ), 1994727.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( state.velocity_( 0 ), 1.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( state.velocity_( 1 ), 0.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( state.velocity_( 2 ), -1.0, 1.0E-15 );
    BOOST_CHECK( state.referenceEpoch_ == state.referenceEpoch_ );

    std::remove( stationInfoFile.c_str( ) );
    std::remove( sinexFile.c_str( ) );
}

BOOST_AUTO_TEST_CASE( testSinexSiteIdMonumentMapping )
{
    const std::string sinexSiteIdFile = createTempPath( ".snx" );
    writeSyntheticSinexSiteIdOnlyFile( sinexSiteIdFile );

    std::map< int, std::string > monumentMap = input_output::readMonumentNumbers( sinexSiteIdFile );
    std::map< std::string, std::string > stationToDomes = input_output::readDomesIdNumbers( sinexSiteIdFile );
    std::map< std::string, std::string > domesToStation = input_output::readGroundStationNames( sinexSiteIdFile );

    BOOST_CHECK_EQUAL( monumentMap.at( 7110 ), "7110" );
    BOOST_CHECK_EQUAL( stationToDomes.at( "7110" ), "12345M001" );
    BOOST_CHECK_EQUAL( domesToStation.at( "12345M001" ), "7110" );

    std::remove( sinexSiteIdFile.c_str( ) );
}

BOOST_AUTO_TEST_CASE( testIlrsStationRegistryFromSinexSiteId )
{
    const std::string sinexSiteIdFile = createTempPath( ".snx" );
    writeSyntheticSinexSiteIdOnlyFile( sinexSiteIdFile );

    const std::map< int, input_output::IlrsStationRegistryEntry > stationRegistry =
            input_output::readIlrsStationRegistryFromSinexSiteId( sinexSiteIdFile );

    BOOST_CHECK_EQUAL( stationRegistry.size( ), 1 );
    BOOST_CHECK( stationRegistry.count( 7110 ) > 0 );
    BOOST_CHECK_EQUAL( stationRegistry.at( 7110 ).stationCode_, 7110 );
    BOOST_CHECK_EQUAL( stationRegistry.at( 7110 ).stationName_, "ARECIBO" );
    BOOST_CHECK_EQUAL( stationRegistry.at( 7110 ).domesId_, "12345M001" );
    BOOST_CHECK_CLOSE_FRACTION( stationRegistry.at( 7110 ).approximateLongitude_, 1.1650072757062149, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( stationRegistry.at( 7110 ).approximateLatitude_, 0.31997702953229373, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( stationRegistry.at( 7110 ).approximateHeight_, 400.0, 1.0E-15 );

    std::remove( sinexSiteIdFile.c_str( ) );
}

BOOST_AUTO_TEST_CASE( testSinexEccentricityParsing )
{
    const std::string eccentricityFile = createTempPath( ".snx" );
    writeSyntheticSinexEccentricityFile( eccentricityFile );

    const std::map< std::string, std::vector< input_output::SinexStationEccentricity > > eccentricityData =
            input_output::readSinexStationEccentricities( eccentricityFile );
    BOOST_CHECK_EQUAL( eccentricityData.size( ), 1 );
    BOOST_CHECK( eccentricityData.count( "12345M001" ) > 0 );
    BOOST_CHECK_EQUAL( eccentricityData.at( "12345M001" ).size( ), 2 );

    const input_output::SinexStationEccentricity firstEntry = eccentricityData.at( "12345M001" ).at( 0 );
    const input_output::SinexStationEccentricity secondEntry = eccentricityData.at( "12345M001" ).at( 1 );

    BOOST_CHECK_EQUAL( firstEntry.stationCode_, 7110 );
    BOOST_CHECK_EQUAL( secondEntry.stationCode_, 7110 );
    BOOST_CHECK_CLOSE_FRACTION( firstEntry.eccentricity_( 0 ), 0.1, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( firstEntry.eccentricity_( 1 ), 0.2, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( firstEntry.eccentricity_( 2 ), 0.3, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( secondEntry.eccentricity_( 0 ), 0.4, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( secondEntry.eccentricity_( 1 ), 0.5, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( secondEntry.eccentricity_( 2 ), 0.6, 1.0E-15 );
    BOOST_CHECK_EQUAL( secondEntry.hasOpenEnd_, true );

    std::remove( eccentricityFile.c_str( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
