/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>
#include <utility>
#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/estimation_setup/observations.h"


// Temporary
#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"

namespace tio = tudat::input_output;

namespace tudat
{
namespace unit_tests
{

template< typename T >
void printArr(const T& arr)
{
  std::cout << "[ ";
  for (const auto& i : arr) {
    std::cout << i << ' ';
  }
  std::cout << "]\n";
}

// Here, we implement a function that specifies a standard format for a file.
std::unique_ptr<tio::TrackingTxtFileContents> readVikingRangeFile(const std::string& fileName)
{
  std::vector<tio::TrackingFileField> columnTypes({
                                                      tio::TrackingFileField::spacecraft_id,
                                                      tio::TrackingFileField::dsn_transmitting_station_id,
                                                      tio::TrackingFileField::dsn_receiving_station_id,
                                                      tio::TrackingFileField::year,
                                                      tio::TrackingFileField::month_three_letter,
                                                      tio::TrackingFileField::day,
                                                      tio::TrackingFileField::hour,
                                                      tio::TrackingFileField::minute,
                                                      tio::TrackingFileField::second,
                                                      tio::TrackingFileField::round_trip_light_time_microseconds,
                                                      tio::TrackingFileField::light_time_measurement_delay_microseconds
                                                  });
  auto vikingFile = createTrackingTxtFileContents(fileName, columnTypes);
  vikingFile->addMetaData(tio::TrackingFileField::file_title, "Viking lander range data");

  return vikingFile;
}

// Setting some path variables for the test files
const std::string vikingRangePath = tudat::paths::getTudatTestDataPath() + "vikingrange.txt";
const std::string marsPathfinderRangePath = tudat::paths::getTudatTestDataPath() + "mars-pathfinder-range.txt";
const std::string junoRangePath = tudat::paths::getTudatTestDataPath() + "juno_range.txt";

BOOST_AUTO_TEST_SUITE(test_generic_txt_file_reader);


BOOST_AUTO_TEST_CASE(JplRangeDataCustomFunction)
{
  std::shared_ptr<tio::TrackingTxtFileContents> vikingFile = readVikingRangeFile(vikingRangePath);

  auto rawDataMap = vikingFile->getRawDataMap();
  auto dataMap = vikingFile->getDoubleDataMap();
  std::vector<tio::TrackingDataType> dataColumnTypes = vikingFile->getDataColumnTypes();

  std::map<tio::TrackingDataType, double> dataBlock3;
  for (auto columnType : dataColumnTypes) {
    dataBlock3[columnType] = dataMap.at(columnType)[3];
  }

  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::spacecraft_id], 1);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::dsn_transmitting_station_id], 43);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::dsn_receiving_station_id], 43);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::year], 1976);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::month], 7);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::day], 22);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::hour], 6);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::minute], 2);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::second], 32);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::two_way_light_time], 2290.150246895);

  // FIXME: this is incorrect!
  observation_models::LinkEnds linkEnds{{observation_models::reflector, observation_models::linkEndId("viking")}};

  std::map<observation_models::ObservableType, observation_models::LinkEnds> observableTypes = {{observation_models::n_way_range, linkEnds}};
//  std::shared_ptr<observation_models::ObservationCollection<double, double>>
//      observationCollection = observation_models::createTrackingTxtFileObservationCollection(vikingFile, observableTypes);
}


//BOOST_AUTO_TEST_CASE( JplRangeDataCustomClass )
//{
//  auto vikingFile = readVikingRangeFile( vikingRangePath );
//  std::map< std::shared_ptr< tio::TxtFieldType >, double > dataBlock3 = vikingFile->dataVector_.at( 3 );
//
//  BOOST_CHECK_EQUAL( dataBlock3[tft::spacecraftIdentifier], 1 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::dsnTransmittingStation], 43 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::dsnReceivingStation], 43 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcYear], 1976 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcMonth], 7 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcDay], 22 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcHour], 6 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcMinute], 2 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcSecond], 32 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::roundTripLightTimeMicroSec], 2290150246.895 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::lightTimeMeasurementAccuracyMicroSec], 0.047 );
//}
//
//
//BOOST_AUTO_TEST_CASE( GenericTxtFile )
//{
//  std::vector< std::shared_ptr< tio::TxtFieldType > > fieldTypeVector { tft::spacecraftIdentifier,
//                                                                        tft::dsnTransmittingStation,
//                                                                        tft::dsnReceivingStation, tft::utcYear,
//                                                                        tft::utcMonth, tft::utcDay, tft::utcHour,
//                                                                        tft::utcMinute, tft::utcSecond,
//                                                                        tft::roundTripLightTimeMicroSec,
//                                                                        tft::lightTimeMeasurementAccuracyMicroSec
//  };
//  std::shared_ptr< tio::TrackingTxtFileContents > fileContents = tio::createTxtFileContents( vikingRangePath,
//                                                                                     fieldTypeVector,
//                                                                                     '#',
//                                                                                     ",: \t" );
//
//  std::map< std::shared_ptr< tio::TxtFieldType >, double > dataBlock3 = fileContents->dataVector_.at( 3 );
//
//  BOOST_CHECK_EQUAL( dataBlock3[tft::spacecraftIdentifier], 1 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::dsnTransmittingStation], 43 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::dsnReceivingStation], 43 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcYear], 1976 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcMonth], 7 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcDay], 22 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcHour], 6 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcMinute], 2 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcSecond], 32 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::roundTripLightTimeMicroSec], 2290150246.895 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::lightTimeMeasurementAccuracyMicroSec], 0.047 );
//
//  BOOST_CHECK_EQUAL( fileContents->getNumColumns( ), 11 );
//  BOOST_CHECK_EQUAL( fileContents->dataVector_.size( ), 1258 );
//}
//
//
//BOOST_AUTO_TEST_CASE( JunoRead )
//{
//  std::vector< std::shared_ptr< tio::TxtFieldType > > fieldTypeVector { tft::spacecraftIdentifier,
//                                                                        tft::dsnTransmittingStation,
//                                                                        tft::dsnReceivingStation, tft::utcYear,
//                                                                        tft::utcMonth, tft::utcDay, tft::utcHour,
//                                                                        tft::utcMinute, tft::utcSecond,
//                                                                        tft::roundTripLightTimeSec,
//                                                                        tft::lightTimeMeasurementDelayMicroSec,
//                                                                        tft::planetNumber, tft::tbdTimeJ2000,
//                                                                        tft::xPlanetFrame, tft::yPlanetFrame,
//                                                                        tft::zPlanetFrame, tft::vXPlanetFrame,
//                                                                        tft::vYPlanetFrame, tft::vZPlanetFrame
//  };
//  std::shared_ptr< tio::TrackingTxtFileContents > fileContents = tio::createTxtFileContents( junoRangePath,
//                                                                                     fieldTypeVector,
//                                                                                     '#',
//                                                                                     ",: \t" );
//
//  std::map< std::shared_ptr< tio::TxtFieldType >, double > dataBlock3 = fileContents->dataVector_.at( 0 );
//
//  BOOST_CHECK_EQUAL( dataBlock3[tft::spacecraftIdentifier], 61 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::dsnTransmittingStation], 55 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::dsnReceivingStation], 55 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcYear], 2016 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcMonth], 8 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcDay], 27 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcHour], 13 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcMinute], 45 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::utcSecond], 6 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::roundTripLightTimeSec], 6355.0487233317 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::lightTimeMeasurementAccuracyMicroSec], 0.0 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::planetNumber], 5 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::tbdTimeJ2000], 525574396.542800 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::xPlanetFrame], 976.985733 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::yPlanetFrame], 68435.520227 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::zPlanetFrame], 32772.692214 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::vXPlanetFrame], 0.727110 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::vYPlanetFrame], 26.571899 );
//  BOOST_CHECK_EQUAL( dataBlock3[tft::vZPlanetFrame], -51.299726 );
//
//  BOOST_CHECK_EQUAL( fileContents->getNumColumns( ), 19 );
//  BOOST_CHECK_EQUAL( fileContents->dataVector_.size( ), 4 );


BOOST_AUTO_TEST_SUITE_END();
}
}
