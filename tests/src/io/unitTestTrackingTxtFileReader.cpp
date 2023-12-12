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

// Some simplifications for shorter lines
namespace tio = tudat::input_output;
namespace tss = tudat::simulation_setup;

namespace tudat
{
namespace unit_tests
{

//! Temporary utility function to print arrays to std::cout
template< typename T >
void printArr(const T& arr)
{
  std::cout << "[ ";
  for (const auto& i : arr) {
    std::cout << i << ' ';
  }
  std::cout << "]\n";
}

// Utility function to get a single block from a datMap that maps keys to vectors
template< typename K, typename V >
std::map<K, V> extractBlockFromVectorMap(const std::map<K, std::vector<V>>& vectorMap, int blockIndex)
{
  std::map<K, V> singleBlock;
  for (const auto& pair : vectorMap) {
    singleBlock[pair.first] = pair.second[blockIndex];
  }
  return singleBlock;
}

//! A function that specifies a standard format for a file. A user can also do this if they often read the same file format
std::unique_ptr<tio::TrackingTxtFileContents> readVikingRangeFile(const std::string& fileName)
{
  std::vector<tio::TrackingFileField> columnTypes({
                                                      tio::TrackingFileField::spacecraft_id,
                                                      tio::TrackingFileField::dsn_transmitting_station_nr,
                                                      tio::TrackingFileField::dsn_receiving_station_nr,
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

BOOST_AUTO_TEST_SUITE(test_tracking_txt_file_reader);

BOOST_AUTO_TEST_CASE(VikingRangeDataCustomFunction)
{

  std::shared_ptr<tio::TrackingTxtFileContents> rawVikingFile = readVikingRangeFile(vikingRangePath);

  std::string spacecraftName = "Viking";

  auto rawDataMap = rawVikingFile->getRawDataMap();
  auto dataMap = rawVikingFile->getDoubleDataMap();
  auto dataBlock3 = extractBlockFromVectorMap(dataMap, 3);

  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::spacecraft_id], 1);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::dsn_transmitting_station_nr], 43);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::dsn_receiving_station_nr], 43);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::year], 1976);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::month], 7);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::day], 22);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::hour], 6);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::minute], 2);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::second], 32);
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::two_way_light_time], 2290.150246895);

  std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents>
      processedVikingFile = std::make_shared<observation_models::ProcessedTrackingTxtFileContents>(rawVikingFile, spacecraftName);
  auto observationCollection = observation_models::createTrackingTxtFileObservationCollection<double, double>(processedVikingFile);

  BOOST_CHECK_EQUAL(observationCollection->getTotalObservableSize(), 1258);
}

BOOST_AUTO_TEST_CASE(marsPathfinderRangeSimpleReading)
{
  std::vector<tio::TrackingFileField> fieldTypeVector{
      tio::TrackingFileField::spacecraft_id,
      tio::TrackingFileField::dsn_transmitting_station_nr,
      tio::TrackingFileField::dsn_receiving_station_nr,
      tio::TrackingFileField::year,
      tio::TrackingFileField::month_three_letter,
      tio::TrackingFileField::day,
      tio::TrackingFileField::hour,
      tio::TrackingFileField::minute,
      tio::TrackingFileField::second,
      tio::TrackingFileField::round_trip_light_time_microseconds,
      tio::TrackingFileField::light_time_measurement_accuracy_microseconds,
  };

  auto rawTrackingFile = createTrackingTxtFileContents(marsPathfinderRangePath, fieldTypeVector, '#', ",: \t");
  auto dataMap = rawTrackingFile->getDoubleDataMap();
  auto dataBlock4 = extractBlockFromVectorMap(dataMap, 4);

  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::spacecraft_id], 3);
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::dsn_transmitting_station_nr], 65);
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::dsn_receiving_station_nr], 65);
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::year], 1997);
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::month], 7);
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::day], 25);
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::hour], 18);
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::minute], 17);
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::second], 02);
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::two_way_light_time], 1420.476556473);
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::light_time_measurement_accuracy], 6.7e-8);

  BOOST_CHECK_EQUAL(rawTrackingFile->getNumColumns(), 11);
}

//
BOOST_AUTO_TEST_CASE(junoSimpleReading)
{
  std::vector<tio::TrackingFileField> fieldTypeVector{
      tio::TrackingFileField::spacecraft_id,
      tio::TrackingFileField::dsn_transmitting_station_nr,
      tio::TrackingFileField::dsn_receiving_station_nr,
      tio::TrackingFileField::year,
      tio::TrackingFileField::month,
      tio::TrackingFileField::day,
      tio::TrackingFileField::hour,
      tio::TrackingFileField::minute,
      tio::TrackingFileField::second,
      tio::TrackingFileField::round_trip_light_time_seconds,
      tio::TrackingFileField::light_time_measurement_delay_microseconds,
      tio::TrackingFileField::planet_nr,
      tio::TrackingFileField::tdb_seconds_j2000,
      tio::TrackingFileField::x_planet_frame_km,
      tio::TrackingFileField::y_planet_frame_km,
      tio::TrackingFileField::z_planet_frame_km,
      tio::TrackingFileField::vx_planet_frame_kms,
      tio::TrackingFileField::vy_planet_frame_kms,
      tio::TrackingFileField::vz_planet_frame_kms
  };

  auto rawTrackingFile = createTrackingTxtFileContents(junoRangePath, fieldTypeVector, '#', ",: \t");
  auto dataMap = rawTrackingFile->getDoubleDataMap();
  auto dataBlock0 = extractBlockFromVectorMap(dataMap, 0);

  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::spacecraft_id], 61);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::dsn_transmitting_station_nr], 55);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::dsn_receiving_station_nr], 55);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::year], 2016);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::month], 8);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::day], 27);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::hour], 13);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::minute], 45);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::second], 6);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::two_way_light_time], 6355.0487233317);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::light_time_measurement_accuracy], 0.0);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::planet_nr], 5);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::tdb_time_j2000], 525574396.542800, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::x_planet_frame], 976985.733, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::y_planet_frame], 68435520.227, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::z_planet_frame], 32772692.214, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::vx_planet_frame], 0727.110, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::vy_planet_frame], 26571.899, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::vz_planet_frame], -51299.726, 1e-4);

  BOOST_CHECK_EQUAL(rawTrackingFile->getNumColumns(), 19);
}

// Todo: Add a rigorous test with the observationcollection
//// Based on the odf file reader test
//BOOST_AUTO_TEST_CASE(testProcessTrackingFile)
//{
//
//  spice_interface::loadStandardSpiceKernels();
//
//  // presets
//  std::string spacecraftName = "Viking";
//
//  // Create system of bodies
//  std::vector<std::string> bodiesToCreate = {"Earth"};
//  tss::BodyListSettings bodySettings = tss::getDefaultBodySettings(bodiesToCreate);
//  bodySettings.at("Earth")->groundStationSettings = tss::getDsnStationSettings();
//
//  tss::SystemOfBodies bodies = createSystemOfBodies(bodySettings);
//
//  // Load Tracking file
//  std::shared_ptr<tio::TrackingTxtFileContents> rawTrackingTxtContents = readVikingRangeFile(vikingRangePath);
//
//  // Process Tracking file
//  std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents> processedTrackingTxtFileContents =
//      std::make_shared<observation_models::ProcessedTrackingTxtFileContents>(rawTrackingTxtContents, spacecraftName);
//
//  // Todo:
//  //  Implement Tests for the observation collection
//
//  // Check Time
//  std::pair<double, double> startAndEndTime = processedTrackingTxtFileContents->getStartAndEndTime();
//  double startTime = startAndEndTime.first;
//  double endTime = startAndEndTime.second;
//
//  // Extract Observation Collection
//  auto observationCollection = createTrackingTxtFileObservationCollection(processedTrackingTxtFileContents);
//}

BOOST_AUTO_TEST_SUITE_END();

}
}
