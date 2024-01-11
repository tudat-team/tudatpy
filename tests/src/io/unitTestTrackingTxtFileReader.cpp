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
using namespace observation_models;

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

//! Utility function to print link Ends
std::ostream& operator<<(std::ostream& os, const observation_models::LinkEnds& linkEnds)
{
  os << "LinkEnds{\n";
  for (auto& linkEnd : linkEnds) {
    os << "\t{type: " << std::to_string(linkEnd.first)
       << ", body: " << linkEnd.second.getBodyName()
       << ", station: " << linkEnd.second.getStationName()
       << "}\n";
  }
  os << "}";
  return os;
}

// Utility function to get a single block from a datMap that maps keys to vectors
template< typename K, typename V >
std::map<K, V> extractBlockFromVectorMap(const std::map<K, std::vector<V>>& vectorMap, int blockIndex)
{
  std::map<K, V> singleBlock;
  for (const auto& pair : vectorMap) {
    if (blockIndex < 0) {
      blockIndex += pair.second.size();
    }
    singleBlock[pair.first] = pair.second.at(blockIndex);
  }
  return singleBlock;
}

//! A function that specifies a standard format for a file. A user can also do this if they often read the same file format
std::unique_ptr<tio::TrackingTxtFileContents> readVikingRangeFile(const std::string& fileName)
{
  std::vector<std::string> columnTypes({
                                           "spacecraft_id",
                                           "dsn_transmitting_station_nr",
                                           "dsn_receiving_station_nr",
                                           "year",
                                           "month_three_letter",
                                           "day",
                                           "hour",
                                           "minute",
                                           "second",
                                           "round_trip_light_time_microseconds",
                                           "light_time_measurement_delay_microseconds"
                                       });
  auto vikingFile = tio::createTrackingTxtFileContents(fileName, columnTypes);
  vikingFile->addMetaData(tio::TrackingDataType::file_name, "Viking lander range data");
  return vikingFile;
}

//! A function that specifies a standard format for a file. A user can also do this if they often read the same file format
std::unique_ptr<tio::TrackingTxtFileContents> readJuiceFdetsFile(const std::string& fileName)
{
  std::vector<std::string>
      columnTypes({"utc_datetime_string", "signal_to_noise_ratio", "normalised_spectral_max", "doppler_measured_frequency_hz", "doppler_noise_hz",});

  auto rawFileContents = tio::createTrackingTxtFileContents(fileName, columnTypes, '#', ", \t");
  rawFileContents->addMetaData(tio::TrackingDataType::file_name, "JUICE Fdets Test File");
  return rawFileContents;
}

// Setting some path variables for the test files
const std::string vikingRangePath = tudat::paths::getTudatTestDataPath() + "vikingrange.txt";
const std::string marsPathfinderRangePath = tudat::paths::getTudatTestDataPath() + "mars-pathfinder-range.txt";
const std::string junoRangePath = tudat::paths::getTudatTestDataPath() + "juno_range.txt";
const std::string marinerRangePath = tudat::paths::getTudatTestDataPath() + "mariner9obs.txt";
const std::string juiceFdetsDopplerPath = tudat::paths::getTudatTestDataPath() + "Fdets.jui2023.04.26.Hb.0006.r2i.txt";

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
  BOOST_CHECK_EQUAL(dataBlock3[tio::TrackingDataType::n_way_light_time], 2290.150246895);

  std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents>
      processedVikingFile = std::make_shared<observation_models::ProcessedTrackingTxtFileContents>(rawVikingFile, spacecraftName);
  auto observationCollection = observation_models::createTrackingTxtFileObservationCollection<double, double>(processedVikingFile);

  BOOST_CHECK_EQUAL(observationCollection->getTotalObservableSize(), 1258);
}

BOOST_AUTO_TEST_CASE(marsPathfinderRangeSimpleReading)
{
  std::vector<std::string> fieldTypeVector{
      "spacecraft_id",
      "dsn_transmitting_station_nr",
      "dsn_receiving_station_nr",
      "year",
      "month_three_letter",
      "day",
      "hour",
      "minute",
      "second",
      "round_trip_light_time_microseconds",
      "light_time_measurement_accuracy_microseconds",
  };

  auto rawTrackingFile = tio::createTrackingTxtFileContents(marsPathfinderRangePath, fieldTypeVector, '#', ",: \t");
  rawTrackingFile->addMetaData(tio::TrackingDataType::spacecraft_transponder_delay, 420.e-6);
  rawTrackingFile->addMetaData(tio::TrackingDataType::uplink_frequency, 7.2e9);
  rawTrackingFile->addMetaData(tio::TrackingDataType::downlink_frequency, 8.4e9);

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
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::n_way_light_time], 1420.476556473);
  BOOST_CHECK_EQUAL(dataBlock4[tio::TrackingDataType::light_time_measurement_accuracy], 6.7e-8);

  BOOST_CHECK_EQUAL(rawTrackingFile->getNumColumns(), 11);
  const auto& metaDataDoubleMap = rawTrackingFile->getMetaDataDoubleMap();
  BOOST_CHECK_EQUAL(metaDataDoubleMap.at(tio::TrackingDataType::spacecraft_transponder_delay), 420.e-6);
}

//
BOOST_AUTO_TEST_CASE(junoSimpleReading)
{
  std::vector<std::string> fieldTypeVector{
      "spacecraft_id",
      "dsn_transmitting_station_nr",
      "dsn_receiving_station_nr",
      "year",
      "month",
      "day",
      "hour",
      "minute",
      "second",
      "round_trip_light_time_seconds",
      "light_time_measurement_delay_microseconds",
      "planet_nr",
      "tdb_spacecraft_seconds_j2000",
      "x_planet_frame_km",
      "y_planet_frame_km",
      "z_planet_frame_km",
      "vx_planet_frame_kms",
      "vy_planet_frame_kms",
      "vz_planet_frame_kms"
  };

  auto rawTrackingFile = tio::createTrackingTxtFileContents(junoRangePath, fieldTypeVector, '#', ",: \t");
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
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::n_way_light_time], 6355.0487233317);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::light_time_measurement_accuracy], 0.0);
  BOOST_CHECK_EQUAL(dataBlock0[tio::TrackingDataType::planet_nr], 5);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::tdb_spacecraft_j2000], 525574396.542800, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::x_planet_frame], 976985.733, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::y_planet_frame], 68435520.227, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::z_planet_frame], 32772692.214, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::vx_planet_frame], 0727.110, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::vy_planet_frame], 26571.899, 1e-4);
  BOOST_CHECK_CLOSE(dataBlock0[tio::TrackingDataType::vz_planet_frame], -51299.726, 1e-4);

  BOOST_CHECK_EQUAL(rawTrackingFile->getNumColumns(), 19);
}

BOOST_AUTO_TEST_CASE(marinerSimpleReading)
{
  std::vector<std::string> fieldTypeVector{
      "year",
      "month",
      "day",
      "hour",
      "minute",
      "second",
      "round_trip_light_time_microseconds",
      "light_time_measurement_accuracy_microseconds",
      "residual_de405_microseconds"
//      "test"
  };

  auto rawTrackingFile = tio::createTrackingTxtFileContents(marinerRangePath, fieldTypeVector, '#', ",: \t");
  rawTrackingFile->addMetaData(input_output::TrackingDataType::observation_body, "Earth");
  rawTrackingFile->addMetaData(input_output::TrackingDataType::observed_body, "Mars");

  auto dataMap = rawTrackingFile->getDoubleDataMap();
  auto dataBlockLast = extractBlockFromVectorMap(dataMap, -1);


//  1972 10 12 00:06:02  2610383946.989   0.475  -0.226

  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::year], 1972);
  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::month], 10);
  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::day], 12);
  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::hour], 0);
  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::minute], 6);
  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::second], 2);
  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::n_way_light_time], 2610.383946989);
  BOOST_CHECK_CLOSE(dataBlockLast[tio::TrackingDataType::light_time_measurement_accuracy], 0.475e-6, 1e-10);
  BOOST_CHECK_CLOSE(dataBlockLast[tio::TrackingDataType::residual_de405], -0.226e-6, 1e-10);

  BOOST_CHECK_EQUAL(rawTrackingFile->getNumColumns(), 9);

  const auto& metaDataStrMap = rawTrackingFile->getMetaDataStrMap();
  BOOST_CHECK_EQUAL(metaDataStrMap.at(input_output::TrackingDataType::observation_body), "Earth");
  BOOST_CHECK_EQUAL(metaDataStrMap.at(input_output::TrackingDataType::observed_body), "Mars");
}

BOOST_AUTO_TEST_CASE(TestVikingRangeDataObservationCollection)
{

  // Load the observations from the Viking file
  std::shared_ptr<tio::TrackingTxtFileContents> rawVikingFile = readVikingRangeFile(vikingRangePath);
  auto observationCollection = observation_models::createTrackingTxtFileObservationCollection<double, double>(rawVikingFile, "Viking", {n_way_range});

  // Check size of observations
  BOOST_CHECK_EQUAL(observationCollection->getTotalObservableSize(), 1258);

  // Check if n-way-range is present in the collection
  const auto& observationTypeStartSize = observationCollection->getObservationTypeStartAndSize();
  BOOST_CHECK(observationTypeStartSize.find(n_way_range) != observationTypeStartSize.end());

  //  Checking the first element with station 63 - 63 in the Viking file
  const LinkDefinition linkDefDsn63({
                                        {transmitter, LinkEndId("Earth", "DSS-63")},
                                        {reflector, LinkEndId("Viking", "Antenna")},
                                        {receiver, LinkEndId("Earth", "DSS-63")},
                                    });

  auto observationsAndTimesDsn63 = observationCollection->getSingleLinkObservationsAndTimes(n_way_range, linkDefDsn63);
  auto observationsDsn63 = observationsAndTimesDsn63.first;
  auto timesDsn63 = observationsAndTimesDsn63.second;

  BOOST_CHECK_CLOSE(observationsDsn63(0, 0), 2371.564782809, 1e-12);
  BOOST_CHECK_CLOSE(timesDsn63[0], 2442999.2528009 - basic_astrodynamics::JULIAN_DAY_ON_J2000, 1e-8);


//  const auto& concatenatedObservations = observationCollection->getObservationVectorReference();
//  std::vector<double> observationTimes = observationCollection->getConcatenatedTimeVector();
//  std::cout << observationTimes[0] << "\t" << concatenatedObservations(0, 0) << "\n";

//  std::vector<int> linkEndIds = observationCollection->getConcatenatedLinkEndIds();
//  std::map<int, observation_models::LinkEnds> linkEndIdMap = observationCollection->getInverseLinkEndIdentifierMap();
//  for (auto& pair : linkEndIdMap) {
//    std::cout << pair.first << " - " << pair.second;
//  }
}

BOOST_AUTO_TEST_CASE(TestJuiceFile)
{
  spice_interface::loadStandardSpiceKernels();

  std::shared_ptr<tio::TrackingTxtFileContents> rawFdetsDopplerFile = readJuiceFdetsFile(juiceFdetsDopplerPath);
  rawFdetsDopplerFile->addMetaData(tio::TrackingDataType::doppler_base_frequency, 8420.0e6);
  rawFdetsDopplerFile->addMetaData(tio::TrackingDataType::doppler_bandwidth, 2.0e3);
  rawFdetsDopplerFile->addMetaData(tio::TrackingDataType::vlbi_station_name, "Hb");


  // CHECK THE RAW FILE
  BOOST_CHECK_EQUAL(rawFdetsDopplerFile->getNumColumns(), 5);

  auto dataMap = rawFdetsDopplerFile->getDoubleDataMap();
  auto dataBlockLast = extractBlockFromVectorMap(dataMap, -1);

  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::tdb_time_j2000], 735687970.0 + 32.184 + 37);
  //  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::tdb_time_j2000], basic_astrodynamics::convertTTtoTAI(735687970.0));

  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::signal_to_noise], 6.766652540970647242e+05);
  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::spectral_max], 5.754946258897545704e+03);
  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::doppler_measured_frequency], 5977954.253958693705);
  BOOST_CHECK_EQUAL(dataBlockLast[tio::TrackingDataType::doppler_noise], -1.3803018149815216e-02);


  // Process file
  auto processedFdetsDopplerFile = std::make_shared<ProcessedTrackingTxtFileContents>(rawFdetsDopplerFile, "JUICE");
  auto observationCollection = observation_models::createTrackingTxtFileObservationCollection<double, double>(processedFdetsDopplerFile);

}

BOOST_AUTO_TEST_SUITE_END();

}
}
