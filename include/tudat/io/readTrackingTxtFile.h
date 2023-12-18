/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_READ_GENERIC_TXT_FILE_H
#define TUDAT_READ_GENERIC_TXT_FILE_H

#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>
#include <string>
#include <cstdarg>
#include <fstream>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "tudat/io/fieldType.h"
#include "tudat/astro/basic_astro.h"
#include "tudat/astro/observation_models/observableTypes.h"

// TODO: Split into header and source file

// A TrackingDataType is a unique form of data that can be used by Tudat to derive the observables
// A TrackingFileField is an identifier that specifies a specific type of column data in an input file (with a specific format)
// A TrackingFileFieldConverter is an interface that can convert a raw string input field to its associated double value as expected for the TrackingDataType
// For example:
// - You are reading a file with a column that shows the round trip light time in microseconds
// - The corresponding TrackingFileField will be TrackingFileField::round_trip_light_time_microseconds
// - Tudat wants to refer to the data by TrackingDataType::n_way_light_time in seconds
// - The converter

namespace tudat
{
namespace input_output
{

//! Utility function to get a value from a string map where the keys are all uppercase
template< typename T >
T upperCaseFromMap(const std::string& strValue, const std::map<std::string, T>& upperCaseMapping)
{
  std::string upperCaseStrValue = boost::to_upper_copy(strValue);
  const auto iter = upperCaseMapping.find(upperCaseStrValue);
  if (iter != upperCaseMapping.cend()) {
    return iter->second;
  }
  throw std::runtime_error("Invalid key found in map while converting tracking data");
}

//! Enum describing a unique data type that can later be used to process the information. These are in SI units
enum class TrackingDataType
{
  year,
  month,
  day,
  hour,
  minute,
  second,
  observation_time_scale,
  file_name,
  n_way_light_time,
  light_time_measurement_delay,
  light_time_measurement_accuracy,
  dsn_transmitting_station_nr,
  dsn_receiving_station_nr,
  spacecraft_id,
  planet_nr,
  tdb_time_j2000,
  x_planet_frame,
  y_planet_frame,
  z_planet_frame,
  vx_planet_frame,
  vy_planet_frame,
  vz_planet_frame,
  residual_de405,
};

//! Enum describing a unique data type and format that can be present in a column of a file. Note that multiple `TrackingFileField`
//! types can represent the same type of `TrackingDataType`
enum class TrackingFileField
{
  ignore = -1,
  year,
  month,
  month_three_letter,
  day,
  hour,
  minute,
  second,
  round_trip_light_time_microseconds,
  round_trip_light_time_seconds,
  time_scale,
  file_title,
  light_time_measurement_delay_microseconds,
  light_time_measurement_delay_seconds,
  light_time_measurement_accuracy_microseconds,
  spacecraft_id,
  dsn_transmitting_station_nr,
  dsn_receiving_station_nr,
  planet_nr,
  tdb_seconds_j2000,
  x_planet_frame_km,
  y_planet_frame_km,
  z_planet_frame_km,
  vx_planet_frame_kms,
  vy_planet_frame_kms,
  vz_planet_frame_kms,
  residual_de405_microseconds,
};

//! Simple converter class that can convert a string data field to a double. One can inherit from this and overload the
//! `toDouble()` method to extend the supported formats
class TrackingFileFieldConverter
{
public:
  explicit TrackingFileFieldConverter(TrackingDataType trackingDataType)
      : doubleDataType_(trackingDataType) {}

  virtual ~TrackingFileFieldConverter() = default;

  virtual double toDouble(std::string& rawField) const
  {
    try {
      return std::stod(rawField);
    } catch (std::invalid_argument&) {
      throw std::runtime_error(
          "The tracking file field cannot be converted correctly. Check your columnTypes. Raw field was \"" + rawField
              + "\".\n");
    }
  }
  const TrackingDataType& getTrackingDataType() { return doubleDataType_; }

private:
  TrackingDataType doubleDataType_;
};

//! A converter specifically for month fields in three-letter representation ("jan", "Jan", ...)
class TrackingFileMonthFieldConverter : public TrackingFileFieldConverter
{
public:
  TrackingFileMonthFieldConverter(TrackingDataType trackingDataType) : TrackingFileFieldConverter(trackingDataType) {}
  double toDouble(std::string& rawField) const
  {
    std::map<std::string, double> monthsMap{{"JAN", 1.},
                                            {"FEB", 2.},
                                            {"MAR", 3.},
                                            {"APR", 4.},
                                            {"MAY", 5.},
                                            {"JUN", 6.},
                                            {"JUL", 7.},
                                            {"AUG", 8.},
                                            {"SEP", 9.},
                                            {"OCT", 10.},
                                            {"NOV", 11.},
                                            {"DEC", 12.}};

    return upperCaseFromMap(rawField, monthsMap);
  }
};

//! Converter that will convert the raw string to double and then apply a scalar multiplier as specified
class TrackingFileFieldMultiplyingConverter : public TrackingFileFieldConverter
{
public:
  TrackingFileFieldMultiplyingConverter(TrackingDataType trackingDataType, double multiplier)
      : TrackingFileFieldConverter(trackingDataType), multiplier_(multiplier) {}
  double toDouble(std::string& rawField) const
  {
    return multiplier_ * TrackingFileFieldConverter::toDouble(rawField);
  }
private:
  double multiplier_;
};

//! Mapping the `TrackingFileField` to the correct converter, including the `TrackingDataType` it will represent
static const std::map<TrackingFileField, std::shared_ptr<TrackingFileFieldConverter>> trackingFileFieldConverterMap = {
    {TrackingFileField::spacecraft_id, std::make_shared<TrackingFileFieldConverter>(TrackingDataType::spacecraft_id)},
    {
        TrackingFileField::dsn_transmitting_station_nr,
        std::make_shared<TrackingFileFieldConverter>(TrackingDataType::dsn_transmitting_station_nr)
    },
    {
        TrackingFileField::dsn_receiving_station_nr,
        std::make_shared<TrackingFileFieldConverter>(TrackingDataType::dsn_receiving_station_nr)
    },
    {TrackingFileField::year, std::make_shared<TrackingFileFieldConverter>(TrackingDataType::year)},
    {TrackingFileField::month, std::make_shared<TrackingFileFieldConverter>(TrackingDataType::month)},
    {TrackingFileField::month_three_letter, std::make_shared<TrackingFileMonthFieldConverter>(TrackingDataType::month)},
    {TrackingFileField::day, std::make_shared<TrackingFileFieldConverter>(TrackingDataType::day)},
    {TrackingFileField::hour, std::make_shared<TrackingFileFieldConverter>(TrackingDataType::hour)},
    {TrackingFileField::minute, std::make_shared<TrackingFileFieldConverter>(TrackingDataType::minute)},
    {TrackingFileField::second, std::make_shared<TrackingFileFieldConverter>(TrackingDataType::second)},
    {TrackingFileField::round_trip_light_time_seconds, std::make_shared<TrackingFileFieldConverter>(TrackingDataType::n_way_light_time)},
    {
        TrackingFileField::round_trip_light_time_microseconds,
        std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::n_way_light_time, 1.e-6)
    },
    {
        TrackingFileField::light_time_measurement_delay_microseconds,
        std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::light_time_measurement_delay, 1.e-6)
    },
    {
        TrackingFileField::light_time_measurement_accuracy_microseconds,
        std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::light_time_measurement_accuracy, 1.e-6)
    },
    {TrackingFileField::planet_nr, std::make_shared<TrackingFileFieldConverter>(TrackingDataType::planet_nr)},
    {TrackingFileField::tdb_seconds_j2000, std::make_shared<TrackingFileFieldConverter>(TrackingDataType::tdb_time_j2000)},
    {TrackingFileField::x_planet_frame_km, std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::x_planet_frame, 1.e3)},
    {TrackingFileField::y_planet_frame_km, std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::y_planet_frame, 1.e3)},
    {TrackingFileField::z_planet_frame_km, std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::z_planet_frame, 1.e3)},
    {TrackingFileField::vx_planet_frame_kms, std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::vx_planet_frame, 1.e3)},
    {TrackingFileField::vy_planet_frame_kms, std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::vy_planet_frame, 1.e3)},
    {TrackingFileField::vz_planet_frame_kms, std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::vz_planet_frame, 1.e3)},
    {TrackingFileField::residual_de405_microseconds, std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::residual_de405, 1.e-6)},
};

//! Class to extract the raw data from a file with the appropriate conversion to doubles
class TrackingTxtFileContents
{
public:
  TrackingTxtFileContents(std::string fileName,
                          std::vector<TrackingFileField> columnTypes,
                          char commentSymbol = '#',
                          std::string valueSeparators = ",: \t")
      : fileName_(std::move(fileName)), columnFieldTypes_(std::move(columnTypes)), commentSymbol_(commentSymbol),
        valueSeparators_(std::move(valueSeparators))
  {
    parseData();
  }

  void parseData();

  void readRawDataMap(std::ifstream& dataFile);

  void addLineToRawDataMap(std::string& rawLine);

  void convertDataMap();

  void addMetaData(TrackingFileField fieldType, const std::string& value);

// Getters
public:
  size_t getNumColumns() const { return columnFieldTypes_.size(); }
  size_t getNumRows() const { return rawDataMap_.at(columnFieldTypes_[0]).size(); }
  const std::vector<TrackingFileField>& getRawColumnTypes() { return columnFieldTypes_; }

  const std::vector<TrackingDataType>& getDataColumnTypes()
  {
    columnDataTypes_.clear();
    for (auto& pair : doubleDataMap_) {
      columnDataTypes_.push_back(pair.first);
    }
    return columnDataTypes_;
  }

  const auto& getRawDataMap() { return rawDataMap_; }
  const auto& getDoubleDataMap() { return doubleDataMap_; }
  const auto& getMetaDataMap() { return metaDataMap_; }

private:
  std::string fileName_ = "None";
  std::string separators_ = ":, \t";
  std::vector<TrackingFileField> columnFieldTypes_;
  std::vector<TrackingDataType> columnDataTypes_;
  char commentSymbol_;
  std::string valueSeparators_;

  std::map<TrackingFileField, std::vector<std::string>> rawDataMap_;
  std::map<TrackingFileField, std::string> metaDataMap_;

  std::map<TrackingDataType, std::vector<int>> intDataMap_;
  std::map<TrackingDataType, std::vector<double>> doubleDataMap_;

// Utility variables
private:
  std::vector<std::string> currentSplitRawLine_;
};

static inline std::unique_ptr<TrackingTxtFileContents> createTrackingTxtFileContents(const std::string& fileName,
                                                                                     std::vector<TrackingFileField>& columnTypes,
                                                                                     char commentSymbol = '#',
                                                                                     const std::string& valueSeparators = ",: \t")
{
  return std::make_unique<TrackingTxtFileContents>(fileName, columnTypes, commentSymbol, valueSeparators);
}

} // namespace input_output
} // namespace tudat

#endif // TUDAT_READ_GENERIC_TXT_FILE_H
