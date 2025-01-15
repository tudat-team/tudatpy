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
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/interface/spice/spiceInterface.h"

/*!
 * Tool to read out files that are structured with tracking measurements in rows and data types in columns
 * A TrackingDataType is a unique form of data that can be used by Tudat to derive the observables
 * A TrackingFileFieldConverter is an interface that can convert a raw string input field to its associated double value as expected for the TrackingDataType
 * Reading a trackingFile can be done by specifying a list of trackingFileFields, each corresponding to such a converter.
 */

namespace tudat
{
namespace input_output
{

/*!
 * Enum describing a unique data type that can later be used to process the information.
 * These are in SI units
 */
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
  observation_body, // In case observations corrected for body center.
  observed_body, // In case observations corrected for body center.
  spacecraft_id,
  spacecraft_name,
  planet_nr,
  tdb_reception_time_j2000,
  utc_reception_time_j2000,
  utc_ramp_referencee_j2000,
  tdb_spacecraft_j2000,
  x_planet_frame,
  y_planet_frame,
  z_planet_frame,
  vx_planet_frame,
  vy_planet_frame,
  vz_planet_frame,
  residual_de405,
  spacecraft_transponder_delay,
  uplink_frequency,
  downlink_frequency,
  signal_to_noise,
  spectral_max,
  doppler_measured_frequency,
  doppler_averaged_frequency,
  doppler_base_frequency,
  doppler_noise,
  doppler_bandwidth,
  receiving_station_name,
  transmitting_station_name,
  time_tag_delay,
  sample_number,
  utc_day_of_year,
  reference_body_distance,
  transmission_frequency_constant_term,
  transmission_frequency_linear_term,
  doppler_predicted_frequency_hz,
  doppler_troposphere_correction,
  scan_nr
};
/*!
 * Simple converter class that can convert a string data field to a double.
 * One can inherit from this and overload the `toDouble()` method to extend the supported formats
 */
class TrackingFileFieldConverter
{
public:

  /*!
   * Constructor
   * @param trackingDataType the SI data type of the double representation
   */
  explicit TrackingFileFieldConverter(TrackingDataType trackingDataType)
      : doubleDataType_(trackingDataType) {}

  //! Destructor
  virtual ~TrackingFileFieldConverter() = default;

  /*!
   * Default toDouble implementation.
   * @param rawField string input as read from the file
   * @return value converted to double
   */
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

  //! Getter for doubleDataType_
  const TrackingDataType& getTrackingDataType() { return doubleDataType_; }

private:
  //! Data type of the double that this converter will return
  TrackingDataType doubleDataType_;
};

//! A converter specifically for month fields in three-letter representation ("jan", "Jan", ...)
class TrackingFileMonthFieldConverter : public TrackingFileFieldConverter
{
public:

  /*!
   * Constructor
   * @param trackingDataType Data type of the double representation (SI)
   */
  TrackingFileMonthFieldConverter(TrackingDataType trackingDataType) : TrackingFileFieldConverter(trackingDataType) {}

  /*!
   * Implementation to convert the months in three letter string to their double representation
   * @param rawField string representing a month (3 letter in lower or upper case)
   * @return double representing a whole number between 1 and 12
   */
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

    return utilities::upperCaseFromMap(rawField, monthsMap);
  }
};

//! Converter that will convert the raw string to double and then apply a scalar multiplier as specified
class TrackingFileFieldMultiplyingConverter : public TrackingFileFieldConverter
{
public:

  /*!
   * Constructor
   * @param trackingDataType Data type of the double representation (SI)
   * @param multiplier The multiplier applied to the direct conversion from string to double. e.g. this can be a scaling factor to convert units
   */
  TrackingFileFieldMultiplyingConverter(TrackingDataType trackingDataType, double multiplier)
      : TrackingFileFieldConverter(trackingDataType), multiplier_(multiplier) {}

  /*!
   * Convert string to double and apply given multiplication factor.
   * @param rawField string
   * @return double representation
   */
  double toDouble(std::string& rawField) const
  {
    return multiplier_ * TrackingFileFieldConverter::toDouble(rawField);
  }
private:

  double multiplier_;
};

//! Converter that will convert the raw string to a ephemeris time representation (TDB) since J2000.
class TrackingFileFieldUTCTimeConverter : public TrackingFileFieldConverter
{
public:
  TrackingFileFieldUTCTimeConverter(TrackingDataType trackingDataType = TrackingDataType::utc_reception_time_j2000)
      : TrackingFileFieldConverter(trackingDataType) {}
  double toDouble(std::string& rawField) const
  {
      return basic_astrodynamics::dateTimeFromIsoString( rawField ).epoch< double >( );
  }
};

//! Mapping the `TrackingFileField` to the correct converter, including the `TrackingDataType` it will represent
static const std::map<std::string, std::shared_ptr<TrackingFileFieldConverter>> trackingFileFieldConverterMap = {
    {"spacecraft_id", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::spacecraft_id)},
    {"dsn_transmitting_station_nr", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::dsn_transmitting_station_nr)},
    {"dsn_receiving_station_nr", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::dsn_receiving_station_nr)},
    {"year", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::year)},
    {"month", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::month)},
    {"month_three_letter", std::make_shared<TrackingFileMonthFieldConverter>(TrackingDataType::month)},
    {"day", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::day)},
    {"hour", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::hour)},
    {"minute", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::minute)},
    {"second", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::second)},
    {
        "time_tag_delay_microseconds",
        std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::time_tag_delay, 1.e-6)
    },
    {"round_trip_light_time_seconds", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::n_way_light_time)},
    {"round_trip_light_time_microseconds", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::n_way_light_time, 1.e-6)},
    {
        "light_time_measurement_delay_microseconds",
        std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::light_time_measurement_delay, 1.e-6)
    },
    {
        "light_time_measurement_accuracy_microseconds",
        std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::light_time_measurement_accuracy, 1.e-6)
    },
    {"planet_nr", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::planet_nr)},
    {"tdb_spacecraft_seconds_j2000", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::tdb_spacecraft_j2000)},
    {"x_planet_frame_km", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::x_planet_frame, 1.e3)},
    {"y_planet_frame_km", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::y_planet_frame, 1.e3)},
    {"z_planet_frame_km", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::z_planet_frame, 1.e3)},
    {"vx_planet_frame_kms", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::vx_planet_frame, 1.e3)},
    {"vy_planet_frame_kms", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::vy_planet_frame, 1.e3)},
    {"vz_planet_frame_kms", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::vz_planet_frame, 1.e3)},
    {"residual_de405_microseconds", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::residual_de405, 1.e-6)},
    {"signal_to_noise_ratio", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::signal_to_noise)},
    {"normalised_spectral_max", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::spectral_max)},
    {"doppler_measured_frequency_hz", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::doppler_measured_frequency)},
    {"doppler_averaged_frequency_hz", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::doppler_averaged_frequency)},
    {"doppler_noise_hz", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::doppler_noise)},
    {"utc_datetime_string", std::make_shared<TrackingFileFieldUTCTimeConverter>(TrackingDataType::utc_reception_time_j2000)},
    {"utc_day_of_year", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::utc_day_of_year)},
    {"tdb_seconds_since_j2000", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::tdb_reception_time_j2000)},
    {"sample_number", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::sample_number)},
    {"reference_body_distance", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::reference_body_distance,1.0E3)},
    {"ramp_reference_time", std::make_shared<TrackingFileFieldUTCTimeConverter>(TrackingDataType::utc_ramp_referencee_j2000)},
    {"transmission_frequency_constant_term", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::transmission_frequency_constant_term)},
    {"transmission_frequency_linear_term", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::transmission_frequency_linear_term)},
    {"doppler_predicted_frequency_hz", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::doppler_predicted_frequency_hz)},
    {"doppler_troposphere_correction", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::doppler_troposphere_correction)},
    {"scan_number", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::scan_nr)}

};

enum TrackingTxtFileReadFilterType
{
    no_tracking_txt_file_filter,
    ifms_tracking_txt_file_filter
};

/*!
 * Class to extract the raw data from a file with the appropriate conversion to doubles. Data fields that do not have an
 * appropriate converter are simply stored as raw strings.
 */
class TrackingTxtFileContents
{
public:

  /*!
   * Constructor
   * @param fileName Path to file that should be read
   * @param columnTypes vector of strings representing the column types. These should be part of `trackingFileFieldConverterMap` to be usefule for tudat
   * Strings that aren't part of that will only be stored as raw fields in this object.
   * @param commentSymbol Lines starting with this symbol are ignored
   * @param valueSeparators string with characters representing separation between columns. ",: " means space , or : mark a new column
   */
  TrackingTxtFileContents(const std::string fileName,
                          const std::vector<std::string> columnTypes,
                          const char commentSymbol = '#',
                          const std::string valueSeparators = ",: \t",
                          const bool ignoreOmittedColumns = false,
                          const TrackingTxtFileReadFilterType dataFilterMethod = no_tracking_txt_file_filter )
      : fileName_(fileName), columnFieldTypes_(columnTypes), commentSymbol_(commentSymbol),
        valueSeparators_(valueSeparators), ignoreOmittedColumns_( ignoreOmittedColumns )
  {
    parseData( dataFilterMethod );
  }

private:

  //! Main parsing sequence to read and process the file
  void parseData( const TrackingTxtFileReadFilterType dataFilterMethod );

  /*!
   * Read out the raw data map from a filestream
   * @param dataFile filestream
   */
  void readRawDataMap(std::ifstream& dataFile);

  /*!
   * Process a single raw line from the file and add it to the data maps
   * @param rawLine single line from the file as a string
   */
  void addLineToRawDataMap(std::string& rawLine);

  bool validateCurrentLineProcessing( const TrackingTxtFileReadFilterType dataFilterMethod, const std::vector<std::string>& rawVector );

    /*!
   * Use the appropriate converters to convert the raw data to the correct double representations
   */
  void convertDataMap( const TrackingTxtFileReadFilterType dataFilterMethod );

// Getters
public:

  /*!
   * Add a trackingdatatype that is constant for the entire file as metadata
   * @param dataType TrackingDataType of interest
   * @param value double or string value that belongs to the trackingDataType.
   */
  void addMetaData(TrackingDataType dataType, double value) { metaDataMapDouble_[dataType] = value; }

  void addMetaData(TrackingDataType dataType, const std::string& value) { metaDataMapStr_[dataType] = value; }

  //! Number of columns expected in the file
  size_t getNumColumns() const { return columnFieldTypes_.size(); }

  //! Number of rows read out in the raw data map
  size_t getNumRows() const { return rawDataMap_.at(columnFieldTypes_[0]).size(); }

  //! Getter for the field types defined by the user
  const std::vector<std::string>& getRawColumnTypes() { return columnFieldTypes_; }

  //! Get the TrackingDataTypes present in the converted data map
  const std::vector<TrackingDataType>& getDataColumnTypes()
  {
    columnDataTypes_.clear();
    for (auto& pair : doubleDataMap_)
    {
      columnDataTypes_.push_back(pair.first);
    }
    return columnDataTypes_;
  }

  //! Getter for `rawDataMap_`
  const auto& getRawDataMap() { return rawDataMap_; }

  //! Getter for `doubleDataMap_`
  const auto& getDoubleDataMap() { return doubleDataMap_; }

  /*!
   * Get the entire converted column for a specific TrackingDataType, or a vector of the correct size with a default val
   * @param dataType TrackingDataType of interest that should be extracted from the datamap
   * @param defaultVal Default double value in case the dataType is not found
   * @return column of data points for the requested data type
   */
  const std::vector<double> getDoubleDataColumn(TrackingDataType dataType, double defaultVal);

  /*!
   * Get the entire converted column for a specific TrackingDataType
   * @param dataType TrackingDataType of interest that should be extracted from the datamap
   * @return column of data points for the requested data type
   */
  const std::vector<double> getDoubleDataColumn(TrackingDataType dataType);

  //! Getter for `metaDataMapDouble_`
  const auto& getMetaDataDoubleMap() { return metaDataMapDouble_; }

  //! Getter for `metaDataMapStr_`
  const auto& getMetaDataStrMap() { return metaDataMapStr_; }

  //! Get a vector of all data types that are present as meta data (valid for the entire file)
  const std::vector<TrackingDataType> getMetaDataTypes();

  //! Get all available data types (either as metadata or individually per row)
  const std::vector<TrackingDataType> getAllAvailableDataTypes();

  void subtractColumnType( const TrackingDataType& columnToSubtractFrom, const TrackingDataType& columnToSubtract );

private:
  //! Path of the file name of interest
  std::string fileName_ = "None";

  //! Vector of strings representing the field types of columns. If known by tudat (`trackingFileFieldConverterMap`), they indicate a converter, otherwise values are stored raw
  std::vector<std::string> columnFieldTypes_;

  //! Vector of data types for the columns that have a double converter
  std::vector<TrackingDataType> columnDataTypes_;

  //! Symbol that marks a comment line in the data file. These lines are ignored
  char commentSymbol_;

  //! String of separator characters that mark a gap between two columns
  std::string valueSeparators_ = ":, \t";

  bool ignoreOmittedColumns_ = false;

  //! Map to link a columnfieldtype (as provided by the user) to a vector of values (read from file)
  std::map<std::string, std::vector<std::string>> rawDataMap_;

  //! Map of meta data (valid for the entire file) with values as doubles
  std::map<TrackingDataType, double> metaDataMapDouble_;

  //! Map of meta data (valid for the entire file) with values as strings
  std::map<TrackingDataType, std::string> metaDataMapStr_;

  //! Map containing the converted double values for the columns where a valid converter was known to tudat
  std::map<TrackingDataType, std::vector<double>> doubleDataMap_;
};

/*!
 * Function to create a read out a tracking data file to raw contents
 * @param fileName
 * @param columnTypes column types (string). If known to Tudat, this will define a tracking data type, otherwise, it is not processed and kept as raw data.
 * @param commentSymbol lines that start with this symbol are ignored
 * @param valueSeparators String of characters that separate columns. E.g. ",:" means that every , and : in the file will create a new column
 * @return TrackingFileContents
 */
static inline std::shared_ptr<TrackingTxtFileContents> createTrackingTxtFileContents(const std::string& fileName,
                                                                                     const std::vector<std::string>& columnTypes,
                                                                                     char commentSymbol = '#',
                                                                                     const std::string& valueSeparators = ",: \t",
                                                                                     const bool ignoreOmittedColumns = false )
{
  return std::make_shared<TrackingTxtFileContents>(fileName, columnTypes, commentSymbol, valueSeparators, ignoreOmittedColumns);
}


inline std::shared_ptr< TrackingTxtFileContents> readIfmsFile(const std::string& fileName, const bool applyTroposphereCorrection = true )
{
    std::vector<std::string>
        columnTypes({"sample_number",
                     "utc_datetime_string",
                     "utc_day_of_year",
                     "tdb_seconds_since_j2000",
                     "reference_body_distance",
                     "ramp_reference_time",
                     "transmission_frequency_constant_term",
                     "transmission_frequency_linear_term",
                     "doppler_averaged_frequency_hz",
                     "doppler_predicted_frequency_hz",
                     "doppler_troposphere_correction",
                     "doppler_noise_hz"});

    auto rawFileContents = createTrackingTxtFileContents(fileName, columnTypes, '#', ", \t",true);
    rawFileContents->addMetaData( TrackingDataType::file_name, fileName );
    if( applyTroposphereCorrection )
    {
        rawFileContents->subtractColumnType( input_output::TrackingDataType::doppler_averaged_frequency, input_output::TrackingDataType::doppler_troposphere_correction );

    }
    return rawFileContents;
}

inline std::shared_ptr< TrackingTxtFileContents> readFdetsFile(
    const std::string& fileName,
    const std::vector<std::string>& columnTypes ={ "utc_datetime_string", "signal_to_noise_ratio", "normalised_spectral_max", "doppler_measured_frequency_hz", "doppler_noise_hz" } )
{
    auto rawFileContents = createTrackingTxtFileContents(fileName, columnTypes, '#', ", \t");
    rawFileContents->addMetaData(TrackingDataType::file_name, fileName );
    return rawFileContents;
}


} // namespace input_output
} // namespace tudat

#endif // TUDAT_READ_GENERIC_TXT_FILE_H
