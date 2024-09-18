/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readTrackingTxtFile.h"

namespace tudat
{
namespace input_output
{

void TrackingTxtFileContents::parseData()
{
  std::ifstream dataFile(fileName_);
  if (!dataFile.good()) {
    throw std::runtime_error("Error when opening Tracking txt file: file " + fileName_ + " could not be opened.");
  }
  readRawDataMap(dataFile);
  convertDataMap();
}

void TrackingTxtFileContents::readRawDataMap(std::ifstream& dataFile)
{
  std::string currentLine;

  while (std::getline(dataFile, currentLine)) {
    if (!currentLine.empty() && currentLine.at(0) != commentSymbol_) {
      addLineToRawDataMap(currentLine);
    }
  }
}

void TrackingTxtFileContents::addLineToRawDataMap(std::string& rawLine)
{

  size_t numColumns = getNumColumns();
  std::vector<std::string> currentSplitRawLine_;

  // Trim the line and split based on the separators
  boost::algorithm::trim(rawLine);
  boost::algorithm::split(currentSplitRawLine_,
                          rawLine,
                          boost::is_any_of(valueSeparators_),
                          boost::algorithm::token_compress_on);

  // Check if the expected number of columns is present in this line
  if (currentSplitRawLine_.size() != numColumns) {
    unsigned int columnsFound = currentSplitRawLine_.size();
    for (auto a : currentSplitRawLine_) {
      std::cout << a << "\n";
    }
    throw std::runtime_error(
        "The current line in file " + fileName_ + " has " + std::to_string(columnsFound) + " columns but "
            + std::to_string(numColumns) + " columns were expected.\nRaw line:" + rawLine);

  }

  // Populate the dataMap_ with a new row on each of the vectors
  for (std::size_t i = 0; i < numColumns; ++i) {
    std::string currentFieldType = columnFieldTypes_.at(i);
    std::string currentValue = currentSplitRawLine_.at(i);
    rawDataMap_[currentFieldType].push_back(currentValue);
  }
}

void TrackingTxtFileContents::convertDataMap()
{
  // Loop over all the column types to convert them to doubles
  for (std::string columnType : columnFieldTypes_) {

    // If the columnType (requested by the user) does not have a known converter in tudat, it is skipped here and only stored as raw string
    if (!trackingFileFieldConverterMap.count(columnType)) {
      std::cout << "Warning: '" << columnType << "' is not recognised as a column type by Tudat. The data is available in the raw format.\n";
      continue;
    }

    // Get the raw vector, correct converter and target data type represented by the doubles
    std::shared_ptr<TrackingFileFieldConverter> converter = trackingFileFieldConverterMap.at(columnType);
    const std::vector<std::string>& rawVector = rawDataMap_.at(columnType);
    const TrackingDataType& dataType = converter->getTrackingDataType();

    // Convert and store double vector
    std::vector<double> dataVector;
    for (std::string rawValue : rawVector) {
      dataVector.push_back(converter->toDouble(rawValue));
    }
    doubleDataMap_[dataType] = dataVector;
  }
}

const std::vector<double> TrackingTxtFileContents::getDoubleDataColumn(TrackingDataType dataType, double defaultVal)
{
  if (doubleDataMap_.find(dataType) != doubleDataMap_.end()) {
    return doubleDataMap_.at(dataType);
  } else if (metaDataMapDouble_.find(dataType) != metaDataMapDouble_.end()) {
    return std::vector<double>(getNumRows(), metaDataMapDouble_[dataType]);
  } else if (!std::isnan(defaultVal)) {
    return std::vector<double>(getNumRows(), defaultVal);
  } else {
    throw std::runtime_error("Error while working with trackinig data. Required datatype not found");
  }
}

const std::vector<double> TrackingTxtFileContents::getDoubleDataColumn(TrackingDataType dataType)
{
  // Call the version with default value, using a default value that triggers an error
  return getDoubleDataColumn(dataType, std::numeric_limits<double>::quiet_NaN());
}

const std::vector<TrackingDataType> TrackingTxtFileContents::getMetaDataTypes()
{
  std::vector<TrackingDataType> metaDataTypes;

  for (const auto& pair : metaDataMapDouble_) {
    metaDataTypes.push_back(pair.first);
  }

  for (const auto& pair : metaDataMapStr_) {
    metaDataTypes.push_back(pair.first);
  }

  return metaDataTypes;
}

const std::vector<TrackingDataType> TrackingTxtFileContents::getAllAvailableDataTypes()
{
  auto columnTypes = getDataColumnTypes();
  auto metaTypes = getMetaDataTypes();
  columnTypes.insert(columnTypes.end(), metaTypes.begin(), metaTypes.end());
  return columnTypes;
}

} // namespace input_output

} // namespace tudat
