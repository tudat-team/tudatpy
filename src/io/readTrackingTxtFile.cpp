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
    throw std::runtime_error("Error when opening Jpl Range file, file " + fileName_ + " could not be opened.");
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

  // Trim the line and split based on the separators
  boost::algorithm::trim(rawLine);
  boost::algorithm::split(currentSplitRawLine_,
                          rawLine,
                          boost::is_any_of(valueSeparators_),
                          boost::algorithm::token_compress_on);

  // Check if the expected number of columns is present in this line
  if (currentSplitRawLine_.size() != numColumns) {
    unsigned int columnsFound = currentSplitRawLine_.size();
    throw std::runtime_error(
        "The current line in file " + fileName_ + " has " + std::to_string(columnsFound) + " columns but "
            + std::to_string(numColumns) + " columns were expected.\nRaw line:" + rawLine);
  }

  // Populate the dataMap_ with a new row on each of the vectors
  for (std::size_t i = 0; i < numColumns; ++i) {
    TrackingFileField currentFieldType = columnFieldTypes_.at(i);
    std::string currentValue = currentSplitRawLine_.at(i);
    rawDataMap_[currentFieldType].push_back(currentValue);
  }
}

void TrackingTxtFileContents::convertDataMap()
{
  for (TrackingFileField columnType : columnFieldTypes_) {
    const std::vector<std::string>& rawVector = rawDataMap_.at(columnType);
    std::shared_ptr<TrackingFileFieldConverter> converter = trackingFileFieldConverterMap.at(columnType);

    const TrackingDataType& dataType = converter->getTrackingDataType();

    std::vector<double> dataVector;
    for (std::string rawValue : rawVector) {
      dataVector.push_back(converter->toDouble(rawValue));
    }
    doubleDataMap_[dataType] = dataVector;
  }
}

void TrackingTxtFileContents::addMetaData(TrackingFileField fieldType, const std::string& value)
{
  metaDataMap_[fieldType] = value;
}

} // namespace input_output

} // namespace tudat
