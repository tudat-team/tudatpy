/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>

#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"

namespace tudat
{
namespace observation_models
{

std::vector<ObservableType> findAvailableObservableTypes(std::vector<input_output::TrackingDataType> availableDataTypes)
{
  std::vector<ObservableType> availableObservableTypes;
  for (const auto& pair : dataTypeToObservableMap) {
    std::vector<input_output::TrackingDataType> requiredDataTypeSet = pair.first;
    if (containsAll(availableDataTypes, requiredDataTypeSet))
      availableObservableTypes.push_back(pair.second);
  }
  return availableObservableTypes;
}

void ProcessedTrackingTxtFileContents::updateObservationTimes()
{
  // TODO: The timescale is currently not considered. Need to find a solution to make that consistent

  observationTimes_.clear();
  const auto& dataMap = rawTrackingTxtFileContents_->getDoubleDataMap();
  const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();
  TimeRepresentation timeRepresentation = getTimeRepresentation();

  switch (timeRepresentation) {
    case day_time: {
      observationTimes_ = convertVectors(convertCalendarDateToJulianDaySinceJ2000<double>,
                                         dataMap.at(input_output::TrackingDataType::year),
                                         dataMap.at(input_output::TrackingDataType::month),
                                         dataMap.at(input_output::TrackingDataType::day),
                                         dataMap.at(input_output::TrackingDataType::hour),
                                         dataMap.at(input_output::TrackingDataType::minute),
                                         dataMap.at(input_output::TrackingDataType::second));
      break;
    }
    default: {
      throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
    }
  }
}

void ProcessedTrackingTxtFileContents::updateLinkEnds()
{
  linkEndsVector_.clear();
  const auto& dataMap = rawTrackingTxtFileContents_->getDoubleDataMap();
  const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();
  LinkEndsRepresentation linkEndsRepresentation = getLinkEndsRepresentation();

  switch (linkEndsRepresentation) {

    // TODO: make a cleaner implementation to allow adding different ways of providing the link ends easily
    case dsn_transmitting_receiving_station_nr: {
      const auto& dsnTransmitterIds = dataMap.at(input_output::TrackingDataType::dsn_transmitting_station_nr);
      const auto& dsnReceiverIds = dataMap.at(input_output::TrackingDataType::dsn_receiving_station_nr);

      for (size_t i = 0; i < numDataRows; ++i) {
        LinkEnds currentLinkEnds{
            {transmitter, LinkEndId("Earth", getStationNameFromStationId("DSS-", static_cast<int>(dsnTransmitterIds[i])))},
            {reflector, LinkEndId(spacecraftName_, "Antenna")},
            {receiver, LinkEndId("Earth", getStationNameFromStationId("DSS-", static_cast<int>(dsnReceiverIds[i])))},
        };
        linkEndsVector_.push_back(currentLinkEnds);
      }
      break;
    }
    default: {
      throw std::runtime_error("Error while processing tracking txt file: LinkEnds representation not recognised or implemented.");
    }
  }

  // Creating a set with all the various linkEnds
  linkEndsSet_ = vectorToSet(linkEndsVector_);
}

void ProcessedTrackingTxtFileContents::updateObservations()
{
  const auto& observableTypes = getObservableTypes();
}

ProcessedTrackingTxtFileContents::TimeRepresentation ProcessedTrackingTxtFileContents::getTimeRepresentation()
{
  auto const& availableDataTypes = rawTrackingTxtFileContents_->getDataColumnTypes();
  if (containsAll(availableDataTypes,
                  std::vector<input_output::TrackingDataType>{
                      input_output::TrackingDataType::year,
                      input_output::TrackingDataType::month,
                      input_output::TrackingDataType::day,
                      input_output::TrackingDataType::hour,
                      input_output::TrackingDataType::minute,
                      input_output::TrackingDataType::second
                  })) {
    return day_time;
  }
  throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
}

ProcessedTrackingTxtFileContents::LinkEndsRepresentation ProcessedTrackingTxtFileContents::getLinkEndsRepresentation()
{
  auto const& availableDataTypes = rawTrackingTxtFileContents_->getDataColumnTypes();
  if (containsAll(availableDataTypes,
                  std::vector<input_output::TrackingDataType>{
                      input_output::TrackingDataType::dsn_transmitting_station_nr,
                      input_output::TrackingDataType::dsn_receiving_station_nr
                  })) {
    return dsn_transmitting_receiving_station_nr;
  }
  throw std::runtime_error("Error while processing tracking txt file: Link Ends representation not recognised or implemented.");
}

} // namespace observation_models
} // namespace tudat
