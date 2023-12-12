/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References: 820-013, TRK-2-18 Tracking System Interfaces Orbit Data File Interface, Revision E, 2008, JPL/DSN
 */


#ifndef TUDAT_PROCESSTRACKINGTXTFILE_H
#define TUDAT_PROCESSTRACKINGTXTFILE_H

// TEMPORARY -- REMOVE THIS
#define TODO_IMPLEMENT std::cerr << __FILE__ << "(" << __LINE__ << "): " << "ERROR! Function not yet implemented\n"

#include "tudat/basics/utilities.h"
#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/io/basicInputOutput.h"
//#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/basic_astro/timeConversions.h"
//#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/simulation/estimation_setup/observations.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/math/interpolators/lookupScheme.h"
//#include "tudat/math/quadrature/trapezoidQuadrature.h"

namespace tudat
{
namespace observation_models
{

// This map specifies what fields are expected in the TrackingTxtFileContents to extract a specific observable
std::map<ObservableType, std::vector<input_output::TrackingDataType>> requiredDataTypesForObservable = {
    {one_way_range, {}},
    {angular_position, {}},
    {position_observable, {}},
    {one_way_doppler, {}},
    {one_way_differenced_range, {}},
    {n_way_range, {}},
    {two_way_doppler, {}},
    {euler_angle_313_observable, {}},
    {velocity_observable, {}},
    {relative_angular_position, {}},
    {n_way_differenced_range, {}},
    {relative_position_observable, {}},
    {dsn_one_way_averaged_doppler, {}},
    {dsn_n_way_averaged_doppler, {}},
};

// Function to check if a container contains all elements of another container
template< typename T, typename U >
bool containsAll(const T& referenceArray, const U searchArray)
{
  return std::includes(referenceArray.begin(), referenceArray.end(), searchArray.begin(), searchArray.end());
}

// Function to apply the convert function to vectors of arguments
template< typename ConvertFunc, typename FirstArg, typename... Args >
std::vector<double> convertVectors(ConvertFunc convertFunc, const std::vector<FirstArg>& firstArg, const std::vector<Args>& ... args)
{
  std::vector<double> result;

  // Assuming all vectors have the same size
  size_t N = firstArg.size();

  for (size_t i = 0; i < N; ++i) {
    result.push_back(convertFunc(firstArg[i], args[i]...));
  }

  return result;
}

// Create ancillary settings
template< typename TimeType = double >
ObservationAncilliarySimulationSettings createOdfAncillarySettings(
    std::shared_ptr<input_output::TrackingTxtFileContents> trackingTxtFileContents)
{
  // Currently not implemented
  return ObservationAncilliarySimulationSettings();
}

// Additional Helper functions
std::string getStationNameFromStationId(const std::string networkPrefix, const int stationId)
{
  return networkPrefix + std::to_string(stationId);
}

// Class containing the processed file contents for a Tracking txt data file
class ProcessedTrackingTxtFileContents
{
  using ttdt = tudat::input_output::TrackingDataType;

public:
  ProcessedTrackingTxtFileContents(std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtContents, std::string spacecraftName)
      : rawTrackingTxtFileContents_(rawTrackingTxtContents), spacecraftName_(spacecraftName)
  {

  }

  void initialise()
  {
    updateObservationTimes();
    updateObservations();
    updateLinkEnds();
    initialised_=true;
  }

  std::pair<double, double> getStartAndEndTime()
  {
    const std::vector<double>& times = getObservationTimes();
    std::pair<double, double> startEndTime({times.front(), times.back()});
    return startEndTime;
  }

//  Update  information from the raw file data
public:

  // TODO: Make sure it is TDB!
  // This will (re)calculate the times from the file and return the vector
  void updateObservationTimes()
  {
    observationTimes_.clear();
    const auto& dataMap = rawTrackingTxtFileContents_->getDoubleDataMap();
    const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();
    TimeRepresentation timeRepresentation = getTimeRepresentation();

    switch (timeRepresentation) {
      case day_time: {
        observationTimes_ = convertVectors(tudat::basic_astrodynamics::convertCalendarDateToJulianDay<double>,
                                           dataMap.at(ttdt::year),
                                           dataMap.at(ttdt::month),
                                           dataMap.at(ttdt::day),
                                           dataMap.at(ttdt::hour),
                                           dataMap.at(ttdt::minute),
                                           dataMap.at(ttdt::second));
        break;
      }
      default: {
        throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
      }
    }
  }

  void updateLinkEnds()
  {
    linkEndsVector_.clear();
    const auto& dataMap = rawTrackingTxtFileContents_->getDoubleDataMap();
    const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();
    LinkEndsRepresentation linkEndsRepresentation = getLinkEndsRepresentation();

    switch (linkEndsRepresentation) {

      // TODO: make a cleaner implementation to allow adding different ways of providing the link ends easily
      case dsn_transmitting_receiving_station_nr: {
        const auto& dsnTransmitterIds = dataMap.at(ttdt::dsn_transmitting_station_nr);
        const auto& dsnReceiverIds = dataMap.at(ttdt::dsn_receiving_station_nr);

        for (size_t i = 0; i < numDataRows; ++i) {
          LinkEnds currentLinkEnds{
              {observation_models::transmitter, observation_models::LinkEndId("Earth", getStationNameFromStationId("DSS", dsnTransmitterIds[i]))},
              {observation_models::reflector, observation_models::LinkEndId(spacecraftName, "Antenna")},
              {observation_models::receiver, observation_models::LinkEndId("Earth", getStationNameFromStationId("DSS", dsnReceiverIds[i]))},
          };
          linkEndsVector_[i] = currentLinkEnds;
        }
        break;
      }
      default: {
        throw std::runtime_error("Error while processing tracking txt file: LinkEnds representation not recognised or implemented.");
      }
    }
  }

  void updateObservations()
  {
    TODO_IMPLEMENT;
    const auto& observableTypes = getObservableTypes();
  }

  void updateObservableTypes()
  {
    observableTypes_.clear();

    TODO_IMPLEMENT;
    observableTypes_ = std::vector<ObservableType>();
  }

// Settings interpreting the file format.
private:

  enum TimeRepresentation
  {
    day_time,
  };

  enum LinkEndsRepresentation
  {
    dsn_transmitting_receiving_station_nr,
  };

  TimeRepresentation getTimeRepresentation()
  {
    auto const& availableDataTypes = rawTrackingTxtFileContents_->getDataColumnTypes();
    if (containsAll(availableDataTypes, std::vector<ttdt>{ttdt::year, ttdt::month, ttdt::day, ttdt::hour, ttdt::minute, ttdt::second})) {
      return day_time;
    }
    throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
  }

  LinkEndsRepresentation getLinkEndsRepresentation()
  {
    auto const& availableDataTypes = rawTrackingTxtFileContents_->getDataColumnTypes();
    if (containsAll(availableDataTypes, std::vector<ttdt>{ttdt::dsn_transmitting_station_nr, ttdt::dsn_receiving_station_nr})) {
      return dsn_transmitting_receiving_station_nr;
    }
  }


// Getters
public:

  const std::vector<ObservableType>& getObservableTypes()
  {
    return observableTypes_;
  }

  const std::vector<double>& getObservationTimes()
  {
    return observationTimes_;
  }

  const std::vector<LinkEnds>& getLinkEndsVector()
  {
    return linkEndsVector_;
  }

private:
  std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtFileContents_;
  std::string spacecraftName_;
  std::vector<double> observationTimes_;
  std::vector<ObservableType> observableTypes_;
  std::vector<LinkEnds> linkEndsVector_;
  bool initialised_;



};

// FIXME: THIS ONE
// Create observation collection
// TODO: This assumes the ancillary settings are the same for all the observables. For txt files this will usually be the case
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFileObservationCollection(
//    std::shared_ptr<input_output::TrackingTxtFileContents> trackingTxtFileContents,
    std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents> processedTrackingTxtFileContents,
    std::map<ObservableType, LinkEnds> observableTypesToProcess = std::map<ObservableType, LinkEnds>(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings(),
    std::pair<TimeType, TimeType> startAndEndTimesToProcess = std::make_pair<TimeType, TimeType>(TUDAT_NAN, TUDAT_NAN))
{

  std::map<ObservableType,
           std::map<LinkEnds,
                    std::vector<std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> > > > >
      observationSets = {};

  //
//
//
////  const auto& availableDataTypes = processedTrackingTxtFileContents->getDataColumnTypes();
////  const auto& dataMap = processedTrackingTxtFileContents->getDoubleDataMap();
////  const auto& numDataRows = processedTrackingTxtFileContents->getNumRows();
//
//
//
//
//  // Check expected sizes
//  if (observableTypesToProcess.empty()) {
//    throw std::runtime_error("Error while processing tracking txt file: Empty list of observableTypes to process");
//  }
//
//
//
//
//
//  // C++ 17 allows looping over two arrays at once.
//  for (auto& pair : observableTypesToProcess) {
//    ObservableType currentObservableType = pair.first;
//    LinkEnds currentLinkEnds = pair.second;
//
//    // Check if the file contains the necessary data types.
//    const auto& requiredDataTypes = requiredDataTypesForObservable[currentObservableType];
//    if (!containsAll(availableDataTypes, requiredDataTypes)) {
//      throw std::runtime_error(
//          "Error while processing tracking txt file: Not all requiredDataTypes are present in the provided Tracking file contents");
//    }
//
//    // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
//    std::vector<TimeType> observationTimes;
//    std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> > observables;
//
////    // TODO: Do this conversion in a separate function (including more options for how the time is stored)
////    auto years = dataMap.at(input_output::TrackingDataType::year);
////    auto months = dataMap.at(input_output::TrackingDataType::month);
////    auto days = dataMap.at(input_output::TrackingDataType::day);
////    auto hours = dataMap.at(input_output::TrackingDataType::hour);
////    auto minutes = dataMap.at(input_output::TrackingDataType::minute);
////    auto seconds = dataMap.at(input_output::TrackingDataType::second);
////    for (size_t i = 0; i < numDataRows; ++i) {
////      observationTimes.push_back(tudat::basic_astrodynamics::convertCalendarDateToJulianDay<TimeType>(years[i],
////                                                                                                      months[i],
////                                                                                                      days[i],
////                                                                                                      hours[i],
////                                                                                                      minutes[i],
////                                                                                                      seconds[i]));
////    }
//
//    switch (currentObservableType) {
//    case n_way_range:
//      observables = std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> >(
//          dataMap.at(input_output::TrackingDataType::two_way_light_time)
//      );
//      break;
//    default:throw std::runtime_error("Error while processing tracking txt file: The requested observable type is not implemented");
//    }
//
//    // FIXME: Fill the values and convert the times
//
//
//    // Look at separateSingleLinkOdfData;
//
//
//    // Todo: Maybe, truncating the observations might be useful. Currently not implemented, because most files until now relatively small
//
//
//    // Create the single observation sets and save them
//    observationSets[currentObservableType][currentLinkEnds].push_back(std::make_shared<SingleObservationSet<ObservationScalarType, TimeType> >(
//        currentObservableType,
//        currentLinkEnds,
//        observables,
//        observationTimes,
//        receiver,
//        std::vector<Eigen::VectorXd>(),
//        nullptr,
//        std::make_shared<ObservationAncilliarySimulationSettings>(ancillarySettings)));
//  }
  return std::make_shared<observation_models::ObservationCollection<ObservationScalarType, TimeType> >(observationSets);
}
//
//// TODO: Add a variation to the function where two vectors are provided instead of a map!
//
//
//
//



























//
//  for (auto observableTypeIterator = trackingTxtFileContents->getProcessedDataBlocks().begin();
//       observableTypeIterator != trackingTxtFileContents->getProcessedDataBlocks().end(); ++observableTypeIterator) {
//    ObservableType currentObservableType = observableTypeIterator->first;
//
//
//
//    // Check if an observation set should be created for the current observable type
//    if (std::count(observableTypesToProcess.begin(), observableTypesToProcess.end(), currentObservableType) == 0) {
//      continue;
//    }
//
//    for (auto linkEndsIterator = observableTypeIterator->second.begin();
//         linkEndsIterator != observableTypeIterator->second.end(); ++linkEndsIterator) {
//      LinkEnds currentLinkEnds = linkEndsIterator->first;
//      std::shared_ptr<ProcessedOdfFileSingleLinkData> currentOdfSingleLinkData = linkEndsIterator->second;
//
//      // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
//      std::vector<std::vector<TimeType> > observationTimes;
//      std::vector<std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> > > observables;
//      std::vector<ObservationAncillarySimulationSettings> ancillarySettings;
//
//      // Fill vectors
//      separateSingleLinkOdfData(
//          currentObservableType, currentOdfSingleLinkData, observationTimes, observables, ancillarySettings);
//
//      std::vector<std::vector<TimeType> > truncatedObservationTimes;
//      std::vector<std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> > > truncatedObservables;
//
//      if (std::isnan(static_cast< double >( startAndEndTimesToProcess.first )) &&
//          std::isnan(static_cast< double >( startAndEndTimesToProcess.second ))) {
//        truncatedObservationTimes = observationTimes;
//        truncatedObservables = observables;
//      } else {
//        for (unsigned int i = 0; i < observationTimes.size(); ++i) {
//          std::vector<TimeType> singleTruncatedObservationTimes;
//          std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> > singleTruncatedObservables;
//
//          for (unsigned int j = 0; j < observationTimes.at(i).size(); ++j) {
//            if (((!std::isnan(static_cast< double >( startAndEndTimesToProcess.first ))
//                && observationTimes.at(i).at(j) >= startAndEndTimesToProcess.first) ||
//                std::isnan(static_cast< double >( startAndEndTimesToProcess.first ))) &&
//                ((!std::isnan(static_cast< double >( startAndEndTimesToProcess.second ))
//                    && observationTimes.at(i).at(j) <= startAndEndTimesToProcess.second) ||
//                    std::isnan(static_cast< double >( startAndEndTimesToProcess.second )))) {
//              singleTruncatedObservationTimes.push_back(observationTimes.at(i).at(j));
//              singleTruncatedObservables.push_back(observables.at(i).at(j));
//            }
//          }
//          truncatedObservationTimes.push_back(singleTruncatedObservationTimes);
//          truncatedObservables.push_back(singleTruncatedObservables);
//        }
//      }
//
//      // Create the single observation sets and save them
//      for (unsigned int i = 0; i < observationTimes.size(); ++i) {
//        sortedObservationSets[currentObservableType][currentLinkEnds].push_back(
//            std::make_shared<SingleObservationSet<ObservationScalarType, TimeType> >(
//                currentObservableType, currentLinkEnds, truncatedObservables.at(i), truncatedObservationTimes.at(i),
//                receiver,
//                std::vector<Eigen::VectorXd>(),
//                nullptr, std::make_shared<ObservationAncillarySimulationSettings>(
//                    ancillarySettings.at(i))));
//      }
//    }
//  }
//
//  return std::make_shared<ObservationCollection<ObservationScalarType, TimeType> >(
//      sortedObservationSets);
//}

} // namespace observation_models
} // namespace tudat

#endif // TUDAT_PROCESSTRACKINGTXTFILE_H
