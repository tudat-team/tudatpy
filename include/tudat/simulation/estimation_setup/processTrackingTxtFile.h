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

using ttdt = tudat::input_output::TrackingDataType;

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

// Additional Helper functions
std::string getStationNameFromStationId(const std::string networkPrefix, const int stationId)
{
  return networkPrefix + std::to_string(stationId);
}

// This map contains all the tracking data types that correspond to an observable
std::map<std::vector<input_output::TrackingDataType>, ObservableType> dataTypeToObservableMap = {
    {{ttdt::two_way_light_time}, n_way_range},
    // Todo: More to be implemented
};

// Get a vector of all the observableTypes that the provided dataTypes contain
std::vector<ObservableType> findAvailableObservableTypes(std::vector<ttdt> availableDataTypes)
{

  std::vector<ObservableType> availableObservableTypes;
  for (const auto& pair : dataTypeToObservableMap) {
    std::vector<ttdt> requiredDataTypeSet = pair.first;
    if (containsAll(availableDataTypes, requiredDataTypeSet))
      availableObservableTypes.push_back(pair.second);
  }
  return availableObservableTypes;
}

// Create ancillary settings
template< typename TimeType = double >
ObservationAncilliarySimulationSettings createOdfAncillarySettings(
    std::shared_ptr<input_output::TrackingTxtFileContents> trackingTxtFileContents)
{
  // Currently not implemented
  return ObservationAncilliarySimulationSettings();
}

// Class containing the processed file contents for a Tracking txt data file
class ProcessedTrackingTxtFileContents
{

public:
  ProcessedTrackingTxtFileContents(std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtContents, std::string spacecraftName)
      : rawTrackingTxtFileContents_(rawTrackingTxtContents), spacecraftName_(spacecraftName)
  {
    initialise();
  }

  void initialise()
  {
    updateObservationTimes();
    updateObservations();
    updateLinkEnds();
    updateObservableTypes();
    initialised_ = true;
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
    linkEndsSet_.clear();
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
              {observation_models::reflector, observation_models::LinkEndId(spacecraftName_, "Antenna")},
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

    // Creating a set with all the various linkEnds
    for (LinkEnds& linkEnds : linkEndsVector_) {
      linkEndsSet_.insert(linkEnds);
    }
  }

  void updateObservations()
  {
    TODO_IMPLEMENT;
    const auto& observableTypes = getObservableTypes();

  }

  void updateObservableTypes()
  {
    observableTypes_ = findAvailableObservableTypes(rawTrackingTxtFileContents_->getDataColumnTypes());
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
    throw std::runtime_error("Error while processing tracking txt file: Link Ends representation not recognised or implemented.");
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
  const std::set<LinkEnds>& getLinkEndsSet()
  {
    return linkEndsSet_;
  }

  // TODO:
  //  Instead of this, process the data in the updateObservations functions!
  const std::map<ttdt, std::vector<double>>& getDoubleDataMap()
  {
    return rawTrackingTxtFileContents_->getDoubleDataMap();
  }

private:
  std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtFileContents_;
  std::string spacecraftName_;
  std::vector<double> observationTimes_;
  std::vector<ObservableType> observableTypes_;
  std::vector<LinkEnds> linkEndsVector_;
  std::set<LinkEnds> linkEndsSet_;
  bool initialised_;

};

// Create observation collection
// TODO: This assumes the ancillary settings are the same for all the observables. For txt files this will usually be the case
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<observation_models::ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFileObservationCollection(
    std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents> processedTrackingTxtFileContents,
    std::vector<ObservableType> observableTypesToProcess = std::vector<ObservableType>(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings(),
    std::pair<TimeType, TimeType> startAndEndTimesToProcess = std::make_pair<TimeType, TimeType>(TUDAT_NAN, TUDAT_NAN))
{

  const auto& dataMap = processedTrackingTxtFileContents->getDoubleDataMap();

  // Check expected sizes
  if (observableTypesToProcess.empty()) {
    observableTypesToProcess = processedTrackingTxtFileContents->getObservableTypes();
  }
  std::set<LinkEnds> linkEndsToProcess = processedTrackingTxtFileContents->getLinkEndsSet();

  std::map<ObservableType, std::map<LinkEnds, std::vector<std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> > > > >
      observationSets;

  std::map<ObservableType, std::map<LinkEnds, std::vector<TimeType>>> observationTimesMap;
  std::map<ObservableType, std::map<LinkEnds, std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> >>> observablesMap;

  // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
  std::vector<TimeType> allObservationTimes = processedTrackingTxtFileContents->getObservationTimes();
  std::vector<LinkEnds> linkEndsVector = processedTrackingTxtFileContents->getLinkEndsVector();
  std::set<LinkEnds> linkEndsSet = processedTrackingTxtFileContents->getLinkEndsSet();

  for (size_t i = 0; i < allObservationTimes.size(); ++i) {
    LinkEnds& currentLinkEnds = linkEndsVector[i];

    for (auto& pair : dataTypeToObservableMap) {
      std::vector<ttdt> currentTrackingDataTypes = pair.first;
      ObservableType& currentObservableType = pair.second;

      if (containsAll(observableTypesToProcess, std::vector<ObservableType>{pair.second})) {
        Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> currentObservable;

        for (size_t j = 0; j < currentTrackingDataTypes.size(); ++j) {
          currentObservable[j] = dataMap.at(currentTrackingDataTypes[j])[i];
        }
        observationTimesMap[currentObservableType][currentLinkEnds].push_back(allObservationTimes[i]);
        observablesMap[currentObservableType][currentLinkEnds].push_back(currentObservable);
      }
    }
  }

  for (ObservableType& currentObservableType : observableTypesToProcess) {
    for (const LinkEnds& currentLinkEnds : linkEndsSet) {
      observationSets[currentObservableType][currentLinkEnds].push_back(std::make_shared<SingleObservationSet<ObservationScalarType, TimeType> >(
          currentObservableType,
          currentLinkEnds,
          observablesMap[currentObservableType][currentLinkEnds],
          observationTimesMap[currentObservableType][currentLinkEnds],
          observation_models::receiver, // Todo: Not sure why this is (copied from odfFiles)
          std::vector<Eigen::VectorXd>(),
          nullptr,
          std::make_shared<ObservationAncilliarySimulationSettings>(ancillarySettings)));
    }
  }

  return std::make_shared<observation_models::ObservationCollection<ObservationScalarType, TimeType> >(observationSets);
}

// TODO:
//  Ancillary settings (frequency, time scale, ...)
//  Format checks
//  Unit Checks
//  Add time conversions



} // namespace observation_models
} // namespace tudat

#endif // TUDAT_PROCESSTRACKINGTXTFILE_H
