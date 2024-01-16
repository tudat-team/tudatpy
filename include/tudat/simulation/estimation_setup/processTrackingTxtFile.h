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

#include "tudat/basics/utilities.h"
#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/simulation/estimation_setup/observations.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/math/interpolators/lookupScheme.h"

namespace tudat
{
namespace observation_models
{

//! Map containing all the tracking data types that correspond to an observable. If a specific set of TrackingDataTypes (keys)
//! Is present in the file, the corresponding ObservableType (value) can later be added to the  ObservationCollection
static const std::map<std::vector<input_output::TrackingDataType>, ObservableType> dataTypeToObservableMap = {
    {{input_output::TrackingDataType::n_way_light_time}, n_way_range},
    {{input_output::TrackingDataType::doppler_measured_frequency}, one_way_doppler},
    // Todo: More to be implemented
};

//! Function to extract a vector of all the observableTypes that can be deduced from a list of DataTypes
std::vector<ObservableType> findAvailableObservableTypes(std::vector<input_output::TrackingDataType> availableDataTypes);

//! Function to create ancillary settings
template< typename TimeType = double >
ObservationAncilliarySimulationSettings createTrackingAncillarySettings(
    std::shared_ptr<input_output::TrackingTxtFileContents> trackingTxtFileContents)
{
  // TODO: Currently not used
  return ObservationAncilliarySimulationSettings();
}

//! Class containing the processed file contents for a Tracking txt data file
class ProcessedTrackingTxtFileContents
{
//  Constructors
public:
  ProcessedTrackingTxtFileContents(std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtContents,
                                   std::string spacecraftName,
                                   const std::map<std::string, Eigen::Vector3d>& earthFixedGroundStationPositions)
      : rawTrackingTxtFileContents_(rawTrackingTxtContents), spacecraftName_(spacecraftName),
        earthFixedGroundStationPositions_(earthFixedGroundStationPositions), initialised_(false)
  {
    initialise();
  }

  void initialise()
  {
    updateLinkEnds(); // Needs to be first
    updateObservationTimes();
    updateObservations();
    updateObservableTypes();
    initialised_ = true;
  }

// Update  information from the raw file data
public:

  //! Method to extract the first and last time of the observations
  std::pair<double, double> getStartAndEndTime()
  {
    const std::vector<double>& times = getObservationTimes();
    std::pair<double, double> startEndTime({times.front(), times.back()});
    return startEndTime;
  }

  //! Method to create a vector of the linkEnds for each of the lines in the file (different ground stations might be used)
  void updateLinkEnds();
  //! Method to (re)calculate the times from the Tracking file
  void updateObservationTimes();

  void updateObservations();

  void updateObservableTypes()
  {
    observableTypes_ = findAvailableObservableTypes(rawTrackingTxtFileContents_->getDataColumnTypes());
  }

  std::vector<double> computeObservationTimesTdbFromJ2000(std::vector<double> observationTimesUtc);

  //! Utility function to get the ground station id
  static std::string getStationNameFromStationId(const std::string networkPrefix, const int stationId)
  {
    return networkPrefix + std::to_string(stationId);
  }

// Settings interpreting the file format.
private:

  enum TimeRepresentation
  {
    calendar_day_time,
    tdb_seconds_j2000,
  };

  enum LinkEndsRepresentation
  {
    dsn_transmitting_receiving_station_nr,
    vlbi_station,
  };

  TimeRepresentation getTimeRepresentation();

  LinkEndsRepresentation getLinkEndsRepresentation();


// Getters
public:
  bool is_initialised() const
  {
    return initialised_;
  }

  const std::vector<ObservableType>& getObservableTypes() const
  {
    return observableTypes_;
  }

  const std::vector<double>& getObservationTimes() const
  {
    return observationTimes_;
  }

  const std::vector<LinkEnds>& getLinkEndsVector() const
  {
    return linkEndsVector_;
  }
  const std::set<LinkEnds>& getLinkEndsSet() const
  {
    return linkEndsSet_;
  }
  int getNumRows() const
  {
    return rawTrackingTxtFileContents_->getNumRows();
  }

  // TODO:
  //  Instead of this, process the data in the updateObservations functions!
  const std::map<input_output::TrackingDataType, std::vector<double>>& getDoubleDataMap()
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
  std::map<std::string, Eigen::Vector3d> earthFixedGroundStationPositions_;
  bool initialised_ = false;
};

//! Function to create an observation collection from the processed Tracking file data
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<observation_models::ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFileObservationCollection(
    std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents> processedTrackingTxtFileContents,
    std::vector<ObservableType> observableTypesToProcess = std::vector<ObservableType>(),
    std::map<std::string, Eigen::Vector3d> earthFixedGroundStationPositions = simulation_setup::getApproximateGroundStationPositionsFromFile(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings(),
    std::pair<TimeType, TimeType> startAndEndTimesToProcess = std::make_pair<TimeType, TimeType>(TUDAT_NAN, TUDAT_NAN))
{

  if (!processedTrackingTxtFileContents->is_initialised()) {
    throw std::runtime_error("Error while processing tracking txt file: processedTrackingTxtFileContents was never initialised.");
  }

  const auto& dataMap = processedTrackingTxtFileContents->getDoubleDataMap();
  std::set<LinkEnds> linkEndsToProcess = processedTrackingTxtFileContents->getLinkEndsSet();

  // Check observable types to process. If empty, process all available, if impossible, throw error
  std::vector<ObservableType> availableObservableTypes = processedTrackingTxtFileContents->getObservableTypes();
  if (observableTypesToProcess.empty()) {
    observableTypesToProcess = availableObservableTypes;
  }
  if (!utilities::containsAll(availableObservableTypes, observableTypesToProcess)) {
    throw std::runtime_error("Error while processing Tracking txt file. Not enough information to extract requested observables");
  }

  std::map<ObservableType, std::map<LinkEnds, std::vector<std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> > > > >
      observationSets;
  std::map<ObservableType, std::map<LinkEnds, std::vector<TimeType>>> observationTimesMap;
  std::map<ObservableType, std::map<LinkEnds, std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> >>> observablesMap;

  // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
  std::vector<TimeType> allObservationTimes = processedTrackingTxtFileContents->getObservationTimes();
  std::vector<LinkEnds> linkEndsVector = processedTrackingTxtFileContents->getLinkEndsVector();
  std::set<LinkEnds> linkEndsSet = processedTrackingTxtFileContents->getLinkEndsSet();

  // Prepare maps that order all observations per observable type and link ends
  // This is necessary for files where the linkends are not always the same
  for (size_t i = 0; i < allObservationTimes.size(); ++i) {
    LinkEnds& currentLinkEnds = linkEndsVector[i];

    for (auto& pair : dataTypeToObservableMap) {
      std::vector<input_output::TrackingDataType> currentTrackingDataTypes = pair.first;
      const ObservableType& currentObservableType = pair.second;

      if (utilities::containsAll(observableTypesToProcess, std::vector<ObservableType>{pair.second})) {
        size_t currentSize = currentTrackingDataTypes.size();
        Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> currentObservable(currentSize);

        for (size_t j = 0; j < currentSize; ++j) {
          currentObservable[j] = dataMap.at(currentTrackingDataTypes[j])[i];
        }
        observationTimesMap[currentObservableType][currentLinkEnds].push_back(allObservationTimes[i]);
        observablesMap[currentObservableType][currentLinkEnds].push_back(currentObservable);

      }
    }
  }

  // Fill the observation collection
  for (ObservableType& currentObservableType : observableTypesToProcess) {
    for (const LinkEnds& currentLinkEnds : linkEndsSet) {
      observationSets[currentObservableType][currentLinkEnds].push_back(std::make_shared<SingleObservationSet<ObservationScalarType, TimeType> >(
          currentObservableType,
          currentLinkEnds,
          observablesMap[currentObservableType][currentLinkEnds],
          observationTimesMap[currentObservableType][currentLinkEnds],
          receiver, // Todo: Not sure how this is used (copied from odfFiles)
          std::vector<Eigen::VectorXd>(),
          nullptr,
          std::make_shared<ObservationAncilliarySimulationSettings>(ancillarySettings)));
    }
  }

  return std::make_shared<ObservationCollection<ObservationScalarType, TimeType> >(observationSets);
}

//! Function to create an observation collection from the raw Tracking file data
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<observation_models::ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFileObservationCollection(
    std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtFileContents,
    std::string spacecraftName,
    std::vector<ObservableType> observableTypesToProcess = std::vector<ObservableType>(),
    std::map<std::string, Eigen::Vector3d> earthFixedGroundStationPositions = simulation_setup::getApproximateGroundStationPositionsFromFile(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings(),
    std::pair<TimeType, TimeType> startAndEndTimesToProcess = std::make_pair<TimeType, TimeType>(TUDAT_NAN, TUDAT_NAN))
{
  auto processedTrackingTxtFileContents = std::make_shared<observation_models::ProcessedTrackingTxtFileContents>(rawTrackingTxtFileContents,
                                                                                                                 spacecraftName,
                                                                                                                 earthFixedGroundStationPositions);
  return createTrackingTxtFileObservationCollection(processedTrackingTxtFileContents,
                                                    observableTypesToProcess,
                                                    earthFixedGroundStationPositions,
                                                    ancillarySettings,
                                                    startAndEndTimesToProcess);
}

} // namespace observation_models
} // namespace tudat

#endif // TUDAT_PROCESSTRACKINGTXTFILE_H
