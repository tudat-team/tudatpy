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

/*!
 * Map containing all the tracking data types that correspond to an observable. If a specific set of TrackingDataTypes (keys)
 * Is present in the file, the corresponding ObservableType (value) can later be added to the  ObservationCollection
 */
//static const std::map<std::vector<input_output::TrackingDataType>, ObservableType> dataTypeToObservableMap = {
//    {{input_output::TrackingDataType::n_way_light_time}, n_way_range},
//    // {{input_output::TrackingDataType::doppler_measured_frequency}, one_way_doppler},
//    // Todo: More to be implemented
//};
static const std::map<ObservableType, std::vector<input_output::TrackingDataType>> observableRequiredDataTypesMap = {
    {n_way_range, {input_output::TrackingDataType::n_way_light_time}},
    // {one_way_doppler, {input_output::TrackingDataType::doppler_measured_frequency}},
    // Todo: More to be implemented
};


/*!
 * Function to extract a vector of all the observableTypes that can be deduced from a list of DataTypes
 * @param availableDataTypes Vector of available TrackingDataTypes
 * @return vector of observable types that can be derived from the tracking data types
 */
std::vector<ObservableType> findAvailableObservableTypes(std::vector<input_output::TrackingDataType> availableDataTypes);

// TODO: Ancillary settings not used right now
////! Function to create ancillary settings
//template< typename TimeType = double >
//ObservationAncilliarySimulationSettings createTrackingAncillarySettings(
//    std::shared_ptr<input_output::TrackingTxtFileContents> trackingTxtFileContents)
//{
//  return ObservationAncilliarySimulationSettings();
//}

//! Class containing the processed file contents for a Tracking txt data file
class ProcessedTrackingTxtFileContents
{
public:
  /*!
   * Constructor
   * @param rawTrackingTxtContents Raw contents object
   * @param spacecraftName Name used to refer to the spacecraft of interest
   * @param earthFixedGroundStationPositions map of earth fixed positions of the ground stations.
   */
  ProcessedTrackingTxtFileContents(const std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtContents,
                                   const std::string spacecraftName,
                                   const std::map<std::string, Eigen::Vector3d>& earthFixedGroundStationPositions)
      : rawTrackingTxtFileContents_(rawTrackingTxtContents), spacecraftName_(spacecraftName),
        earthFixedGroundStationPositions_(earthFixedGroundStationPositions), initialised_(false)
  {
    initialise();
  }

  //! Main initialisation sequence to create the processed file contents.
  void initialise()
  {
    updateLinkEnds();
    updateObservationTimes();
    updateObservableTypes();
    updateObservations();
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

  //! Method to update the observations
  // TODO: Currently not implemented. This should contain some of the things that are currently done in `createTrackingTxtFileObservationCollection`
  void updateObservations();

  //! Update the available observableTypes described by the processed data file
  void updateObservableTypes()
  {
    observableTypes_ = findAvailableObservableTypes(rawTrackingTxtFileContents_->getDataColumnTypes());
  }

  /*!
   * Convert the UTC observation times to TDB
   * @param observationTimesUtc
   * @return vector of TDB observation times
   */
  std::vector<double> computeObservationTimesTdbFromJ2000(std::vector<double> observationTimesUtc);

  //! Utility function to get the ground station id
  //! TODO: This is temporary until the functionality to read from file is implemented
  static std::string getStationNameFromStationId(const int stationId, const std::string networkPrefix="DSS-")
  {
    return networkPrefix + std::to_string(stationId);
  }

// Settings interpreting the file format.
private:

  //! Enum of implemented time representations. Combination of columns in the processed file contents
  enum TimeRepresentation
  {
    calendar_day_time,
    tdb_seconds_j2000,
  };

  //! Enum of implemented link ends representations. Combination of columns in the processed file contents
  enum LinkEndsRepresentation
  {
    dsn_transmitting_receiving_station_nr,
    vlbi_station,
  };

  //! Deduce the time representation from available columns
  TimeRepresentation getTimeRepresentation();

  //! Deduce the link ends representation from available columns
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

  const std::map<ObservableType, std::vector<double>> getObservationMap() const
  {
    return observationMap_;
  }

private:
  //! TrackingTxtFileContents raw file contents
  std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtFileContents_;

  //! Name of the spacecraft (as chosen by the user)
  std::string spacecraftName_;

  //! Vector of TDB observation times
  std::vector<double> observationTimes_;

  //! Map of observations
  std::map<ObservableType, std::vector<double>> observationMap_;

  //! Vector of observable types that can be derived from the tracking data
  std::vector<ObservableType> observableTypes_;

  //! Vector of linkEnds for every data type
  std::vector<LinkEnds> linkEndsVector_;

  //! Set of distinct linkEnds
  std::set<LinkEnds> linkEndsSet_;

  //! Map with earth fixed positions for all ground stations
  std::map<std::string, Eigen::Vector3d> earthFixedGroundStationPositions_;

  //! Boolean verifying that the initialisation is complete
  bool initialised_ = false;
};

/*!
 * Function to create an observation collection from the processed Tracking file data
 * @param processedTrackingTxtFileContents
 * @param observableTypesToProcess
 * @param earthFixedGroundStationPositions
 * @param ancillarySettings
 * @param startAndEndTimesToProcess
 * @return observation collection
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<observation_models::ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFileObservationCollection(
    const std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents> processedTrackingTxtFileContents,
    std::vector<ObservableType> observableTypesToProcess = std::vector<ObservableType>(),
    const std::map<std::string, Eigen::Vector3d> earthFixedGroundStationPositions = simulation_setup::getApproximateDsnGroundStationPositions(),
    // const std::map<std::string, Eigen::Vector3d> earthFixedGroundStationPositions = simulation_setup::getApproximateGroundStationPositionsFromFile(), // TODO: get ground stations from file
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings(),
    const std::pair<TimeType, TimeType> startAndEndTimesToProcess = std::make_pair<TimeType, TimeType>(TUDAT_NAN, TUDAT_NAN))
{

  // Make sure processing the tracking file was successful
  if (!processedTrackingTxtFileContents->is_initialised()) {
    throw std::runtime_error("Error while processing tracking txt file: processedTrackingTxtFileContents was never initialised.");
  }

  // Get the double map and set of distinct LinkEnds
//  const auto& dataMap = processedTrackingTxtFileContents->getDoubleDataMap(); // FIXME: REMOVE

  const auto& allObservationsMap = processedTrackingTxtFileContents->getObservationMap();
  std::set<LinkEnds> linkEndsToProcess = processedTrackingTxtFileContents->getLinkEndsSet();


  // Check observable types to process. If empty, process all available, if impossible, throw error
  std::vector<ObservableType> availableObservableTypes = processedTrackingTxtFileContents->getObservableTypes();
  if (observableTypesToProcess.empty()) {
    observableTypesToProcess = availableObservableTypes;
  }
  if (!utilities::containsAll(availableObservableTypes, observableTypesToProcess)) {
    throw std::runtime_error("Error while processing Tracking txt file. Not enough information to extract requested observables");
  }

  // Initialise necessary maps
  std::map<ObservableType, std::map<LinkEnds, std::vector<std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> > > > > observationSets;
  std::map<ObservableType, std::map<LinkEnds, std::vector<TimeType>>> observationTimesMap;
  std::map<ObservableType, std::map<LinkEnds, std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> >>> observablesMap;

  // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
  std::vector<TimeType> allObservationTimes = utilities::staticCastVector< TimeType, double >( processedTrackingTxtFileContents->getObservationTimes() );
  std::vector<LinkEnds> linkEndsVector = processedTrackingTxtFileContents->getLinkEndsVector();
  std::set<LinkEnds> linkEndsSet = processedTrackingTxtFileContents->getLinkEndsSet();

  // Prepare maps that order all observations per observable type and link ends
  // This is necessary for files where the linkends are not always the same
  // Loop over all observation times (rows in the file)
  for (size_t i = 0; i < allObservationTimes.size(); ++i) {
    // Get link ends of this observation
    LinkEnds& currentLinkEnds = linkEndsVector[i];

   // Loop over the available observable types
    for (const ObservableType& currentObservableType : observableTypesToProcess)
    {
      // Get vector of appropriate tracking data columns for this observable
      std::vector<input_output::TrackingDataType> currentTrackingDataTypes = observableRequiredDataTypesMap.at(currentObservableType);

      // Prepare a container to store the data columns
      Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> currentObservable(1);
      currentObservable(0) = allObservationsMap.at(currentObservableType)[i];

      // Store in correct maps
      observationTimesMap[currentObservableType][currentLinkEnds].push_back(allObservationTimes[i]);
      observablesMap[currentObservableType][currentLinkEnds].push_back(currentObservable);
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
          receiver, // TODO: Not sure how this is used (copied from odfFiles)
          std::vector<Eigen::VectorXd>(),
          nullptr,
          std::make_shared<ObservationAncilliarySimulationSettings>(ancillarySettings)));
    }
  }

  // Return as shared pointer
  return std::make_shared<ObservationCollection<ObservationScalarType, TimeType> >(observationSets);
}

/*!
 * Function to create an observation collection from the raw Tracking file data
 * @param rawTrackingTxtFileContents The raw tracking file contents
 * @param spacecraftName Name of the spacecraft
 * @param observableTypesToProcess Vector of observable types that need to be in the collection
 * @param earthFixedGroundStationPositions map of ground station positions
 * @param ancillarySettings
 * @param startAndEndTimesToProcess
 * @return observation collection
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<observation_models::ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFileObservationCollection(
    const std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtFileContents,
    const std::string spacecraftName,
    const std::vector<ObservableType> observableTypesToProcess = std::vector<ObservableType>(),
    const std::map<std::string, Eigen::Vector3d> earthFixedGroundStationPositions = simulation_setup::getApproximateDsnGroundStationPositions(),
//    std::map<std::string, Eigen::Vector3d> earthFixedGroundStationPositions = simulation_setup::getApproximateGroundStationPositionsFromFile(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings(),
    std::pair<TimeType, TimeType> startAndEndTimesToProcess = std::make_pair<TimeType, TimeType>(TUDAT_NAN, TUDAT_NAN))
{
  // Create processed tracking file contents
  auto processedTrackingTxtFileContents = std::make_shared<observation_models::ProcessedTrackingTxtFileContents>(rawTrackingTxtFileContents,
                                                                                                                 spacecraftName,
                                                                                                                 earthFixedGroundStationPositions);

  // Create observation collection and return
  return createTrackingTxtFileObservationCollection(processedTrackingTxtFileContents,
                                                    observableTypesToProcess,
                                                    earthFixedGroundStationPositions,
                                                    ancillarySettings,
                                                    startAndEndTimesToProcess);
}

} // namespace observation_models
} // namespace tudat

#endif // TUDAT_PROCESSTRACKINGTXTFILE_H
